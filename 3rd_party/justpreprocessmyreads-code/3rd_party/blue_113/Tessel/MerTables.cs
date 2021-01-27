using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;
using System.Threading;
using CommonCode;

namespace MerCollections
{
    public class MerTables
    {
        // allow 2 status bits (1 base) at the end of each singleton mer
        public const ulong singletonRCFlagMask = 0x0000000000000001;        // rc flag is lowest order bit
        public const ulong singletonPromotedFlagMask = 0x0000000000000001;  // same bit re-used to resolve races during promotion
        public const ulong singletonActiveFlagMask = 0x0000000000000002;    // active/promoted flag is 2nd bit (1==active)
        public const ulong singletonMerMask = 0xfffffffffffffffc;           // mask out the last 2 bits (status bits in singleton filter) for comparisons
        public const ulong resetActiveSingletonMask = 0xfffffffffffffffd;   // just clear the active bit 
        public const int int31Mask = 0x7fffffff;                            // all but the top bit of an int 32
        public const ulong fullMerMask = 0xffffffffffffffff;                // all bits in the mer are used in comparisons

        public int noOfPartitions = 0;                              // no. of mer partitions (hash-distributed mers)
        public int noSingletonPartitions = 0;                       // base-prefix ordered entries
        public int singletonPrefixBits = 0;                         // how many LHS bits are used for partitioning singletons
        int[] maxSingletonCapacity;                                 // max no. of entries allowed in each singleton table before flushing
        int maxTableSize = 40000000;                                // ensure table partitions are not too big (used in initial partition sizing)
        int minTableSize = 5000000;                                 // nor too small
        int maxSingletonSize = 40000000;
        int minSingletonSize = 5000000;
        const long minKeepDepth = 3;                                // only keep k-mers of at least this depth in the repeated kmers table at the end-of-file flush

        public MerDictionary[] repeatedMers = null;                 // hash-partitioned (closed) hash tables holding repeat mers
                                                                    // these tables hold canonical (lowest of as-read/RC) mers
                                                                    // shared amongst all threads - insertion should be lock free. no deletions or resizes done
        public bool[] repeatedMersFull = null;                      // repeatedMers partition is full, new mers go in per-thread overflow tables
        public MerDictionary[] overflowMers = null;                 // per-thread overflow repeat mer tables

        public LowRepMerBuffer culledBuffer = null;                 // temporary shared structures used during multi-threaded low-rep mer flush
        public object culledLock = new object();                    // (only ever done if there are more than 2 seq files being tiled - e.g. multi-lane datasets)
        public List<LowRepMerBuffer> filledCulledBuffers = null;    // only used if culledBuffer becomes full

        public MerCollection[] singletonFilters = null;             // partitioned hash sets holding possible singleton mers
                                                                    // these tables contain partitioned mers with an RC status bit at the RHS of the ulong
                                                                    // [0] says whether the mer is the RC form (true) or the as-is (false) form (used to calculate initial counts)
        public List<MerCollection>[] singletonsWaitingFlush = null; // (probably empty) list of singleton buffers waiting to be flushed. Used to consolidate flush files
        public List<int>[] singletonsWaitingFlushCounts;            // (how many singletons are in each of the waiting collections)
        object lockSingletons = new object();                       // ensure clean replacement/flush of singletons array
        public List<string>[] flushedSingletonFNs = null;           // per-partition list of file names of flushed singleton partitions
        public List<ulong>[] firstFlushedSingletonMer = null;       // lowest mer in each flush singleton file
        public List<ulong>[] lastFlushedSingletonMer = null;        // highest mer in each flush singleton file
        public int[] flushSingletonNumber = null;                   // next flush number for each partition (used to construct unique file names)
        public List<string> flushedLowRepsFNs = null;               // list of file names of flushed LowReps partitions
        public List<ulong> firstFlushedLowRepsMer = null;           // lowest mer in each flush LowReps file
        public List<ulong> lastFlushedLowRepsMer = null;            // highest mer in each flush LowReps file
        public int[] flushLowRepsNumber = null;                     // next flush number for each partition (used to construct unique file names)
        string tempDirectory = "";                                  // where the temporary flush files are to be written

        const bool tracing = false;
        public Queue<TraceEntry> traceUpdates = new Queue<TraceEntry>(1000000);
        const int maxTrace = 1000000;

        // constructor
        public MerTables(long dictionarySize, string tempDir, int noThreads)
        {
            // scale genome size to compensate for the number of repeated error k-mers
            dictionarySize = dictionarySize * 2;

            // how many shared mer partitions are needed to safely hold this many distinct k-mers?
            this.noOfPartitions = (int)(dictionarySize / maxTableSize + 1);
            if (this.noOfPartitions < 1)
                this.noOfPartitions = 1;
            // and how big should the partitions be? 
            int partitionSize = (int)(dictionarySize / noOfPartitions);
            if (partitionSize < minTableSize)
                partitionSize = minTableSize;

            // how many singleton partitions are desirable?
            // per-partition singleton tables can use of a considerable amount of memory, but too few will increase the number of concurrently open files during merge
            if (noOfPartitions > 64)
            {
                noSingletonPartitions = 32;
                singletonPrefixBits = 5;
            }
            if (noOfPartitions <= 64)
            {
                noSingletonPartitions = 16;
                singletonPrefixBits = 4;
            }
            if (noOfPartitions <= 16)
            {
                noSingletonPartitions = 4;
                singletonPrefixBits = 2;
            }
            if (noOfPartitions <= 4)
            {
                noSingletonPartitions = 1;
                singletonPrefixBits = 0;
            }

            this.tempDirectory = tempDir;

            repeatedMers = new MerDictionary[noOfPartitions];               // create partitioned dictionaries
            repeatedMersFull = new bool[noOfPartitions];                    // create full flags array (default is false)
            overflowMers = new MerDictionary[noThreads];                    // create per-thread overflow tables
            singletonFilters = new MerCollection[noSingletonPartitions];    // create partitioned singleton filters
            singletonsWaitingFlush = new List<MerCollection>[noSingletonPartitions];
            singletonsWaitingFlushCounts = new List<int>[noSingletonPartitions];
            flushedSingletonFNs = new List<string>[noSingletonPartitions];
            firstFlushedSingletonMer = new List<ulong>[noSingletonPartitions];
            lastFlushedSingletonMer = new List<ulong>[noSingletonPartitions];
            flushSingletonNumber = new int[noSingletonPartitions];
            maxSingletonCapacity = new int[noSingletonPartitions];

            flushedLowRepsFNs = new List<string>();
            firstFlushedLowRepsMer = new List<ulong>();
            lastFlushedLowRepsMer = new List<ulong>();

            // initialise per-partition structures
            for (int i = 0; i < noOfPartitions; i++)
            {
                repeatedMers[i] = new MerDictionary(partitionSize, fullMerMask);
            }

            // initialise per-singleton-partition structures
            for (int i = 0; i < noSingletonPartitions; i++)
            {
                int scaledSingletonSize = 3 * partitionSize / noSingletonPartitions * (noSingletonPartitions - i);
                scaledSingletonSize = Math.Min(scaledSingletonSize, maxSingletonSize);
                scaledSingletonSize = Math.Max(scaledSingletonSize, minSingletonSize);
                singletonFilters[i] = new MerCollection(scaledSingletonSize, singletonMerMask);
                maxSingletonCapacity[i] = singletonFilters[i].length * 9 / 10;
                flushedSingletonFNs[i] = new List<string>();
                flushSingletonNumber[i] = 1;
                firstFlushedSingletonMer[i] = new List<ulong>();
                lastFlushedSingletonMer[i] = new List<ulong>();
            }

            // initialise per-thread structures
            for (int i = 0; i < noThreads; i++)
            {
                // overflowMers[i] = new MerDictionary(singletonSize); // allocated when first used to save on memory space
            }
        }

        public void AddOrIncrement(ulong mer, int threadNo)
        {
            long addingIncrement = 0x0000000100000000;                  // assume we've got the as-read form is the canonical form
            ulong rcFlagToBeSet = 0x0;                                  // and that we don't want to set the RC flag

            // generate canonical k-mer first
            ulong rcMer = MerStrings.ReverseComplement(mer);
            if (rcMer < mer)
            {
                mer = rcMer;
                addingIncrement = 0x0000000000000001;                   // increment the low part of the count pair
                rcFlagToBeSet = singletonRCFlagMask;                             // remember if the canonical k-mer was the RC form
            }

            int absMerHashCode = mer.GetHashCode() & int31Mask;
            int partitionNo = absMerHashCode % noOfPartitions;
            int singletonPartitionNo = singletonPrefixBits == 0 ? 0 : (int)(mer >> (64 - singletonPrefixBits));

            // this mer may have been seen before, so first try updating it in one of the repeated mer tables
            bool updatedRepeat = UpdateRepeatedMer(partitionNo, mer, threadNo, mer, addingIncrement);

            if (updatedRepeat)
                return;

            // handling a k-mer for the first time - try adding it to the singletons table
            // ----------------------------------------------------------------------------

            // get a stable pointer to the current singetons table (in case someone else fills it and initiates a flush while we're still busy with it)
            MerCollection thisSingletonPartition = singletonFilters[singletonPartitionNo];
            Interlocked.Increment(ref thisSingletonPartition.activeCount);

            // try to add this mer to this partition's singletons collection (and fetch the existing singleton+flag if it's already there)
            int filterIdx;
            ulong fMer = mer | rcFlagToBeSet | singletonActiveFlagMask;
            bool added = thisSingletonPartition.TryInsertKey(fMer, out filterIdx);

            if (added)
            {
                // successfully added this mer so we must be seeing it for the first time 

                // if singleton table is already full enough, flush it out and empty the table
                if (thisSingletonPartition.Count >= maxSingletonCapacity[singletonPartitionNo])
                {
                    bool flushNeeded = true;
                    int flushNumberToUse = 0;

                    // lock this section to avoid two threads trying to flush/replace the same singleton buffer concurrently
                    lock (lockSingletons)
                    {
                        // test entry condition now that we have the lock (filter may have been reset while we were waiting)
                        if (!thisSingletonPartition.flushed)
                        {
                            // allocate a replacement table for the other threads to use while we're flushing this one
                            int newSingletonLength = thisSingletonPartition.length + thisSingletonPartition.length / 4;
                            if (newSingletonLength > maxSingletonSize)
                                newSingletonLength = maxSingletonSize;
                            MerCollection emptySingletonFilter = new MerCollection(newSingletonLength, singletonMerMask);     // allocate new local filter for the partition

                            singletonFilters[singletonPartitionNo] = emptySingletonFilter;                  // make it visible to the concurrent threads (single point assignment)
                            maxSingletonCapacity[singletonPartitionNo] = newSingletonLength * 8 / 10;
                            thisSingletonPartition.flushed = true;
                            flushNumberToUse = flushSingletonNumber[singletonPartitionNo];
                            flushSingletonNumber[singletonPartitionNo]++;
                        }
                        else
                            flushNeeded = false;
                    }

                    if (flushNeeded)
                    {
                        while (thisSingletonPartition.activeCount > 1)
                        {
                            // pause briefly to let any inflight updates to this table to complete
                            Thread.Sleep(100);
                        }
                        FlushSingletons(thisSingletonPartition, singletonPartitionNo, flushNumberToUse);
                    }
                    //flushes++;
                }
            }
            else
            {
                // Insert failed, so must be seeing this k-mer for second (or rarely more) time. Mark as inactive in singletons and add to a repeats table with appropriate counts. 
                // There can be a race here with two threads trying to concurrently promote the same singleton. This is resolved by atomically clearing the singleton
                // active flag - and only one of the threads will get the 'active' flag returned from the Exchange. This thread does the promotion - and then sets the 
                // promotion-complete bit for the singleton. The other threads will spin until they find this bit has been set.

                if (tracing)
                    lock (traceUpdates)
                    {
                        traceUpdates.Enqueue(new TraceEntry(threadNo, 1, singletonPartitionNo, filterIdx, (ulong)thisSingletonPartition.entries[filterIdx].key));
                        if (traceUpdates.Count > maxTrace)
                        {
                            traceUpdates.Dequeue();
                        }
                    }

                // get the current value of this singleton entry (safe because the promotion changes are progressive)
                ulong merFromFilter = (ulong)thisSingletonPartition.entries[filterIdx].key;
                // and see if this singleton may have already been promoted
                bool activeSingleton = (merFromFilter & singletonActiveFlagMask) != 0;

                // if this singleton may be 'active', try to promote it
                if (activeSingleton)
                {
                    ulong inactiveMer = mer & singletonMerMask;                      // build what the inactive-but-being-promoted entry should look like
                    // if no-one else has altered the singleton entry, then set it to inactive-but-being-promoted
                    long currentMerFromFilter = Interlocked.CompareExchange(ref thisSingletonPartition.entries[filterIdx].key, (long)inactiveMer, (long)merFromFilter);

                    if (tracing)
                        lock (traceUpdates)
                        {
                            traceUpdates.Enqueue(new TraceEntry(threadNo, 2, singletonPartitionNo, filterIdx, (ulong)currentMerFromFilter));
                            if (traceUpdates.Count > maxTrace)
                            {
                                traceUpdates.Dequeue();
                            }
                        }

                    // if this thread successfully set the singleton to 'inactive', it will take care of the promotion
                    if (currentMerFromFilter == (long)merFromFilter)
                    {
                        ulong rcFlag = merFromFilter & singletonRCFlagMask;          // non-zero --> RC found in singletons 

                        long initialCount = 0;
                        if (rcFlag != 0)	                                // singleton was seen in RC form
                            initialCount = 0x0000000000000001;
                        else	                                            // singleton was seen in as-is form
                            initialCount = 0x0000000100000000;

                        if (repeatedMersFull[partitionNo])
                        {
                            if (overflowMers[threadNo] == null)
                            {
                                overflowMers[threadNo] = new MerDictionary(repeatedMers[partitionNo].lengthEntries / 10, fullMerMask);
                                //Console.WriteLine("added overflow for thread " + threadNo + " for [" + partitionNo + "]");
                            }

                            bool full = overflowMers[threadNo].Add(mer, initialCount);
                            if (full)
                                overflowMers[threadNo].Resize();
                        }
                        else
                        {
                            bool full = repeatedMers[partitionNo].Add(mer, initialCount);
                            if (full)
                                repeatedMersFull[partitionNo] = true;
                        }

                        // now that the mer has been promoted, set the 'promoted' flag
                        inactiveMer = inactiveMer | (long)singletonPromotedFlagMask;
                        thisSingletonPartition.entries[filterIdx].key = (long)inactiveMer;

                        if (tracing)
                            lock (traceUpdates)
                            {
                                traceUpdates.Enqueue(new TraceEntry(threadNo, 3, singletonPartitionNo, filterIdx, (ulong)thisSingletonPartition.entries[filterIdx].key));
                                if (traceUpdates.Count > maxTrace)
                                {
                                    traceUpdates.Dequeue();
                                }
                            }
                    }
                }

                // singleton is now known to be no longer active, so wait (if necessary) for the 'promoted' flag to be set and increment the repeat counter

                merFromFilter = (ulong)thisSingletonPartition.entries[filterIdx].key;

                if (tracing)
                    lock (traceUpdates)
                    {
                        traceUpdates.Enqueue(new TraceEntry(threadNo, 4, singletonPartitionNo, filterIdx, merFromFilter));
                        if (traceUpdates.Count > maxTrace)
                        {
                            traceUpdates.Dequeue();
                        }
                    }

                bool promotionComplete = (merFromFilter & singletonPromotedFlagMask) != 0;
                bool alreadySlept = false;
                while (!promotionComplete)
                {
                    promotionComplete = (((ulong)thisSingletonPartition.entries[filterIdx].key & singletonPromotedFlagMask) != 0);
                    if (alreadySlept && !promotionComplete)
                    {
                        if (tracing)
                        {
                            lock (traceUpdates)
                            {
                                StreamWriter trace = new StreamWriter("trace.txt");
                                foreach (TraceEntry t in traceUpdates)
                                    trace.WriteLine(t.place + "\t" + t.thread + "\t" + t.partition + "\t" + t.index + "\t" + t.value.ToString("x16"));
                                trace.Close();
                            }
                            Console.WriteLine("promotion still not complete after sleep");
                        }
                    }
                    if (!promotionComplete)
                        Thread.Sleep(100);
                    alreadySlept = true;
                }

                UpdateRepeatedMerAfterPromotion(partitionNo, mer, threadNo, mer, addingIncrement);
                //if (!updateSucceeded)
                //{
                //    lock (traceUpdates)
                //    {
                //        StreamWriter trace = new StreamWriter("trace.txt");
                //        foreach (TraceEntry t in traceUpdates)
                //            trace.WriteLine(t.thread + "\t" + t.place + "\t" + t.partition + "\t" + t.index + "\t" + t.value.ToString("x16"));
                //        trace.Close();
                //    }
                //    Console.WriteLine("UpdateRepeatedMerRetry failed after waiting for promotion to complete");
                //}

            }

            Interlocked.Decrement(ref thisSingletonPartition.activeCount);
        }

        private bool UpdateRepeatedMer(int partitionNo, ulong pMer, int threadNo, ulong mer, long addingIncrement)
        {
            bool updated = false;

            // update the count for this mer if it's in the shared repeated set 
            updated = repeatedMers[partitionNo].UpdateIfPresent(pMer, addingIncrement);
            if (updated)
            {
                //repeats++;
                return true;
            }

            // mer could be in the thread-local overflow table
            updated = overflowMers[threadNo] != null && overflowMers[threadNo].UpdateIfPresent(mer, addingIncrement);
            if (updated)
            {
                //overflow++;
                return true;
            }

            return false;
        }

        private void UpdateRepeatedMerAfterPromotion(int partitionNo, ulong pMer, int threadNo, ulong mer, long addingIncrement)
        {
            bool updated = false;

            // update the count for this mer if it's in the shared repeated set 
            updated = repeatedMers[partitionNo].UpdateIfPresent(pMer, addingIncrement);
            if (updated)
            {
                //repeats++;
                return;
            }

            // mer could be in the thread-local overflow table (or it may have been promoted by another thread)
            updated = overflowMers[threadNo] != null && overflowMers[threadNo].UpdateIfPresent(mer, addingIncrement);
            if (updated)
            {
                //overflow++;
                return;
            }

            // if we get here, the promoted singleton must have gone to another thread's overflow table (and the partition must be full)
            // so create a new overflow table if we don't already have one
            if (repeatedMersFull[partitionNo] && overflowMers[threadNo] == null)
            {
                overflowMers[threadNo] = new MerDictionary(repeatedMers[partitionNo].lengthEntries / 10, fullMerMask); // set length properly *****
                Console.WriteLine("added overflow for thread " + threadNo);
            }

            // finally add this promoted singleton to the overflow table for this thread 
            if (overflowMers[threadNo] != null)
            {
                overflowMers[threadNo].Add(mer, addingIncrement);
                return;
            }

            // the shared repeated table isn't full but we didn't find the singleton there, could have been a race on a bucket that resulted in an orphan
            repeatedMers[partitionNo].Add(pMer, addingIncrement);

        }

        public void FlushDeferredSingletons(int partitionNo, int flushNo)
        {
            FlushSingletons(null, partitionNo, flushNo);
        }

        private void FlushSingletons(MerCollection singletons, int partitionNo, int flushNo)
        {
            int singletonCapacity = 0;
            int singletonsCount = 0;
            int merIdx = 0;
            int inactiveSingletons = 0;
            ulong[] sortedSingletons = null;

            if (singletons != null)
            {
                // get this now before the Sort clears away the HashTable structures
                singletonCapacity = singletons.Capacity;
                // extract the actual singleton mers from the MerCollection Entries 
                singletonsCount = singletons.Sort();
                // get a pointer to this singletons array
                sortedSingletons = singletons.sortedKeys;

                // compress the singletons array in-place by skipping over inactive entries (shift active towards the start of the array)

                for (int s = 0; s < singletonsCount; s++)
                {
                    ulong mer = sortedSingletons[s];
                    if ((mer & singletonActiveFlagMask) != 0)
                    {
                        sortedSingletons[merIdx] = mer;
                        merIdx++;
                    }
                    else
                        inactiveSingletons++;
                }
            }

            // if we're writing out few active singletons, save this load for the next flush
            if (inactiveSingletons > singletonsCount / 2)
            {
                if (singletonsWaitingFlush[partitionNo] == null)
                {
                    singletonsWaitingFlush[partitionNo] = new List<MerCollection>();
                    singletonsWaitingFlushCounts[partitionNo] = new List<int>();
                }

                singletonsWaitingFlush[partitionNo].Add(singletons);
                singletonsWaitingFlushCounts[partitionNo].Add(merIdx);

                //Console.WriteLine("deferring singleton write[" + partitionNo + "] " + merIdx + " mers");
            }
            else
            {
                // got a singleton set worth flushing so add any pending singleton sets to it (or we're finishing up)
                if (singletonsWaitingFlush[partitionNo] != null && singletonsWaitingFlush[partitionNo].Count != 0)
                {
                    int startingWaiter = 0;

                    // find out how many singletons are waiting
                    int sumWaiting = 0;
                    foreach (int s in singletonsWaitingFlushCounts[partitionNo])
                        sumWaiting += s;

                    // if we're just flushing at end then turn the first waiting set into the base set for the merge
                    if (singletons == null)
                    {
                        sortedSingletons = singletonsWaitingFlush[partitionNo][0].sortedKeys;
                        singletonCapacity = sortedSingletons.Length;
                        startingWaiter = 1;
                        merIdx = singletonsWaitingFlushCounts[partitionNo][0];
                        sumWaiting -= merIdx;
                    }

                    // ensure the buffer is big enough
                    if (merIdx + sumWaiting > singletonCapacity)
                        Array.Resize<ulong>(ref sortedSingletons, merIdx + sumWaiting);

                    // copy the waiting singletons into this buffer
                    for (int i = startingWaiter; i < singletonsWaitingFlush[partitionNo].Count; i++)
                    {
                        int singlesInArray = singletonsWaitingFlushCounts[partitionNo][i];
                        ulong[] pendingSingles = singletonsWaitingFlush[partitionNo][i].sortedKeys;
                        for (int s = 0; s < singlesInArray; s++)
                        {
                            sortedSingletons[merIdx] = pendingSingles[s];
                            merIdx++;
                        }
                    }

                    singletonsWaitingFlush[partitionNo].Clear();
                    singletonsWaitingFlushCounts[partitionNo].Clear();

                    // and re-sort all the singletons
                    Array.Sort<ulong>(sortedSingletons, 0, merIdx);
                }

                if (merIdx > 0)
                {
                    // and write out the (possibly augmented) set of singletons
                    //Console.WriteLine("wrote  " + merIdx + "/" + singletonsCount + " singletons for flush " + partitionNo + "-" + flushNo);

                    // Add process ID at the beginning of the file name
                    int processID = Process.GetCurrentProcess().Id;
                    string binaryfileName = tempDirectory + processID + "_" + partitionNo + "_" + flushNo + ".bst";

                    // Binary writer
                    BinaryWriter flushWriter = new BinaryWriter(File.Open(binaryfileName, FileMode.Create, FileAccess.Write));

                    // Write the number of k-mers in the flush file at the start
                    flushWriter.Write(merIdx);
                    // Write out the in-use singletons
                    for (int i = 0; i < merIdx; i++)
                        flushWriter.Write(sortedSingletons[i]);

                    flushWriter.Close();

                    flushedSingletonFNs[partitionNo].Add(binaryfileName);                                   // Add name of file to list of flush files for the partition
                    firstFlushedSingletonMer[partitionNo].Add(sortedSingletons[0] & singletonMerMask);          // remember lowest and highest mer in each flush file 
                    lastFlushedSingletonMer[partitionNo].Add(sortedSingletons[merIdx - 1] & singletonMerMask);  // (without the RC bit)

                    if (singletons != null)
                        singletons.sortedKeys = null;
                }
            }
        }

        // flush the low-rep mers from the repeat tables, condense the remaining repeated mers and fold in the per-thread repeats. Can only be called after all the
        // threads have finished for a seq data file. This code is *not* thread-safe.
        public void FlushLowRepMers(MerTables merTable, int fileNo)
        {
            // allocate a buffer to hold the flushed low-rep mers
            //int initialBufferLength = 500000;
            int initialBufferLength = this.repeatedMers[0].Capacity;
            culledBuffer = new LowRepMerBuffer();
            culledBuffer.keys = new ulong[initialBufferLength + noOfPartitions];
            culledBuffer.values = new long[initialBufferLength + noOfPartitions];
            culledBuffer.idx = 0;
            culledBuffer.bufferActive = true;
            culledBuffer.bufferNo = 1;
            culledBuffer.limit = initialBufferLength;
            culledLock = new object();

            FlushingThreadParams[] flushingParams = new FlushingThreadParams[noOfPartitions];
            Thread[] flushingThreads = new Thread[noOfPartitions];

            for (int p = 0; p < noOfPartitions; p++)
            {
                flushingParams[p] = new FlushingThreadParams();
                flushingParams[p].merTable = merTable;
                flushingParams[p].partitionNo = p;
                flushingThreads[p] = new Thread(new ParameterizedThreadStart(MerTables.FlushLowRepMersInPartition));
                flushingThreads[p].Priority = ThreadPriority.BelowNormal;
                flushingThreads[p].Start(flushingParams[p]);
            }

            for (int p = 0; p < noOfPartitions; p++)
            {
                flushingThreads[p].Join();
                flushingThreads[p] = null;
            }

            // write out any filled culled buffers
            int bufferNo = 0;
            if (filledCulledBuffers != null)
            {
                for (int i = 0; i < filledCulledBuffers.Count; i++)
                {
                    WriteLowRepMers(fileNo, bufferNo, filledCulledBuffers[i], filledCulledBuffers[i].keys.Length);
                    bufferNo++;
                    filledCulledBuffers[i] = null;
                }
                filledCulledBuffers = null;
            }
            // finally write out the remaining culled low-rep mers
            WriteLowRepMers(fileNo, bufferNo, culledBuffer, culledBuffer.idx);

            // return the temporary buffers
            culledBuffer = null;

            // finally push the per-thread dictionaries to the shared dictionary
            for (int t = 0; t < overflowMers.Length; t++)
            {
                if (overflowMers[t] == null)
                    continue;

                MerDictionary currentOverflow = overflowMers[t];
                MerDictionary replacementOverflow = new MerDictionary(currentOverflow.Capacity, fullMerMask);

                foreach (KeyValuePair<ulong, long> kvp in currentOverflow)
                {
                    int absMerHashCode = kvp.Key.GetHashCode() & int31Mask;
                    int partitionNo = absMerHashCode % noOfPartitions;

                    if (repeatedMersFull[partitionNo])
                        replacementOverflow.Add(kvp.Key, kvp.Value);
                    else
                    {
                        bool OK = repeatedMers[partitionNo].Add(kvp.Key, kvp.Value);
                        if (!OK)
                            repeatedMersFull[partitionNo] = true;
                    }
                }

                overflowMers[t] = replacementOverflow;
            }

        }

        private static void FlushLowRepMersInPartition(object threadParams)
        {
            FlushingThreadParams flushingParams = (FlushingThreadParams)threadParams;
            MerTables merTable = flushingParams.merTable;
            int partitionNo = flushingParams.partitionNo;
            int culledMers = merTable.repeatedMers[partitionNo].Reduce(minKeepDepth, merTable);
            // and it can't possibly be full if we've just culled some entries for it
            if (culledMers > 0)
                merTable.repeatedMersFull[partitionNo] = false;
            //Console.WriteLine("flushed " + culledMers + " from partition " + partitionNo + ". " + merTable.repeatedMers[partitionNo].Count + " left");
        }

        private void WriteLowRepMers(int fileNo, int bufferNo, LowRepMerBuffer culledBuffer, int noMers)
        {
            ulong[] culledBufferKeys = culledBuffer.keys;
            long[] culledBufferValues = culledBuffer.values;

            //for (int i = 0; i < noMers; i++)
            //    if (culledBufferValues[i] == 0)
            //        Debugger.Break();

            // sort the singletons (only active)
            Array.Sort<ulong, long>(culledBufferKeys, culledBufferValues, 0, noMers);

            // Add process ID at the beginning of the file name
            int processID = Process.GetCurrentProcess().Id;
            string binaryfileName = tempDirectory + processID + "_" + fileNo + "_" + bufferNo + ".bfm";

            // Binary writer
            BinaryWriter flushWriter = new BinaryWriter(File.Open(binaryfileName, FileMode.Create, FileAccess.Write));

            // Write the number of k-mers in the flush file at the start
            flushWriter.Write(noMers);
            // Write out the in-use singletons
            for (int i = 0; i < noMers; i++)
            {
                flushWriter.Write(culledBufferKeys[i]);
                flushWriter.Write(culledBufferValues[i]);
            }

            flushWriter.Close();

            flushedLowRepsFNs.Add(binaryfileName);                       // Add name of file to list of flush files for the partition
            firstFlushedLowRepsMer.Add(culledBufferKeys[0]);             // remember lowest and highest mer in each flush file 
            lastFlushedLowRepsMer.Add(culledBufferKeys[noMers - 1]);

            //Console.WriteLine("flushed " + noMers + " low rep mers (" + fileNo + "_" + bufferNo + ")");
        }

    }

    public class TraceEntry
    {
        public int thread;
        public int place;
        public int partition;
        public int index;
        public ulong value;

        public TraceEntry(int thread, int place, int partition, int index, ulong value)
        {
            this.thread = thread;
            this.place = place;
            this.partition = partition;
            this.index = index;
            this.value = value;
        }
    }

    class FlushingThreadParams
    {
        public MerTables merTable;
        public int partitionNo;
    }

    public class LowRepMerBuffer
    {
        public int idx;
        public int limit;
        public bool bufferActive;
        public int bufferNo;
        public ulong[] keys;
        public long[] values;
    }
}