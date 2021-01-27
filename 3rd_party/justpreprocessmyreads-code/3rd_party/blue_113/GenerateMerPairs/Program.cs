using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading;
using System.Diagnostics;
using CommonCode;

namespace GenerateMerPairs
{
    // Generates a k-mer pairs file from a set of reads and a .cbt file.
    //
    // usage: generateMerPairs readsPattern cbtFN
    //
    // The resulting pairs file will be have the same names but the .cbt will be replaced with .prs. 
    // This file starts with an int32 (gap) and this is followed by a series of (pair, count) values (ulong, int).
    // Each 'pair' is a pair of 16-base mer 'stubs' packed into a ulong. Each stub is the first 16bp of - the first parts of a pair of k-mers from a read. 
    // Both of these k-mers must not be singletons (or close) - checked against the .cbt table.
    // The 'count' is an average of the repetition depth for each of these k-mers. 
    //

    class Program
    {
        const string version = "1.1.2";

        const int merStubSize = 16;         // DO NOT CHANGE THIS VALUE. Later code assumes that there will be two 16-base mers making up a 32-base ulong
        const ulong merStubMask = 0xffffffff00000000;   // used to extract merStubs from full k-mers
        const int endGuard = 16;            // don't get too close to the error-prone end
        const int minGap = 16;              // don't let the gap get too small
        const int maxGap = 64;              // and don't let the gap be too big

        static int gap = 0;                 // (constant) gap between mers
        static int merSize = 0;             // size of the mers in the .cbt file
        static int stepSize = 2;            // gaps between pairs (1 for small reads, 2 for larger)
        static int pairStride = 0;          // how big is the effective k-mer 2*stub + gap

        const int maxReadSize = 1000;       // max length for a read
        const int batchSize = 1000;         // how many reads are read at once by a thread
        const int bufferSize = 10000;       // how many pairs are read to start each writing batch
        const int defaultReadLength = 300;
        const int defaultHeaderLength = 100;

        static long progressReadsProcessed = 0;
        static long totalReadsProcessed = 0;
        static long totalReadsRead = 0;
        static long totalPairsWritten = 0;
        static long totalPairsGenerated = 0;
        static long totalDeepUnbalancedReads = 0;

        //static ulong lastWritten = 0;

        // progress monitors
        const int reportInterval = 60000;               // report back at this interval (ms)
        static bool stopMonitor = false;                // tell monitor thread to finish
        static bool mergingPhase = false;               // either counting or merging

        static MerTable<ulong> uniqueMers = null;       // the collection of unique mers from tiled reads (after discarding low-rep ones)
        static int averageDepth;                        // average depth of the loaded tiled mers

        static MerCollections.MerTables uniquePairs;    // shared Mer dictionary (ulong, long) - thread safe as there are no deletions or resizes

        static EventWaitHandle[] threadFinished;        // set when each worker thread finishes their task
        static StreamReader[] readsFiles;               // the set of open reads files
        static int readsFormat;

        static void Main(string[] args)
        {
            if (args.Length < 2)
            {
                Console.WriteLine("usage: GenerateMerPairs [-m min] [-t threads] cbtFN readsPattern or file names (" + version + ")");
                return;
            }

            List<string> FNParams = new List<string>();     // the .cbt name and the set of file names or patterns
            int noThreads = 1;                              // no. of healing threads to run in parallel (1 thread is default)
            int minLoadReps = 3;                            // min rep count needed before mer will be loaded into uniqueMers table or saved as a pair

            for (int p = 0; p < args.Length; p++)
            {
                if (args[p][0] == '-')
                {
                    args[p] = args[p].ToLower();

                    if (args[p] == "-m" || args[p] == "-min")
                    {
                        if (!CheckForParamValue(p, args.Length, "minReps number expected after -m|-min"))
                            return;
                        try
                        {
                            minLoadReps = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -m|-min parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-t" || args[p] == "-threads")
                    {
                        if (!CheckForParamValue(p, args.Length, "number expected after -t|-threads"))
                            return;
                        try
                        {
                            noThreads = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -t|-threads parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    Console.WriteLine("unrecognised option: " + args[p]);
                    Console.WriteLine("usage: generateMerPairs [-m min] [-t threads] cbtFN readsPattern or file names (" + version + ")");
                    return;
                }

                FNParams.Add(args[p]);
            }

            if (FNParams.Count < 2)
            {
                Console.WriteLine("expected a cbt file name and at least one reads file name or pattern");
                return;
            }

            // take the cbt file name from the start of the non-option list
            string cbtFN = FNParams[0];
            FNParams.RemoveAt(0);

            if (FNParams.Count == 0)
            {
                Console.WriteLine("did not find any reads file names or patterns");
                return;
            }

            string pairsFN = cbtFN.Replace(".cbt", ".prs");

            List<string> readsFileNames = new List<string>(FNParams.Count);
            List<string> readsFilePaths = new List<string>(FNParams.Count);
            foreach (string readsFNP in FNParams)
            {
                string readsFileName;
                string readsFilePath;
                GetPathFN(readsFNP, out readsFilePath, out readsFileName);
                readsFilePaths.Add(readsFilePath);
                readsFileNames.Add(readsFileName);
            }

            List<string> expandedReadsFNs = new List<string>();
            for (int f = 0; f < FNParams.Count; f++)
            {
                string[] matchedReadsFNs = Directory.GetFiles(readsFilePaths[f], readsFileNames[f], SearchOption.TopDirectoryOnly);
                foreach (string matchedReadsFN in matchedReadsFNs)
                    expandedReadsFNs.Add(matchedReadsFN);
            }

            // make sure there aren't any duplicates in the file list (seems to be a bug on the Cherax SGI HPC system and it returns each file name twice)
            List<string> distinctReadsFNs = new List<string>();
            foreach (string fn in expandedReadsFNs)
                if (!distinctReadsFNs.Contains(fn))
                    distinctReadsFNs.Add(fn);

            // finally... the set of fully qualified, distinct reads files
            string[] readsFNs;
            readsFNs = distinctReadsFNs.ToArray();

            if (readsFNs.Length == 0)
            {
                Console.WriteLine("No matching read files found");
                return;
            }

            int noOfReadsFiles = readsFNs.Length;
            readsFiles = new StreamReader[noOfReadsFiles];
            for (int f = 0; f < noOfReadsFiles; f++)
            {
                string readsFN = readsFNs[f];
                readsFiles[f] = new StreamReader(readsFN);
            }

            // look at the first file to determine the file format and possible read length
            StreamReader testReader = new StreamReader(readsFNs[0]);
            char headerChar = (char)testReader.Peek();
            if (headerChar == '>')
                readsFormat = MerStrings.formatFNA;
            if (headerChar == '@')
                readsFormat = MerStrings.formatFASTQ;
            int readLength = 0;
            for (int i = 0; i < 20; i++)
            {
                string nextRead = MerStrings.ReadRead(testReader, readsFormat);
                if (nextRead == null)
                    break;
                int nextLength = nextRead.Length;
                if (nextLength > readLength)
                    readLength = nextLength;
            }
            testReader.Close();

            // have to able to fit at least two full mers into the read (no overlaps)
            if (readLength < 2 * merSize)
            {
                Console.WriteLine("reads too short to generate pairs: " + readLength);
                return;
            }

            if (!File.Exists(cbtFN))
            {
                Console.WriteLine(".cbt file not found: " + cbtFN);
                return;
            }

            //string knownPairsFN = "C.sporogenesRaw_25_Copy_1.prs";
            //BinaryReader knownPairs = new BinaryReader(File.Open(knownPairsFN, FileMode.Open, FileAccess.Read));

            //knownPairs.ReadInt32();

            //while (true)
            //{
            //    ulong mer = 0;
            //    int count = 0;

            //    try
            //    {
            //        mer = knownPairs.ReadUInt64();
            //        count = knownPairs.ReadInt32();

            //        goodPairs.Add(mer, count);
            //    }
            //    catch
            //    {
            //        break;
            //    }
            //}

            //knownPairs.Close();
            //Console.WriteLine("loaded " + goodPairs.Count + " good mers from " + knownPairsFN);


            long loadedUniqueMers = 0;
            long loadedTotalMers = 0;

            // load the .cbt file into a merTable (either a hash table (small) or a sorted array (large))
            MerStrings.LoadCBTFile(cbtFN, minLoadReps, 0, 0, minLoadReps,
                                   out uniqueMers, out merSize, out averageDepth, out loadedUniqueMers, out loadedTotalMers);

            if (merSize < merStubSize)
            {
                Console.WriteLine("mers in .cbt file are shorter than merStub size: " + merSize + " < " + merStubSize);
                return;
            }

            uniquePairs = new MerCollections.MerTables(loadedUniqueMers, noThreads);

            // calculate a gap size based on the first read
            gap = (readLength - endGuard)/2 - (merStubSize * 2);
            if (gap < minGap)
                gap = minGap;
            if (gap > maxGap)
                gap = maxGap;

            pairStride = merStubSize + gap + merStubSize;

            // start the monitor/synchronising thread
            Thread monitorProgress = new Thread(RateReporter);
            monitorProgress.Priority = ThreadPriority.AboveNormal;
            monitorProgress.Start();

            DateTime pairingStart = DateTime.Now;

            foreach (string readsFN in readsFNs)
            {
                Console.WriteLine("Generating pairs from " + readsFN);
                StreamReader reads = new StreamReader(readsFN, Encoding.ASCII, false, 1000000);
                BufferedReader bufferedReads = new BufferedReader(readsFormat, reads, null);

                threadFinished = new EventWaitHandle[noThreads];
                int threadNo = 0;
                for (int i = 0; i < noThreads; i++)
                    threadFinished[i] = new EventWaitHandle(false, EventResetMode.AutoReset);

                for (int t = 0; t < noThreads; t++)
                {
                    threadParams workerParam = new threadParams();
                    workerParam.threadNo = threadNo;
                    workerParam.bufferedReadsFile = bufferedReads;
                    ThreadPool.QueueUserWorkItem(new WaitCallback(PairWorker), workerParam);
                    threadNo++;
                }
                //  and wait for them all to finish
                for (int t = 0; t < noThreads; t++)
                {
                    threadFinished[t].WaitOne();
                }
            }

            BinaryWriter pairsFile = new BinaryWriter(File.Open(pairsFN, FileMode.Create, FileAccess.Write));
            pairsFile.Write(gap);

            for (int pi = 0; pi < uniquePairs.noOfPartitions; pi++)
                totalPairsGenerated += uniquePairs.repeatedMers[pi].Sort();

            for (int ti = 0; ti < noThreads; ti++)
                if (uniquePairs.overflowMers[ti] != null)
                    totalPairsGenerated += uniquePairs.overflowMers[ti].Sort();

            MergeAndWrite(pairsFile, uniquePairs.repeatedMers, uniquePairs.overflowMers);

            pairsFile.Close();

            StopMonitorThread(monitorProgress);

            //Console.WriteLine(totalDeepUnbalancedReads + " deep unbalanced reads");
            //Console.WriteLine(totalReadsProcessed + " reads processed");
            Console.WriteLine("wrote " + totalPairsWritten + " pairs from " + totalReadsRead + " reads in " + (DateTime.Now - pairingStart).TotalSeconds.ToString("#.0") + "s");
        }

        private static void GetPathFN(string readsFN, out string readsPath, out string readsFNP)
        {
            char FSC = Path.DirectorySeparatorChar;
            string FSS = new string(FSC, 1);
            readsPath = null;
            if (readsFN.Contains(FSS))
            {
                readsPath = readsFN.Substring(0, readsFN.LastIndexOf(FSC));
                readsFNP = readsFN.Substring(readsFN.LastIndexOf(FSC) + 1);
            }
            else
            {
                readsPath = Directory.GetCurrentDirectory();
                readsFNP = readsFN;
            }
        }

        private static bool CheckForParamValue(int p, int argsLength, string msg)
        {
            if (p == argsLength)
            {
                Console.WriteLine(msg);
                return false;
            }
            return true;
        }

        static void PairWorker(object param)
        {
            threadParams threadParam = (threadParams)param;
            int threadNo = (int)threadParam.threadNo;
            BufferedReader readsFile = threadParam.bufferedReadsFile;
            bool EOF = false;
            Sequence[] readHeaderBatch = new Sequence[batchSize];
            Sequence[] readBatch = new Sequence[batchSize];
            for (int i = 0; i < batchSize; i++)
            {
                readHeaderBatch[i] = new Sequence(defaultHeaderLength);
                readBatch[i] = new Sequence(defaultReadLength);
            }
            int readsInBatch = 0;
            long threadReadsRead = 0;
            long threadReadsProcessed = 0;

            ulong[] mersFromRead = new ulong[1000];
            bool[] merValid = new bool[1000];
            ulong[] canonicalMersFromRead = new ulong[1000];
            int[] plusDepths = new int[1000];
            int[] rcDepths = new int[1000];
            bool deepUnbalanced = false;
            long threadDeepUnbalancedCount = 0;

            int minDepth = averageDepth / 20;

            while (!EOF)
            {
                lock (readsFile)
                {
                    readsInBatch = readsFile.ReadReads(batchSize, readHeaderBatch, readBatch, null, null);

                    if (readsInBatch != batchSize)
                        EOF = true;

                    threadReadsRead += readsInBatch;
                }

                progressReadsProcessed += readsInBatch;

                for (int r = 0; r < readsInBatch; r++)
                {
                    threadReadsProcessed++;

                    Sequence read = readBatch[r];
                    int readLength = read.Length;

                    if (readLength < 2 * merSize)
                        continue;

                    if (readLength < 200)
                        stepSize = 1;
                    else
                        stepSize = 2;

                    //string target = "GTATATAATAAAGTTTTTTATAAAATTTTAAAAGATCATTATAAAAATATAATAACAATTAATATAATATTAATATACTTTAGTTATAGCTATAAATCTTT";
                    //if (read.ToString() == target)
                    //    Debugger.Break();

                    int merCount = MerStrings.GenerateMersFromRead(read, merSize, ref mersFromRead, ref merValid);

                    for (int i = 0; i < merCount; i++)
                        if (merValid[i])
                        {
                            ulong rcMer = MerStrings.ReverseComplement(mersFromRead[i], merSize);
                            if (rcMer < mersFromRead[i])
                                canonicalMersFromRead[i] = rcMer;
                            else
                                canonicalMersFromRead[i] = mersFromRead[i];
                        }

                    GetDepthsForRead(merCount, mersFromRead, canonicalMersFromRead,  merValid, plusDepths, rcDepths, minDepth, out deepUnbalanced);

                    if (deepUnbalanced)
                    {
                        threadDeepUnbalancedCount++;
                        continue;
                    }

                    ulong pair;
                    int pairDepth;
                    bool gotPair;
                    int startingM = 0;
                    int lastM = read.Length - pairStride; // generate pairs up to the end of the read (used to only generate from first part)

                    while (startingM < lastM)
                    {
                        if (merValid[startingM])
                        {
                            gotPair = GeneratePairFromRead(mersFromRead, merValid, plusDepths, rcDepths, startingM, merCount, minDepth, out pair, out pairDepth);

                            if (gotPair)
                            {
                                ulong rcPair = MerStrings.ReverseComplement(pair, 32);
                                if (rcPair < pair)
                                    pair = rcPair;

                                //if (pair == 0x054A0985B90B34D1)
                                //    Debugger.Break();

                                uniquePairs.AddIfNotPresent(pair, pairDepth, threadNo);

                                //lock (pairDictionary)
                                //{
                                //    if (!pairDictionary.ContainsKey(pair))
                                //        pairDictionary.Add(pair, pairDepth);
                                //}

                                //Interlocked.Increment(ref GPTrue);
                                //gotPairFromRead = true;
                            }
                            //else
                                //Interlocked.Increment(ref GPfalse);
                        }

                        startingM += stepSize;
                    }

                    //if (!gotPairFromRead)
                    //    threadReadsWithNoPairs++;
                }
            }

            Interlocked.Add(ref totalReadsProcessed, threadReadsProcessed);
            Interlocked.Add(ref totalReadsRead, threadReadsRead);
            Interlocked.Add(ref totalDeepUnbalancedReads, threadDeepUnbalancedCount);

            threadFinished[threadNo].Set();
        }

        private static void GetDepthsForRead(int merCount, ulong[] mersFromRead, ulong[] canonicalMers, bool[] merValid, int[] plusDepths, int[] rcDepths, int minDepth, out bool deepUnbalanced)
        {
            int sumPlus = 0;
            int sumRC = 0;
            int nonZeroMers = 0;

            for (int m = 0; m < merCount; m++)
            {
                int plusDepth = 0;
                int rcDepth = 0;

                if (merValid[m])
                {
                    ulong mer = mersFromRead[m];
                    ulong canonicalMer = canonicalMers[m];
                    ulong countPair = 0;
                    bool rcMerWasCanonical = (mer != canonicalMer);

                    //Interlocked.Add(ref summedGDMers, (long)canonicalMer);

                    if (uniqueMers.TryGetValue(canonicalMer, out countPair))
                    {
                        // extract the plus, RC and qual values from the packed ulong value
                        if (rcMerWasCanonical)
                        {
                            rcDepth = (int)(countPair >> 32);
                            plusDepth = (int)(countPair & 0xFFFFFFFF);
                        }
                        else
                        {
                            plusDepth = (int)(countPair >> 32);
                            rcDepth = (int)(countPair & 0xFFFFFFFF);
                        }
                    }
                    //Interlocked.Add(ref summedCountPairs, (long)countPair);

                    if (uniqueMers.ContainsKey(canonicalMer))
                    {
                        countPair = uniqueMers[canonicalMer];
                        //Interlocked.Add(ref summedCountPairsViaKey, (long)countPair);
                    }

                }

                plusDepths[m] = plusDepth;
                rcDepths[m] = rcDepth;

                sumPlus += plusDepth;
                sumRC += rcDepth;
                if (plusDepth > minDepth || rcDepth > minDepth)
                    nonZeroMers++;
            }

            int averageNonZeroDepth = nonZeroMers > 0 ? (sumPlus + sumRC) / nonZeroMers : 0;
            deepUnbalanced = (averageNonZeroDepth > averageDepth) && (sumPlus > 1000 * sumRC) | (sumRC > 1000 * sumPlus);
        }

        private static bool GeneratePairFromRead(ulong[] mersFromRead, bool[] merValid, int[] plusDepths, int[] rcDepths, int m, int merCount, int minDepth, out ulong pair, out int pairDepth)
        {
            pair = 0;
            pairDepth = 0;
            ulong firstMer;
            ulong secondMer;
            int firstMerDepthFwd = 0;
            int secondMerDepthFwd = 0;
            int firstMerDepthRC = 0;
            int secondMerDepthRC = 0;

            // generate the two k-mers (separated by a 'gap' gap) -- mer was checked for validity before calling this method
            firstMer = mersFromRead[m];
            int secondM = m + merStubSize + gap + merStubSize - merSize;
            if (!merValid[secondM])
                return false;
            if (secondM >= merCount)
                return false;
            secondMer = mersFromRead[secondM];
            firstMerDepthFwd = plusDepths[m];
            firstMerDepthRC = rcDepths[m];
            secondMerDepthFwd = plusDepths[secondM];
            secondMerDepthRC = rcDepths[secondM];

            // look up both k-mers in the mers table. Very low repetition mers are not loaded into the table, so they will be regarded as 'not found'
            if (firstMerDepthFwd == 0 && firstMerDepthRC == 0)
            {
                //if (tracingNow)
                //    trace.WriteLine("couldn't find first mer: " + MerStrings.ExpandMer(firstMer));
                return false;
            }
            if (secondMerDepthFwd == 0 && secondMerDepthRC == 0)
            {
                //if (tracingNow)
                //    trace.WriteLine("couldn't find second mer: " + MerStrings.ExpandMer(secondMer));
                return false;
            }

            // Don't create a pair based on an effectively zero mer (a low-rep error mer derived from a high rep good mer)
            int firstTotalDepth = firstMerDepthFwd + firstMerDepthRC;
            int secondTotalDepth = secondMerDepthFwd + secondMerDepthRC;

            //// first check to reduce number of unnecessary alternative checks
            //if (firstTotalDepth * 10 < secondTotalDepth || secondTotalDepth * 10 < firstTotalDepth)
            //{
            //    ulong merToCheck = firstMer;
            //    int depthToCheck = firstTotalDepth;
            //    if (secondTotalDepth < firstTotalDepth)
            //    {
            //        merToCheck = secondMer;
            //        depthToCheck = secondTotalDepth;
            //    }

            //    int merFill = 64 - merSize * 2;
            //    ulong rshiftedMer = merToCheck >> merFill & 0xFFFFFFFFFFFFFFFC;
            //    int sumVariantCounts = 0;
            //    // generate all variants of this mer, look them up in the table and sum their counts
            //    for (ulong b = 0; b <= 3; b++)
            //    {
            //        ulong variantMer = (rshiftedMer | b) << merFill;
            //        ulong rcVariantMer = MerStrings.ReverseComplement(variantMer, merSize);
            //        if (rcVariantMer < variantMer)
            //            variantMer = rcVariantMer;
            //        ulong variantCountPair = 0;
            //        if (uniqueMers.TryGetValue(variantMer, out variantCountPair))
            //            sumVariantCounts += (int)(variantCountPair >> 32) + (int)(variantCountPair & 0xFFFFFFFF);
            //    }
            //    // and reject the current mer if it's counts are tiny relative to the sum for the variants
            //    if (((depthToCheck * 100) / sumVariantCounts) <= 5)
            //        return false;
            //}

            if ((firstTotalDepth * 20 < secondTotalDepth && firstTotalDepth < minDepth) || (secondTotalDepth * 20 < firstTotalDepth && secondTotalDepth < minDepth))
                return false;

            // construct the pair ulong
            pair = (firstMer & merStubMask) | ((secondMer >> (64 - merSize * 2)) & ~merStubMask);
            //if (pair == 0x054A0985B90B34D1)
            //    Debugger.Break();
            // calculate the average depth for the pair of k-mers
            pairDepth = firstMerDepthFwd + firstMerDepthRC + secondMerDepthFwd + secondMerDepthRC;
            pairDepth = pairDepth / 2;

            //if (firstTotalDepth > (secondTotalDepth * 20))
            //    Debugger.Break();

            return true;
        }

        private static void MergeAndWrite(BinaryWriter pairsFile,  MerCollections.MerDictionary[] repeatedMers, MerCollections.MerDictionary[] overflowMers)
        {
            mergingPhase = true;

            int noOfOverflows = 0;
            for (int p = 0; p < overflowMers.Length; p++)
                if (overflowMers[p] != null)
                    noOfOverflows++;

            //                 shared mers           overflow        
            int noMerSources = repeatedMers.Length + noOfOverflows; 
            MerSource[] merSources = new MerSource[noMerSources];
            int sourceCounts = 0;

            int nextSource = 0;
            // shared repeated mers partitions
            for (int i = 0; i < repeatedMers.Length; i++)
            {
                merSources[nextSource] = new MerDictionarySource(repeatedMers[i]);
                nextSource++;
                sourceCounts += repeatedMers[i].Count;
                //Console.WriteLine("repeatedMers[" + i + "]=" + repeatedMers[i].Count);
            }
            // all the overflow mer tables
            for (int i = 0; i < overflowMers.Length; i++)
            {
                if (overflowMers[i] != null)
                {
                    merSources[nextSource] = new MerDictionarySource(overflowMers[i]);
                    nextSource++;
                    sourceCounts += overflowMers[i].Count;
                    //Console.WriteLine("overflowMers[" + i + "]=" + overflowMers[i].Count);
                }
            }

            //Console.WriteLine("Total mers=" + sourceCounts);
            //Console.WriteLine("Dictionary=" + pairDictionary.Count);

            WriteBufferDelegate wbd = new WriteBufferDelegate(WriteBuffer);

            // now just merge and write until all mers have been written
            bool mersLeft = true;

            ulong[][] bufferMers = new ulong[2][];
            bufferMers[0] = new ulong[noMerSources * bufferSize];
            bufferMers[1] = new ulong[noMerSources * bufferSize];
            ulong[][] bufferValues = new ulong[2][];
            bufferValues[0] = new ulong[noMerSources * bufferSize];
            bufferValues[1] = new ulong[noMerSources * bufferSize];
            int[] bufferCount = new int[2];
            IAsyncResult[] iarWriteBuffer = new IAsyncResult[2];
            int currentBuffer = 0;
            int previousBuffer = 1;
            ulong highestMerInBuffer = 0;

            merSources[0].Open();               // just being polite

            while (mersLeft)
            {
                if (iarWriteBuffer[currentBuffer] != null)
                    wbd.EndInvoke(iarWriteBuffer[currentBuffer]);

                mersLeft = FillBuffer(merSources, ref bufferMers[currentBuffer], ref bufferValues[currentBuffer], out bufferCount[currentBuffer], highestMerInBuffer, out highestMerInBuffer);

                if (!mersLeft)
                    break;

                if (iarWriteBuffer[previousBuffer] != null && !iarWriteBuffer[previousBuffer].IsCompleted)
                    iarWriteBuffer[previousBuffer].AsyncWaitHandle.WaitOne();
                iarWriteBuffer[currentBuffer] = wbd.BeginInvoke(pairsFile, bufferMers[currentBuffer], bufferValues[currentBuffer], bufferCount[currentBuffer], null, null);

                previousBuffer = currentBuffer;
                if (currentBuffer == 0)
                    currentBuffer = 1;
                else
                    currentBuffer = 0;
            }

            for (int i = 0; i < 2; i++)
                if (iarWriteBuffer[i] != null && !iarWriteBuffer[i].IsCompleted)
                    wbd.EndInvoke(iarWriteBuffer[i]);

            //for (int s = 0; s < merSources.Length; s++)
            //{
            //    Console.WriteLine("skipped[" + s + "]=" + merSources[s].repeatsSkipped);
            //}
        }

        private static bool FillBuffer(MerSource[] merSources, ref ulong[] bufferMers, ref ulong[] bufferValues, out int bufferCount, ulong startingMer, out ulong endingMer)
        {
            bufferCount = 0;

            bool foundSomeMers = false;
            ulong highestFromFirstSource = ulong.MaxValue;

            // initial buffer loading from the first repeated mer table (it's as good a one as any)
            MerSource repeats = merSources[0];
            if (repeats.valid)
            {
                foundSomeMers = true;

                // get a buffer load of mers from this source
                for (int i = 0; i < bufferSize; i++)
                {
                    bufferMers[bufferCount] = repeats.lowestMer;
                    bufferValues[bufferCount] = repeats.value;
                    highestFromFirstSource = repeats.lowestMer;
                    bufferCount++;

                    repeats.MoveToNextDistinctMer();

                    if (!repeats.valid)
                    {
                        highestFromFirstSource = ulong.MaxValue;
                        break;
                    }
                }
            }

            endingMer = highestFromFirstSource;
            //Console.WriteLine("buffer: " + startingMer.ToString("X16") + " " + endingMer.ToString("X16"));

            // now add the matching mers/counts from all the other sources
            for (int s = 1; s < merSources.Length; s++)
            {
                MerSource ms = merSources[s];

                // if this source has not yet been opened but overlaps the current buffer range, then open it
                if (!ms.opened)
                {
                    if (startingMer <= ms.firstMer && ms.firstMer <= endingMer)
                    {
                        ms.Open();
                        //Console.WriteLine("opened " + ms.sourceType + "-" + s + ": " + startingMer.ToString("X16") + " " + endingMer.ToString("X16") + " " + ms.firstMer.ToString("X16") + " " + ms.lastMer.ToString("X16"));
                    }
                }

                if (!ms.valid)
                    continue;

                foundSomeMers = true;

                while (highestFromFirstSource >= ms.lowestMer)
                {
                    if (bufferCount == bufferMers.Length)
                    {
                        Array.Resize<ulong>(ref bufferMers, (bufferMers.Length + bufferMers.Length / 4));
                        Array.Resize<ulong>(ref bufferValues, (bufferValues.Length + bufferValues.Length / 4));
                    }

                    bufferMers[bufferCount] = ms.lowestMer;
                    bufferValues[bufferCount] = ms.value;
                    bufferCount++;

                    ms.MoveToNextDistinctMer();
                    if (!ms.valid)
                        break;
                }
            }

            if (foundSomeMers)
                Array.Sort<ulong, ulong>(bufferMers, bufferValues, 0, bufferCount);

            return foundSomeMers;
        }

        private static void GetLowestMer(ulong[] bufferMers, ulong[] bufferValues, ref int lowestIdx, int bufferCount, out ulong lowestMer, out ulong lowestCount)
        {
            lowestMer = bufferMers[lowestIdx];
            lowestCount = bufferValues[lowestIdx];

            lowestIdx++;
            // could be repeats (coming from different sources)
            while (lowestIdx < bufferCount && bufferMers[lowestIdx] == lowestMer)
            {
                if (bufferValues[lowestIdx] > lowestCount)
                    lowestCount = bufferValues[lowestIdx];
                lowestIdx++;
            }
        }

        private delegate void WriteBufferDelegate(BinaryWriter pairsFile, ulong[] bufferMers, ulong[] bufferCountPairs, int bufferCount);

        private static void WriteBuffer(BinaryWriter pairsFile, ulong[] bufferMers, ulong[] bufferCountPairs, int bufferCount)
        {
            int lowestIdx = 0;
            ulong lowestMer;
            ulong lowestCountPair;

            while (lowestIdx < bufferCount)
            {
                GetLowestMer(bufferMers, bufferCountPairs, ref lowestIdx, bufferCount, out lowestMer, out lowestCountPair);
                WriteMerToFile(lowestMer, lowestCountPair, pairsFile);
            }
        }

        private static void WriteMerToFile(ulong mer, ulong value, BinaryWriter pairsFile)
        {
            int intValue = (int)(value & 0x00000000FFFFFFFF);
            pairsFile.Write(mer);
            pairsFile.Write(intValue);
            totalPairsWritten++;

            //if (mer <= lastWritten && mer != 0)
            //    Debugger.Break();
        }

        static void RateReporter()
        {
            long currentReadsRead = 0;
            long lastReadsRead = 0;

            while (!stopMonitor)
            {
                Thread.Sleep(reportInterval);
                currentReadsRead = progressReadsProcessed;

                if (!mergingPhase)
                {
                    Console.WriteLine("tiled " + currentReadsRead + " reads (+" + (currentReadsRead - lastReadsRead) + ")");
                    lastReadsRead = currentReadsRead;
                    //Console.WriteLine("tiled " + currentValidMers + "/" + currentMersTiled + " " + kLength + "-mers from " + currentReadsRead + " reads");
                }
                else
                    Console.WriteLine("wrote " + totalPairsWritten + " pairs");
            }
        }

        private static void StopMonitorThread(Thread monitorProgress)
        {
            stopMonitor = true;
            monitorProgress.Abort();
            monitorProgress.Join();
        }
    }

    public abstract class MerSource
    {
        public ulong lowestMer;             // current mer available from this source (next to be consumed)
        public ulong value;                 // and corresponding value
        public bool valid;                  // is there a waiting mer from this source or is it finished/not open
        public bool opened;                 // has this source been opened yet?
        public ulong firstMer;              // lowest mer in this source
        public ulong lastMer;               // highest mer in this source
        public string sourceType;           // used for tracing
        public int repeatsSkipped = 0;

        public abstract void MoveToNextDistinctMer();
        public abstract void Open();
    }

    public class MerDictionarySource : MerSource
    {
        ulong[] mers;
        ulong[] values;
        int idx = 0;
        int maxIdx;

        public MerDictionarySource(MerCollections.MerDictionary source)
        {
            if (source != null)
                InitialiseMerDictionarySource(source);
            else
            {
                idx = 0;
                maxIdx = -1;
                valid = false;
            }
            opened = false;
            sourceType = "MerDictionary";
        }

        private void InitialiseMerDictionarySource(MerCollections.MerDictionary source)
        {
            int merCount = source.Sort();

            mers = source.sortedKeys;
            values = source.sortedValues;
            idx = -1;
            maxIdx = mers.Length - 1;
            valid = false;

            if (merCount > 0)
            {
                firstMer = mers[0];
                lastMer = mers[maxIdx];
                valid = true;
                MoveToNextDistinctMer();
            }
        }

        public override void Open()
        {
            opened = true;
        }

        public override void MoveToNextDistinctMer()
        {
            if (!valid)
                return;

            // no more available?
            if (idx == maxIdx)
            {
                valid = false;
                lowestMer = ulong.MaxValue;
                value = 0;
                return;
            }

            // must be some left so move on to the next one
            idx++;
            lowestMer = mers[idx];
            value = values[idx];

            // and skip any replicates of the 'next' mer
            while ((idx + 1) < maxIdx && mers[idx + 1] == lowestMer)
            {
                idx++;
                if (values[idx] > value)
                    value = values[idx];
                repeatsSkipped++;
            }
        }
    }

    public class threadParams
    {
        public int threadNo;
        public BufferedReader bufferedReadsFile;
    }
}
