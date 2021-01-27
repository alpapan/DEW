﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;
using System.Threading;
using CommonCode;

namespace MerCollections
{
    // A read-only (once loaded) mer table built on top of Tessel's MerDictionary. These tables are lock-free - and there can be no deletions or resizes while
    // they are being loaded. If any of the partitions get full, any subsequent entries will be put into a per-thread overflow table (these can be resized). 
    // The primary partitions and any overflow tables will be merged by the calling code. 

    public class MerTables
    {
        public int noOfPartitions = 0;                              // no. of mer partitions (hash-distributed mers)
        int maxTableSize = 40000000;                                // ensure table partitions are not too big (used in initial partition sizing)
        int minTableSize = 5000000;                                 // nor too small

        public const int int31Mask = 0x7fffffff;                    // all but the top bit of an int 32
        public const ulong fullMerMask = 0xffffffffffffffff;        // all bits in the mer are used in comparisons

        public MerDictionary[] repeatedMers = null;                 // hash-partitioned (closed) hash tables holding repeat mers
        // these tables hold canonical (lowest of as-read/RC) mers
        // shared amongst all threads - insertion should be lock free. no deletions or resizes done
        public bool[] repeatedMersFull = null;                      // repeatedMers partition is full, new mers go in per-thread overflow tables
        public MerDictionary[] overflowMers = null;                 // per-thread overflow repeat mer tables

        // compatibility only
        public LowRepMerBuffer culledBuffer = null;                 // temporary shared structures used during multi-threaded low-rep mer flush
        public object culledLock = new object();                    // (only ever done if there are more than 2 seq files being tiled - e.g. multi-lane datasets)
        public List<LowRepMerBuffer> filledCulledBuffers = null;    // only used if culledBuffer becomes full

        // constructor
        public MerTables(long dictionarySize, int noThreads)
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

            repeatedMers = new MerDictionary[noOfPartitions];               // create partitioned dictionaries
            repeatedMersFull = new bool[noOfPartitions];                    // create full flags array (default is false)
            overflowMers = new MerDictionary[noThreads];                    // create per-thread overflow tables

            // initialise per-partition structures
            for (int i = 0; i < noOfPartitions; i++)
            {
                repeatedMers[i] = new MerDictionary(partitionSize, fullMerMask);
            }

            // initialise per-thread structures
            for (int i = 0; i < noThreads; i++)
            {
                // overflowMers[i] = new MerDictionary(singletonSize); // allocated when first used to save on memory space
            }
        }

        public bool Contains(ulong mer)
        {
            int absMerHashCode = mer.GetHashCode() & int31Mask;
            int partitionNo = absMerHashCode % noOfPartitions;
            return repeatedMers[partitionNo].FindEntry(mer) >= 0;
        }

        public void AddIfNotPresent(ulong mer, long value, int threadNo)
        {
            int absMerHashCode = mer.GetHashCode() & int31Mask;
            int partitionNo = absMerHashCode % noOfPartitions;

            // this mer may have been seen before
            int idx = repeatedMers[partitionNo].FindEntry(mer);

            if (idx >= 0)
            {
                long storedValue = repeatedMers[partitionNo].entries[idx].value;
                if (value > storedValue)
                    repeatedMers[partitionNo].entries[idx].value = value;

                return;
            }

            // new mer
            if (repeatedMersFull[partitionNo])
            {
                // no space in main table so add it to this thread's overflow table
                if (overflowMers[threadNo] == null)
                {
                    overflowMers[threadNo] = new MerDictionary(repeatedMers[partitionNo].lengthEntries / 10, fullMerMask);
                    Console.WriteLine("added overflow for thread " + threadNo + " for [" + partitionNo + "]");
                }

                bool full = overflowMers[threadNo].Add(mer, value);
                // add will always work but could return 'no more please' status
                if (full)
                    overflowMers[threadNo].Resize();
            }
            else
            {
                // space in main table - so add it. Could be a race here where the mer is added twice. Resolved during writing phase.
                bool full = repeatedMers[partitionNo].Add(mer, value);
                if (full)
                    repeatedMersFull[partitionNo] = true;
            }
        }

        public long Count
        {
            get
            {
                long totalCount = 0;
                for (int p = 0; p < noOfPartitions; p++)
                    totalCount += repeatedMers[p].Count;
                for (int t = 0; t < overflowMers.Length; t++)
                    if (overflowMers[t] != null)
                        totalCount += overflowMers[t].Count;

                return totalCount;
            }
        }

        public long Capacity
        {
            get
            {
                long totalCapacity = 0;
                for (int p = 0; p < noOfPartitions; p++)
                    totalCapacity += repeatedMers[p].Capacity;
                for (int t = 0; t < overflowMers.Length; t++)
                    if (overflowMers[t] != null)
                        totalCapacity += overflowMers[t].Capacity;

                return totalCapacity;
            }
        }

    }

    class TraceEntry
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

    // Tessel compatibility only
    class FlushingThreadParams
    {
        public MerTables merTable = null;
        public int partitionNo = 0;
    }

    // Tessel compatibility only
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