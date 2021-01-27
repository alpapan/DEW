using System;
using System.Collections.Generic;
//using System.Linq;
using System.Text;
using System.Runtime.InteropServices;
using System.Threading;
using System.Diagnostics;

namespace MerCollections
{
    // A specialised Dictionary for holding k-mers for Tessel
    // These dictionaries never have entries deleted and are never automatically resized. These restrictions allow them to be lock-free.
    // the shared k-mer dictionaries are never resized, but the per-thread overflow dictionaries can be manually resized (safely)

    public class MerDictionary
    {
        private int[] buckets;                              // hash index array. Contains indexes into starts of linked lists for each hash bucket
        public Entry[] entries;                             // hash table contents. An array of entries linked into lists
        private int lengthBuckets = 0;                      // current length of buckets array
        public int lengthEntries = 0;                       // and the length of the entries array
        private int count;                                  // no. of allocated (in-use) entries
        private int fullTrigger;                            // table is full when this many entries have been added 
        private IEqualityComparer<ulong> comparer;          // mer comparer (equality and hash)

        public ulong[] sortedKeys = null;                   // sorted, flattened arrays used after mers have all been counted
        public ulong[] sortedValues = null;                 // will be null until Sort() is called, and then dictionary structures will be released

        public struct Entry
        {
            public int hashCode;                            // hash code that took us to this bucket. Optimisation to reduce key comparison costs. Works because folding to buckets.length means multiple hashes will go to the same bucket. Removed entries are marked by -1 in this field.
            public int next;                                // next list pointer for this entry list
            public ulong key;                               // the mer being stored. Mers are partitioned over multiple tables, based on their top bases, so the LHS of each mer will be empty and can be used for status flags. These are masked off for comparisons.
            public long value;                              // the value (count pair) associated with the key. Should never hit sign bit. Interlocked.Add will only accept signed parameters
        }

        // constructor. Must always pass in the capacity. No () constructor
        public MerDictionary(int capacity, ulong mask)
        {
            if (capacity > 0)
                Initialize(capacity);
            else
                throw (new ApplicationException("mer dictionary size <= 0"));
            this.comparer = new MerEqualityComparer(mask);
        }

        // A lock-free add to a mer dictionary (which is guaranteed not to have deletions). The 'lock free' nature of this code will rarely allow
        // the same entry to be present multiple time in the table. Post-race increments will always update the first entry added, and the
        // merge phase will ensure correct counts.
        // Return value says whether this add pushed the table past its 'full' point.
        //
        public bool Add(ulong key, long value)
        {
            //int index = FindEntry(key);
            //if (index != -1)
            //{
            //    Interlocked.Add(ref this.entries[index].value, value);
            //    return true;
            //}

            int hashCode = comparer.GetHashCode(key);
            int targetBucket = hashCode % this.buckets.Length;

            // no need to search the free list as there isn't one (no deletions). Just get the next unused entry.
            int index = Interlocked.Increment(ref this.count) - 1;

            this.entries[index].hashCode = hashCode;
            this.entries[index].next = this.buckets[targetBucket];  // point this new entry's next to the last entry pointed by buckets[index]
            this.entries[index].key = key;
            this.entries[index].value = value;
            this.buckets[targetBucket] = index;     // if collision occurs, update value in buckets[index] to point to new slot in entries[]

            return this.count >= this.fullTrigger;
        }

        public bool UpdateIfPresent(ulong key, long value)
        {
            int index = this.FindEntry(key);
            if (index >= 0)
            {
                Interlocked.Add(ref this.entries[index].value, value);
                return true;
            }

            return false;
        }

        public bool ContainsKey(ulong key)
        {
            return (this.FindEntry(key) >= 0);
        }

        public bool TryGetValue(ulong key, out long value)
        {
            int index = this.FindEntry(key);
            if (index >= 0)
            {
                value = this.entries[index].value;
                return true;
            }
            value = default(long);
            return false;
        }

        // called from setter (never used)
        private void Set(ulong key, long value)
        {
            int hashCode = comparer.GetHashCode(key);
            int targetBucket = hashCode % this.buckets.Length;

            // search per-bucket linked list for entry (and index)
            for (int i = this.buckets[targetBucket]; i >= 0; i = this.entries[i].next)
            {
                if ((this.entries[i].hashCode == hashCode) && comparer.Equals(this.entries[i].key, key))
                {
                    // found the entry so just replace its value
                    this.entries[i].value = value;
                    return;
                }
            }

            // not present so add it. As there are no deletions, there is no need for a free list
            int index;

            if (this.count >= this.fullTrigger)
                throw (new ApplicationException("dictionary full"));

            index = Interlocked.Increment(ref this.count);

            this.entries[index].hashCode = hashCode;
            this.entries[index].next = this.buckets[targetBucket];  // point this new entry's next to the last entry pointed by buckets[index]
            this.entries[index].key = key;
            this.entries[index].value = value;
            this.buckets[targetBucket] = index;     // if collision occurs, update value in buckets[index] to point to new slot in entries[]
        }


        public int FindEntry(ulong key)
        {
            int hashCode = comparer.GetHashCode(key);
            for (int i = this.buckets[hashCode % this.lengthBuckets]; i >= 0; i = this.entries[i].next)
            {
                if ((this.entries[i].hashCode == hashCode) && comparer.Equals(this.entries[i].key, key))
                {
                    return i;
                }
            }

            return -1;
        }

        private void Initialize(int capacity)
        {
            int prime = HashHelpers.GetPrime(capacity * 8 / 10);
            this.buckets = new int[prime];
            for (int i = 0; i < prime; i++)
            {
                this.buckets[i] = -1;
            }
            this.entries = new Entry[capacity];
            this.lengthBuckets = prime;
            this.lengthEntries = capacity;
            this.fullTrigger = (int)(((long)capacity * 90) / 100);              // 90% is enough to trigger a resize or 'full' condition
        }

        // can be called manually but is not thread-safe. Never called for shared repeatedMers but called for per-thread overflow tables
        public void Resize()
        {
            int prime = HashHelpers.GetPrime(this.lengthBuckets * 5 / 4);
            int[] newBuckets = new int[prime];
            for (int i = 0; i < newBuckets.Length; i++)
            {
                newBuckets[i] = -1;
            }
            Entry[] newEntries = new Entry[this.lengthEntries * 5 / 4];
            Array.Copy(this.entries, 0, newEntries, 0, this.count);
            for (int j = 0; j < this.count; j++)
            {
                int bucket = newEntries[j].hashCode % prime;
                newEntries[j].next = newBuckets[bucket];
                newBuckets[bucket] = j;
            }
            this.buckets = newBuckets;
            this.entries = newEntries;
            this.lengthBuckets = prime;
            this.lengthEntries = newEntries.Length;
            this.fullTrigger = (int)(((long)this.lengthEntries * 90) / 100);
        }

        public void Clear()
        {
            if (count > 0)
            {
                for (int i = 0; i < buckets.Length; i++)
                    buckets[i] = -1;
                Array.Clear(entries, 0, count);
                count = 0;
            }
        }

        // An in-place reduction of a MerDictionary. All entries with a depth < minReps are copied to the culled arrays 
        // and entries to be kept are moved towards the beginning of the table. Assumes that the entries are inserted 
        // sequentially and that there are no empty slots. Culled entries are placed in the shared 'culled' arrays. If these arrays
        // are filled, they are moved to the 'filled' lists.
        //
        public int Reduce(long minReps, MerTables merTable)
        {
            int culledCount = 0;
            LowRepMerBuffer culledBuffer = merTable.culledBuffer;

            // clear out all the buckets
            for (int i = 0; i < buckets.Length; i++)
                buckets[i] = -1; 

            // and either re-add or cull each entry
            int originalCount = count;
            count = 0;
            for (int i = 0; i < originalCount; i++)
            {
                ulong key = this.entries[i].key;
                long value = this.entries[i].value;
                long foldedValue = (value >> 32) + (value & 0xffffffff);

                if (foldedValue > minReps)
                    Add(key, value);                            // must always be room in the table - worst case is rewriting all entries in-place
                else
                {
                    int ci = Interlocked.Increment(ref culledBuffer.idx) - 1;

                    if (ci >= culledBuffer.limit)
                    {
                        lock (merTable.culledLock)
                        {
                            // re-test now that we have the lock (some other thread could have flushed the array while we were waiting for the lock)
                            if (culledBuffer.bufferActive)
                            {
                                // this buffer is being saved so it's no longer active
                                culledBuffer.bufferActive = false;
                                // add the full buffer to the list of filled buffers
                                if (merTable.filledCulledBuffers == null)
                                    merTable.filledCulledBuffers = new List<LowRepMerBuffer>();
                                merTable.filledCulledBuffers.Add(merTable.culledBuffer);

                                // prepare a new culledBuffer
                                LowRepMerBuffer newCulledBuffer = new LowRepMerBuffer();
                                newCulledBuffer.keys = new ulong[merTable.culledBuffer.keys.Length];
                                newCulledBuffer.values = new long[merTable.culledBuffer.values.Length];
                                newCulledBuffer.idx = 0;
                                newCulledBuffer.limit = culledBuffer.limit;
                                newCulledBuffer.bufferActive = true;
                                newCulledBuffer.bufferNo = culledBuffer.bufferNo + 1;

                                // and (atomically) make it available to any concurrent threads
                                merTable.culledBuffer = newCulledBuffer;
                            }
                        }

                        // remember to use the new buffer for this thread    
                        culledBuffer = merTable.culledBuffer;   
                        // and get a new index for this new buffer
                        ci = Interlocked.Increment(ref culledBuffer.idx) - 1;
                    }

                    culledBuffer.keys[ci] = key;
                    culledBuffer.values[ci] = value;
                    culledCount++;
                }
            } // for each entry in the original entries table

            return culledCount;
        }

        public int Sort()
        {
            // just return if we've already sorted this table
            if (sortedKeys != null)
                return count;

            sortedKeys = new ulong[count];
            sortedValues = new ulong[count];

            for (int i = 0; i < count; i++)
            {
                sortedKeys[i] = entries[i].key;
                sortedValues[i] = (ulong)entries[i].value;
            }

            entries = null;
            buckets = null;

            Array.Sort<ulong, ulong>(sortedKeys, sortedValues);

            return count;
        }

        public IEnumerator<KeyValuePair<ulong, long>> GetEnumerator()
        {
            for (int i = 0; i < count; i++)
            {
                if (entries[i].hashCode >= 0)
                {
                    KeyValuePair<ulong, long> kvp = new KeyValuePair<ulong, long>(entries[i].key, entries[i].value);
                    yield return kvp;
                }
            }
        }

        public int Count
        {
            get
            {
                return (this.count);
            }
        }

        public int Capacity
        {
            get
            {
                if (this.entries == null)
                    return this.count;
                else
                    return this.entries.Length;
            }
        }

        public long this[ulong key]
        {
            get
            {
                int index = this.FindEntry(key);
                if (index >= 0)
                {
                    return this.entries[index].value;
                }
                throw (new ApplicationException("not found"));
            }
            set
            {
                this.Set(key, value);
            }
        }
    }

    // A hash-based collection of unique mers. Effectively a Dictionary without the values, or a HashSet without the 'set' operations.
    // Adding entries to this collection is thread-safe if there are no deletions.
    public class MerCollection
    {
        private int[] buckets;                              // hash index array. Contains indexes into starts of linked lists for each hash bucket
        public Entry[] entries;                             // hash table contents. An array of entries linked into lists
        public int length = 0;                              // current length of buckets and entries arrays
        private int count = 0;                              // count of allocated entries (either still in use or free)
        private int freeCount;                              // no. of allocated entries available for re-use
        private int freeList;                               // start of list of allocated but now available entries
        public IEqualityComparer<ulong> comparer;           // mer comparer (equality and hash)

        public bool flushed = false;                        // this collection has already been flushed
        public int activeCount = 0;                         // how many threads are actively updating this collection?
        public ulong[] sortedKeys = null;                   // sorted array form of the collection. Used during merge. HashSet data structures released after copying.

        public struct Entry
        {
            public int hashCode;                            // hash code that took us to this bucket. Optimisation to reduce key comparison costs. Works because folding to buckets.length means multiple hashes will go to the same bucket.
            public int next;                                // next list pointer for this entry list
            public long key;                                // the mer being stored. Mers are partitioned over multiple tables, based on their top bases, so the LHS of each mer will be empty and can be used for status flags. These are masked off for comparisons.
            // this is really a ulong but the Interlocked primitive only accepts longs... so much casting needed
        }

        // constructor. Must always pass in the capacity. No () constructor
        public MerCollection(int capacity, ulong mask)
        {
            if (capacity > 0)
                Initialize(capacity);
            else
                throw (new ApplicationException("mer collection size <= 0"));
            this.comparer = new MerEqualityComparer(mask);
        }

        private void Initialize(int capacity)
        {
            int prime = HashHelpers.GetPrime(capacity * 8 / 10);
            this.buckets = new int[prime];
            for (int i = 0; i < prime; i++)
            {
                this.buckets[i] = -1;
            }
            this.entries = new Entry[capacity];
            this.freeList = -1;
            this.length = capacity;
            this.count = 0;
            this.freeCount = 0;
        }

        public void Add(ulong key)
        {
            this.Insert(key);
        }

        public bool ContainsKey(ulong key)
        {
            return (this.FindEntry(key) >= 0);
        }

        public bool TryGetKey(ulong mer, out ulong key)
        {
            int index = this.FindEntry(mer);
            if (index >= 0)
            {
                key = (ulong)this.entries[index].key;
                return true;
            }
            key = default(ulong);
            return false;
        }

        private void Insert(ulong key)
        {
            int freeList;
            int num = comparer.GetHashCode(key);
            int index = num % this.buckets.Length;
            for (int i = this.buckets[index]; i >= 0; i = this.entries[i].next)
            {
                if ((this.entries[i].hashCode == num) && comparer.Equals((ulong)this.entries[i].key, key))
                    throw (new ApplicationException("mer already present"));
            }

            if (this.freeCount > 0)
            {
                freeList = this.freeList;
                this.freeList = this.entries[freeList].next;
                this.freeCount--;
            }
            else
            {
                if (this.count == this.entries.Length)
                {
                    this.Resize();
                    index = num % this.buckets.Length;
                }
                freeList = Interlocked.Increment(ref this.count) - 1;
            }

            this.entries[freeList].hashCode = num;
            this.entries[freeList].next = this.buckets[index];      // point this new entry's next to the last entry pointed by buckets[index]
            this.entries[freeList].key = (long)key;
            this.buckets[index] = freeList;                         // if collision occurs, update value in buckets[index] to point to new slot in entries[]
        }

        // either inserts a key into the hash table (and returns true) or, if the key is already present, returns the existing key and false.
        // the mask is used to delete the status bits at the end of the mer key.
        public bool TryInsertKey(ulong key, out int currentIdx)
        {
            int freeList;
            int hash = comparer.GetHashCode(key);
            int index = hash % this.buckets.Length;
            currentIdx = -1;

            for (int i = this.buckets[index]; i >= 0; i = this.entries[i].next)
            {
                if ((this.entries[i].hashCode == hash) && comparer.Equals((ulong)this.entries[i].key, key))
                {
                    currentIdx = i;
                    return false;                           // insert attempt failed - return index of entry
                }
            }

            if (this.freeCount > 0)
            {
                freeList = this.freeList;
                this.freeList = this.entries[freeList].next;
                this.freeCount--;
            }

            else
            {
                if (this.count == this.entries.Length)
                {
                    this.Resize();
                    index = hash % this.buckets.Length;
                }
                freeList = Interlocked.Increment(ref this.count) - 1;
            }

            this.entries[freeList].hashCode = hash;
            this.entries[freeList].next = this.buckets[index];  // point this new entry's next to the last entry pointed by buckets[index]
            this.entries[freeList].key = (long)key;
            this.buckets[index] = freeList;                     // if collision occurs, update value in buckets[index] to point to new slot in entries[]
            currentIdx = index;
            return true;                                        // return true if Insert succeeded
        }

        public int FindEntry(ulong key)
        {
            if (this.buckets != null)
            {
                int num = comparer.GetHashCode(key);
                for (int i = this.buckets[num % this.buckets.Length]; i >= 0; i = this.entries[i].next)
                {
                    if ((this.entries[i].hashCode == num) && comparer.Equals((ulong)this.entries[i].key, key))
                    {
                        return i;
                    }
                }
            }
            return -1;
        }

        private void Resize()
        {
            int prime = HashHelpers.GetPrime(this.count + this.count / 4);
            int[] numArray = new int[prime];
            for (int i = 0; i < numArray.Length; i++)
            {
                numArray[i] = -1;
            }
            Entry[] destinationArray = new Entry[prime];
            Array.Copy(this.entries, 0, destinationArray, 0, this.count);
            for (int j = 0; j < this.count; j++)
            {
                int index = destinationArray[j].hashCode % prime;
                destinationArray[j].next = numArray[index];
                numArray[index] = j;
            }
            this.buckets = numArray;
            this.entries = destinationArray;
            this.length = prime;
        }

        public int Sort()
        {
            int mersInCollection = this.count - freeCount;

            if (sortedKeys != null)
                return mersInCollection;

            sortedKeys = new ulong[this.Capacity];
            int keyIdx = 0;

            for (int i = 0; i < this.count; i++)
            {
                if (entries[i].hashCode >= 0)
                {
                    sortedKeys[keyIdx] = (ulong)entries[i].key;
                    keyIdx++;
                }
            }

            entries = null;
            buckets = null;
            this.count = keyIdx;

            Array.Sort<ulong>(sortedKeys, 0, keyIdx);

            return mersInCollection;
        }

        public void Clear()
        {
            if (this.count >= 0)
            {
                for (int i = 0; i < buckets.Length; i++)
                    buckets[i] = -1;
                Array.Clear(entries, 0, this.count);
                this.freeList = -1;
                this.count = 0;
                this.freeCount = 0;
            }
        }

        public bool Remove(ulong key)
        {
            int hashCode = comparer.GetHashCode(key);
            int bucket = hashCode % buckets.Length;
            int last = -1;
            for (int i = buckets[bucket]; i >= 0; last = i, i = entries[i].next)
            {
                if (entries[i].hashCode == hashCode && comparer.Equals((ulong)entries[i].key, key))
                {
                    if (last < 0)
                    {
                        buckets[bucket] = entries[i].next;
                    }
                    else
                    {
                        entries[last].next = entries[i].next;
                    }
                    entries[i].hashCode = -1;
                    entries[i].next = freeList;
                    entries[i].key = (long)default(ulong);
                    freeList = i;
                    freeCount++;
                    return true;
                }
            }
            return false;
        }

        public IEnumerator<ulong> GetEnumerator()
        {
            for (int i = 0; i < this.count; i++)
            {
                if (entries[i].hashCode >= 0)
                {
                    yield return (ulong)entries[i].key;
                }
            }
        }

        public int Count
        {
            get
            {
                return (this.count - this.freeCount);
            }
        }

        public int Capacity
        {
            get
            {
                return (this.entries.Length);
            }
        }

        public ulong this[ulong key]
        {
            set
            {
                this.Insert(key);
            }
        }
    }


    public class MerEqualityComparer : IEqualityComparer<ulong>
    {
        const ulong positiveInt32Mask = 0x7fffffff;

        ulong mask;

        public MerEqualityComparer(ulong mask)
        {
            this.mask = mask;
        }

        public bool Equals(ulong x, ulong y)
        {
            return (x & mask) == (y & mask);
        }

        public int GetHashCode(ulong mer)
        {
            return (int)(HashMer(mer) & positiveInt32Mask);
        }

        public UInt32 HashMer(ulong data)
        {
            data = data & mask;
            UInt32 seed = 0xc58f1a7b;
            const UInt32 m = 0x5bd1e995;
            const Int32 r = 24;

            UInt32 h = seed ^ (UInt32)8;

            // this is unrolled loop that originally handled 32 bits at a time
            UInt32 k = (UInt32)(data >> 32);
            k *= m;
            k ^= k >> r;
            k *= m;

            h *= m;
            h ^= k;

            // second interation of loop
            k = (UInt32)(data & 0x00000000ffffffff);
            k *= m;
            k ^= k >> r;
            k *= m;

            h *= m;
            h ^= k;
            // Do a few final mixes of the hash to ensure the last few
            // bytes are well-incorporated.

            h ^= h >> 13;
            h *= m;
            h ^= h >> 15;

            return h;
        }

    }

    internal static class HashHelpers
    {
        // Fields
        internal static readonly int[] primes = new int[] { 
        3, 7, 11, 0x11, 0x17, 0x1d, 0x25, 0x2f, 0x3b, 0x47, 0x59, 0x6b, 0x83, 0xa3, 0xc5, 0xef, 
        0x125, 0x161, 0x1af, 0x209, 0x277, 0x2f9, 0x397, 0x44f, 0x52f, 0x63d, 0x78b, 0x91d, 0xaf1, 0xd2b, 0xfd1, 0x12fd, 
        0x16cf, 0x1b65, 0x20e3, 0x2777, 0x2f6f, 0x38ff, 0x446f, 0x521f, 0x628d, 0x7655, 0x8e01, 0xaa6b, 0xcc89, 0xf583, 0x126a7, 0x1619b, 
        0x1a857, 0x1fd3b, 0x26315, 0x2dd67, 0x3701b, 0x42023, 0x4f361, 0x5f0ed, 0x72125, 0x88e31, 0xa443b, 0xc51eb, 0xec8c1, 0x11bdbf, 0x154a3f, 0x198c4f, 
        0x1ea867, 0x24ca19, 0x2c25c1, 0x34fa1b, 0x3f928f, 0x4c4987, 0x5b8b6f, 0x6dda89
     };

        // Methods
        internal static int GetPrime(int min)
        {
            for (int i = 0; i < primes.Length; i++)
            {
                int num2 = primes[i];
                if (num2 >= min)
                {
                    return num2;
                }
            }
            for (int j = min | 1; j < 0x7fffffff; j += 2)
            {
                if (IsPrime(j))
                {
                    return j;
                }
            }
            return min;
        }

        internal static bool IsPrime(int candidate)
        {
            if ((candidate & 1) == 0)
            {
                return (candidate == 2);
            }
            int num = (int)Math.Sqrt((double)candidate);
            for (int i = 3; i <= num; i += 2)
            {
                if ((candidate % i) == 0)
                {
                    return false;
                }
            }
            return true;
        }
    }


}