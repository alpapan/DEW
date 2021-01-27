using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Diagnostics;

namespace MergeCBTFiles
{
    class Program
    {
        // Merges any number of .cbt files into a single .cbt file.
        //
        // usage: MergeCBTFiles minReps <cbt file patterns or names> mergedFN
        //

        const int bufferSize = 10000;
        static Dictionary<int, long> sumReps = new Dictionary<int, long>();             // for generating histogram of sum counts
        static long mersWritten = 0;                                                    // no. of mers written to merged .cbt file 
        static long mersDropped = 0;                                                    // no. of mers dropped (too few reps) 
        static long uniqueMersWritten = 0;                                              // no. of unique mers written to merged .cbt file 
        static long uniqueMersDropped = 0;                                              // no. of unique mers dropped (too few reps) 

        static void Main(string[] args)
        {
            if (args.Length < 3)
            {
                Console.WriteLine("usage: MergeCBTFiles minReps <cbt file patterns or names> mergedFN");
                return;
            }

            int minReps = Convert.ToInt32(args[0]);
            List<string> fnPatterns = new List<string>();
            for (int i = 1; i < args.Length-1; i++)
                fnPatterns.Add(args[i]);
            string mergedFN = args[args.Length-1];

            List<string> cbtFNs = new List<string>();
            foreach (string fnPattern in fnPatterns)
            {
                string[] fns = Directory.GetFiles(Directory.GetCurrentDirectory(), fnPattern);
                if (fns.Length == 0)
                {
                    Console.WriteLine(fnPattern + " did not match any files");
                    return;
                }
                foreach (string fn in fns)
                    cbtFNs.Add(fn);
            }

            int noCBTFiles = cbtFNs.Count;
            BinaryReader[] cbtFiles = new BinaryReader[noCBTFiles];
            int[] cbtMerLengths = new int[noCBTFiles];

            // open all the .cbt files
            for (int i = 0; i < noCBTFiles; i++)
            {
                cbtFiles[i] = new BinaryReader(File.Open(cbtFNs[i], FileMode.Open, FileAccess.Read));
                cbtMerLengths[i] = cbtFiles[i].ReadInt32();
            }

            // check that they all use the same k-mer length
            int merSize = cbtMerLengths[0];
            for (int i = 1; i < noCBTFiles; i++)
            {
                if (cbtMerLengths[i] != merSize)
                {
                    Console.WriteLine("inconsistent k-mer sizes in .cbt files - expected " + merSize + " but found a " + cbtMerLengths[i]);
                    return;
                }
            }

            // open the merged .cbt file
            string normalEnding = "_" + merSize + ".cbt";
            if (mergedFN.EndsWith(normalEnding))
                mergedFN.Substring(0, mergedFN.Length - normalEnding.Length);
            string histoFN = mergedFN + "_" + merSize + "_histo.txt";
            mergedFN = mergedFN + normalEnding;
            BinaryWriter cbtFile = new BinaryWriter(File.Open(mergedFN, FileMode.Create, FileAccess.Write));
            StreamWriter histo = new StreamWriter(File.Open(histoFN, FileMode.Create, FileAccess.Write));

            // write out the k-mer length
            cbtFile.Write(merSize);

            // now just merge and write until all mers have been written 
            bool mersLeft = true;
            Stopwatch mergingTimer = new Stopwatch();
            mergingTimer.Start();

            CBTSource[] merSources = new CBTSource[noCBTFiles];
            for (int i = 0; i < noCBTFiles; i++)
                merSources[i] = new CBTSource(cbtFiles[i]);

            WriteBufferDelegate wbd = new WriteBufferDelegate(WriteBuffer);

            ulong[][] bufferMers = new ulong[2][];
            bufferMers[0] = new ulong[noCBTFiles * bufferSize];
            bufferMers[1] = new ulong[noCBTFiles * bufferSize];
            ulong[][] bufferCountPairs = new ulong[2][];
            bufferCountPairs[0] = new ulong[noCBTFiles * bufferSize];
            bufferCountPairs[1] = new ulong[noCBTFiles * bufferSize];
            int[] bufferCount = new int[2];
            IAsyncResult[] iarWriteBuffer = new IAsyncResult[2];
            int currentBuffer = 0;
            int previousBuffer = 1;
            ulong highestMerInBuffer = 0;

            // fill merged buffers from all of the .cbt files and write them out
            while (mersLeft)
            {
                if (iarWriteBuffer[currentBuffer] != null)
                    wbd.EndInvoke(iarWriteBuffer[currentBuffer]);

                mersLeft = FillBuffer(merSources, bufferSize, ref bufferMers[currentBuffer], ref bufferCountPairs[currentBuffer], out bufferCount[currentBuffer], highestMerInBuffer, out highestMerInBuffer);

                if (!mersLeft)
                    break;

                if (iarWriteBuffer[previousBuffer] != null && !iarWriteBuffer[previousBuffer].IsCompleted)
                    iarWriteBuffer[previousBuffer].AsyncWaitHandle.WaitOne();
                iarWriteBuffer[currentBuffer] = wbd.BeginInvoke(cbtFile, minReps, bufferMers[currentBuffer], bufferCountPairs[currentBuffer], bufferCount[currentBuffer], null, null);

                previousBuffer = currentBuffer;
                if (currentBuffer == 0)
                    currentBuffer = 1;
                else
                    currentBuffer = 0;
            }

            // and flush the remaining buffered k-mers
            for (int i = 0; i < 2; i++)
                if (iarWriteBuffer[i] != null && !iarWriteBuffer[i].IsCompleted)
                    wbd.EndInvoke(iarWriteBuffer[i]);

            mergingTimer.Stop();

            Console.WriteLine("Merged " + uniqueMersWritten + "/" + mersWritten + " " + merSize + "-mers from " + noCBTFiles + " files. " +
                               uniqueMersDropped + "/" + mersDropped + " " + merSize + "-mers dropped (depth < " + minReps + ") in " + mergingTimer.Elapsed.TotalSeconds.ToString("#.0") + "s");

            Console.WriteLine("Generating histogram and stats...");
            int[] sums = new int[sumReps.Count];
            long[] repsReps = new long[sumReps.Count];

            sumReps.Keys.CopyTo(sums, 0);
            sumReps.Values.CopyTo(repsReps, 0);     

            Array.Sort(sums, repsReps);

            Process myProcess = Process.GetCurrentProcess();
            string myProcessNameAndArgs = myProcess.ProcessName;
            foreach (string a in args)
                myProcessNameAndArgs = myProcessNameAndArgs + " " + a;
            histo.WriteLine(">" + myProcessNameAndArgs);

            histo.WriteLine(">sums");
            long totalMers = mersWritten + mersDropped;
            histo.WriteLine(">copies\tcounts\t" + totalMers);
            for (int i = 0; i < sums.Length; i++)
            {
                histo.Write(sums[i]);
                histo.Write('\t');
                histo.Write(repsReps[i]);
                histo.Write('\t');
                long mersInBucket = sums[i] * repsReps[i];
                histo.Write(mersInBucket);
                histo.Write('\t');
                histo.Write((((float)mersInBucket / (float)totalMers) * 100.0).ToString("F2"));
                histo.WriteLine();
            }

            histo.WriteLine();
            histo.WriteLine(uniqueMersWritten + "\tdistinct mers written to cbt file");
            histo.WriteLine(mersWritten + "\ttotal mers written to cbt file");
            histo.WriteLine(mersDropped + "\tmers dropped (too few reps)");
            histo.WriteLine(mergingTimer.Elapsed.TotalSeconds.ToString("#.0") + "\ts merging");

            histo.Close();

        }

        private static bool FillBuffer(CBTSource[] merSources, int initialBufferSize, ref ulong[] bufferMers, ref ulong[] bufferCountPairs, out int bufferCount, ulong startingMer, out ulong endingMer)
        {
            bufferCount = 0;

            bool foundSomeMers = false;
            ulong highestMer = ulong.MaxValue;

            // initial buffer loading from the first .cbt file
            CBTSource first = merSources[0];
            if (first.valid)
            {
                foundSomeMers = true;

                // get a buffer load of mers from this source
                for (int i = 0; i < initialBufferSize; i++)
                {
                    bufferMers[bufferCount] = first.lowestMer;
                    bufferCountPairs[bufferCount] = first.countPair;
                    highestMer = first.lowestMer;
                    bufferCount++;

                    first.MoveToNextMer();

                    if (!first.valid)
                    {
                        highestMer = ulong.MaxValue;
                        break;
                    }
                }
            }

            endingMer = highestMer;
            
            // now add the matching mers/counts from all the other sources
            for (int s = 1; s < merSources.Length; s++)
            {
                CBTSource ms = merSources[s];

                if (!ms.valid)
                    continue;

                foundSomeMers = true;

                while (highestMer >= ms.lowestMer)
                {
                    if (bufferCount == bufferMers.Length)
                    {
                        Array.Resize<ulong>(ref bufferMers, (bufferMers.Length + bufferMers.Length / 4));
                        Array.Resize<ulong>(ref bufferCountPairs, (bufferCountPairs.Length + bufferCountPairs.Length / 4));
                    }

                    bufferMers[bufferCount] = ms.lowestMer;
                    bufferCountPairs[bufferCount] = ms.countPair;
                    bufferCount++;

                    ms.MoveToNextMer();
                    if (!ms.valid)
                        break;
                }
            }

            if (foundSomeMers)
                Array.Sort<ulong, ulong>(bufferMers, bufferCountPairs, 0, bufferCount);

            return foundSomeMers;
        }

        private static void GetLowestMer(ulong[] bufferMers, ulong[] bufferCountPairs, ref int lowestIdx, int bufferCount, out ulong lowestMer, out ulong lowestCountpair)
        {
            lowestMer = bufferMers[lowestIdx];
            lowestCountpair = bufferCountPairs[lowestIdx];

            lowestIdx++;
            while (lowestIdx < bufferCount && bufferMers[lowestIdx] == lowestMer)
            {
                lowestCountpair += bufferCountPairs[lowestIdx];
                lowestIdx++;
            }
        }

        private delegate void WriteBufferDelegate(BinaryWriter cbtFile, int minCount, ulong[] bufferMers, ulong[] bufferCountPairs, int bufferCount);

        private static void WriteBuffer(BinaryWriter cbtFile, int minCount, ulong[] bufferMers, ulong[] bufferCountPairs, int bufferCount)
        {
            int lowestIdx = 0;
            ulong lowestMer;
            ulong lowestCountPair;

            while (lowestIdx < bufferCount)
            {
                GetLowestMer(bufferMers, bufferCountPairs, ref lowestIdx, bufferCount, out lowestMer, out lowestCountPair);
                WriteMerToFile(lowestMer, lowestCountPair, minCount, cbtFile);
            }
        }

        private static void WriteMerToFile(ulong mer, ulong count, int minCount, BinaryWriter cbtFile)
        {
            int countAsRead = (int)(count >> 32);
            int countRC = (int)(count & 0x00000000FFFFFFFF);
            int sumCount = countAsRead + countRC;

            if (sumCount >= minCount)
            {
                cbtFile.Write(mer);
                cbtFile.Write(countAsRead);
                cbtFile.Write(countRC);
                mersWritten += sumCount;
                uniqueMersWritten++;
            }
            else
            {
                mersDropped += sumCount;
                uniqueMersDropped++;
            }

            if (sumReps.ContainsKey(sumCount))
                sumReps[sumCount]++;
            else
                sumReps.Add(sumCount, 1);
        }

        public class CBTSource
        {
            BinaryReader cbtFile;
            public bool valid = true;
            public ulong lowestMer = 0;
            public ulong countPair = 0;

            public CBTSource(BinaryReader cbtFile)
            {
                this.cbtFile = cbtFile;
            }

            public void MoveToNextMer()
            {
                if (!valid)
                    return;

                try
                {
                    lowestMer = cbtFile.ReadUInt64();
                    countPair = cbtFile.ReadUInt64();
                }
                catch
                {
                    valid = false;
                    lowestMer = ulong.MaxValue;
                    countPair = 0;
                }
            }
        }
    }


}
