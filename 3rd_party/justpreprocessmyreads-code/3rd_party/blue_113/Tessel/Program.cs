using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.IO;
using System.Diagnostics;
using CommonCode;
using MerCollections;

namespace Tessel
{
    class Program
    {
        // Tiles a set of read files into canonical k-mers and writes out a file containing these k-mers and how many times each of them appears.
        //
        // usage: Tessel -k k-merLength -g genomeLength [-t #threads] [-f fasta|fastq] [-tmp tempDir] [-s] cbtFN readsFN... or readsPattern
        //
        // Tessel  -k 25 -g 6000000 -t 4 CsporParamTest s_1_sequence_merged.fastq
        // -k 25 -g 60000000 -t 16 -min 1  -tmp g:\temp\ TesselTest-16 s_1_?_sequence.fastq

        const string version = "1.1.2";             // version - update on every significant change

        // progress monitors
        const int reportInterval = 60000;           // report back at this interval (ms)
        static bool stopMonitor = false;            // tell monitor thread to finish

        // current statistics
        static bool mergingPhase = false;           // to report on statistics in tiling/counting or merging/writing phase
        static long progressReadsRead = 0;          // no. of reads
        static long progressMersWritten = 0;        // no. of mers written to cbt file
        static long progressMersDropped = 0;        // no. of mers dropped (too few reps) 
        static int kLength = 0;                     // k-mer length - only for reporting

        // final statistics
        static object statsLock = new object();
        static long totalReadsRead = 0;             // no. of reads read from all files
        static long totalMersTiled = 0;             // no. of mers tiled from these reads
        static long totalMersWithNs = 0;            // no. of mers not written because they contained Ns
        static long totalValidMers = 0;             // no. of mers that are valid (do not contain Ns)
        static long totalUniqueMers = 0;            // no. of unique mers found (including generated RC forms)
        static long totalMersWritten = 0;           // no. of mers counted
        static long totalUniqueMersDropped = 0;     // no. of unique mers dropped (too few reps)
        static long totalMersDropped = 0;           // no. of mers dropped

        static MerTables merTables = null;
        static Dictionary<int, long> sumReps = new Dictionary<int, long>();         // for generating histogram of sum counts

        const int maxReadSize = 1000;
        const int batchSize = 1000;
        const int bufferSize = 10000;
        const int defaultReadLength = 300;
        const int defaultHeaderLength = 100;

        static void Main(string[] args)
        {
            if (args.Length == 0)
            {
                // Tessel -k k-merLength -g genomeLength [-t #threads] [-f fasta|fastq] [-tmp tempDir] [-min minCount] [-s] cbtName readsFN... or readsPattern
                WriteUsage();
                return;
            }

            int readsFormat = MerStrings.formatNone; // format e.g. fasta, fasta, sfa
            string cbtName = null;               // output binary file
            bool textMers = false;               // write out k-mers to a text file
            int merSize = 0;                     // tiling length
            long genomeSize = 0;                 // estimated genome size
            int noThreads = 1;                   // no. of parallel threads used in couting and writing
            int minCount = 1;                    // only write k-mers with counts >= to .cbt file
            string tempDir = "";                 // temporary directory for saving flushed singletons
            System.IO.SearchOption readsSearchOption  = SearchOption.TopDirectoryOnly;  // default is top-directory only; -s option forces recursive search

            string myProcessNameAndArgs = null;  // extract the command arguments from the program
            List<string> FNParams = new List<string>();    // the set of file names or patterns to be (expanded and) tiled. The first will be the .cbt name.
            string fnSeparator = Path.DirectorySeparatorChar.ToString();  // \ for Windows; / for Unix/Linux

            for (int p = 0; p < args.Length; p++)
            {
                if (args[p][0] == '-')
                {
                    args[p] = args[p].ToLower();

                    if (args[p] == "-k")
                    {
                        if (!CheckForParamValue(p, args.Length, "k-mer length number expected after -k"))
                            return;
                        try
                        {
                            merSize = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -k parameter: " + args[p + 1]);
                            return;
                        }
                        if (merSize > 31)
                        {
                            Console.WriteLine("k-mer length must be <= 31");
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-g" || args[p] == "-genome")
                    {
                        if (!CheckForParamValue(p, args.Length, "genome length number expected after -g|-genome"))
                            return;
                        try
                        {
                            genomeSize = Convert.ToInt64(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -g|-genome parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-f" || args[p] == "-format")
                    {
                        if (!CheckForParamValue(p, args.Length, "reads format expected after -f|-format"))
                            return;
                        string readsFormatParam = args[p + 1].ToLower();
                        if (readsFormatParam == "fna")
                            readsFormat = MerStrings.formatFNA;
                        else if (readsFormatParam == "fasta")
                            readsFormat = MerStrings.formatFNA;
                        else if (readsFormatParam == "fa")
                            readsFormat = MerStrings.formatFNA;
                        else if (readsFormatParam == "fastq")
                            readsFormat = MerStrings.formatFASTQ;
                        else if (readsFormatParam == "fq")
                            readsFormat = MerStrings.formatFASTQ;
                        else
                        {
                            Console.WriteLine("reads format must be fasta or fastq: " + args[p + 1]);
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

                    if (args[p] == "-m" || args[p] == "-min")
                    {
                        if (!CheckForParamValue(p, args.Length, "number expected after -m|-min"))
                            return;
                        try
                        {
                            minCount = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -m|-min parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-tmp")
                    {
                        if (!CheckForParamValue(p, args.Length, "temp directory name expected after -tmp"))
                            return;
                        tempDir = args[p + 1];
                        p++;
                        continue;
                    }

                    if (args[p] == "-s")
                    {
                        readsSearchOption = SearchOption.AllDirectories;
                    }

                    Console.WriteLine("unrecognised option: " + args[p]);
                    WriteUsage();
                    return;
                }

                FNParams.Add(args[p]);
            }


            // track what program and args produced the result files
            Process myProcess = Process.GetCurrentProcess();
            myProcessNameAndArgs = myProcess.ProcessName;
            foreach (string a in args)
                myProcessNameAndArgs = myProcessNameAndArgs + " " + a;

            // validate the temp directory
            if (tempDir != "")
            try
            {
                // add a trailing \ if the temp directory name doesn't already have one
                if (!tempDir.EndsWith(fnSeparator))
                    tempDir += fnSeparator;
                string testTempFN = tempDir + "43EDD23F-5F68-47f0-B7B9-66AE9EE3BF0B.txt";
                StreamWriter testTemp = new StreamWriter(testTempFN);
                testTemp.Close();
                File.Delete(testTempFN);
            }
            catch
            {
                Console.WriteLine("Temp directory: " + tempDir + " was invalid");
                return;
            }

            kLength = merSize;
            MerStrings.Initialise(merSize);

            string seqFilesDir = Directory.GetCurrentDirectory();
            string singletonsDir = tempDir == "" ? (seqFilesDir+fnSeparator) : tempDir;
            merTables = new MerTables(genomeSize, singletonsDir, noThreads);

            if (FNParams.Count < 2)
            {
                Console.WriteLine("Expected a cbt file name and at least one reads file name or pattern");
                return;
            }

            // the .cbt name is expected to be the first of the non-option-like parameters
            cbtName = FNParams[0];
            FNParams.RemoveAt(0);
            // and remove the .cbt suffix if it was present
            if (cbtName.ToLower().EndsWith(".cbt"))
                cbtName = cbtName.Substring(0, cbtName.Length - ".cbt".Length);

            // find out if we're to write the .cbt to a directory - and make sure the directory is present
            if (cbtName.Contains(fnSeparator))
            {
                string cbtDir = cbtName.Substring(0, cbtName.LastIndexOf(fnSeparator));
                if (!Directory.Exists(cbtDir))
                    Directory.CreateDirectory(cbtDir);
            }

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
                string[] matchedReadsFNs = Directory.GetFiles(readsFilePaths[f], readsFileNames[f], readsSearchOption);
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

            // if we weren't told what format the reads are in, use the format from the first file 
            if (readsFormat == MerStrings.formatNone)
            {
                StreamReader formatTester = new StreamReader(readsFNs[0]);
                string firstLine = formatTester.ReadLine();
                if (firstLine[0] == '>')
                    readsFormat = MerStrings.formatFNA;
                if (firstLine[0] == '@')
                    readsFormat = MerStrings.formatFASTQ;
                formatTester.Close();
            }

            // start the monitor/synchronising thread
            Thread monitorProgress = new Thread(RateReporter);
            monitorProgress.Priority = ThreadPriority.AboveNormal;
            monitorProgress.Start();

            Stopwatch countingTimer = new Stopwatch();
            Stopwatch mergingTimer = new Stopwatch();
            Stopwatch totalTimer = new Stopwatch();
            int fileNo = 0;
            countingTimer.Start();
            totalTimer.Start();

            /* Tiling and Counting phase */
            foreach (string readsFN in readsFNs)
            {
                Console.WriteLine("Reading and tiling " + readsFN);
                StreamReader reads = new StreamReader(readsFN, Encoding.ASCII, false, 1000000);
                BufferedReader bufferedReads = new BufferedReader(readsFormat, reads, null);
                fileNo++;

                //// check that the file appears to be in the expected format
                //char firstChar = (char)reads.Peek();
                //if (readsFormat == MerStrings.formatFASTQ && firstChar != '@')
                //{
                //    Console.WriteLine(readsFN + " does not appear to be in FASTQ format");
                //    StopMonitorThread(monitorProgress);
                //    return;
                //}
                //if (readsFormat == MerStrings.formatFNA && firstChar != '>')
                //{
                //    Console.WriteLine(readsFN + " does not appear to be in FASTA format");
                //    StopMonitorThread(monitorProgress);
                //    return;
                //}

                countingThreadParams[] countingParams = new countingThreadParams[noThreads];
                Thread[] countingThreads = new Thread[noThreads];

                for (int t = 0; t < noThreads; t++)
                {
                    countingParams[t] = new countingThreadParams();
                    countingParams[t].threadNumber = t;
                    countingParams[t].readsFile = bufferedReads;
                    countingParams[t].merSize = merSize;
                    countingThreads[t] = new Thread(new ParameterizedThreadStart(Program.CountingThread));
                    countingThreads[t].Priority = ThreadPriority.BelowNormal;
                    countingThreads[t].Start(countingParams[t]);
                }

                for (int t = 0; t < noThreads; t++)
                {
                    countingThreads[t].Join();
                    countingThreads[t] = null;  
                    //Console.WriteLine("finished counting thread " + (t + 1));
                }

                //long repeatedUniqueMers = 0;
                //long repeatedMerCount = 0;
                //long totalRepeatedUniqueMers = 0;
                //for (int p = 0; p < merTables.noOfPartitions; p++)
                //{
                //    repeatedUniqueMers += merTables.repeatedMers[p].Count;
                //    totalRepeatedUniqueMers += merTables.repeatedMers[p].entries.Length;
                //    foreach (KeyValuePair<ulong, long> kvp in merTables.repeatedMers[p])
                //    {
                //        long count = kvp.Value;
                //        long rcCount = (long)count & 0xffffffff;
                //        long plusCount = (long)count >> 32;
                //        repeatedMerCount += rcCount + plusCount;
                //    }
                //}
                //Console.WriteLine("repeats: " + repeatedUniqueMers + "/" + totalRepeatedUniqueMers + " unique mers with " + repeatedMerCount + " repeats");
                //for (int p = 0; p < merTables.noOfPartitions; p++)
                //    Console.WriteLine("p" + p + ": " + merTables.repeatedMers[p].Count + "/" + merTables.repeatedMers[p].entries.Length + " repeats");
                //for (int s = 0; s < merTables.noSingletonPartitions; s++)
                //    Console.WriteLine("s" + s + ": " + merTables.singletonFilters[s].Count + "/" + merTables.singletonFilters[s].entries.Length + " singles, " + merTables.flushSingletonNumber[s] + " flushes");
                //for (int t = 0; t < noThreads; t++)
                //{
                //    if (merTables.overflowMers[t] != null)
                //        Console.WriteLine("o" + t + ": " + "overflow=" + merTables.overflowMers[t].Count);
                //}

                // don't clean up for either the first or last files
                if (fileNo < readsFNs.Length && fileNo != 1)
                {
                    //Console.WriteLine("memory before condensing tables: " + GC.GetTotalMemory(false));
                    merTables.FlushLowRepMers(merTables, fileNo);
                    // and tidy up the heap after all this table clean-up
                    //GC.Collect();
                    //Console.WriteLine("memory after condensing tables: " + GC.GetTotalMemory(false));
                }

            } // for each file

            countingTimer.Stop();

            // convert the per-thread/partition tables to simple sorted arrays (not done in parallel to avoid excess memory usage)
            for (int p = 0; p < merTables.noSingletonPartitions; p++)
            {
                merTables.singletonFilters[p].Sort();
                merTables.FlushDeferredSingletons(p, 0);
            }

            //GC.Collect();   
            for (int t = 0; t < noThreads; t++)
            {
                if (merTables.overflowMers[t] != null)
                    merTables.overflowMers[t].Sort();
            }
            //GC.Collect();

            // Merging and Writing stage 
            mergingTimer.Start();

            //for (int t = 0; t < noThreads; t++)
            //    merTables.MergeSingletonFiles(t);

            BinaryWriter cbtFile = new BinaryWriter(File.Open(cbtName + "_" + merSize + ".cbt", FileMode.Create, FileAccess.Write));
            StreamWriter histo = new StreamWriter(File.Open(cbtName + "_" + merSize + "_histo.txt", FileMode.Create, FileAccess.Write));

            //Console.WriteLine(GC.GetTotalMemory(false) + " memory at end of counting");
            //for (int p = 0; p < merTables.noOfPartitions; p++)
            //    Console.WriteLine("merTable[" + p + "].length = " + merTables.repeatedMers[p].lengthEntries + "\tcount= " + merTables.repeatedMers[p].Count);

            //histo.WriteLine(GC.GetTotalMemory(false) + " memory at end of counting");
            //for (int p = 0; p < merTables.noOfPartitions; p++)
            //    histo.WriteLine("merTable[" + p + "].length = " + merTables.repeatedMers[p].length + "\t" + merTables.repeatedMers[p].Count);

            // first 4 bytes of the .cbt file is the k-mer size
            cbtFile.Write(merSize);

            MergeAndWrite(cbtFile, minCount, merTables.repeatedMers, merTables.overflowMers, merTables.singletonFilters, 
                            merTables.flushedSingletonFNs, merTables.firstFlushedSingletonMer, merTables.lastFlushedSingletonMer,
                            merTables.flushedLowRepsFNs, merTables.firstFlushedLowRepsMer, merTables.lastFlushedLowRepsMer);

            mergingTimer.Stop();

            cbtFile.Close();

            // Delete temporary ".bst" flushed files

            Console.WriteLine("Deleting temporary flush files...");
            for (int p = 0; p < merTables.noSingletonPartitions; p++)
                foreach (string f in merTables.flushedSingletonFNs[p])
                {
                    File.Delete(f);
                }
            foreach (string f in merTables.flushedLowRepsFNs)
            {
                File.Delete(f);
            }

            StopMonitorThread(monitorProgress);

            // Dump expanded k-mers and counts to a text file
            if (textMers)
            {
                Console.WriteLine("Writing expanded k-mers to file...");
                StreamWriter results = new StreamWriter(File.Open(cbtName + "_" + merSize + "_dump.txt", FileMode.Create, FileAccess.Write));
                BinaryReader cbtReader = new BinaryReader(File.Open(cbtName + "_" + merSize + ".cbt", FileMode.Open, FileAccess.Read));
                WriteExpandedMers(cbtReader, results);
                results.Close();
                cbtReader.Close();
            }

            // Extract, sort and write histogram for original and sum counts to a text file
            Console.WriteLine("Generating histogram and stats...");
            int[] sums = new int[sumReps.Count];
            long[] repsReps = new long[sumReps.Count];
            int i = 0;

            foreach (KeyValuePair<int, long> kvp in sumReps)
            {
                sums[i] = kvp.Key;
                repsReps[i] = kvp.Value;
                i++;
            }
            Array.Sort(sums, repsReps);

            histo.WriteLine(">" + myProcessNameAndArgs);
            histo.WriteLine(">sums");
            long totalMers = totalMersWritten + totalMersDropped;
            histo.WriteLine(">copies\tcounts\t" + totalMers);
            for (i = 0; i < sums.Length; i++)
            {
                histo.Write(sums[i]);
                histo.Write('\t');
                histo.Write(repsReps[i]);
                histo.Write('\t');
                long mersInBucket = sums[i] * repsReps[i];
                histo.Write(mersInBucket);
                histo.Write('\t');
                histo.Write((((float)mersInBucket/(float)totalMers)*100.0).ToString("F2"));
                histo.WriteLine();
            }

            histo.WriteLine();
            histo.WriteLine(totalReadsRead + "\treads");
            histo.WriteLine(totalMersTiled + "\tmers tiled from reads");
            histo.WriteLine(totalMersWithNs + "\tN-mers");
            histo.WriteLine(totalValidMers + "\tvalid mers");
            histo.WriteLine(totalUniqueMers + "\tdistinct mers written to cbt file");
            histo.WriteLine(totalMersWritten + "\ttotal mers written to cbt file");
            histo.WriteLine(totalMersDropped + "\tmers dropped (too few reps)");

            histo.WriteLine(countingTimer.Elapsed.TotalSeconds.ToString("#.0") + "\ts counting");
            histo.WriteLine(mergingTimer.Elapsed.TotalSeconds.ToString("#.0") + "\ts merging");

            //histo.WriteLine(merTables.singles + "\tsingles\t" + merTables.promoted + "\tpromoted\t" + merTables.repeatsFull + "\trepeatsFull\t" + merTables.repeatsQuick + "\trepeatsQuick\t" + merTables.flushes + "\tflushes");
            //for (int t = 0; t < noThreads; t++)
            //{
            //    histo.WriteLine(threadLockTimers[t].Elapsed.TotalSeconds.ToString("#.000") + "\t secs waiting by thread " + t);
            //}

            //histo.WriteLine("Partition stats");
            //histo.WriteLine("p#\tsingles\tdictionary\trepeats\tpromoted\tflushes\tresizes\tAs\tdict0\twaiting");
            //for (int p = 0; p < merSets.noOfPartitions; p++)
            //{
            //    histo.Write(p + "\t");
            //    histo.Write(merSets.pSingles[p] + "\t");
            //    histo.Write(merSets.pDictionary[p] + "\t");
            //    histo.Write(merSets.repeatedMers[p].Count + "\t");
            //    histo.Write(merSets.pPromoted[p] + "\t");
            //    histo.Write(merSets.pFlushes[p] + "\t");
            //    histo.Write(merSets.repeatedMers[p].resizeCount + "\t");
            //    histo.Write(merSets.As[p] + "\t");
            //    if (merSets.repeatedMers[p].ContainsKey(0))
            //        histo.Write((merSets.repeatedMers[p][0]>>32) + "\t" + (merSets.repeatedMers[p][0] & 0xffffffff) + "\t");
            //    else
            //        histo.Write("\t\t");
            //    histo.WriteLine(merSets.pWaitingTimer[p].Elapsed.TotalSeconds.ToString("#.000"));
            //}

            histo.Close();

            totalTimer.Stop();
            Console.WriteLine("Tiling took " + totalTimer.Elapsed.TotalSeconds.ToString("#.0") + "s");
        }

        private static void StopMonitorThread(Thread monitorProgress)
        {
            stopMonitor = true;
            monitorProgress.Abort();
            monitorProgress.Join();
        }
        
        private static void CountingThread(object threadParams)
        {
            countingThreadParams theseParams = (countingThreadParams)threadParams;
            int threadNo = theseParams.threadNumber;
            BufferedReader readsFile = theseParams.readsFile;
            int merSize = theseParams.merSize;

            ulong[] merSet = new ulong[maxReadSize];
            bool[] merValid = new bool[maxReadSize];

            Sequence[] readHeaderBatch = new Sequence[batchSize];
            Sequence[] readBatch = new Sequence[batchSize];
            for (int i = 0; i < batchSize; i++)
            {
                readHeaderBatch[i] = new Sequence(defaultHeaderLength);
                readBatch[i] = new Sequence(defaultReadLength);
            }

            countingStats threadStats = new countingStats();

            //ulong testMer;
            //MerStrings.CondenseMer("CCTAAAAAAAAAAAAAAAAAAAAAA", out testMer);
            //ulong testMerRC = MerStrings.ReverseComplement(testMer);

            bool EOF = false;
            int readsInBatch = 0;

            while (!EOF)
            {
                lock (readsFile)
                {
                    readsInBatch = readsFile.ReadReads(batchSize, readHeaderBatch, readBatch, null, null);

                    if (readsInBatch != batchSize)
                        EOF = true;

                    progressReadsRead += readsInBatch;
                }

                for (int r = 0; r < readsInBatch; r++)
                {
                    threadStats.reads++;

                    if (readBatch[r].Length < merSize)
                        continue;

                    threadStats.noOfTiledMers++;
                    //progressMersTiled++;

                    int mersInRead = MerStrings.GenerateMersFromRead(readBatch[r], merSize, ref merSet, ref merValid);

                    for (int i = 0; i < mersInRead; i++)
                    {
                        threadStats.noOfTiledMers++;
                        //progressMersTiled++;
                        if (merValid[i])
                        {
                            merTables.AddOrIncrement(merSet[i], threadNo);
                            threadStats.noOfValidMers++;
                            //progressValidMers++;
                        }
                        else
                            threadStats.noOfNMers++;
                    } // for each mer in the read

                } // for each read in the batch

            } // while there are still reads in the file

            lock (statsLock)
            {
                totalReadsRead += threadStats.reads;
                totalMersTiled += threadStats.noOfTiledMers;
                totalMersWithNs += threadStats.noOfNMers;
                totalValidMers += threadStats.noOfValidMers;
            }
        }

        private static void MergeAndWrite(BinaryWriter cbtFile, int minCount, MerDictionary[] repeatedMers, MerDictionary[] overflowMers, 
                                            MerCollection[] singletons, List<string>[] singletonFNs, List<ulong>[] lowestFlushedMer, List<ulong>[] highestFlushedMer,
                                            List<string> flushedLowRepsFNs, List<ulong> firstFlushedLowRepsMer, List<ulong> lastFlushedLowRepsMer)
        {
            mergingPhase = true;

            int noOfFlushFiles = 0;
            int maxFlushFiles = 0;
            for (int p = 0; p < singletonFNs.Length; p++)
            {
                int flushFilesInPartition = singletonFNs[p].Count;
                noOfFlushFiles += flushFilesInPartition;
                if (flushFilesInPartition > maxFlushFiles)
                    maxFlushFiles = flushFilesInPartition;
            }
            int noOfOverflows = 0;
            for (int p = 0; p < overflowMers.Length; p++)
                if (overflowMers[p] != null)
                    noOfOverflows++;

            //                 shared mers           overflow        unflushed           flushed          flushed lowrep mers
            int noMerSources = repeatedMers.Length + noOfOverflows + singletons.Length + noOfFlushFiles + flushedLowRepsFNs.Count;
            int maxOpenMerSources = repeatedMers.Length + noOfOverflows + singletons.Length + maxFlushFiles + flushedLowRepsFNs.Count;
            MerSource[] merSources = new MerSource[noMerSources];

            int nextSource = 0;
            // shared repeated mers partitions
            for (int i = 0; i < repeatedMers.Length; i++)
            {
                merSources[nextSource] = new MerDictionarySource(repeatedMers[i]);
                nextSource++;
            }
            //Console.WriteLine("created " + repeatedMers.Length + " MerDictionarySources (primary)");
            // all the overflow mer tables
            for (int i = 0; i < overflowMers.Length; i++)
            {
                if (overflowMers[i] != null)
                {
                    merSources[nextSource] = new MerDictionarySource(overflowMers[i]);
                    nextSource++;
                }
            }
            //Console.WriteLine("created " + noOfOverflows + " MerDictionarySources (overflows)");
            // all the unflushed singleton tables
            for (int i = 0; i < singletons.Length; i++)
            {
                merSources[nextSource] = new MerSingletonsSource(singletons[i]);
                nextSource++;
            }
            // all the flushed low-rep mers
            for (int i = 0; i < flushedLowRepsFNs.Count; i++)
            {
                merSources[nextSource] = new MerFlushedLowRepSource(flushedLowRepsFNs[i], firstFlushedLowRepsMer[i], lastFlushedLowRepsMer[i]);
                nextSource++;
            }
            //Console.WriteLine("created " + flushedLowRepsFNs.Count + " MerFlushedLowRepSources");
            // and the flushed singleton files (all start being closed, and opened only when needed)
            int totalSingletonFiles = 0;
            for (int p = 0; p < singletonFNs.Length; p++)
            {
                for (int f = 0; f < singletonFNs[p].Count; f++)
                {
                    string singletonFN = singletonFNs[p][f];
                    MerFlushedSingletonSource flushedSource = new MerFlushedSingletonSource(singletonFN, lowestFlushedMer[p][f], highestFlushedMer[p][f]);        
                    merSources[nextSource] = flushedSource;
                    nextSource++;
                    totalSingletonFiles++;
                }
            }
            //Console.WriteLine("created " + totalSingletonFiles + " MerFlushedSingletonSources");
            //for (int i = 0; i < nextSource; i++)
            //    Console.WriteLine(i + " " + merSources[i].valid + "  " + merSources[i].opened + " " + merSources[i].lowestMer.ToString("x16") + " " + merSources[i].firstMer.ToString("x16") + " " + merSources[i].lastMer.ToString("x16"));

            WriteBufferDelegate wbd = new WriteBufferDelegate(WriteBuffer);

            // now just merge and write until all mers have been written
            bool mersLeft = true;

            ulong[][] bufferMers = new ulong[2][];
            bufferMers[0] = new ulong[maxOpenMerSources * bufferSize];
            bufferMers[1] = new ulong[maxOpenMerSources * bufferSize];
            ulong[][] bufferCountPairs = new ulong[2][];
            bufferCountPairs[0] = new ulong[maxOpenMerSources * bufferSize];
            bufferCountPairs[1] = new ulong[maxOpenMerSources * bufferSize];
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

                mersLeft = FillBuffer(merSources, ref bufferMers[currentBuffer], ref bufferCountPairs[currentBuffer], out bufferCount[currentBuffer], highestMerInBuffer, out highestMerInBuffer);

                if (!mersLeft)
                    break;

                if (iarWriteBuffer[previousBuffer] != null && !iarWriteBuffer[previousBuffer].IsCompleted)
                    iarWriteBuffer[previousBuffer].AsyncWaitHandle.WaitOne();
                iarWriteBuffer[currentBuffer] = wbd.BeginInvoke(cbtFile, minCount, bufferMers[currentBuffer], bufferCountPairs[currentBuffer], bufferCount[currentBuffer], null, null);

                previousBuffer = currentBuffer;
                if (currentBuffer == 0)
                    currentBuffer = 1;
                else
                    currentBuffer = 0;
            }

            for (int i = 0; i < 2; i++)
                if (iarWriteBuffer[i] != null && !iarWriteBuffer[i].IsCompleted)
                    wbd.EndInvoke(iarWriteBuffer[i]);
        }

        private static bool FillBuffer(MerSource[] merSources, ref ulong[] bufferMers, ref ulong[] bufferCountPairs, out int bufferCount, ulong startingMer, out ulong endingMer)
        {
            bufferCount = 0;

            bool foundSomeMers = false;
            ulong highestRepeat = ulong.MaxValue;

            // initial buffer loading from the first repeated mer table (it's as good a one as any)
            MerSource repeats = merSources[0];
            if (repeats.valid)
            {
                foundSomeMers = true;

                // get a buffer load of mers from this source
                for (int i = 0; i < bufferSize; i++)
                {
                    bufferMers[bufferCount] = repeats.lowestMer;
                    bufferCountPairs[bufferCount] = repeats.countPair;
                    highestRepeat = repeats.lowestMer;
                    bufferCount++;

                    repeats.MoveToNextMer();

                    if (!repeats.valid)
                    {
                        highestRepeat = ulong.MaxValue;
                        break;
                    }
                }

                // and load any duplicates of the highest value as well - only possible if there was an insert race at just this mer
                // (can never overflow the buffer as we're just loading from the first repeated partition and there will be some singletons)
                while (repeats.valid && repeats.lowestMer == highestRepeat)
                {
                    bufferMers[bufferCount] = repeats.lowestMer;
                    bufferCountPairs[bufferCount] = repeats.countPair;
                    bufferCount++;
                    repeats.MoveToNextMer();
                }
            }

            endingMer = highestRepeat;
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

                while (highestRepeat >= ms.lowestMer)
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
                progressMersWritten += sumCount;
                totalUniqueMers++;
                totalMersWritten += sumCount;
            }
            else
            {
                progressMersDropped += sumCount;
                totalUniqueMersDropped++;
                totalMersDropped += sumCount;
            }

            if (sumReps.ContainsKey(sumCount))
                sumReps[sumCount]++;
            else
                sumReps.Add(sumCount, 1);
        }

        public abstract class MerSource
        {
            public ulong lowestMer;             // current mer available from this source (next to be consumed)
            public ulong countPair;             // and corresponding count pair
            public bool valid;                  // is there a waiting mer from this source or is it finished/not open
            public bool opened;                 // has this source been opened yet?
            public ulong firstMer;              // lowest mer in this source
            public ulong lastMer;               // highest mer in this source
            public string sourceType;           // used for tracing

            protected const ulong asReadSingleCount = 0x0000000100000000;
            protected const ulong rcSingleCount =     0x0000000000000001;

            public abstract void MoveToNextMer();
            public abstract void Open();
        }

        public class MerDictionarySource : MerSource
        {
            ulong[] mers;
            ulong[] counts;
            int idx = 0;
            int maxIdx;

            public MerDictionarySource(MerDictionary source)
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

            private void InitialiseMerDictionarySource(MerDictionary source)
            {
                int merCount = source.Sort();

                mers = source.sortedKeys;
                counts = source.sortedValues;
                idx = 0;
                maxIdx = mers.Length - 1;
                valid = false;

                if (merCount > 0)
                {
                    lowestMer = mers[0];
                    countPair = counts[0];
                    firstMer = mers[0];
                    lastMer = mers[maxIdx];
                    valid = true;
                }
            }

            public override void Open()
            {
                opened = true;
            }

            public override void MoveToNextMer()
            {
                if (!valid)
                    return;

                if (idx == maxIdx)
                {
                    valid = false;
                    lowestMer = ulong.MaxValue;
                    countPair = 0;
                    return;
                }

                idx++;
                lowestMer = mers[idx];
                countPair = counts[idx];
            }
        }

        public class MerSingletonsSource : MerSource
        {
            ulong[] mers;
            int idx = -1;
            int maxIdx;

            public MerSingletonsSource(MerCollection source)
            {
                int merCount = source.Count;

                mers = source.sortedKeys;
                maxIdx = mers.Length - 1;

                if (merCount > 0)
                {
                    valid = true;
                    firstMer = 0;
                    lastMer = MerTables.singletonMerMask;
                    //firstMer = mers[0] & MerTables.merMask;
                    //lastMer = mers[maxIdx] & MerTables.merMask;
                    MoveToNextMer();
                }
                else
                    valid = false;

                opened = false;
                sourceType = "Singletons";
            }

            public override void Open()
            {
                opened = true;
            }

            public override void MoveToNextMer()
            {
                if (!valid)
                    return;

                bool foundSingleton = false;
              
                // skip until we find the next active singleton
                while (idx+1 <= maxIdx)
                {
                    idx++;
                    lowestMer = mers[idx];
                    if ((lowestMer & MerTables.singletonActiveFlagMask) == 0)
                        continue;
                    if ((lowestMer & MerTables.singletonRCFlagMask) == 0)
                        countPair = asReadSingleCount;
                    else
                        countPair = rcSingleCount;
                    lowestMer &= MerTables.singletonMerMask;
                    foundSingleton = true;
                    break;
                }
                
                if (!foundSingleton)
                {
                    valid = false;
                    lowestMer = ulong.MaxValue;
                    countPair = 0;
                }
            }
        }

        public class MerFlushedSingletonSource : MerSource
        {
            BinaryReader flushedSingletons;
            int length;
            int count = 0;
            string flushedFN;

            public MerFlushedSingletonSource(string flushedFN, ulong firstMer, ulong lastMer)
            {
                this.flushedFN = flushedFN;
                this.firstMer = firstMer;
                this.lastMer = lastMer;
                valid = false;
                opened = false;
                sourceType = "FlushedSingletons";
            }

            public override void Open()
            {
                flushedSingletons = new BinaryReader(File.Open(flushedFN, FileMode.Open, FileAccess.Read));
                length = flushedSingletons.ReadInt32();

                lowestMer = flushedSingletons.ReadUInt64();
                if ((lowestMer & MerTables.singletonRCFlagMask) == 0)
                    countPair = asReadSingleCount;
                else
                    countPair = rcSingleCount;
                lowestMer &= MerTables.singletonMerMask;
                valid = true;
                count = 1;
                opened = true;
                //Console.WriteLine("opened " + flushedFN + " with " + length + " singletons. lowest=" + lowestMer.ToString("x16") + " " + firstMer.ToString("x16") + " " + lastMer.ToString("x16"));
            }

            public override void MoveToNextMer()
            {
                if (!valid)
                    return;

                if (count == length)
                {
                    valid = false;
                    lowestMer = ulong.MaxValue;
                    countPair = 0;
                    flushedSingletons.Close();
                    return;
                }
                else
                {
                    lowestMer = flushedSingletons.ReadUInt64();
                    count++;
                    if ((lowestMer & MerTables.singletonRCFlagMask) == 0)
                        countPair = asReadSingleCount;
                    else
                        countPair = rcSingleCount;
                    lowestMer &= MerTables.singletonMerMask;
                    return;
                }
            }
        }

        public class MerFlushedLowRepSource : MerSource
        {
            BinaryReader flushedLowReps;
            int length;
            int count = 0;
            string flushedFN;

            public MerFlushedLowRepSource(string flushedFN, ulong firstMer, ulong lastMer)
            {
                this.flushedFN = flushedFN;
                this.firstMer = firstMer;
                this.lastMer = lastMer;
                valid = false;
                opened = false;
                sourceType = "FlushedLowReps";
            }

            public override void Open()
            {
                flushedLowReps = new BinaryReader(File.Open(flushedFN, FileMode.Open, FileAccess.Read));
                length = flushedLowReps.ReadInt32();

                lowestMer = flushedLowReps.ReadUInt64();
                countPair = flushedLowReps.ReadUInt64();
                valid = true;
                count = 1;
                opened = true;
                //Console.WriteLine("opened " + flushedFN + " with " + length + " lowreps. lowest=" + lowestMer.ToString("x16") + " " + firstMer.ToString("x16") + " " + lastMer.ToString("x16"));
            }

            public override void MoveToNextMer()
            {
                if (!valid)
                    return;

                if (count == length)
                {
                    valid = false;
                    lowestMer = ulong.MaxValue;
                    countPair = 0;
                    flushedLowReps.Close();
                    return;
                }
                else
                {
                    lowestMer = flushedLowReps.ReadUInt64();
                    countPair = flushedLowReps.ReadUInt64();
                    count++;    ;
                    return;
                }
            }
        }

        public static void WriteExpandedMers(BinaryReader cbtFile, StreamWriter results)
        {
            bool EOF = false;
            ulong packedMer = 0;
            int merAsReadCount = 0;
            int merRCCount = 0;

            int merSize = cbtFile.ReadInt32();

            while (!EOF)
            {
                try
                {
                    packedMer = cbtFile.ReadUInt64();
                    merAsReadCount = cbtFile.ReadInt32();
                    merRCCount = cbtFile.ReadInt32();
                }
                catch
                {
                    EOF = true;
                    cbtFile.Close();
                }
                if (EOF)
                    break;

                string key = MerStrings.ExpandMer(packedMer);
                results.WriteLine(key + "\t" + merAsReadCount + "\t" + merRCCount);
            }
        }

        private static void WriteUsage()
        {
            Console.WriteLine("usage: Tessel -k k-merLength -g genomeLength [-t #threads] [-f fasta|fastq] [-tmp tempDir] [-min minCount] [-s] cbtName readsFN... or readsPattern (" + version + ")");
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

        static void RateReporter()
        {
            long currentReadsRead = 0;
            long currentMersWritten = 0;
            long currentMersDropped = 0;
            long lastReadsRead = 0;

            while (!stopMonitor)
            {
                Thread.Sleep(reportInterval);
                currentReadsRead = progressReadsRead;
                currentMersWritten = progressMersWritten;
                currentMersDropped = progressMersDropped;
                if (!mergingPhase)
                {
                    Console.WriteLine("tiled " + currentReadsRead + " reads (+" + (currentReadsRead - lastReadsRead) + ")");
                    lastReadsRead = currentReadsRead;
                }
                else
                    Console.WriteLine("wrote " + (currentMersWritten+currentMersDropped) + "/" + totalValidMers + " " + kLength +
                                        "-mers (kept " + currentMersWritten + ", dropped " + currentMersDropped + ")");
            }
        }
    }

    public class countingThreadParams
    {
        public int threadNumber;
        public BufferedReader readsFile;
        public int merSize;
    }

    public class countingStats
    {
        public long reads = 0;
        public long totalNoOfReads = 0;
        public long noOfTiledMers = 0;
        public long noOfValidMers = 0;
        public long noOfNMers = 0;
    }
}