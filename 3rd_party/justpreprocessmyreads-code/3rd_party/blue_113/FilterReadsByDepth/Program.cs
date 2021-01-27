using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.IO;
using System.Diagnostics;
using CommonCode;

namespace FilterReadsByDepth
{
    // FilterReadsByDepth -min minDepth -max maxDepth [-reduce maxCopies] [-histoOnly] [-stats statsFN] [-f format] [-t #threads] cbtFN fileNames

    // -histoOnly DML_121g_GTGAAA_L007_R_001_h_25.cbt  DML_121g_GTGAAA_L007_R?_001_h_0_10000.fastq -t 24
    // -max 10000 -r 100 DML_121g_25.cbt DML_121g_GTGAAA_L007_R?_001_h.fastq -t 32
    // -min 1500 -t 16 -s MCb4.5MiSeqH_fhisto.txt MCb4.5_25.cbt H:\SeqData\Xiaoyi\MCb4-5_S8_L001_R?_001_H.fastq
    // -min 1000 -t 16 -s Brodie_FilteringStats.txt h:\temp\Brodie_25.cbt f:\temp\Brodie\A1_R?.fastq

    class Program
    {
        const int maxReadLength = 1000;
        static MerTable<ulong> uniqueMers = null;   // the collection of tiled reads 
        static int merSize = 0;                     // how big the mers are in the tiled files
        static int averageDepth = 0;                // and the average depth of coverage for all loaded mers
        static int readsFormat;                     // fastq or fasta

        const int batchSize = 1000;                 // how many reads are fetched by each thread for a batch

        static int minDepth = 0;                    // reads must be at least this depth
        static int maxDepth = int.MaxValue;         // and less than or equal to this value
        static int reducedDepth = int.MaxValue;     // reduce over-represented reads down to this level
        static bool reducingReads = false;          // are we thinning reads?
        static string statsFN = null;               // names used for stats/histo file
        static bool histoOnly = false;              // only generate histogram - don't actually write any filtered reads

        static StreamReader[] readsFiles;
        static StreamWriter[] filteredReads;
        static int qualBase;                        // qual base used for Fastq char-int conversion - not used here but kept for compatibility with Blue code
        static bool fullQualHeaders = true;         // does the fastq file use full headers for the qual lines or just "+"
        static string outputDir = null;             // where to write the output files

        const int monitorInterval = 100;            // monitor thread progress at this interval
        const int reportInterval = 60000;           // report back at this interval (ms)
        static bool stopMonitor = false;            // tell monitor thread to finish
        // rough (non-mutexed) stats for progress reporting
        static long progressReads = 0;              // how many reads have been processed so far
        static long progressWantedReads = 0;        // how many reads were written to filtered set
        static long reducedReads = 0;               // how many reads were rejected as duplicates of other high rep reads
        static long totalReads = 0;                 // and how many reads did we examine?
        static long discardedReads = 0;             // how many reads were dropped (reduced or pair of a reduced read)

        static string myProcessNameAndArgs;

        static Dictionary<string, int> highRepSeqs; // long seqs + counts for high depth reads 

        static void Main(string[] args)
        {
            if (args.Length == 0)
            {
                Console.WriteLine("FilterReadsByDepth -min minDepth -max maxDepth [-reduce maxCopies] [-histoOnly] [-stats statsFN] [-f format] [-t #threads] cbtFN fileNames");
                return;
            }
            List<string> FNParams = new List<string>();    // the .cbt name and the set of file names or patterns
            int noThreads = 2;

            // find out who we are so we can track what program & args produced the result files
            Process myProcess = Process.GetCurrentProcess();
            myProcessNameAndArgs = myProcess.ProcessName;
            foreach (string a in args)
                myProcessNameAndArgs = myProcessNameAndArgs + " " + a;

            for (int p = 0; p < args.Length; p++)
            {
                if (args[p][0] == '-')
                {
                    args[p] = args[p].ToLower();

                    if (args[p] == "-s" || args[p] == "-stats")
                    {
                        if (!CheckForParamValue(p, args.Length, "stats file name string expected after -s|-stats"))
                            return;
                        statsFN = args[p + 1];
                        p++;
                        continue;
                    }

                    if (args[p] == "-min")
                    {
                        if (!CheckForParamValue(p, args.Length, "minDepth number expected after -min"))
                            return;
                        try
                        {
                            minDepth = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -min parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-max")
                    {
                        if (!CheckForParamValue(p, args.Length, "maxDepth number expected after -max"))
                            return;
                        try
                        {
                            maxDepth = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -max parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-r" || args[p] == "-reduce")
                    {
                        if (!CheckForParamValue(p, args.Length, "reduced depth number expected after -reduce"))
                            return;
                        try
                        {
                            reducedDepth = Convert.ToInt32(args[p + 1]);
                            reducingReads = true;
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -reduce parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-histoonly" || args[p] == "-ho")
                    {
                        histoOnly = true;
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

                    if (args[p] == "-f" || args[p] == "-format")
                    {
                        if (!CheckForParamValue(p, args.Length, "reads format expected after -f|-format"))
                            return;
                        string readsFormatParam = args[p+1].ToLower();
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
                            Console.WriteLine("reads format must be fasta or fastq: " + args[p+1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-o" || args[p] == "-output")
                    {
                        if (!CheckForParamValue(p, args.Length, "directory name expected after -o|-output"))
                            return;
                        outputDir = args[p + 1];
                        p++;
                        continue;
                    }

                }

                FNParams.Add(args[p]);
            }

            if (FNParams.Count < 2)
            {
                Console.WriteLine("expected a cbt file name and at least one reads file name or pattern");
                return;
            }

            // validate the output directory & set the output prefix string
            string fnSeparator = Path.DirectorySeparatorChar.ToString();  // \ for Windows; / for Unix/Linux
            if (outputDir != null)
            {
                try
                {
                    // add a trailing \ if the output directory name doesn't already have one
                    if (!outputDir.EndsWith(fnSeparator))
                        outputDir += fnSeparator;
                    string testOutputFN = outputDir + "43EDD23F-5F68-47f0-B7B9-66AE9EE3BF0B.txt";
                    StreamWriter testTemp = new StreamWriter(testOutputFN);
                    testTemp.Close();
                    File.Delete(testOutputFN);
                }
                catch
                {
                    Console.WriteLine("Output directory: " + args[6] + " was invalid");
                    return;
                }
            }

            // take the cbt file name from the start of the non-option list
            string cbtFN = FNParams[0];
            FNParams.RemoveAt(0);

            if (FNParams.Count == 0)
            {
                Console.WriteLine("did not find any reads file names or patterns");
                return;
            }

            if (!File.Exists(cbtFN))
            {
                Console.WriteLine("k-mer consensus (.cbt) file not found: " + cbtFN);
                return;
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
            Array.Sort(readsFNs);

            int noOfReadsFiles = distinctReadsFNs.Count;
            if (noOfReadsFiles == 0)
            {
                Console.WriteLine("No matching reads files found");
                return;
            }

            StreamReader formatTester = new StreamReader(readsFNs[0]);
            string firstLine = formatTester.ReadLine();
            if (firstLine[0] == '>')
                readsFormat = MerStrings.formatFNA;
            if (firstLine[0] == '@')
                readsFormat = MerStrings.formatFASTQ;
            formatTester.Close();
            formatTester = null;

            if (statsFN == null)
            {
                // construct a stats
                statsFN = readsFileNames[0].Substring(0, readsFileNames[0].LastIndexOf('.'));
                statsFN = statsFN.Replace('?', '_');
                statsFN = statsFN.Replace('*', '_');
                statsFN = statsFN.Replace('/', '_');
                statsFN = statsFN.Replace('\\', '_');
                statsFN = statsFN.Replace("__", "_");
                statsFN = statsFN.Replace("__", "_");
                statsFN = statsFN + "_fstats.txt";
                statsFN = statsFN.Replace("__", "_");
            }

            // calculate the min load depth from the min reps depth - don't need to load all of the singletons and other errors into memory
            //int minLoadDepth = minDepth / 2;
            //if (minLoadDepth <= 1)
            //    minLoadDepth = 2;
            int minLoadDepth = minDepth;

            long loadedUniqueMers = 0;
            long loadedTotalMers = 0;

            // load the .cbt file into a merTable (either a hash table (small) or a sorted array (large))
            long mersLoaded = MerStrings.LoadCBTFile(cbtFN, minLoadDepth, 0, 0, minDepth,
                                                      out uniqueMers, out merSize, out averageDepth, out loadedUniqueMers, out loadedTotalMers);

            if (merSize < 1 || merSize > 32)
            {
                Console.WriteLine("bad k-mer size found at start of .cbt file");
                return;
            }

            MerStrings.Initialise(merSize);

            highRepSeqs = new Dictionary<string, int>(10000000);
            highRepSeqs.Add(new string('A', 40), 0);
            highRepSeqs.Add(new string('C', 40), 0);
            highRepSeqs.Add(new string('G', 40), 0);
            highRepSeqs.Add(new string('T', 40), 0);

            // resolve the FASTQ qual ambiguity by reading through quals until one is encountered that can only come from either of the alternative sets
            if (readsFormat == MerStrings.formatFASTQ)
                qualBase = MerStrings.ResolveFastqQualAmbiguity(readsFNs[0], out fullQualHeaders);
            // and check whether we've got Unix data so we can write out the corrected files in the same format
            string lfConvention = MerStrings.LFConvention(readsFNs[0]);

            // start the monitor/synchronising thread
            Thread monitorProgress = new Thread(RateReporter);
            monitorProgress.Priority = ThreadPriority.AboveNormal;
            monitorProgress.Start();

            readsFiles = new StreamReader[2];
            filteredReads = new StreamWriter[2];
            Dictionary<int, int> readDepths = new Dictionary<int, int>(1000);

            // filter a pair of files at a time (allowing us to filter many files in a single run while keeping pairedness) 
            for (int f = 0; f < noOfReadsFiles; f += 2)
            {
                // for each file in the pair
                for (int p = 0; p < 2; p++)
                {
                    if (f + p < noOfReadsFiles)
                    {
                        string fullReadsFN = readsFNs[f + p];
                        string readsPath;
                        string readsFN;
                        GetPathFN(fullReadsFN, out readsPath, out readsFN);
                        string fileSuffix = readsFN.Substring(readsFN.LastIndexOf('.'));
                        string fileWithoutSuffix = readsFN.Substring(0, readsFN.LastIndexOf('.'));

                        readsFiles[p] = new StreamReader(fullReadsFN, Encoding.ASCII, false, 1000000);
                        Console.WriteLine("filtering " + readsFN);

                        // check that the file appears to be in the expected format
                        char firstChar = (char)readsFiles[p].Peek();
                        if (readsFormat == MerStrings.formatFASTQ && firstChar != '@')
                        {
                            Console.WriteLine(readsFN + " does not appear to be in FASTQ format");
                            return;
                        }
                        if (readsFormat == MerStrings.formatFNA && firstChar != '>')
                        {
                            Console.WriteLine(readsFN + " does not appear to be in FASTA format");
                            return;
                        }

                        string outputPath = outputDir == null ? readsPath + fnSeparator : outputDir;
                        if (!histoOnly)
                        {
                            string maxDepthString = maxDepth.ToString();
                            if (maxDepth == int.MaxValue)
                                maxDepthString = "max";
                            filteredReads[p] = new StreamWriter(outputPath + fileWithoutSuffix + "_" + minDepth + "_" + maxDepthString + fileSuffix,
                                                                false, readsFiles[p].CurrentEncoding, 1000000);
                            filteredReads[p].NewLine = lfConvention;
                        }
                    }
                    else
                    {
                        readsFiles[p] = null;
                        filteredReads[p] = null;
                    }
                }

                filterThreadParams[] filterParams = new filterThreadParams[noThreads];
                Thread[] filteringThreads = new Thread[noThreads];

                // ready a new thread for each parallel healer
                for (int b = 0; b < noThreads; b++)
                {
                    filterParams[b] = new filterThreadParams();
                    filterParams[b].threadNumber = b + 1;
                    filterParams[b].readsFiles = readsFiles;
                    filterParams[b].filteredReads = filteredReads;
                    filteringThreads[b] = new Thread(new ParameterizedThreadStart(Program.FilteringThread));
                    filteringThreads[b].Priority = ThreadPriority.BelowNormal;
                    filteringThreads[b].Name = b.ToString();
                    filteringThreads[b].Start(filterParams[b]);
                }

                // and wait for all threads to finish
                for (int b = 0; b < noThreads; b++)
                {
                    filteringThreads[b].Join();
                    filteringThreads[b] = null;
                    //Console.WriteLine("finished healing thread " + b);
                }

                foreach (StreamWriter r in filteredReads)
                    if (r != null)
                        r.Close();

                // merge the per-thread histograms
                for (int b = 0; b < noThreads; b++)
                {
                    Dictionary<int, int> threadReadDepths = filterParams[b].depthHisto;
                    foreach (KeyValuePair<int, int> kvp in threadReadDepths)
                    {
                        if (readDepths.ContainsKey(kvp.Key))
                            readDepths[kvp.Key] += kvp.Value;
                        else
                            readDepths.Add(kvp.Key, kvp.Value);
                    }
                }
            } // for a pair of files

            StreamWriter histo = new StreamWriter(statsFN);
            histo.WriteLine(myProcessNameAndArgs);
            histo.WriteLine();
            histo.WriteLine("depth\tcount");
            int[] depths = readDepths.Keys.ToArray<int>();
            int[] counts = readDepths.Values.ToArray<int>();
            Array.Sort<int, int>(depths, counts);
            for (int i = 0; i < readDepths.Count; i++)
            {
                histo.WriteLine(depths[i] + "\t" + counts[i]);
            }

            Console.WriteLine("discarded " + reducedReads + "/" + discardedReads + " of " + totalReads + " reads");
            histo.WriteLine("discarded " + reducedReads + "/" + discardedReads + " of " + totalReads + " reads");
            histo.Close();

            stopMonitor = true;
            monitorProgress.Join();
        }

        private static void FilteringThread(object threadParams)
        {
            filterThreadParams theseParams = (filterThreadParams)threadParams;
            int filterNumber = theseParams.threadNumber;                // which healing thread is this one?
            StreamReader[] readsFiles = theseParams.readsFiles;         // the (shared) read files to be processed
            StreamWriter[] filteredReads = theseParams.filteredReads;   // corresponding (shared) streams for filtered reads

            int noReadFNs = readsFiles.Length;
            bool[] fileActive = new bool[noReadFNs];                    // have not yet reached EOF on this reads file
            bool[,] readValid = new bool[batchSize, noReadFNs];         // did the last read from this file produce a read?
            int filesStillActive = 0;                                   // how many active reads files are still around

            string[,] readHeaderSet = new string[batchSize, noReadFNs]; // a batch of sets of read headers
            string[,] readSet = new string[batchSize, noReadFNs];       // a batch of sets of reads, possibly one from each file
            string[,] qualHeaderSet = new string[batchSize, noReadFNs]; // 
            string[,] qualsSet = new string[batchSize, noReadFNs];      // text form of the quals
            int[] merDepths = new int[maxReadLength];

            int[] depths = new int[noReadFNs];                          // depths for each read in the set
            bool[] rightDepth = new bool[noReadFNs];                    // are depths within the requested bounds?
            bool[] keepThisReadSet = new bool[batchSize];               // at least one of the set is of the desired depth, so keep the lot

            Dictionary<int, int> readDepths = new Dictionary<int, int>(1000);   // depth histogram for this thread

            for (int f = 0; f < noReadFNs; f++)
                if (readsFiles[f] != null)
                {
                    fileActive[f] = true;                               // stays true until EOF
                    filesStillActive++;                                 
                }

            // get the next set of reads and check their depths  
            while (filesStillActive > 0)
            {
                lock (readsFiles)
                {
                    // try getting the next batch of reads 
                    for (int b = 0; b < batchSize; b++)
                        for (int f = 0; f < noReadFNs; f++)
                            if (fileActive[f])                              // only if we haven't already reached EOF on this file
                            {
                                readSet[b, f] = MerStrings.ReadRead(readsFiles[f], null, readsFormat, out readHeaderSet[b, f], out qualHeaderSet[b,f], out qualsSet[b, f]);
                                if (readSet[b, f] == null)                            // this read failed - now at EOF for the file
                                {
                                    fileActive[f] = false;
                                    readValid[b, f] = false;
                                    filesStillActive--;
                                }
                                else
                                {
                                    readValid[b, f] = true;
                                    Interlocked.Increment(ref totalReads);
                                    progressReads++;
                                }
                            }
                            else
                                readValid[b, f] = false;

                } // lock to ensure synchronised reading from all reads files

                // now have a set of reads (the n'th read from each file, if they exist. So filter each one in turn. 


                for (int b = 0; b < batchSize; b++)
                {
                    keepThisReadSet[b] = true;
                    for (int f = 0; f < noReadFNs; f++)
                    {
                        if (readValid[b, f])
                        {
                            depths[f] = CalculateReadDepth(readSet[b, f], merDepths);

                            //if (depths[f] > 100000)
                            //    Debugger.Break();

                            if (reducingReads && !histoOnly)
                            {
                                if (depths[f] >= minDepth)                                  // possibly in the allowable range
                                {
                                    if (depths[f] >= maxDepth)                              // above the max level, so a candidate for thinning
                                    {
                                        // extract and test all the long read keys
                                        int keyReps = 0;
                                        for (int i = 0; i < readSet[b, f].Length - 40; i++)
                                        {
                                            string readKey = readSet[b, f].Substring(i, 40);

                                            // ignore them if they contain an N
                                            if (readKey.Contains('N'))
                                                continue;

                                            // look the next seq in the table
                                            if (highRepSeqs.ContainsKey(readKey))
                                            {
                                                highRepSeqs[readKey]++;
                                                keyReps = highRepSeqs[readKey];
                                            }

                                            // and break if we found it
                                            if (keyReps > 0)
                                                break;
                                        }

                                        if (keyReps > reducedDepth)
                                        {
                                            rightDepth[f] = false;                      // we already have enough of these reads, so mark it to be discarded
                                            Interlocked.Increment(ref reducedReads);
                                        }

                                        if (keyReps == 0)                               // didn't find this read already, so remember it for the future
                                        {
                                            string readKey = readSet[b, f].Substring(0, 40);
                                            if (!readKey.Contains('N'))
                                                lock (highRepSeqs)
                                                {
                                                    if (!highRepSeqs.ContainsKey(readKey))
                                                        highRepSeqs.Add(readKey, 1);
                                                }
                                            rightDepth[f] = true;                       // and let the read through
                                        }

                                    }
                                    else
                                        rightDepth[f] = true;                               // reducing but read between min and max so let it through
                                }
                                else
                                    rightDepth[f] = false;                                  // reducing, but below the requested min depth
                            }
                            else
                                rightDepth[f] = depths[f] >= minDepth && depths[f] <= maxDepth; // not reducing, so must be between min and max
                        }
                        else
                        {
                            depths[f] = 0;
                            rightDepth[f] = false;
                        }

                        // keep the read only if all members of the set should be kept (if paired)
                        keepThisReadSet[b] = keepThisReadSet[b] & rightDepth[f];

                        if (readDepths.ContainsKey(depths[f]))
                            readDepths[depths[f]]++;
                        else
                            readDepths.Add(depths[f], 1);
                    }
                } // end of checking a batch

                for (int b = 0; b < batchSize; b++)
                {
                    if (filesStillActive > 0 && !histoOnly)
                    {
                        lock (filteredReads)
                        {
                            for (int f = 0; f < noReadFNs; f++)
                                if (readValid[b, f])
                                {
                                    if (keepThisReadSet[b])
                                    {
                                        SaveFilteredReadAndQual(filteredReads[f], readHeaderSet[b, f], readSet[b, f], qualsSet[b, f]);
                                        progressWantedReads++;
                                    }
                                    else
                                        Interlocked.Increment(ref discardedReads);

                                }
                        } // writing out a set of healed reads
                    }
                }
            } // end of file reading/healing loop

            theseParams.depthHisto = readDepths;
        }

        private static void SaveFilteredReadAndQual(StreamWriter filteredReads, string header, string filteredRead, string quals)
        {
            MerStrings.WriteRead(filteredReads, header, filteredRead, readsFormat);
            if (readsFormat == MerStrings.formatFASTQ)
            {
                if (fullQualHeaders)
                    filteredReads.WriteLine("+" + header.Substring(1));
                else
                    filteredReads.WriteLine("+");
                filteredReads.WriteLine(quals);
            }
        }

        // Generate a set of mers from a read and calculate their read depths. 
        private static int GetMerDepths(string read, int[] merDepths)
        {
            int readLength = read.Length;
            int mersInRead = readLength - merSize + 1;

            bool merIsValid = false;
            ulong lastMer = 0;

            // read too short to tile for mers
            if (readLength < merSize)
            {
                mersInRead = 0;
                return mersInRead;
            }

            for (int i = 0; i < mersInRead; i++)
            {
                if (merIsValid)
                    merIsValid = MerStrings.CondenseMerIncremental(merSize, lastMer, read, i, out lastMer);
                else
                    merIsValid = MerStrings.CondenseMer(read.Substring(i, merSize), out lastMer);
                if (merIsValid)
                {
                    int plusCount = 0;
                    int rcCount = 0;
                    if (FindMerInUniqueMers(lastMer, out plusCount, out rcCount))
                    {
                        int sumCount = plusCount + rcCount;
                        // don't count huge unbalanced counts
                        if ((plusCount > 100 * rcCount || rcCount > 100 * plusCount))
                            merDepths[i] = 0;
                        else
                            merDepths[i] = sumCount;
                    }
                    else
                        merDepths[i] = 0;
                }
                else
                    merDepths[i] = 0;
            }

            return mersInRead;
        }

        private static bool FindMerInUniqueMers(ulong mer, out int plusCount, out int rcCount)
        {
            bool foundMer = true;
            ulong rcMer = MerStrings.ReverseComplement(mer);
            ulong countPair;
            bool rcMerWasCanonical = false;

            rcMerWasCanonical = rcMer < mer;
            if (rcMerWasCanonical)
                mer = rcMer;

            if (!uniqueMers.TryGetValue(mer, out countPair))
            {
                //string missingMer = MerStrings.ExpandMer(packedMer);
                countPair = 0;                                  // not in the table
                foundMer = false;
            }

            // extract the plus, RC and qual values from the packed ulong value
            if (rcMerWasCanonical)
            {
                rcCount = (int)(countPair >> 32);
                plusCount = (int)(countPair & 0xFFFFFFFF);
            }
            else
            {
                plusCount = (int)(countPair >> 32);
                rcCount = (int)(countPair & 0xFFFFFFFF);
            }

            return foundMer;
        }

        private static int CalculateReadDepth(string read, int[] merDepths)                                                
        {
            int noOKMers = 0;                       // no. of mers that are better than OK depth
            int OKMerAvg = 0;                       // average of counts for these mers

            int merCount = GetMerDepths(read, merDepths);

            // calculate the average depth for the passable mers
            double invSum = 0.0f;

            for (int m = 0; m < merCount; m++)
            {
                int depth = merDepths[m];
                if (depth > 0)
                {
                    noOKMers++;                                         // working out an average of the non-zero mers 
                    invSum += 1.0f / (double)depth;
                }
            }

            if (noOKMers > 0)                                           // set 'OK' depth to be 50% of average 'good' for this read
                OKMerAvg = (int)((double)noOKMers / invSum);            // harmonic mean

            //if (OKMerAvg > 10000)
            //    Debugger.Break();

            return OKMerAvg;
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
            int monitorCounter = 0;
            DateTime lastTimeAwake = DateTime.Now;
            long lastReadsCount = 0;

            while (!stopMonitor)
            {
                Thread.Sleep(monitorInterval);

                monitorCounter += monitorInterval;
                if (monitorCounter >= reportInterval)
                {
                    monitorCounter = 0;
                    DateTime timeNow = DateTime.Now;
                    double timeTaken = (timeNow - lastTimeAwake).TotalSeconds;
                    lastTimeAwake = timeNow;
                    long currentReadsCount = progressReads;
                    long readsInTime = currentReadsCount - lastReadsCount;
                    lastReadsCount = currentReadsCount;
                    int readsRate = (int)((double)readsInTime / timeTaken);
                    string reducedMsg = "";
                    if (reducingReads)
                        reducedMsg = reducedReads + " high-rep reads discarded";

                    Console.WriteLine("kept " + progressWantedReads + " from " + currentReadsCount + " reads at " + readsRate + " rps. " + reducedMsg);
                }
            }
        }
    }

    public class filterThreadParams
    {
        public int threadNumber;
        public StreamReader[] readsFiles;
        public StreamWriter[] filteredReads;
        public Dictionary<int, int> depthHisto;
    }
}

