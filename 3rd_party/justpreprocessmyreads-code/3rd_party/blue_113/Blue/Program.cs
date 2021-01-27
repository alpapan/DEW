using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Diagnostics;
using System.Threading;
using CommonCode;

namespace Blue
{
    // Blue corrects a set of reads using a k-mer consensus table derived from a set of reads (a .cbt file created by Tessel).
    // 
    // usage: Blue <options> <k-mer file> <list of file names or patterns for reads files to be corrected>
    //                  [-r run] : tag inserted into corrected reads files, and name used for stats file. Default is "corrected_<depth>"
    //                   -m minReps : min k-mer depth for OK when scanning reads (used only when dynamic calculation fails)
    //                  [-f file format] (fasta or fastq) : Blue will try to work the format out from the file so this parameter is optional
    //                  [-hp] : If this parameter is set, Blue will check for possible errors at the end of every homopolymer string - intended for 454 and IonTorrent
    //                  [-good [good%] : Just save reads that look good.
    //                  [-t #threads] : No of threads to use for the correction. Default is 1 thread.
    //                  [-fixed] : fixed length reads, corrected reads are always the same length as raw reads (default)
    //                  [-v] : variable length reads, may be trimmed or extended during healing
    //                  [-paired]: force paired reads even if single file (e.g. for merged Velvet-ready paired read sets)
    //                  [-unpaired] : don't treat a pair of files as a pair set
    //                  [-o outputDir] : directory for corrected reads etc  
    //
    // The first filename-like parameter (not preceded by a '-') is the k-mer consensus table (.cbt file).
    // This is followed by a list of file names or patterns for the sequence data files to be corrected. 
    // Blue accepts either fastq or fasta input files and writes corrected reads in the same format. 
    // Fastq files include qual scores. Fasta files may have their
    // quals in a separate .qual file (and these are expected to be lists of comma-separated integers)
    //
    // -r hg -m 50 -good 80 -t 8 Cspor_25.cbt  s_1_?_sequence.fastq

    // types of corrections made to a k-mer
    public enum FixTypes
    {
        fixNone = 0,                      // no change (default)
        fixSub = 1,                       // substitution 
        fixDel = 2,                       // fix deletion by inserting
        fixIns = 3,                       // fix insertion by deleting
        fixN = 4,                         // fix by replacing N
        fixAbandon = 5,                   // abandoned attempt at fixing mer
        maxFixTypes = 6
    }

    // results of checking whether read should be healed (and what happened when we tried heaing it)
    public enum ReadState
    {
        ReadOK = 0,
        ReadHasErrors = 1,
        ReadSkipDeep = 2,
        ReadSkipUnbalanced = 3,
        ReadSkipLowComplexity = 4,
        ReadAbandonedRewriting = 5,
        ReadAbandonedTreeSize = 6,
        ReadAbandonedNs = 7
    }

    // reasons for calling TryHealingMer
    enum CheckResult
    {
        checkNone = 0,                    // not checking this mer at all
        checkDepth = 1,                   // sudden drop in depth
        checkHP = 2,                      // end of 454 HP run
        checkAlt = 3,                     // viable alternate k-mers
        checkPoor = 4,                    // below OK or asymmetric counts
        checkBad = 5                      // below minDepth
    }

    class Program
    {
        const string version = "1.1.2";             // version - update on every significant change
        const bool usingAbandonedCache = false;     // are we using the cache of indeterminate results?
        const bool saveSlowReads = false;           // save slow-correcting reads for analysis
        const bool perfTrace = false;               // save performance stats for each read

        const int monitorInterval = 100;            // monitor thread progress at this interval
        const int reportInterval = 60000;           // report back at this interval (ms)
        static bool stopMonitor = false;            // tell monitor thread to finish
        // rough (non-mutexed) stats for progress reporting
        static long progressReads = 0;              // how many reads have been processed so far
        static long progressHealedReads = 0;        // how many reads were healed
        // states for monitor thread reporting
        static bool inHealingPhase = false;
        static bool inLoadingPhase = false;

        const int maxMerSize = 32;                  // max bases that can fit in a single ulong
        const int batchSize = 1000;                 // how many reads are fetched by each healing thread for a batch
        const int defaultReadLength = 200;          // default char[] length (resized as needed)
        const int defaultHeaderLength = 50;         // default header length


        // trace string consts - always keep in sync with corresponding enums
        static char[] fixTypes = new char[] { '-', 'S', 'D', 'I', 'N', 'A' };
        static string[] fixNames = new string[] { "None", "Sub", "Del", "Ins", "N", "Abn" };
        static string[] checkNames = new string[] { "", "(checkDepth)", "(checkHP)", "(checkAlt)", "(poor)", "(bad)" };

        const int maxNs = 3;                        // max no. of N's in a single read (avoid combinatoric explosion)
        const int maxFollowerRepairs = 3;           // max no. of repairs that can be made while trying to find followers
        const int maxCacheAllowed = 10000;          // ** fix ** adjust to read length?
        const int maxTHMAllowed = 5000;             // ** fix ** adjust to read length?
        const int almostEnd = 2;                    // how many bases from the end is 'close'
        const int maxConsecutivePoor = 100;         // max consecutive 'poor' k-mers allowed in CountFollowers recursive loop
        const int maxConsecutiveNFixesAllowed = 10; // max no.of consecutive Ns that can be fixed
        const int highDepthFactor = 10;             // high depth reads are this many times the average
        const int veryHighDepthFactor = 100;        // and very deep reads are... 

        const bool forwardPass = false;             // which way through the read (parameter to TryHealingRead)
        const bool reversePass = true;

        // consts for types of variants to generate
        const int varyLast = 1;                     // vary last base only
        const int varySingle = 2;                   // allow any single base to vary
        const int varyTwo = 3;                      // allow any two bases to vary

        const int modelMostlySubs = 1;              // largely sub errors and the occasional random indel
        const int modelIndelsCommon = 2;            // indels common after homopolymer runs so do extra checks       
        static int errorModel = modelMostlySubs;    // what type of data are we dealing with

        static bool readsFixedLength = true;        // reads are fixed length and must be the same length when written out (false --> reads lengths can change)
        static bool trimReads = false;              // trim reads back after correction (remove trailing bases with poor coverage or lacking pair support)

        static int readsFormat = MerStrings.formatNone; // what format are the reads files (see MerString.formatXXX enumeration)
        // this will be set automatically by looking at the file names/content if not explicitly set by a parameter
        static bool pairingSpecified = false;       // pairing parameter set - so don't set by default
        static bool pairedReads = false;            // do we correct reads in sets or singly? (true by default if we have multiple input files)
        static bool fullQualHeaders = false;        // do the reads use full FASTQ qual headers or just "+"?
        static int qualBase = 0;                    // qual offset for fastq reads. Set to 33 or 64 once format is known. Set to 0 for Fasta files
        static char replacementQual = (char)35;     // qual value used for fixed/inserted bases. (0-40)

        const bool findBiggest = true;              // look for biggest in list (for FindBestValue)
        const bool findSmallest = false;            // look for smallest in list (for FindBestValue)

        static int saveReads = saveAll;             // save healed reads option (default is save all)
        const int saveAll = 1;                      // all reads - faulty or not
        const int saveGoodOnly = 2;                 // only 'good' reads
        static int saveGoodMaxPoor = 5;             // % of k-mers that must be poor to qualify the whole read as 'not good' 
        static bool saveNonGoodReads = false;       // (param) save any non-good reads to 'problems' files

        static StreamWriter perf = null;
        //static StreamWriter trace = new StreamWriter("trace.txt");

        static string myProcessNameAndArgs;

        // precise stats (accumulated locally & merged under lock protection)
        static healingStats stats = new healingStats();
        static Stopwatch timeSpentReading = new Stopwatch();
        static Stopwatch timeSpentWriting = new Stopwatch();
        static Stopwatch[] threadWaitingRead;
        static Stopwatch[] threadWaitingWrite;

        //[ThreadStatic]
        //static bool previousMerGood = false;
        //static Dictionary<int, long> subCounts = new Dictionary<int, long>(4);

        const int traceReadLimit = 10000;           // stop trace after this many reads have been healed
        static StreamWriter traces;
        static string tracesFN;
        static Dictionary<string, string> refSequences;
        const int traceOff = 0;                     // no tracing
        const int traceChanges = 1;                 // tracing healing at a high level (default trace)
        const int traceChoices = 2;                 // trace choices made (and variants to choose from)
        const int traceRead = 3;                    // detailed tracing for a specified read (set by check for a read)
        const int traceFollowers = 4;               // trace out the recursive followers code (set manually if needed)
        static int tracing = traceOff;              // level of tracing wanted
        static List<string> traceClump = new List<string>(100); // traces are written to this list before being written to the trace file.
        static bool tracingThisRead = false;        // flag to say we're in the trapped read

        static Object traceLock = new Object();     // lock object used to single-thread writing out traces
        static Object statsLock = new Object();     // lock object used when merging local stats with global ones

        // shared application wide data. These are either read-ony after they've been initialised or accessed appropriately
        // this done to avoid passing these variables as (unchanging) parameters to a number of highly-called methods.
        static MerTable<ulong> uniqueMers = null;   // the collection of tiled reads 
        static int merSize = 0;                     // how big the mers are in the tiled files
        static int averageDepth = 0;                // and the average depth of coverage for all loaded mers
        static int minReps = 0;                     // (param) default reps count between bad & poor
        static MerTable<int> merPairs = null;       // pairs of short (16-mer) pairs derived from reads - used for longer-than-k-mer checks (MerIshealingCandidate & FindPlausibleVariant)
        //static MerDictionary<int> merPairs = null;  // MerTable is more memory efficient (resized down after loading) 
        const int merPairFragmentLength = 16;       // pairs of 16-mers
        static int merPairGap = 0;                  // gap between the two 16-mers in a pair
        static int merPairLength = 0;               // length to go backwards to the start of a read
        const bool merPairBack = true;              // m points to the start of the k-mer containing the second fragment
        const bool merPairForward = false;          // m points to the start of the k-mer containing the first fragment
        const bool merPairBackOnly = true;          // only try for 'back' pairs
        const bool merPairEither = false;
        const int pairNotPossible = 0;              // too close to the start of a read for a pair
        const int pairFound = 1;                    // suitable pair found
        const int pairNotFound = 2;                 // matching pair not found or reps too low
        static HashSet<ulong> lowComplexityTrap;    // small hash set of low-complexity mers. These are detected and avoided 
        static ulong lowComplexityTrapMask = 0xffffffffff000000; // just the first 16 bases

        [ThreadStatic]
        static Queue<merProperties> freeMerProperties;      // a cache of recyclable merProperties objects - one per thread to avoid locking overhead
        [ThreadStatic]
        static Queue<Depth[]> freeDepths;                   // a cache of recyclable depth arrays objects - one per thread to avoid locking overhead
        [ThreadStatic]
        static Queue<List<merProperties>> freeVariantSets;  // a cache of variantSets (THM)
        [ThreadStatic]
        static Queue<ulong[]> freeMerVariants;              // a cache of merVariant arrays (THM)
        [ThreadStatic]
        static int[] merDepths;                             // and save allocating this for every read as well
        [ThreadStatic]
        static Queue<Sequence> freeSequences;               // recycling sequences (reads and quals)

        public static void Main(string[] args)
        {
            if (args.Length == 0)
            {
                WriteUsage();
                return;
            }

            string runName = null;               // tag used for stats/trace/corrected files created by this run
            string cbtFN = null;                 // name of .cbt file (and used for create .prs pairs file name)
            string statsFN = null;               // name of stats file (default is cbtFN + "stats.txt")
            string outputDir = null;             // output directory (default is current directory)
            List<string> FNParams = new List<string>();    // the .cbt name and the set of file names or patterns
            int noHealers = 1;                   // no. of healing threads to run in parallel (1 thread is default)
            int minLoadReps = 0;                 // min rep count needed before mer will be loaded into uniqueMers table (derived from minReps)

            for (int p = 0; p < args.Length; p++)
            {
                if (args[p][0] == '-')
                {
                    args[p] = args[p].ToLower();

                    if (args[p] == "-h" || args[p] == "-help")
                    {
                        WriteUsageFull();
                        return;
                    }

                    if (args[p] == "-r" || args[p] == "-run")
                    {
                        if (!CheckForParamValue(p, args.Length, "run name string expected after -r|-run"))
                            return;
                        runName = args[p + 1];
                        p++;
                        continue;
                    }

                    if (args[p] == "-s" || args[p] == "-stats")
                    {
                        if (!CheckForParamValue(p, args.Length, "stats file name string expected after -s|-stats"))
                            return;
                        statsFN = args[p + 1];
                        p++;
                        continue;
                    }

                    if (args[p] == "-m" || args[p] == "-min")
                    {
                        if (!CheckForParamValue(p, args.Length, "minReps number expected after -m|-min"))
                            return;
                        try
                        {
                            minReps = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -m|-min parameter: " + args[p + 1]);
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

                    if (args[p] == "-hp")
                    {
                        errorModel = modelIndelsCommon;
                        continue;
                    }

                    if (args[p] == "-paired")
                    {
                        pairingSpecified = true;
                        pairedReads = true;
                        continue;
                    }

                    if (args[p] == "-unpaired")
                    {
                        pairingSpecified = true;
                        pairedReads = false;
                        continue;
                    }

                    if (args[p] == "-fixed")
                    {
                        readsFixedLength = true;
                        continue;
                    }

                    if (args[p] == "-v" || args[p] == "-variable")
                    {
                        readsFixedLength = false;
                        continue;
                    }

                    if (args[p] == "-g" || args[p] == "-good")
                    {
                        if (!CheckForParamValue(p, args.Length, "%good expected after -g|-good"))
                            return;
                        int saveGoodParam = 0;
                        try
                        {
                            saveGoodParam = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -g|-good parameter: " + args[p + 1]);
                            return;
                        }

                        if (saveGoodParam > 100)
                            saveGoodParam = 100;
                        if (saveGoodParam < 0)
                            saveGoodParam = 0;
                        saveGoodMaxPoor = 100 - saveGoodParam;
                        saveReads = saveGoodOnly;

                        p++;
                        continue;
                    }

                    if (args[p] == "-trace")
                    {
                        tracing = traceChanges;
                        continue;
                    }
                    if (args[p] == "-tracechanges")
                    {
                        tracing = traceChanges;
                        continue;
                    }
                    if (args[p] == "-tracechoices")
                    {
                        tracing = traceChoices;
                        continue;
                    }

                    if (args[p] == "-problems")
                    {
                        saveNonGoodReads = true;
                        continue;
                    }

                    if (args[p] == "-t" || args[p] == "-threads")
                    {
                        if (!CheckForParamValue(p, args.Length, "number expected after -t|-threads"))
                            return;
                        try
                        {
                            noHealers = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -t|-threads parameter: " + args[p + 1]);
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

                    Console.WriteLine("unrecognised option: " + args[p]);
                    WriteUsage();
                    return;
                }

                FNParams.Add(args[p]);
            }

            if (runName == null)
                runName = "corrected_" + minReps;

            if (minReps == 0)
            {
                Console.WriteLine("no minimum k-mer depth specified (-m)");
                return;
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
                    if (!Directory.Exists(outputDir))
                        Directory.CreateDirectory(outputDir);
                    if (!outputDir.EndsWith(fnSeparator))
                        outputDir += fnSeparator;
                    string testOutputFN = outputDir + "43EDD23F-5F68-47f0-B7B9-66AE9EE3BF0B.txt";
                    StreamWriter testTemp = new StreamWriter(testOutputFN);
                    testTemp.Close();
                    File.Delete(testOutputFN);
                }
                catch
                {
                    Console.WriteLine("Output directory: " + outputDir + " was invalid");
                    return;
                }
            }

            if (FNParams.Count == 0)
            {
                Console.WriteLine("did not find any reads file names or patterns");
                return;
            }

            // take the cbt file name from the start of the non-option list
            cbtFN = FNParams[0];
            FNParams.RemoveAt(0);

            if (FNParams.Count == 0)
            {
                Console.WriteLine("did not find any reads file names");
                return;
            }

            // calculate the min load depth from the min reps depth - don't need to load all of the singletons and other errors into memory
            minLoadReps = minReps / 4;
            if (minLoadReps <= 1)
                minLoadReps = 2;

            // find out who we are so we can track what program & args produced the result files
            Process myProcess = Process.GetCurrentProcess();
            myProcessNameAndArgs = myProcess.ProcessName;
            foreach (string a in args)
                myProcessNameAndArgs = myProcessNameAndArgs + " " + a;
            myProcessNameAndArgs += " (" + version + ") " + DateTime.Now.ToString();

            //for (int i = 0; i < 4; i++)
            //    subCounts.Add(i, 0);


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

            int noOfReadsFiles = distinctReadsFNs.Count;
            if (noOfReadsFiles == 0)
            {
                Console.WriteLine("No matching reads files found");
                return;
            }

            foreach (string distinctFN in distinctReadsFNs)
            {
                if (!File.Exists(distinctFN))
                {
                    Console.WriteLine("Reads file " + distinctFN + " not found");
                    return;
                }
            }

            if (statsFN == null)
            {
                // construct a stats file name if one wasn't given
                int dotIdx = readsFileNames[0].LastIndexOf('.');
                if (dotIdx >= 0)
                {
                    statsFN = readsFileNames[0].Substring(0, readsFileNames[0].LastIndexOf('.'));
                    statsFN = statsFN.Replace('?', '_');
                    statsFN = statsFN.Replace('*', '_');
                    statsFN = statsFN.Replace('/', '_');
                    statsFN = statsFN.Replace('\\', '_');
                    statsFN = statsFN.Replace("__", "_");
                    statsFN = statsFN.Replace("__", "_");
                }
                else
                    statsFN = readsFileNames[0];
                statsFN = statsFN + "_" + runName + "_stats.txt";
                statsFN = statsFN.Replace("__", "_");
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

            // resolve the FASTQ qual ambiguity by reading through quals until one is encountered that can only come from either of the alternative sets
            if (readsFormat == MerStrings.formatFASTQ)
                qualBase = MerStrings.ResolveFastqQualAmbiguity(readsFNs[0], out fullQualHeaders);

            // and check whether we've got Unix data so we can write out the corrected files in the same format
            string lfConvention = MerStrings.LFConvention(readsFNs[0]);

            int fileArraySize = noOfReadsFiles;
            if (noOfReadsFiles == 2 & !pairingSpecified)
                pairedReads = true;

            StreamReader[] readsFiles = new StreamReader[fileArraySize];
            BufferedReader[] bufferedReadsFiles = new BufferedReader[fileArraySize];
            StreamReader[] qualFiles = new StreamReader[fileArraySize];         // will be null for all but 454 reads in FASTA format
            //string[] qualFNs = new string[fileArraySize];
            StreamWriter[] healedReads = new StreamWriter[fileArraySize];
            StreamWriter[] healedQuals = new StreamWriter[fileArraySize];       // will be null for all but 454 reads in FASTA format
            StreamWriter[] problemReads = null;
            StreamWriter[] problemQuals = null;                                 // will be null for all but 454 reads in FASTA format
            StreamWriter[] singleReads = null;
            StreamWriter[] singleQuals = null;                                  // will be null for all but 454 reads in FASTA format
            if (saveNonGoodReads || saveSlowReads)
            {
                problemReads = new StreamWriter[fileArraySize];
                problemQuals = new StreamWriter[fileArraySize];
            }
            if (saveReads == saveGoodOnly && pairedReads)
            {
                singleReads = new StreamWriter[fileArraySize];
                singleQuals = new StreamWriter[fileArraySize];
            }
            bool qualsPresent = false;                                          // quals found - so read and write 

            // open original and healed streams for each read file (and qual if necessary)
            for (int f = 0; f < noOfReadsFiles; f++)
            {
                string fullReadsFN = readsFNs[f];
                string readsPath;
                string readsFN;
                GetPathFN(fullReadsFN, out readsPath, out readsFN);
                string fileSuffix = readsFN.Substring(readsFN.LastIndexOf('.'));
                string fileWithoutSuffix = readsFN.Substring(0, readsFN.LastIndexOf('.'));

                readsFiles[f] = new StreamReader(fullReadsFN, Encoding.ASCII, false, 1000000);
                Console.WriteLine("correcting " + fullReadsFN);

                // check that the file appears to be in the expected format (**to do** combine with earlier file status checks) - breaks bufferedReads if left here
                //char firstChar = (char)readsFiles[f].Peek();
                //if (readsFormat == MerStrings.formatFASTQ && firstChar != '@')
                //{
                //    Console.WriteLine(readsFN + " does not appear to be in FASTQ format");
                //    return;
                //}
                //if (readsFormat == MerStrings.formatFNA && firstChar != '>')
                //{
                //    Console.WriteLine(readsFN + " does not appear to be in FASTA format");
                //    return;
                //}

                if (readsFormat == MerStrings.formatFNA)
                {
                    string qualFN = fileWithoutSuffix + ".qual";
                    if (File.Exists(qualFN))
                        qualsPresent = true;
                }

                if (readsFormat == MerStrings.formatFASTQ)
                    qualsPresent = true;

                if (saveReads == saveAll || saveReads == saveGoodOnly)
                {
                    string outputPath = outputDir == null ? readsPath + fnSeparator : outputDir;

                    healedReads[f] = new StreamWriter(outputPath + fileWithoutSuffix + "_" + runName + fileSuffix, false, readsFiles[f].CurrentEncoding, 1000000);
                    healedReads[f].NewLine = lfConvention;
                    if (saveNonGoodReads || saveSlowReads)
                    {
                        problemReads[f] = new StreamWriter(outputPath + fileWithoutSuffix + "_" + runName + "_problems" + fileSuffix);
                        problemReads[f].NewLine = lfConvention;
                    }
                    if (saveReads == saveGoodOnly && pairedReads)
                    {
                        singleReads[f] = new StreamWriter(outputPath + fileWithoutSuffix + "_" + runName + "_singles" + fileSuffix);
                        singleReads[f].NewLine = lfConvention;
                    }

                    if (qualsPresent && readsFormat == MerStrings.formatFNA)
                    {
                        string qualFN = fileWithoutSuffix + ".qual";
                        qualFiles[f] = new StreamReader(qualFN);
                        healedQuals[f] = new StreamWriter(outputPath + fileWithoutSuffix + "_" + runName + ".qual");
                        healedQuals[f].NewLine = lfConvention;
                        if (saveNonGoodReads || saveSlowReads)
                        {
                            problemQuals[f] = new StreamWriter(outputPath + fileWithoutSuffix + "_" + runName + "_problems" + ".qual");
                            problemQuals[f].NewLine = lfConvention;
                        }
                        if (saveReads == saveGoodOnly)
                        {
                            singleQuals[f] = new StreamWriter(outputPath + fileWithoutSuffix + "_" + runName + "_singles" + ".qual");
                            singleQuals[f].NewLine = lfConvention;
                        }
                    }
                    if (readsFormat == MerStrings.formatFASTQ)
                    {
                        healedQuals[f] = healedReads[f];
                        if (saveNonGoodReads || saveSlowReads)
                            problemQuals[f] = problemReads[f];
                        if (saveReads == saveGoodOnly && pairedReads)
                            singleQuals[f] = singleReads[f];
                    }
                }

                bufferedReadsFiles[f] = new BufferedReader(readsFormat, readsFiles[f], qualFiles[f]);
            }

            if (perfTrace)
            {
                perf = new StreamWriter(outputDir + runName + "_perf.txt");
                perf.WriteLine("read#\tfaHeader\ttime\tcacheSize\treverse\tresult\tunhealedRead\tthmCalls\tpoor");
            }

            // load reference sequences corresponding to reads if these exist
            refSequences = new Dictionary<string, string>();
            foreach (string readFN in readsFNs)
            {
                string suffix = readFN.Substring(readFN.LastIndexOf('.'));
                string refFN = readFN.Replace(suffix, ".ref");
                if (File.Exists(refFN))
                    LoadRefSequences(refFN, refSequences);
            }
            tracesFN = runName + "_trace_" + merSize + "_" + minReps + ".txt";
            if (tracing >= traceChanges)
            {
                traces = new StreamWriter(outputDir + tracesFN);
                traces.WriteLine(myProcessNameAndArgs);
            }

            StreamWriter statsFile = new StreamWriter(outputDir + statsFN);

            long loadedUniqueMers = 0;
            long loadedTotalMers = 0;

            //Console.WriteLine(GC.GetTotalMemory(false) + " memory before loadCBT");

            // load the .cbt file into a merTable (either a hash table (small) or a sorted array (large))
            long mersLoaded = MerStrings.LoadCBTFile(cbtFN, minLoadReps, 0, 0, minReps,
                                                      out uniqueMers, out merSize, out averageDepth, out loadedUniqueMers, out loadedTotalMers);

            GC.Collect();
            //Console.WriteLine(GC.GetTotalMemory(false) + " memory after loadCBT");

            if (merSize == 0)
            {
                Console.WriteLine("bad k-mer size found at start of .cbt file");
                return;
            }

            MerStrings.Initialise(merSize);

            if (tracing >= traceChanges)
            {
                traces.WriteLine("loaded " + loadedUniqueMers + "/" + loadedTotalMers + " avg depth = " + averageDepth);
                traces.WriteLine();
            }

            // and load the pairs file if it exists
            string pairsFN = cbtFN.Replace(".cbt", ".prs");

            if (File.Exists(pairsFN))
            {
                FileInfo pairsFI = new FileInfo(pairsFN);
                const int bytesPerPair = (64 + 32) / 8;
                long pairsFileLength = pairsFI.Length;
                long pairsArrayLength = pairsFileLength / bytesPerPair;

                // allocate the pairs table 
                merPairs = new MerTable<int>(pairsArrayLength);

                BinaryReader pairsFile = new BinaryReader(File.Open(pairsFN, FileMode.Open, FileAccess.Read, FileShare.Read));
                merPairGap = pairsFile.ReadInt32();
                merPairLength = merPairFragmentLength + merPairGap + merPairFragmentLength;

                bool EOF = false;
                long merPairsLoaded = 0;
                long totalPairsRead = 0;
                while (!EOF)
                {
                    try
                    {
                        ulong merPair = pairsFile.ReadUInt64();
                        int merPairDepth = pairsFile.ReadInt32();
                        totalPairsRead++;

                        if (merPairDepth > minLoadReps)
                        {
                            merPairs.Add(merPair, merPairDepth);
                            merPairsLoaded++;
                        }

                    }
                    catch (EndOfStreamException e)
                    {
                        EOF = true;
                    }
                    if (EOF)
                        break;
                }

                pairsFile.Close();
                bool prsLoaded = merPairs.LoadFinished();

                if (prsLoaded)
                    Console.WriteLine(merPairsLoaded + " mer pairs loaded from " + pairsFN);
                else
                {
                    Console.WriteLine(".prs load failed - check file is sorted");
                    return;
                }

                GC.Collect();
                //Console.WriteLine(GC.GetTotalMemory(false) + " memory after loadPairs");
            }

            lowComplexityTrap = new HashSet<ulong>();
            AddToLowComplexityMask("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", lowComplexityTrap, lowComplexityTrapMask);
            AddToLowComplexityMask("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", lowComplexityTrap, lowComplexityTrapMask);
            AddToLowComplexityMask("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", lowComplexityTrap, lowComplexityTrapMask);
            AddToLowComplexityMask("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", lowComplexityTrap, lowComplexityTrapMask);
            AddToLowComplexityMask("ACACACACACACACACACACACACACACACAC", lowComplexityTrap, lowComplexityTrapMask);
            AddToLowComplexityMask("AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG", lowComplexityTrap, lowComplexityTrapMask);
            AddToLowComplexityMask("ATATATATATATATATATATATATATATATAT", lowComplexityTrap, lowComplexityTrapMask);
            AddToLowComplexityMask("CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG", lowComplexityTrap, lowComplexityTrapMask);
            AddToLowComplexityMask("CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT", lowComplexityTrap, lowComplexityTrapMask);
            AddToLowComplexityMask("GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT", lowComplexityTrap, lowComplexityTrapMask);

            // start the monitor/synchronising thread
            Thread monitorProgress = new Thread(RateReporter);
            monitorProgress.Priority = ThreadPriority.AboveNormal;
            monitorProgress.Start();
            inHealingPhase = true; inLoadingPhase = false;

            // process all the reads by tiling them and replacing low-rep mers with more likely variants.
            DateTime startHealing = DateTime.Now;

            healingThreadParams[] healingParams = new healingThreadParams[noHealers];
            Thread[] healingThreads = new Thread[noHealers];

            threadWaitingRead = new Stopwatch[noHealers];
            threadWaitingWrite = new Stopwatch[noHealers];

            // ready a new thread for each parallel healer
            for (int b = 0; b < noHealers; b++)
            {
                healingParams[b] = new healingThreadParams();
                healingParams[b].threadNumber = b + 1;
                healingParams[b].bufferedReadsFiles = bufferedReadsFiles;
                healingParams[b].healedReads = healedReads;
                healingParams[b].healedQuals = healedQuals;
                healingParams[b].problemReads = problemReads;
                healingParams[b].problemQuals = problemQuals;
                healingParams[b].singleReads = singleReads;
                healingParams[b].singleQuals = singleQuals;
                healingThreads[b] = new Thread(new ParameterizedThreadStart(Program.HealingThread));
                healingThreads[b].Priority = ThreadPriority.BelowNormal;
                healingThreads[b].Name = b.ToString();
                healingThreads[b].Start(healingParams[b]);
                Thread.Sleep(100);
                if (tracing > traceOff)
                {
                    // single-thread when tracing to give reproducible ordering in trace files
                    healingThreads[b].Join();
                    //Console.WriteLine("finished healing thread " + b);
                }
            }

            // and wait for all threads to finish
            if (tracing == traceOff)
                for (int b = 0; b < noHealers; b++)
                {
                    healingThreads[b].Join();
                    healingThreads[b] = null;
                    //Console.WriteLine("finished healing thread " + b);
                }

            foreach (StreamWriter r in healedReads)
                if (r != null)
                    r.Close();
            if (healedQuals != null)
                foreach (StreamWriter q in healedQuals)
                    if (q != null)
                        q.Close();
            if (problemReads != null)
                foreach (StreamWriter r in problemReads)
                    if (r != null)
                        r.Close();
            if (problemQuals != null)
                foreach (StreamWriter q in problemQuals)
                    if (q != null)
                        q.Close();
            if (singleReads != null)
                foreach (StreamWriter s in singleReads)
                    if (s != null)
                        s.Close();
            if (singleQuals != null)
                foreach (StreamWriter q in singleQuals)
                    if (q != null)
                        q.Close();
            if (perf != null)
                perf.Close();

            long totalReadingTime = 0;
            long totalWritingTime = 0;
            for (int b = 0; b < noHealers; b++)
            {
                totalReadingTime += threadWaitingRead[b].ElapsedMilliseconds;
                totalWritingTime += threadWaitingWrite[b].ElapsedMilliseconds;
            }

            DateTime endHealing = DateTime.Now;
            statsFile.WriteLine(myProcessNameAndArgs);
            statsFile.WriteLine();
            statsFile.WriteLine((endHealing - startHealing).TotalSeconds.ToString("#.0") + "\tsecs");
            statsFile.WriteLine(stats.reads + "\treads");
            statsFile.WriteLine(stats.OKReads + "\treads OK");
            statsFile.WriteLine(stats.OKReadsChecked + "\treads OK but checked");
            statsFile.WriteLine(stats.shortReads + "\tshort reads passed");
            statsFile.WriteLine(stats.longReads + "\tlong reads passed");
            statsFile.WriteLine(stats.healedReads + "\thealed");
            statsFile.WriteLine(stats.RCHealedReads + "\thealed in reverse direction");
            statsFile.WriteLine(stats.fwdHealedReads + "\thealed in third pass");
            statsFile.WriteLine(stats.readsNotHealedFully + "\tpartially healed");
            statsFile.WriteLine(stats.readsNotHealedAtAll + "\tnot healed at all");
            statsFile.WriteLine(stats.abandonedNs + "\tabandoned reads (too many Ns)");
            statsFile.WriteLine(stats.abandonedRewriting + "\tabandoned reads (rewriting)");
            statsFile.WriteLine(stats.abandonedTree + "\tabandoned reads (tree size)");
            statsFile.WriteLine(stats.abandonedTime.ToString("F1") + "\tsecs spent on abandoned reads");
            statsFile.WriteLine(stats.unbalancedReads + "\treads skipped (deep unbalanced)");
            statsFile.WriteLine(stats.tooHighDepthReads + "\treads skipped (too deep)");
            statsFile.WriteLine(stats.lowComplexityReads + "\treads skipped (low complexity)");
            statsFile.WriteLine(stats.singletonReads + "\tsingleton reads");
            statsFile.WriteLine(stats.slowHealings + "\tslow healings");
            statsFile.WriteLine(stats.discarded + "\tdiscarded 'not good' reads");
            statsFile.WriteLine(stats.mers + "\tmers");
            statsFile.WriteLine(stats.replacedMers + "\tfixed mers");
            statsFile.WriteLine(stats.fixesByType[(int)FixTypes.fixSub] + "\tsubs");
            statsFile.WriteLine(stats.fixesByType[(int)FixTypes.fixDel] + "\tdels");
            statsFile.WriteLine(stats.fixesByType[(int)FixTypes.fixIns] + "\tins");
            statsFile.WriteLine(stats.cacheResizes + "\tcache resizes");
            statsFile.WriteLine(stats.extended + "\tIllumina reads extended");
            statsFile.WriteLine(timeSpentReading.Elapsed.TotalSeconds.ToString("#.0") + "\tsecs reading reads (exc locks)");
            statsFile.WriteLine(timeSpentWriting.Elapsed.TotalSeconds.ToString("#.0") + "\tsecs writing reads (exc locks)");
            statsFile.WriteLine(((double)totalReadingTime / 1000.0).ToString("#.0") + "\tsecs reading reads (inc locks)");
            statsFile.WriteLine(((double)totalWritingTime / 1000.0).ToString("#.0") + "\tsecs writing reads (inc locks)");

            //for (int i = 0; i < 4; i++)
            //    statsFile.WriteLine("subs " + "\t" + i + "\t" + subCounts[i]);

            if (tracing > traceOff)
                traces.Close();
            statsFile.Close();

            stopMonitor = true;
            monitorProgress.Join();

            Console.WriteLine("finished healing reads in " + (endHealing - startHealing).TotalSeconds.ToString("#.0") + "s");
            Console.WriteLine(stats.reads + " reads");
            Console.WriteLine(stats.OKReads + " reads OK");
            Console.WriteLine(stats.OKReadsChecked + " reads OK but checked");
            Console.WriteLine(stats.shortReads + " short reads passed");
            Console.WriteLine(stats.longReads + " long reads passed");
            Console.WriteLine(stats.healedReads + " healed");
            Console.WriteLine(stats.RCHealedReads + " healed in reverse direction");
            Console.WriteLine(stats.fwdHealedReads + " healed in third pass");
            Console.WriteLine(stats.readsNotHealedFully + " partially healed");
            Console.WriteLine(stats.readsNotHealedAtAll + " not healed at all");
            Console.WriteLine(stats.abandonedNs + " abandoned reads (too many Ns)");
            Console.WriteLine(stats.abandonedRewriting + " abandoned reads (rewriting)");
            Console.WriteLine(stats.abandonedTree + " abandoned reads (tree size)");
            Console.WriteLine(stats.abandonedTime.ToString("F1") + " secs spent on abandoned reads");
            Console.WriteLine(stats.unbalancedReads + " reads skipped (deep unbalanced)");
            Console.WriteLine(stats.tooHighDepthReads + " reads skipped (too deep)");
            Console.WriteLine(stats.lowComplexityReads + " reads skipped (low complexity)");
            Console.WriteLine(stats.singletonReads + " singleton reads");
            Console.WriteLine(stats.slowHealings + " slow healings");
            Console.WriteLine(stats.discarded + " discarded 'not good' reads");
            Console.WriteLine(stats.mers + " mers");
            Console.WriteLine(stats.replacedMers + " fixed mers");
            Console.WriteLine(stats.fixesByType[(int)FixTypes.fixSub] + " subs");
            Console.WriteLine(stats.fixesByType[(int)FixTypes.fixDel] + " dels");
            Console.WriteLine(stats.fixesByType[(int)FixTypes.fixIns] + " ins");
            Console.WriteLine(stats.cacheResizes + " cache resizes");
            Console.WriteLine(stats.extended + " Illumina reads extended");
            //for (int i = 0; i < 4; i++)
            //    Console.WriteLine("subs " + "\t" + i + "\t" + subCounts[i]);
            //Console.ReadLine();

            //trace.Close();
        }

        private static void AddToLowComplexityMask(string lowComplexitySeq, HashSet<ulong> lowComplexityTrap, ulong lowComplexityTrapMask)
        {
            ulong[] lowComplexityVariants = new ulong[merSize * 4];
            ulong lowComplexityMer = MerStrings.CondenseMer(lowComplexitySeq);
            int variantCount = GenerateMerSubVariants(lowComplexityMer, lowComplexityVariants, 0, varySingle);
            for (int v = 0; v < variantCount; v++)
                if (!lowComplexityTrap.Contains(lowComplexityVariants[v] & lowComplexityTrapMask))
                    lowComplexityTrap.Add(lowComplexityVariants[v] & lowComplexityTrapMask);
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

        private static void WriteUsage()
        {
            Console.WriteLine("usage: Blue [-help] [-r run] -m minReps [-f fasta|fastq] [-hp] [-good %] [-t threads] [-fixed] [-variable]" +
                                    " [-output dir] [-paired] [-unpaired] k-merFN readsFNs or patterns (" + version + ")");
        }

        private static void WriteUsageFull()
        {
            Console.WriteLine("usage: Blue <options> <k-mer file> <list of file names or patterns for reads files to be corrected>\n" +
                                "\t[-h|help] : display this usage help and exit\n" +
                                "\t[-r|run name] : tag inserted into corrected reads files, and name used for stats file. Default is corrected_<depth>\n" +
                                "\t -m|min minReps : min k-mer depth for OK when scanning reads (used only when dynamic calculation fails)\n" +
                                "\t[-f|format fasta or fastq] : Blue will try to work the format out from the file so this parameter is optional\n" +
                                "\t[-hp] : If this parameter is set, Blue will check for possible errors at the end of every homopolymer string - intended for 454 and IonTorrent\n" +
                                "\t[-g|good good%] : Just save reads that look good (min % of good k-mers in read).\n" +
                                "\t[-t|threads #threads] : No of threads to use for the correction. Default is 1 thread.\n" +
                                "\t[-fixed] : fixed length reads, corrected reads are always the same length as raw reads (default).\n" +
                                "\t[-v|variable] : variable length reads, may be trimmed or extended during healing (default if -hp set).\n" +
                                "\t[-paired]: force paired reads even if single reads specified (for merged paired reads).\n" +
                                "\t[-unpaired] : don't treat a pair of files as a pair set.\n" +
                                "\t[-s|stats name] : name of stats file (constructed from first seq file name if not provided).\n" +
                                "\t[-o|output outputDir] : directory for corrected reads etc" +
                                "\t k-mer file : .cbt file created by Tessel. The same name is used to find the pairs file (.prs) if it exists\n");
        }

        private static bool CheckForParamValue(int p, int argsLength, string msg)
        {
            if (p == argsLength - 1)
            {
                Console.WriteLine(msg);
                return false;
            }
            return true;
        }


        private static void LoadRefSequences(string refFN, Dictionary<string, string> refSequences)
        {
            //>FUOCO2J01D9UB4 rank=0000217 x=1635.0 y=834.0 length=120
            //=TGCCGCACCTGCTGCGACAGGGTCGGCTGCGAAACGTGCAGCGCCTCGGCGGCGCGGGTGAAGTTGCGATGCTCGGCGACAGCGAGCAGGTAGCGCAGGTGGCGAAGCAGCATGGAAACATCC 

            //>FUOCO2J01EDV9U rank=0000223 x=1681.0 y=1552.0 length=105
            //=CCGCCAGCAACTTGGCTTCCAGCACTTCGTTGCTGTCGTAGACGTCGTAGACGACCTTGATCCCGGTTTCCTTGGTGAACTTCTCCAGGGTGTCCGGCGCGATGTAG 

            StreamReader refs = new StreamReader(refFN);
            bool EOF = false;
            string line;
            string[] parsedLine;
            char[] blankSeparator = new char[] { ' ' };

            while (!EOF)
            {
                line = refs.ReadLine();
                if (line == null)
                    break;

                if (line.Length > 1 && line[0] == '>')
                {
                    // found next header line
                    parsedLine = line.Split(blankSeparator);
                    string readID = parsedLine[0].Substring(1);     // strip off leading > or @
                    string refSequence = "";
                    bool within = true;

                    while (within)
                    {
                        // skip any lines other than the reference (=), until we hit it or are about to hit the next group
                        if (refs.Peek() == '>')
                            break;
                        refSequence = refs.ReadLine();
                        if (refSequence.Length > 1 && refSequence[0] == '=')
                        {
                            // found the first reference sequence for this read
                            refSequences.Add(readID, refSequence);
                            break;
                        }
                    }
                }
            }
        }

        //private static int CalculateMedian(Dictionary<ulong, ulong> uniqueMers)
        //{
        //    // build reps histogram (for sums)
        //    Dictionary<int, int> sumReps = new Dictionary<int, int>(100);

        //    foreach (KeyValuePair<ulong, ulong> kvp in uniqueMers)
        //    {
        //        ulong sumPair = kvp.Value;
        //        int sum = (int)(sumPair >> 32) + (int)(sumPair & 0xFFFFFFFF);

        //        if (sumReps.ContainsKey(sum))
        //            sumReps[sum]++;
        //        else
        //            sumReps.Add(sum, 1);
        //    }

        //    // extract, sort rep histogram
        //    int[] reps = new int[sumReps.Count];
        //    int[] repsReps = new int[sumReps.Count];
        //    int i = 0;

        //    foreach (KeyValuePair<int, int> kvp in sumReps)
        //    {
        //        reps[i] = kvp.Key;
        //        repsReps[i] = kvp.Value;
        //        i++;
        //    }
        //    Array.Sort(reps, repsReps);

        //    return 0;
        //}

        private static void HealingThread(object threadParams)
        {
            //Console.WriteLine(GC.GetTotalMemory(false) + " memory on entry to HealingThread");
            healingThreadParams theseParams = (healingThreadParams)threadParams;
            int healerNumber = theseParams.threadNumber;            // which healing thread is this one?
            BufferedReader[] bufferedReadsFiles = theseParams.bufferedReadsFiles; // the buffered read layer over the (shared) reads files (and qual files if they exist)
            StreamWriter[] healedReads = theseParams.healedReads;   // (if saving) corresponding (shared) streams for healed reads
            StreamWriter[] healedQuals = theseParams.healedQuals;   // and 454 quals if they exist
            StreamWriter[] problemReads = theseParams.problemReads; // (if saving) corresponding (shared) streams for problem reads
            StreamWriter[] problemQuals = theseParams.problemQuals; // and 454 quals if they exist
            StreamWriter[] singleReads = theseParams.singleReads;   // good but unpaired reads (good only)
            StreamWriter[] singleQuals = theseParams.singleQuals;

            // per-thread stats counters - added to global counters under lock control at end
            healingStats threadStats = new healingStats();
            // follower cache used for all reads corrected in this thread
            //followerCache cachedFollowersFinal = new followerCache();
            followerCache cachedFollowersFinal = null;
            // and the per-thread set of available (recycled) merProperties and depths objects
            freeMerProperties = new Queue<merProperties>(100);
            freeDepths = new Queue<Depth[]>(100);
            freeMerVariants = new Queue<ulong[]>(100);
            freeVariantSets = new Queue<List<merProperties>>(100);
            freeSequences = new Queue<Sequence>(100);
            merDepths = new int[500];

            threadWaitingRead[healerNumber - 1] = new Stopwatch();
            threadWaitingWrite[healerNumber - 1] = new Stopwatch();

            int noReadFNs = bufferedReadsFiles.Length;
            bool[] fileActive = new bool[noReadFNs];                        // have not yet reached EOF on this reads file
            bool[][] readValid = new bool[noReadFNs][];                     // is this read valid or were we past EOF?
            int filesStillActive = noReadFNs;                               // how many active reads files are still around

            Sequence[][] headerSet = new Sequence[noReadFNs][];             // a set of read headers
            Sequence[][] readSet = new Sequence[noReadFNs][];               // a set of reads, possibly one from each file (length+char[])
            Sequence[][] qualHeaderSet = new Sequence[noReadFNs][];         // set of headers for the quals
            Sequence[][] qualsSet = new Sequence[noReadFNs][];              // set of quals (in canonical form)
            Sequence[][] healedReadSet = new Sequence[noReadFNs][];         // and the returned corrected reads
            Sequence[][] healedQualSet = new Sequence[noReadFNs][];         // the corresponding set of quals in canonical form (returned from healing)
            bool[][] readPossiblyChanged = new bool[noReadFNs][];           // read was possibly altered - use list form of quals just in case
            bool[][] slowHealing = new bool[noReadFNs][];                   // did this read take too long to heal (used to track performance problems)
            bool[][] readHasProblem = new bool[noReadFNs][];                // does this read still look like an error (used in saveGood)
            bool[][] readWasAbandoned = new bool[noReadFNs][];              // read was abandoned - and not changed at all
            bool[] someReadHasProblem = new bool[batchSize];                // union of readHasProblem for the 'pair'
            bool problemsNeedsFlushing = false;                             // wrote to Problems - need to flush to ensure tail of file is not lost on premature termination

            for (int f = 0; f < noReadFNs; f++)
            {
                fileActive[f] = true;                                   // stays true until EOF
                headerSet[f] = new Sequence[batchSize];
                readSet[f] = new Sequence[batchSize];
                qualHeaderSet[f] = new Sequence[batchSize];
                qualsSet[f] = new Sequence[batchSize];
                healedReadSet[f] = new Sequence[batchSize];
                healedQualSet[f] = new Sequence[batchSize];
                readPossiblyChanged[f] = new bool[batchSize];
                slowHealing[f] = new bool[batchSize];
                readHasProblem[f] = new bool[batchSize];
                readWasAbandoned[f] = new bool[batchSize];
                readValid[f] = new bool[batchSize];

                for (int b = 0; b < batchSize; b++)
                {
                    headerSet[f][b] = new Sequence(defaultHeaderLength);
                    readSet[f][b] = new Sequence(defaultReadLength);
                    qualHeaderSet[f][b] = new Sequence(defaultHeaderLength);
                    qualsSet[f][b] = new Sequence(defaultReadLength);
                    healedReadSet[f][b] = new Sequence(defaultReadLength);
                    healedQualSet[f][b] = new Sequence(defaultReadLength);
                    readValid[f][b] = true;
                }
            }

            //Console.WriteLine(GC.GetTotalMemory(false) + " memory after HealingThread initialisation " + healerNumber);

            // get the next set of reads and correct them 
            while (filesStillActive > 0)
            {
                threadWaitingRead[healerNumber - 1].Start();

                // read the next read from each of the files together - giving a pair of reads if we have paired reads - and do this for a batch to reduce locking overhead
                lock (bufferedReadsFiles)
                {
                    // fetch the next batch of reads 
                    timeSpentReading.Start();

                    for (int f = 0; f < noReadFNs; f++)
                    {
                        if (!fileActive[f])
                        {
                            for (int b = 0; b < batchSize; b++)
                                readValid[f][b] = false;
                            continue;
                        }

                        int readsRead = bufferedReadsFiles[f].ReadReads(batchSize, headerSet[f], readSet[f], qualHeaderSet[f], qualsSet[f]);
                        if (readsRead != batchSize)
                        {
                            fileActive[f] = false;
                            for (int b = readsRead; b < batchSize; b++)
                                readValid[f][b] = false;
                            filesStillActive--;
                        }
                        threadStats.reads += readsRead;
                        progressReads += readsRead;
                    }

                    timeSpentReading.Stop();

                } // lock to ensure synchronised reading from all reads files

                threadWaitingRead[healerNumber - 1].Stop();

                // process the just-acquired batch of reads
                for (int b = 0; b < batchSize; b++)
                {
                    // now have a set of reads (the n'th read from each file, if they exist. So try healing each one in turn. 
                    someReadHasProblem[b] = false;
                    for (int f = 0; f < noReadFNs; f++)
                        if (readValid[f][b])
                        {
                            if (readsFormat == MerStrings.formatFASTQ)
                                MerStrings.ConvertQualsToCanonicalForm(qualsSet[f][b], readsFormat, qualBase);      // FASTA quals were converted when they were read in
                            HealARead(headerSet[f][b], readSet[f][b], qualsSet[f][b], cachedFollowersFinal, threadStats,
                                      out readPossiblyChanged[f][b], healedReadSet[f][b], healedQualSet[f][b],
                                      out slowHealing[f][b], out readHasProblem[f][b], out readWasAbandoned[f][b]);
                            someReadHasProblem[b] = someReadHasProblem[b] | readHasProblem[f][b];
                        }
                }

                // have now healed this batch, so write out each set (unless one of them is 'bad' and we're only writing good read sets

                problemsNeedsFlushing = false;
                threadWaitingWrite[healerNumber - 1].Start();

                lock (healedReads)
                {
                    timeSpentWriting.Start();
                    for (int b = 0; b < batchSize; b++)
                    {
                        for (int f = 0; f < noReadFNs; f++)
                            if (readValid[f][b])
                            {
                                bool discardThisRead = pairedReads ? someReadHasProblem[b] : readHasProblem[f][b];
                                // write read set to healed reads file (always or if 'good')
                                if (saveReads == saveAll || (saveReads == saveGoodOnly && !discardThisRead))
                                    SaveHealedReadAndQual(healedReads[f], healedQuals[f], headerSet[f][b], healedReadSet[f][b], qualHeaderSet[f][b], healedQualSet[f][b]);
                                else
                                    threadStats.discarded++;

                                // this read set/pair was poor, so write all the original reads to 'problems' if requested and any good reads to 'singles'
                                if ((saveSlowReads && slowHealing[f][b]) || (discardThisRead && saveReads == saveGoodOnly))
                                {
                                    if (saveNonGoodReads | saveSlowReads)
                                        SaveHealedReadAndQual(problemReads[f], problemQuals[f], headerSet[f][b], readSet[f][b], qualHeaderSet[f][b], healedQualSet[f][b]);
                                    if (!readHasProblem[f][b])
                                        SaveHealedReadAndQual(singleReads[f], singleQuals[f], headerSet[f][b], healedReadSet[f][b], qualHeaderSet[f][b], healedQualSet[f][b]);
                                    problemsNeedsFlushing = saveSlowReads;
                                }
                            }
                    }

                    timeSpentWriting.Stop();

                    if (saveSlowReads && problemsNeedsFlushing)
                        foreach (StreamWriter s in problemReads)
                            s.Flush();

                } // writing out a batch of healed reads under lock protection

                threadWaitingWrite[healerNumber - 1].Stop();

                if (tracing > traceOff && threadStats.reads > traceReadLimit)                  // only trace this much 
                    break;

            } // end of reading/healing loop

            //Console.WriteLine("read and healed " + threadStats.reads + " on thread " + healerNumber);
            //Console.WriteLine(freeReadContexts.Count + " pooled read contexts");

            lock (statsLock)
            {
                stats.reads += threadStats.reads;
                stats.shortReads += threadStats.shortReads;
                stats.longReads += threadStats.longReads;
                stats.OKReads += threadStats.OKReads;
                stats.OKReadsChecked += threadStats.OKReadsChecked;
                stats.healedReads += threadStats.healedReads;
                stats.readsNotHealedFully += threadStats.readsNotHealedFully;
                stats.readsNotHealedAtAll += threadStats.readsNotHealedAtAll;
                stats.singletonReads += threadStats.singletonReads;
                stats.RCHealedReads += threadStats.RCHealedReads;
                stats.fwdHealedReads += threadStats.fwdHealedReads;
                stats.abandonedReads += threadStats.abandonedReads;
                stats.abandonedTime += threadStats.abandonedTime;
                stats.abandonedNs += threadStats.abandonedNs;
                stats.abandonedRewriting += threadStats.abandonedRewriting;
                stats.abandonedTree += threadStats.abandonedTree;
                stats.unbalancedReads += threadStats.unbalancedReads;
                stats.lowComplexityReads += threadStats.lowComplexityReads;
                stats.tooHighDepthReads += threadStats.tooHighDepthReads;
                stats.mers += threadStats.mers;
                stats.replacedMers += threadStats.replacedMers;
                stats.slowHealings += threadStats.slowHealings;
                stats.discarded += threadStats.discarded;
                stats.extended += threadStats.extended;
                for (int i = 0; i < (int)FixTypes.maxFixTypes; i++)
                    stats.fixesByType[i] += threadStats.fixesByType[i];
            }
        }

        private static void HealARead(Sequence readHeader, Sequence startingRead, Sequence startingQuals, followerCache cachedFollowersFinal, healingStats threadStats,
                                      out bool changesMadeToRead, Sequence healedRead, Sequence healedQuals, out bool slowHealingRead, out bool readHasProblem, out bool readWasAbandoned)
        {
            int mersInHealedRead = 0;                           // how many mers are in the healed read
            int zeroMers = 0;                                   // how many mers in the healed have 'zero' counts (<loadReps)
            int poorMers = 0;                                   // how many mers in the read are still poor after healing
            int readMersReplaced = 0;                           // how many mers were replaced in this read
            ReadState readBeforeHealing;                        // state of read before first (forward) pass
            ReadState readBeforeHealingRC = ReadState.ReadOK;   // state of read before second (reverse) pass 
            ReadState readBeforeHealingFwd = ReadState.ReadOK;  // state of read before third (forward) pass 
            ReadState readAfterHealing;                         // state of read after forwards attempt
            ReadState readAfterHealingRC = ReadState.ReadOK;    // state of read after trying in reverse direction
            ReadState readAfterHealingFwd = ReadState.ReadOK;   // state of read after second forward pass
            bool readWasChanged = false;                        // last attempt at healing changed something
            bool readWasChangedRC = false;                      // last attempt at RC healing changed something
            bool readWasChangedFwd = false;                     // last attempt at Fwd healing changed something
            bool abandoned = false;                             // was last attempt at healing abandoned?
            int firstGoodMer = -1;                              // first good mer found in read - either original or repaired
            int wasTracing = traceOff;                          // remember previous tracing state so it can be reset (used when breaking on a read)

            Stopwatch healReadTimer = new Stopwatch();
            healReadTimer.Start();

            slowHealingRead = false;
            readHasProblem = false;
            readWasAbandoned = false;
            changesMadeToRead = false;

            const bool breakOnRead = false;
            if (breakOnRead && (startingRead.Matches("GCTCAACAACGCCATGAGCCGCGGCCAGGCCAAGGGCGCGGCGGGCGCCCAGGGGATCGCCGACGCGCGCACCAACCGCCCGATCATCCTCCGCCCGCCCCCCGCGCGCGCCCAGCC") ||
                                readHeader.Matches("@SRR617721.7393796 SOLEXA4:20:D0V34ACXX:1:2114:5208:49282 length=101")))
            {
                wasTracing = tracing;
                //tracing = traceRead; // or manually set 4 to for tracing follower code
                tracingThisRead = true;
                //tracing = traceFollowers; 
                Debugger.Break();
                // force tracing for this read if we're not already tracing
                if (wasTracing == traceOff)
                {
                    traces = new StreamWriter(tracesFN); // bug here
                    traces.WriteLine(myProcessNameAndArgs);
                    traces.WriteLine();
                }
            }

            if (startingRead.Length < merSize)                          // ignore reads shorter than a single k-mer
            {
                threadStats.shortReads++;
                readHasProblem = true;
                startingRead.CopyTo(healedRead);
                startingQuals.CopyTo(healedQuals);
                return;
            }

            threadStats.mers += (startingRead.Length - merSize + 1);

            // try healing the read in a forward direction
            TryHealingRead(forwardPass, 0, startingRead, readHeader, startingQuals, healedRead, healedQuals, cachedFollowersFinal,
                                        out readBeforeHealing, out readAfterHealing, out readWasChanged, out abandoned, out firstGoodMer,
                                        ref mersInHealedRead, ref zeroMers, ref poorMers, ref readMersReplaced, threadStats);

            // This attempt at healing wasn't (completely) successful, perhaps the problem can be fixed
            // by healing in the reverse direction (can fix multiple breaks in the first mer).
            // This reverse healing pass will only attempt to fix errors in the first not-good k-mers (the last after reversal). 
            if ((readAfterHealing == ReadState.ReadHasErrors) && firstGoodMer >= 0 && !abandoned)
            {
                int startingPoint = healedRead.Length - (firstGoodMer - 1) - merSize;
                if (firstGoodMer < 0)
                    startingPoint = 0;

                Sequence partiallyHealedRead = GetSequence();
                healedRead.CopyTo(partiallyHealedRead);
                Sequence partiallyHealedQuals = GetSequence();
                healedQuals.CopyTo(partiallyHealedQuals);

                TryHealingRead(reversePass, startingPoint, partiallyHealedRead, readHeader, partiallyHealedQuals, healedRead, healedQuals, cachedFollowersFinal,
                                            out readBeforeHealingRC, out readAfterHealingRC, out readWasChangedRC, out abandoned, out firstGoodMer,
                                            ref mersInHealedRead, ref zeroMers, ref poorMers, ref readMersReplaced, threadStats);

                ReturnSequence(partiallyHealedRead);
                ReturnSequence(partiallyHealedQuals);
            } // second (reverse) pass

            // nothing fixed on first pass, and reverse pass made some fixes but read still broken. 
            // Try another forward pass as the reverse pass was effectively the first pass in this case
            if (!readWasChanged && readWasChangedRC && (readAfterHealingRC == ReadState.ReadHasErrors) && !abandoned)
            {
                if (tracing >= traceChanges)
                    TraceClumpAdd("Try healing again in forward direction");

                zeroMers = 0;
                poorMers = 0;

                Sequence partiallyHealedRead = GetSequence();
                healedRead.CopyTo(partiallyHealedRead);
                Sequence partiallyHealedQuals = GetSequence();
                healedQuals.CopyTo(partiallyHealedQuals);

                TryHealingRead(forwardPass, 0, partiallyHealedRead, readHeader, partiallyHealedQuals, healedRead, healedQuals, cachedFollowersFinal,
                                            out readBeforeHealingFwd, out readAfterHealingFwd, out readWasChangedFwd, out abandoned, out firstGoodMer,
                                            ref mersInHealedRead, ref zeroMers, ref poorMers, ref readMersReplaced, threadStats);

                ReturnSequence(partiallyHealedRead);
                ReturnSequence(partiallyHealedQuals);
            } // third (forward) pass

            // if abandoned, restore read to its starting state
            if (abandoned)
            {
                startingRead.CopyTo(healedRead);
                startingQuals.CopyTo(healedQuals);
            }

            changesMadeToRead = readWasChanged || readWasChangedRC || readWasChangedFwd;
            bool triedToFixRead = (readBeforeHealing == ReadState.ReadHasErrors) | (readBeforeHealingRC == ReadState.ReadHasErrors) | (readBeforeHealingFwd == ReadState.ReadHasErrors);
            readHasProblem = (poorMers > (healedRead.Length - merSize + 1) * saveGoodMaxPoor / 100);

            if (readsFixedLength && changesMadeToRead && healedRead.Length < startingRead.Length)
            {
                ExtendShortFixedRead(healedRead, startingRead.Length, healedQuals);
                threadStats.extended++;
            }



            threadStats.replacedMers += readMersReplaced;

            if (readBeforeHealing == ReadState.ReadOK)
            {
                // reads that did not need healing at all - determined by first TryHealingRead call
                threadStats.OKReads++;
                if (tracing > traceOff)
                    TraceClumpAdd("Read OK");
            }

            if (readBeforeHealing == ReadState.ReadSkipDeep)
            {
                // reads that did not need healing at all - determined by first TryHealingRead call
                threadStats.tooHighDepthReads++;
                if (tracing > traceOff)
                    TraceClumpAdd("Read skipped - too deep to safely correct");
            }

            if (readBeforeHealing == ReadState.ReadSkipLowComplexity)
            {
                // reads that did not need healing at all - determined by first TryHealingRead call
                threadStats.lowComplexityReads++;
                if (tracing > traceOff)
                    TraceClumpAdd("Read skipped - low complexity");
            }

            if (readBeforeHealing == ReadState.ReadSkipUnbalanced)
            {
                // reads that did not need healing at all - determined by first TryHealingRead call
                threadStats.unbalancedReads++;
                if (tracing > traceOff)
                    TraceClumpAdd("Read skipped - deep, unbalanced");
            }

            if (triedToFixRead && !changesMadeToRead)
                if (poorMers == 0)
                {
                    // reads that we tried to heal but turned out OK after all
                    threadStats.OKReadsChecked++;
                    if (tracing > traceOff)
                        TraceClumpAdd("Read OK (but checked over)");
                }
                else
                {
                    threadStats.readsNotHealedAtAll++;
                    if (tracing > traceOff)
                        TraceClumpAdd("Not healed at all");
                }

            if ((readBeforeHealing == ReadState.ReadHasErrors) && (readAfterHealing == ReadState.ReadOK) && readWasChanged)
            {
                // broken reads that were completely healed in the forward pass
                threadStats.healedReads++;
                progressHealedReads++;
                if (tracing > traceOff)
                    TraceClumpAdd("Read healed (forward)");
            }

            if ((readBeforeHealingRC == ReadState.ReadHasErrors) && (readAfterHealingRC == ReadState.ReadOK) && readWasChangedRC)
            {
                // broken reads that were completely healed after trying in the reverse direction
                threadStats.healedReads++;
                threadStats.RCHealedReads++;
                progressHealedReads++;
                if (tracing > traceOff)
                    TraceClumpAdd("Read healed (reverse)");
            }

            if ((readBeforeHealingFwd == ReadState.ReadHasErrors) && (readAfterHealingFwd == ReadState.ReadOK) && readWasChangedFwd)
            {
                // broken reads that were completely healed after the third pass
                threadStats.healedReads++;
                threadStats.fwdHealedReads++;
                progressHealedReads++;
                if (tracing > traceOff)
                    TraceClumpAdd("Read healed on third pass");
            }

            if (abandoned)
            {
                threadStats.abandonedReads++;
                if (readAfterHealing == ReadState.ReadAbandonedRewriting || readAfterHealingRC == ReadState.ReadAbandonedRewriting || readAfterHealingFwd == ReadState.ReadAbandonedRewriting)
                    threadStats.abandonedRewriting++;
                if (readAfterHealing == ReadState.ReadAbandonedTreeSize || readAfterHealingRC == ReadState.ReadAbandonedTreeSize || readAfterHealingFwd == ReadState.ReadAbandonedTreeSize)
                    threadStats.abandonedTree++;
                if (readAfterHealing == ReadState.ReadAbandonedNs || readAfterHealingRC == ReadState.ReadAbandonedNs || readAfterHealingFwd == ReadState.ReadAbandonedNs)
                    threadStats.abandonedNs++;
                readWasAbandoned = true;
            }

            if (readHasProblem && changesMadeToRead && !abandoned)
            {
                threadStats.readsNotHealedFully++;
                if (tracing > traceOff)
                    TraceClumpAdd("Partially healed read");
            }

            // all k-mer depths are below the 'zero' cutoff - so say this read is a singleton
            if (zeroMers == mersInHealedRead)
            {
                threadStats.singletonReads++;
                if (tracing > traceOff)
                    TraceClumpAdd("Read singleton");
            }

            healReadTimer.Stop();
            float timeTakenToHeal = (float)(healReadTimer.ElapsedMilliseconds) / 1000;
            slowHealingRead = timeTakenToHeal > 1.0;
            if (slowHealingRead)
                threadStats.slowHealings++;
            if (abandoned)
                threadStats.abandonedTime += timeTakenToHeal;

            if (tracing > traceOff)
            {
                TraceClumpAdd("healing took " + timeTakenToHeal.ToString("F3") + (slowHealingRead ? " (slow)" : ""));
                WriteTraceClump(traceClump);
                traces.WriteLine(Thread.CurrentThread.Name);
            }

            if (tracing >= traceRead)                       // reset one-off trace-in-depth toggle  
            {
                tracing = wasTracing;
                tracingThisRead = false;
                if (wasTracing == traceOff && traces != null)
                    traces.Close();
            }
        }  // for a single read 

        private static void SaveHealedReadAndQual(StreamWriter healedReads, StreamWriter healedQuals,
                                                  Sequence readHeader, Sequence read, Sequence qualHeader, Sequence quals)
        {
            SaveHealedRead(healedReads, readHeader, read);
            SaveHealedQual(healedQuals, qualHeader, quals);
        }

        private static void SaveHealedRead(StreamWriter healedFile, Sequence readHeader, Sequence healedRead)
        {
            // @1:1:0:686#0/1
            // NTGGAGAATTCTGGATCCTCGGACTAAAACAATAGCAGTTGATTCGCTCACAGTTCTGGAGGCTAGAGGTATGAAA
            // +1:1:0:686#0/1
            // @hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh\hhhhhhhhhhhhhhhhVg[hhhU^hhWhfgVhc^Dhh`_V
            //
            // >FUOCO2J02IPX2K rank=0000260 x=3458.0 y=3626.0 length=70
            // GAGCAGCTCCTCAAGCAACTGCAACTGGATTGAGCAAGGAATGTTCGCAGCTACCCGACT
            // GACCCGTCTT

            if (healedFile == null)
                return;

            // write out the header if it exists
            if (readHeader != null)
                healedFile.WriteLine(readHeader.Bases, 0, readHeader.Length);

            if (readsFormat == MerStrings.formatFNA)
            {
                // if we concatenated shorter source lines (e.g. from 454 reads) when the data was read,
                // we'll now split it again as it's being written
                int m = 0;
                int hrLen = healedRead.Length;
                while (m < hrLen)
                {
                    int wLen = Math.Min(60, hrLen - m);
                    healedFile.WriteLine(healedRead.Bases, m, wLen);
                    m += 60;
                }
            }
            else
                healedFile.WriteLine(healedRead.Bases, 0, healedRead.Length);
        }

        private static void SaveHealedQual(StreamWriter healedQualFile, Sequence qualHeader, Sequence healedQuals)
        {
            // The quals will come in text form (healedQualsString) if they are unchanged.
            // Altered quals will always come in List form. 

            if (healedQualFile == null)
                return;

            // no current qual, but there were some at some stage... happens with regression test 454 data but should otherwise never happen
            if (healedQuals.Length == 0)
                return;

            if (readsFormat == MerStrings.formatFNA)
            {
                if (qualHeader.Length != 0)
                    healedQualFile.WriteLine(qualHeader.Bases, 0, qualHeader.Length);

                // easier to just convert to list and then write out 60 quals/line than to parse the concatenated qual string 
                // if there were known to be no changes to the text-form quals
                List<int> healedQualsList = new List<int>(500);
                for (int i = 0; i < healedQuals.Length; i++)
                    healedQualsList.Add((int)healedQuals.Bases[i]);

                // if we concatenated shorter qual lines (e.g. from 454 reads) when they were being read,
                // we'll now reformat and split them again as they're being written
                int m = 0;
                int qualCount = healedQualsList.Count;
                while (m < qualCount)
                {
                    int wLen = Math.Min(60, qualCount - m);
                    for (int i = m; i < m + wLen; i++)
                    {
                        int nextQual = healedQualsList[i];
                        healedQualFile.Write(nextQual);
                        healedQualFile.Write(' ');
                    }
                    healedQualFile.WriteLine();
                    m += 60;
                }
            }

            if (readsFormat == MerStrings.formatFASTQ)
            {
                healedQualFile.WriteLine(qualHeader.Bases, 0, qualHeader.Length);

                // convert the canonical quals back to either Sanger or Solexa format
                for (int i = 0; i < healedQuals.Length; i++)
                    healedQuals.Bases[i] += (char)qualBase;

                healedQualFile.WriteLine(healedQuals.Bases, 0, healedQuals.Length);
            }
        }

        private static void WriteTraceClump(List<string> lines)
        {
            if (lines.Count == 0)
                return;

            lock (traceLock)
            {
                foreach (string l in lines)
                    traces.WriteLine(l);
                traces.Flush();
            }
            lines.Clear();
        }

        private static void TraceClumpAdd(string s)
        {
            lock (traceLock)
            {
                traceClump.Add(Thread.CurrentThread.Name + "\t" + s);
            }
        }

        private static void TraceRead(string tag, Sequence read)
        {
            string traceLine = "";
            int merCount = 0;
            int[] merPlusDepths;
            int[] merRCDepths;

            traceLine += tag + new string(read.Bases, 0, read.Length) + " (" + read.Length + ")";
            TraceClumpAdd(traceLine);
            traceLine = "";
            GetMerDepthsForTrace(read, out merPlusDepths, out merRCDepths);
            merCount = merPlusDepths.Length;

            for (int m = 0; m < merCount; m++)
                traceLine += (merPlusDepths[m] + "/" + merRCDepths[m] + "\t");
            TraceClumpAdd(traceLine);
        }

        private static string Indent(int n)
        {
            return new string(' ', n * 5);
        }

        private static void GetMerDepthsForTrace(Sequence read, out int[] merPlusDepths, out int[] merRCDepths)
        {
            int readLength = read.Length;
            int mersInRead = readLength - merSize + 1;
            merPlusDepths = new int[mersInRead];
            merRCDepths = new int[mersInRead];

            // read too short to tile for mers
            if (readLength < merSize)
                return;

            bool merIsValid = false;
            ulong lastMer = 0;

            for (int i = 0; i < mersInRead; i++)
            {
                if (merIsValid)
                    merIsValid = MerStrings.CondenseMerIncremental(merSize, lastMer, read, i, out lastMer);
                else
                    merIsValid = MerStrings.CondenseMer(read, i, merSize, out lastMer);
                if (merIsValid)
                {
                    int plusCount = 0;
                    int rcCount = 0;
                    if (FindMerInUniqueMers(lastMer, out plusCount, out rcCount))
                    {
                        merPlusDepths[i] = plusCount;
                        merRCDepths[i] = rcCount;
                    }
                }
                // depth arrays are zero by default so don't bother setting zero counts in all other cases
            }
        }

        // Makes a single pass through a read, trying to heal it whenever it finds that the mer depth
        // has fallen below the 'minReps' target or the +/RC values are unbalanced.
        //
        //      readNeedsHealing is set on return if there are still below-minDepth regions in the read
        //      readWasHealed is set if at least one successful healing was done
        //
        // If the healing is being done in reverse, the string is RCed and only the start (now end) of the read is corrected,
        // and only if the first k-mer is broken
        // 
        private static void TryHealingRead(bool forwardOrReverse, int startingPoint, Sequence startingRead, Sequence readHeader, Sequence startingQuals,
                                             Sequence healedRead, Sequence healedQuals, followerCache cachedFollowersFinal,
                                             out ReadState readBeforeHealing, out ReadState readAfterHealing, out bool readWasChanged, out bool readWasAbandoned,
                                             out int firstGoodMer, ref int merCount, ref int zeroMers, ref int poorMers,
                                             ref int replacedMers, healingStats threadStats)
        {
            Stopwatch healReadTimer = null;
            if (perfTrace)
            {
                healReadTimer = new Stopwatch();
                healReadTimer.Start();
            }

            string referenceRead = null;                            // reference read (for tracing)
            int minDepth = 0;                                       // minimum depth for an acceptable mer
            int OKDepth = 0;                                        // target depth for a 'good' mer in this read
            CheckResult checkReason = CheckResult.checkNone;        // why are we looking for a replacement mer?
            int belowOKCount = 0;                                   // how many mers fall below the OK mark
            ReadState readChecked;                                  // do we need to try healing this read?

            int followerRepairs = maxFollowerRepairs;               // starting value for the number of follower repairs allowed
            char repBase = '-';                                     // last base of replacement mer (packed form)
            bool cacheFinal = false;                                // not used here (used in CountFollowers to control caching result from TryHealingMer)
            char exceptionTag = ' ';                                // read failed with an exception - record what type for tracing

            char healedMark = ' ';                                  // char form of repT ype for tracing
            FixTypes fixType = FixTypes.fixNone;                    // final type of change to a mer

            followerCache cachedFollowersAbandoned = null;          // cache of abandoned following attempts (re-built for each starting mer)

            int startingReadLength = startingRead.Length;           // starting length - used to ensure fixed-length reads don't grow
            int startingMerCount = startingReadLength - merSize + 1;    // starting length in mers
            int maxConsecutiveFixesAllowed = 0;                     // how many 'consecutive' fixes allowed before we abandon the read because we're just rewriting it (set dynamically)
            int tailOfRead = Math.Max((startingMerCount - 10), (startingMerCount * 8) / 10);                // tail of read is last 20% or 10 (allow consecutive fixes here)
            bool highDepthRead;                                     // read contains high depth kmers 
            bool unbalancedRead;                                    // read is unbalanced (one strand 100x other strand)
            bool lowComplexityRead;                                 // read contains a low-complexity region (cause performance problems)
            bool zeroDepthsPresent;                                 // initial scan found some zero depth k-mers (on one strand or both)
            bool subFixesOnly = false;                              // set whenever it seems sensible to restrict fix types for performance reasons (e.g. high depth reads)
            int startOfFixRun = -1;                                 // start of run of 'consecutive' fixes

            readWasChanged = false;                                 // assume that the read is not altered (until we alter it)
            firstGoodMer = -1;                                      // no good mers to start
            readWasAbandoned = false;                               // assume read is not abandoned
            int thmCalls = 0;                                       // start counting again for each THR call

            // reinitialise cache for next read
            if (cachedFollowersFinal != null)
            {
                cachedFollowersFinal.cachedKeys.Clear();
                cachedFollowersFinal.nextCachedFollower = 0;
            }

            // copy the starting read across as the 'healed' read
            startingRead.CopyTo(healedRead);
            startingQuals.CopyTo(healedQuals);

            if (startingMerCount > merDepths.Length)
            {
                Array.Resize<int>(ref merDepths, startingMerCount + 50);
            }

            // get the depths from this read and see if we need to have a closer look at it
            // get the mer depths for this read
            merCount = GetMerDepths(startingRead, merDepths, out unbalancedRead, out lowComplexityRead, out zeroDepthsPresent);
            // and dynamically calculate the various depth stats 
            readChecked = AnalyseMerDepths(merCount, merDepths, averageDepth, unbalancedRead, lowComplexityRead, zeroDepthsPresent, out minDepth, out OKDepth, out belowOKCount, out highDepthRead);

            // and remember if it looked to be broken (for reporting back up to the caller)
            readBeforeHealing = readChecked;
            readAfterHealing = readChecked;

            // revert to sub-only fixes if we look like we might get into some performance problems with too many variants spawned from high-rep k-mers
            subFixesOnly = highDepthRead | unbalancedRead;

            // nothing to heal in this read - all mers are adequately covered already
            if (readChecked == ReadState.ReadOK)
            {
                if (tracing > traceOff && forwardOrReverse == forwardPass)
                {
                    if (readHeader != null && readHeader.Length != 0)
                        TraceClumpAdd(readHeader.ToString() + " OK=" + OKDepth + " min=" + minDepth);
                    TraceRead("*", startingRead);
                }
                readBeforeHealing = ReadState.ReadOK;
                readAfterHealing = ReadState.ReadOK;
                readWasChanged = false;

                firstGoodMer = 0;
                return;
            }

            //trace.WriteLine(faHeader + (forwardOrReverse==reversePass ? " reverse" : " forward"));

            // skipping because the read contained low-complexity regions (and these will only result in abandoning the read later if we do try to heal it)
            // these will often have generate excessive depth reads as well - but we'll report this cause rather than 'too deep'
            if (readChecked == ReadState.ReadSkipLowComplexity)
            {
                if (tracing > traceOff)
                {
                    if (readHeader != null && readHeader.Length != 0)
                        TraceClumpAdd(readHeader.ToString() + " OK=" + OKDepth + " min=" + minDepth);
                    TraceRead("~", startingRead);
                }
                readBeforeHealing = ReadState.ReadSkipLowComplexity;
                readAfterHealing = ReadState.ReadSkipLowComplexity;
                readWasChanged = false;

                CountPoorZeroMers(startingRead, merDepths, merCount, minDepth, zeroDepthsPresent, out zeroMers, out poorMers);

                firstGoodMer = 0;
                return;
            }

            // skipping because the read has huge unbalanced counts (like Illumina primers)
            if (readChecked == ReadState.ReadSkipUnbalanced)
            {
                if (tracing > traceOff)
                {
                    if (readHeader != null && readHeader.Length != 0)
                        TraceClumpAdd(readHeader + " OK=" + OKDepth + " min=" + minDepth);
                    TraceRead("~", startingRead);
                }
                readBeforeHealing = ReadState.ReadSkipUnbalanced;
                readAfterHealing = ReadState.ReadSkipUnbalanced;
                readWasChanged = false;

                CountPoorZeroMers(startingRead, merDepths, merCount, minDepth, zeroDepthsPresent, out zeroMers, out poorMers);

                firstGoodMer = 0;
                return;
            }

            // skipping because the read was just too deep in places to safely try to correct
            if (readChecked == ReadState.ReadSkipDeep)
            {
                if (tracing > traceOff)
                {
                    if (readHeader != null && readHeader.Length != 0)
                        TraceClumpAdd(readHeader.ToString() + " OK=" + OKDepth + " min=" + minDepth);
                    TraceRead("~", startingRead);
                }
                readBeforeHealing = ReadState.ReadSkipDeep;
                readAfterHealing = ReadState.ReadSkipDeep;
                readWasChanged = false;

                CountPoorZeroMers(startingRead, merDepths, merCount, minDepth, zeroDepthsPresent, out zeroMers, out poorMers);

                firstGoodMer = 0;
                return;
            }

            if (forwardOrReverse == reversePass)
            {
                // we only do a reverse pass if the first k-mer is broken
                bool lowRepVariant = false;
                ulong firstMer;
                bool firstMerIsPoor = false;
                int firstMerDepth = 0;
                bool firstMerUnbalanced = false;
                bool firstMerHasPair = false;

                // generate the first k-mer (if possible) and see if if may need to be corrected still
                bool firstMerValid = MerStrings.CondenseMer(startingRead, 0, merSize, out firstMer);
                if (firstMerValid)
                {
                    LookupMer(firstMer, highDepthRead, OKDepth, out firstMerDepth, out lowRepVariant, out firstMerUnbalanced);
                    firstMerIsPoor = MerIsHealingCandidate(null, 0, firstMer, startingRead.Bases[merSize], firstMerDepth, lowRepVariant,
                                                           firstMerUnbalanced, minDepth, OKDepth, 0, out checkReason, out firstMerHasPair);
                }
                else
                    checkReason = CheckResult.checkPoor;

                // need to correct the read if either the first mer is unable to be converted to binary form or if if fails the normal goodness tests
                firstMerIsPoor = firstMerIsPoor | !firstMerValid;

                // don't do a reverse healing if the mer looks OK
                if (firstMerIsPoor && checkReason == CheckResult.checkAlt)
                    firstMerIsPoor = false;

                if (!firstMerIsPoor)
                {
                    readBeforeHealing = ReadState.ReadHasErrors;               // only get to do a second pass if the read needed repairing
                    readAfterHealing = ReadState.ReadHasErrors;                // and we're not doing anything about it because we're only trying to fix the first k-mer
                    readWasChanged = false;
                    firstGoodMer = 0;
                    return;
                }

                MerStrings.ReverseComplement(healedRead);
                if (healedQuals.Length != 0)
                    healedQuals.Reverse();

                if (tracing > traceOff)
                {
                    TraceClumpAdd("Trying to reverse heal poor first k-mer");
                }
            }

            // trace read before we start healing it
            if (tracing > traceOff)
            {
                if (readHeader != null && readHeader.Length != 0)
                {
                    string faHeaderString = readHeader.ToString();
                    TraceClumpAdd(faHeaderString + " OK=" + OKDepth + " min=" + minDepth);
                    int readIDEnd = faHeaderString.IndexOf(' ');
                    if (readIDEnd < 0)
                        readIDEnd = faHeaderString.Length;
                    string readID = faHeaderString.Substring(1, readIDEnd - 1);
                    if (refSequences.ContainsKey(readID))
                    {
                        referenceRead = refSequences[readID];
                        if (forwardOrReverse == reversePass)
                            referenceRead = "=" + MerStrings.ReverseComplement(referenceRead.Substring(1)).Trim();
                        TraceClumpAdd(referenceRead);
                    }
                }
                TraceRead(">", healedRead);
            }

            int mersInHealedRead = startingPoint;                       // how many mers are in the healed read
            int previousMerDepth = 0;                                   // depth of the previous mer examined/fixed in the scanning loop
            int previousFixM = -1;                                      // location (m) of previous fix - used for counting consecutive fixes
            int consecutiveFixes = 0;                                   // count of consecutive fixes (used to catch rewriting)
            int consecutiveNs = 0;                                      // count of consecutive N fixes - stop rewriting long N regions

            // something may be broken in this read so move along it progressively, trying to fix each mer and 
            // rebuilding the read string as we make fixes. 

            for (int m = startingPoint; m < Math.Min(merCount, startingMerCount); m++)
            {
                merProperties replacementMer = null;                    // stats etc from TryHealingMer call
                ulong packedMer;                                        // the current packed mer
                int merDepth;                                           // the depth sum from the (+/-) pair
                bool lowRepVariant;                                     // looks like an error variant of a high-rep mer?
                bool unbalancedMer;                                     // or suspiciously unbalanced
                fixType = FixTypes.fixNone;                             // assume no change is made to this mer
                string merWithNs = null;                                // saved mer including Ns (for tracing only)
                bool merWasChanged = false;                             // did TryHealingMer make a change?
                char followingBase = ' ';                               // base following this mer
                bool pairWasFound = false;                              // did this starting k-mer have a pair?

                thmCalls = 0;                                           // restart the thmCalls counter afresh

                if (m + merSize < healedRead.Length)
                    followingBase = healedRead.Bases[m + merSize];
                else
                    followingBase = '\0';

                if (!MerStrings.CondenseMer(healedRead, m, merSize, out packedMer))        // next mer must have contained Ns if this failed
                {
                    if (FindBestReplacementForNs(healedRead, m, OKDepth, out packedMer))   // so try to find the best possible replacement
                    {
                        if (tracing > traceOff)
                            merWithNs = healedRead.ToString(m, merSize);
                        string replacedNsMer = MerStrings.ExpandMer(packedMer);
                        for (int r = 0; r < merSize; r++)                           // and replace the N-containing mer
                            healedRead.Bases[m + r] = replacedNsMer[r];
                        readWasChanged = true;
                        fixType = FixTypes.fixN;
                        consecutiveNs++;
                        if (consecutiveNs > maxConsecutiveNFixesAllowed && m < tailOfRead)
                        {
                            readWasAbandoned = true;                           // looks like we're just rewriting this patch of Ns
                            readAfterHealing = ReadState.ReadAbandonedNs;
                            exceptionTag = 'w';
                            if (tracing > traceOff)
                                TraceClumpAdd("read abandoned - rewriting Ns");
                            break;                                          // so we'll give up on it
                        }
                    }
                    else
                    {
                        // can't find a suitable replacement mer, so move on - and it could be fixed in a later reverse pass
                        // (could be a combination of Ns and errors)
                        zeroMers++;
                        poorMers++;
                        mersInHealedRead++;
                        continue;
                    }
                }
                else
                    consecutiveNs = 0;

                if (tracing >= traceRead)
                {
                    bool breakOnMer = false;
                    string traceMer = MerStrings.ExpandMer(packedMer);
                    if (breakOnMer && traceMer == "TCCGGCCAAGGGCCAGGCCAAGGCA")
                    {
                        TraceClumpAdd("Found target mer " + traceMer + " @" + m);
                        Debugger.Break();
                    }
                }

                // look up next mer and retrieve counts
                LookupMer(packedMer, highDepthRead, OKDepth, out merDepth, out lowRepVariant, out unbalancedMer);

                // see if this mer is a candidate for healing -- either broken or worth checking out
                if (MerIsHealingCandidate(healedRead, m, packedMer, followingBase, merDepth, lowRepVariant, unbalancedMer, minDepth, OKDepth, previousMerDepth, out checkReason, out pairWasFound))
                {
                    // abandoned cache
                    if (usingAbandonedCache)
                        if (cachedFollowersFinal != null)                   // start a new cache for abandoned follower attempts (if we're caching at all)
                            cachedFollowersAbandoned = new followerCache(); // entries in this cache depend on starting location, so it is rebuilt every time we move along the read

                    replacementMer = GetMerProperty();                      // properties assigned by TryHealingMer

                    try         // catch runaway tree exploration!
                    {
                        merWasChanged = TryHealingMer(packedMer, merDepth, lowRepVariant, unbalancedMer, pairWasFound, checkReason, followerRepairs, 0, minDepth, OKDepth, ref subFixesOnly,
                                                      healedRead, startingMerCount, m, m, merCount, cachedFollowersFinal, cachedFollowersAbandoned,
                                                      replacementMer, out cacheFinal, true, 0, ref thmCalls);
                    }
                    catch (ApplicationException e)
                    {
                        // must be cache overflow/tree explosion error as that's the only Application Exception thrown
                        if (tracing > traceOff)
                            TraceClumpAdd("tree explosion - read abandoned - " + e.Message);
                        exceptionTag = 'c';
                        readWasAbandoned = true;
                        readAfterHealing = ReadState.ReadAbandonedTreeSize;
                        break;
                    }

                    if (merWasChanged)
                    {
                        fixType = replacementMer.fixType;
                        UpdateHealingRead(fixType, healedRead, m, replacementMer.variant, startingMerCount, ref merCount, ref repBase);
                        replacedMers++;
                        threadStats.fixesByType[(int)fixType]++;

                        // decrement the consecutive fixes count once we're clear of an error-filled patch
                        if ((m - previousFixM) > 3 && (replacementMer.mersToNextFix > 3 || m >= tailOfRead))
                        //if ((m - previousFixM) > 3 || replacementMer.nextFixType != FixTypes.fixAbandon || replacementMer.mersToNextFix > 3)
                        {
                            consecutiveFixes--;
                            if (consecutiveFixes < 0)
                            {
                                consecutiveFixes = 0;
                                startOfFixRun = -1;
                            }
                        }
                        else
                        {
                            //if (m <= tailOfRead && !replacementMer.pairWasFound)    // don't count pair-validated replacements towards the consecutive count
                            if (m <= tailOfRead)
                                consecutiveFixes++;
                            if (startOfFixRun < 0)
                            {
                                startOfFixRun = m;
                                maxConsecutiveFixesAllowed = Math.Min((startingMerCount - m) / 2, startingMerCount / 5);

                                if (maxConsecutiveFixesAllowed > 10 && checkReason == CheckResult.checkAlt)
                                    maxConsecutiveFixesAllowed = 10;
                                if (maxConsecutiveFixesAllowed < 5)
                                    maxConsecutiveFixesAllowed = 5;
                            }
                        }
                        previousFixM = m;

                        if (consecutiveFixes > maxConsecutiveFixesAllowed && startOfFixRun < tailOfRead)
                        {
                            readWasAbandoned = true;                           // looks like we're just rewriting this read
                            readAfterHealing = ReadState.ReadAbandonedRewriting;
                            exceptionTag = 'w';
                            if (tracing > traceOff)
                                TraceClumpAdd("read abandoned - rewriting");
                            break;                                          // so we'll give up on it
                        }
                    }

                    readWasChanged = readWasChanged | merWasChanged;
                    if (merWasChanged)
                        merDepth = replacementMer.depth;
                }

                // mer was either OK or we tried to fix (and either succeeded or failed)
                mersInHealedRead++;
                healedMark = fixTypes[(int)fixType];
                previousMerDepth = merDepth;
                merDepths[m] = merDepth;
                if (merDepth >= OKDepth && firstGoodMer < 0)
                    firstGoodMer = m;

                if (healedQuals.Length != 0 && (merWasChanged || fixType == FixTypes.fixN || (merDepth >= OKDepth && healedQuals.Bases[m + merSize - 1] < replacementQual)))
                {
                    // fixing first mer in the read, so set all of its quals before making any other possible changes
                    if (m == 0)
                    {
                        for (int i = 0; i < merSize; i++)
                            healedQuals.Bases[i] = replacementQual;
                    }
                    // no change but mismatch between merDepth and qual score
                    if (fixType == FixTypes.fixNone)
                    {
                        healedQuals.Bases[m + merSize - 1] = replacementQual;
                    }
                    // fixing a Sub so just replace the qual
                    if (fixType == FixTypes.fixSub)
                    {
                        healedQuals.Bases[m + merSize - 1] = replacementQual;
                    }
                    // repair added a base so we need to shift remaining qual values right by one place
                    if (fixType == FixTypes.fixDel)
                    {
                        healedQuals.Insert(m + merSize - 1, replacementQual);
                    }
                    // repair deleted a base so we need to shift remaining qual values left by one place 
                    if (fixType == FixTypes.fixIns)
                        healedQuals.Remove(m + merSize - 1);
                }

                if (tracing > traceOff && (merWasChanged || fixType == FixTypes.fixN))
                {
                    string traceLine = "";
                    traceLine += healedMark;
                    for (int p = 0; p < m; p++)
                        traceLine += " ";
                    string originalMer;
                    string healedMer;
                    if (fixType == FixTypes.fixN)
                    {
                        originalMer = merWithNs;
                        healedMer = MerStrings.ExpandMer(packedMer);
                        traceLine += healedMer + " for " + originalMer + " @" + m;
                    }
                    else
                    {
                        originalMer = MerStrings.ExpandMer(packedMer);
                        healedMer = MerStrings.ExpandMer(replacementMer.variant);
                        traceLine += healedMer + " (fo=" + replacementMer.goodFollowers + "/" + replacementMer.allFollowers + "/" +
                                     replacementMer.maxFollowers + "/" + replacementMer.healedMerCount +
                                     ", fx=" + replacementMer.fixes + ", s=" + replacementMer.sum + ") for " + originalMer + " @" + m +
                                     " thm=" + thmCalls + " cf=" + consecutiveFixes + " maxCF=" + maxConsecutiveFixesAllowed;
                    }
                    TraceClumpAdd(traceLine);
                }

                if (replacementMer != null)
                    ReturnMerProperty(replacementMer);

            } // for each mer in the read

            // we abandoned the read rather than fixing it
            if (readWasAbandoned)
            {
                readBeforeHealing = readChecked;
                readWasChanged = false;                              // reset this - as we're throwing away all our changes
                if (tracing > traceOff)
                    TraceClumpAdd(exceptionTag + startingRead.ToString());
                //threadStats.abandonedReads++;
                zeroMers = 0;
                poorMers = startingReadLength - merSize + 1;

                if (perfTrace)
                {
                    string healingResult = "abandoned";
                    TracePerf(forwardOrReverse, startingRead.ToString(), readHeader.ToString(), healingResult, healReadTimer, cachedFollowersFinal, thmCalls, poorMers);
                }

                return;
            }

            // if we're doing a reverse-pass, reverse the read again now to get in the same direction as we started
            if (forwardOrReverse == reversePass)
            {
                if (tracing > traceOff)
                    TraceRead("<", healedRead);
                MerStrings.ReverseComplement(healedRead);
                healedQuals.Reverse();
            }

            // if the read has been healed at all, clean up any uncorrected or really poor bases at the end
            if (readWasChanged)
            {
                merCount = mersInHealedRead;

                // trim the (uncorrected) tail from the healing read if we stopped before getting to its end
                if ((healedRead.Length - merSize + 1) > mersInHealedRead)
                {
                    int wantedLength = mersInHealedRead + merSize - 1;
                    int basesToTrim = healedRead.Length - wantedLength;
                    if (tracing >= traceChoices)
                        TraceClumpAdd("trimming uncorrected tail for " + basesToTrim);
                    healedRead.Length = wantedLength;
                    if (startingQuals.Length > 0)
                        healedQuals.Length = wantedLength;
                    merCount -= basesToTrim;
                }

                // don't allow fixed-length data to grow at all
                if (readsFixedLength && healedRead.Length > startingReadLength)
                {
                    int addedBases = healedRead.Length - startingReadLength;
                    healedRead.Length = startingReadLength;
                    healedQuals.Length = startingQuals.Length;
                    merCount -= addedBases;
                }
            } // trimming if read was changed

            if (tracing > traceOff)
            {
                TraceRead("<", healedRead);
                if (referenceRead != null)
                    TraceClumpAdd(referenceRead);
                //string traceLine = "";
                //for (int i = 0; i < merCount; i++)
                //    traceLine += (healedTypes[i] + "\t");
                //TraceClumpAdd(traceLine);
            }

            // need to re-calculate mer depths properly as we only did a partial pass in this case
            if (forwardOrReverse == reversePass)
                merCount = GetMerDepths(healedRead, merDepths, out unbalancedRead, out lowComplexityRead, out zeroDepthsPresent);
            // and recompute the acceptable depths and faulty k-mers to see if we've done with this read
            AnalyseMerDepths(merCount, merDepths, averageDepth, unbalancedRead, lowComplexityRead, zeroDepthsPresent, out minDepth, out OKDepth, out belowOKCount, out highDepthRead);
            CountPoorZeroMers(healedRead, merDepths, merCount, minDepth, zeroDepthsPresent, out zeroMers, out poorMers);

            readAfterHealing = (belowOKCount > 0) ? ReadState.ReadHasErrors : ReadState.ReadOK;

            if (perfTrace)
            {
                if (healReadTimer.ElapsedMilliseconds >= 50)
                {
                    string healingResult = readAfterHealing != ReadState.ReadOK ? "needy" : "OK";
                    TracePerf(forwardOrReverse, startingRead.ToString(), readHeader.ToString(), healingResult, healReadTimer, cachedFollowersFinal, thmCalls, poorMers);
                }
            }

        }

        private static void CountPoorZeroMers(Sequence read, int[] merDepths, int merCount, int poorDepth, bool zeroDepthsPresent, out int zeroMers, out int poorMers)
        {
            zeroMers = 0;
            poorMers = 0;
            // if we've previously detected zero-depth mers in this read (on either strand), we do a more expensive test, otherwise we'll do a cheap one (as we won't have introduced any)
            if (zeroDepthsPresent)
                CountZeroDepthMers(read, poorDepth, out zeroMers, out poorMers);
            else
            {
                for (int i = 0; i < merCount; i++)
                {
                    if (merDepths[i] == 0)
                        zeroMers++;
                    if (merDepths[i] < poorDepth)
                        poorMers++;
                }
            }
        }

        private static void CountZeroDepthMers(Sequence read, int poorDepth, out int zeroMers, out int poorMers)
        {
            int readLength = read.Length;
            int mersInRead = readLength - merSize + 1;

            bool merIsValid = false;
            ulong lastMer = 0;

            zeroMers = 0;
            poorMers = 0;

            for (int i = 0; i < mersInRead; i++)
            {
                if (merIsValid)
                    merIsValid = MerStrings.CondenseMerIncremental(merSize, lastMer, read, i, out lastMer);
                else
                    merIsValid = MerStrings.CondenseMer(read, i, merSize, out lastMer);

                int plusCount = 0;
                int rcCount = 0;
                int sumCount;
                if (merIsValid)
                    FindMerInUniqueMers(lastMer, out plusCount, out rcCount);

                sumCount = plusCount + rcCount;

                // push effectively zero counts down to zero
                if (sumCount > 100 * plusCount)
                    plusCount = 0;
                if (sumCount > 100 * rcCount)
                    rcCount = 0;

                // below 'OK'
                if (sumCount < poorDepth)
                    poorMers++;
                // both strands zero
                if (plusCount == 0 && rcCount == 0)
                    zeroMers++;
                // above OK but one strand zero - count these as poor as well
                if (sumCount >= poorDepth && (plusCount == 0 || rcCount == 0))
                {
                    zeroMers++;
                    poorMers++;
                }

            }

            return;
        }

        //Depth[] merVariantDepths = new Depth[Math.Max(merSize * 4, noMerVariants)];

        private static Depth[] GetDepths(int noMerVariants)
        {
            // freeMerProperties is thread-local so no locking is needed
            Depth[] thisDepth = null;

            if (freeDepths.Count == 0)
            {
                thisDepth = new Depth[noMerVariants];
            }
            else
            {
                thisDepth = freeDepths.Dequeue();
            }

            if (thisDepth.Length < noMerVariants)
                Array.Resize<Depth>(ref thisDepth, noMerVariants);

            return thisDepth;
        }

        private static void ReturnDepths(Depth[] depths)
        {
            freeDepths.Enqueue(depths);
        }

        // ulong[] merVariants = new ulong[merSize * 4]; 

        private static ulong[] GetMerVariants(int noMerVariants)
        {
            // freeMerVariants is thread-local so no locking is needed
            ulong[] thisMerSet = null;

            if (freeMerVariants.Count == 0)
            {
                thisMerSet = new ulong[noMerVariants];
            }
            else
            {
                thisMerSet = freeMerVariants.Dequeue();
            }

            if (thisMerSet.Length < noMerVariants)
                Array.Resize<ulong>(ref thisMerSet, noMerVariants);

            return thisMerSet;
        }

        private static void ReturnMerVariants(ulong[] merVariants)
        {
            freeMerVariants.Enqueue(merVariants);
        }

        private static List<merProperties> GetVariantSet()
        {
            // freeVariantSets is thread-local so no locking is needed
            List<merProperties> thisVariantSet = null;

            if (freeVariantSets.Count == 0)
            {
                thisVariantSet = new List<merProperties>(10);
            }
            else
            {
                thisVariantSet = freeVariantSets.Dequeue();
            }

            return thisVariantSet;
        }

        private static void ReturnVariantSet(List<merProperties> thisVariantSet)
        {
            thisVariantSet.Clear();
            freeVariantSets.Enqueue(thisVariantSet);
        }

        private static Sequence GetSequence()
        {
            // freeReadContexts is thread-local so no locking is needed
            Sequence thisSequence = null;

            if (freeSequences.Count == 0)
            {
                thisSequence = new Sequence(defaultReadLength);
            }
            else
            {
                thisSequence = freeSequences.Dequeue();
            }

            return thisSequence;
        }

        private static void ReturnSequence(Sequence thisSequence)
        {
            freeSequences.Enqueue(thisSequence);
        }

        private static merProperties GetMerProperty()
        {
            // freeMerProperties is thread-local so no locking is needed
            merProperties thisMp = null;

            if (freeMerProperties.Count == 0)
            {
                thisMp = new merProperties();
            }
            else
            {
                thisMp = freeMerProperties.Dequeue();
            }

            return thisMp;
        }

        private static void ReturnMerProperty(merProperties merProperty)
        {
            merProperty.validVariant = true;
            merProperty.markedVariant = false;
            merProperty.fixType = 0;
            merProperty.nextFixType = 0;
            merProperty.mersToNextFix = 0;
            merProperty.perfectFix = false;
            merProperty.variant = 0;
            merProperty.depth = 0;
            merProperty.unbalanced = false;
            merProperty.pairWasFound = false;
            merProperty.sum = 0;
            merProperty.goodFollowers = 0;
            merProperty.allFollowers = 0;
            merProperty.maxFollowers = 0;
            merProperty.healedMerCount = 0;
            merProperty.fixes = 0;
            merProperty.cached = false;
            merProperty.savedGoodFollowers = 0;
            merProperty.savedAllFollowers = 0;
            merProperty.savedMaxFollowers = 0;

            freeMerProperties.Enqueue(merProperty);

        }

        private static void TracePerf(bool forwardOrReverse, string unhealedRead, string faHeader, string result,
                                      Stopwatch healReadTimer, followerCache cachedFollowersFinal, int thmCalls, int poorMers)
        {
            healReadTimer.Stop();
            int cacheSize = 0;
            if (cachedFollowersFinal != null)
                cacheSize = cachedFollowersFinal.cachedKeys.Count;
            lock (perf)
            {
                perf.WriteLine(progressReads + "\t" + faHeader + "\t" + healReadTimer.ElapsedMilliseconds + "\t" + cacheSize + "\t" +
                               " reverse=" + forwardOrReverse + "\t" + result + "\t" + unhealedRead + "\t" + thmCalls + "\t" + poorMers);
            }
        }

        private static bool LookupMer(ulong mer, bool highDepthRead, int OKDepth, out int depth, out bool lowRepVariant, out bool unbalancedDepth)
        {
            int plusCount = 0;
            int rcCount = 0;
            lowRepVariant = false;

            bool merFound = FindMerInUniqueMers(mer, out plusCount, out rcCount);

            // return sum of both plus & RC counts
            depth = plusCount + rcCount;

            // and see if this mer looks like an error variant spin-off from a higher rep primary or if is just unbalanced
            if (depth > 0)
            {
                if (highDepthRead)
                    lowRepVariant = (depth < OKDepth * 50) && MerIsLowRepVariant(mer, depth);
                unbalancedDepth = (plusCount * 1000 / depth > 950) | (rcCount * 1000 / depth > 950);
            }
            else
            {
                lowRepVariant = false;
                unbalancedDepth = false;
            }

            return merFound;
        }

        private static bool FindMerInUniqueMers(ulong mer, out int plusCount, out int rcCount)
        {
            bool foundMer = true;
            ulong rcMer = MerStrings.ReverseComplementPacked(mer);
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

        // Takes a mer string containing one or more Ns and returns the packed form of the
        // variant that has the highest score. Returns false if no such variant could be found.
        private static bool FindBestReplacementForNs(Sequence read, int m, int OKDepth, out ulong bestPackedMer)
        {
            int bestDepth = -1;                     // best depth of all the alternatives

            // don't try replacing the Ns if there are too many of them in the mer - gets very combinatoric
            int nCount = 0;
            for (int i = 0; i < merSize; i++)
                if (read.Bases[m + i] == 'N')
                    nCount++;
            if (nCount > maxNs)
            {
                bestPackedMer = 0;
                bestDepth = 0;
                return false;
            }

            string mer = read.ToString(m, merSize);
            List<string> merVariants = new List<string>();
            bool stillMoreNs = true;
            merVariants.Add(mer);
            Sequence nextMerB = new Sequence(merSize);

            // generate mers with all replacements for the Ns
            while (stillMoreNs)
            {
                string nextMer = merVariants[0];
                nextMerB.CopyFrom(nextMer);
                int nextN = nextMer.IndexOf('N');
                if (nextN < 0)
                    break;
                foreach (char b in "ACGT")
                {
                    nextMerB.Bases[nextN] = b;
                    merVariants.Add(nextMerB.ToString());
                }
                merVariants.RemoveAt(0);

            }

            nextMerB = null;
            bestDepth = -1;
            bestPackedMer = 0;
            int depth;
            bool lowRepVariant;
            bool unbalancedMer;
            ulong tempPackedMer;

            foreach (string mv in merVariants)
            {
                MerStrings.CondenseMer(mv, out tempPackedMer);

                LookupMer(tempPackedMer, false, OKDepth, out depth, out lowRepVariant, out unbalancedMer);

                if (depth > bestDepth && !lowRepVariant && !unbalancedMer)
                {
                    bestDepth = depth;
                    bestPackedMer = tempPackedMer;
                }
            }

            return bestDepth > 0;
        }

        // Generate a set of mers from a read and calculate their read depths. 
        private static int GetMerDepths(Sequence read, int[] merDepths, out bool unbalancedRead, out bool lowComplexity, out bool zeroDepthsPresent)
        {
            int readLength = read.Length;
            int mersInRead = readLength - merSize + 1;
            unbalancedRead = false;
            lowComplexity = false;
            zeroDepthsPresent = false;
            int sumPlus = 0;
            int sumRC = 0;

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
                    merIsValid = MerStrings.CondenseMer(read, i, merSize, out lastMer);
                if (merIsValid)
                {
                    int plusCount = 0;
                    int rcCount = 0;
                    if (FindMerInUniqueMers(lastMer, out plusCount, out rcCount))
                    {
                        merDepths[i] = plusCount + rcCount;
                        sumPlus += plusCount;
                        sumRC += rcCount;
                        if (plusCount == 0 || rcCount == 0)
                            zeroDepthsPresent = true;
                    }
                    else
                        merDepths[i] = 0;
                    lowComplexity |= lowComplexityTrap.Contains(lastMer & lowComplexityTrapMask);
                }
                else
                    merDepths[i] = 0;
            }

            unbalancedRead = (sumPlus > 100 * sumRC) | (sumRC > 100 * sumPlus);

            return mersInRead;
        }

        // Determines whether the read needs healing and sets the initial 'good' (OKDepth) and 'passable' (minDepth) depth levels.
        private static ReadState AnalyseMerDepths(int merCount, int[] merDepths, int averageDepth, bool unbalancedRead, bool lowComplexityRead, bool zeroDepthsPresent,
                                                           out int minDepth, out int OKDepth, out int belowOKCount, out bool highDepthRead)
        {
            ReadState readNeedsHealing = ReadState.ReadOK;
            int noOKMers = 0;                       // no. of mers that are better than OK depth
            int noOKMersForHM = 0;                  // no. of mers used to calcualte harmonic mean (excludes abnormally high depth mers)
            int OKMerAvg = 0;                       // average of counts for these mers
            int maxDepth = 0;

            minDepth = 0;                           // minimum allowable depth (for 'any' followers)
            OKDepth = 0;                            // depth for 'good' followers and detecting places to try error correction
            belowOKCount = 0;                       // count of mers below OK level

            // calculate the average depth for the passable mers
            double invSum = 0.0f;
            for (int m = 0; m < merCount; m++)
            {
                int depth = merDepths[m];
                if (depth >= minReps)
                {
                    if (depth < averageDepth * highDepthFactor)
                        noOKMersForHM++;                                // working out an average of the 'good' mers 
                    noOKMers++;
                    invSum += 1.0f / (double)depth;
                }
                if (depth > maxDepth)
                    maxDepth = depth;
            }

            if (noOKMers > 0)                                           // set 'OK' depth to be 1/3 of average 'good' for this read
            {
                if (noOKMers - noOKMersForHM > noOKMersForHM)           // more 'deep' k-mers than 'normal' ones
                    OKMerAvg = (int)((double)noOKMers / invSum);        // harmonic mean of all 'non-zero' k-mers
                else
                    OKMerAvg = (int)((double)noOKMersForHM / invSum);   // harmonic mean of just the 'normal' k-mres
                OKDepth = OKMerAvg / 3;
                minDepth = OKMerAvg / 6;                                // and 'min' to 1/6 of this average
            }
            else
            {
                OKDepth = averageDepth / 3;                             // so set 'OK' to 33% of overall average
                minDepth = averageDepth / 6;
            }
            if (minDepth < 1)                                           // min depth always > 0
                minDepth = 1;

            highDepthRead = maxDepth > averageDepth * highDepthFactor; // high-depth read (could be a repeat or artefact) - will revert to sub-only fixes
            bool veryDeepRead = maxDepth > averageDepth * veryHighDepthFactor;  // very deep coverage - read best left alone

            int prevMerDepth = 0;
            for (int m = 0; m < merCount; m++)
            {
                if (merDepths[m] < OKDepth)
                    belowOKCount++;
                if (merDepths[m] < minReps)                             // trigger healing if drop below specified minimum
                    readNeedsHealing = ReadState.ReadHasErrors;
                if (!highDepthRead && merDepths[m] < prevMerDepth / 3)  // or we've hit a sudden drop (only if not high-depth)
                    readNeedsHealing = ReadState.ReadHasErrors;
                prevMerDepth = merDepths[m];
            }

            if (veryDeepRead)
            {
                if (unbalancedRead)
                    readNeedsHealing = ReadState.ReadSkipUnbalanced;
                else
                    if (!zeroDepthsPresent)                             // best to skip this read unless there are obvious errors in it
                        readNeedsHealing = ReadState.ReadSkipDeep;
            }
            if (lowComplexityRead)
                readNeedsHealing = ReadState.ReadSkipLowComplexity;


            return readNeedsHealing;
        }

        private static bool MerIsHealingCandidate(Sequence healingRead, int m, ulong packedMer, char followingBase, int depth, bool lowRepVariant,
                                                  bool unbalancedMer, int minDepth, int OKDepth, int previousDepth, out CheckResult checkReason, out bool pairWasFound)
        {
            checkReason = CheckResult.checkNone;
            pairWasFound = false;
            bool checkingHP = errorModel == modelIndelsCommon && MerIsHomopolymerEnd(packedMer, followingBase);
            bool checkingDepth = (depth < previousDepth / 3) && (depth * 10 > previousDepth); // sudden drop but new depth at least 10% of previous one
            bool merPoor = (depth < OKDepth) || lowRepVariant || unbalancedMer;
            bool merBad = (depth < minDepth);

            if (checkingDepth)
            {
                if (!MerHasAlternatives(packedMer, depth, OKDepth))
                    checkingDepth = false;
            }

            if (checkingHP)
                checkReason = CheckResult.checkHP;
            if (checkingDepth)
                checkReason = CheckResult.checkDepth;
            if (merPoor)
                checkReason = CheckResult.checkPoor;
            if (merBad)
                checkReason = CheckResult.checkBad;

            if (checkReason == CheckResult.checkNone)
                if (merPairs != null && healingRead != null)
                {
                    int pairPresent = IsPairPresent(healingRead, m, minDepth, merPairBackOnly);
                    if (pairPresent == pairNotFound)
                        checkReason = CheckResult.checkAlt;
                    if (pairPresent == pairFound)
                        pairWasFound = true;
                }

            return checkReason > CheckResult.checkNone;
        }

        private static bool MerHasAlternatives(ulong packedMer, int depth, int OKDepth)
        {
            int rhsFillSize = 64 - merSize * 2;
            ulong maskedMer = packedMer & (0xffffffffffffffff << (rhsFillSize + 2));
            int depthPlus = 0;
            int depthRC = 0;

            //int noAlternatives = 0;
            int sumDepth = 0;

            // get the summed depths of all the variants
            for (ulong b = 0; b < 4; b++)
            {
                ulong alternateMer = maskedMer | (b << rhsFillSize);
                if (alternateMer == packedMer)
                    sumDepth += depth;
                else
                {
                    FindMerInUniqueMers(packedMer, out depthPlus, out depthRC);
                    sumDepth += depthPlus + depthRC;
                }
            }

            // and say there are viable looking alternatives if the depth of this mer is < 90% of all depths

            return (sumDepth > 0) && ((depth * 100) / sumDepth) < 90;

            //// if we've hit on a huge depth mer, recalculate the number of alternatives, throwing away any low-depth ones
            //if (noAlternatives > 1 && maxDepth > averageDepth * 20)
            //{
            //    noAlternatives = 0;
            //    for (ulong b = 0; b < 4; b++)
            //    {
            //        ulong alternateMer = maskedMer | (b << rhsFillSize);
            //        if (alternateMer == packedMer)
            //            alternateDepth = depth;
            //        else
            //            LookupMer(alternateMer, out alternateDepth, out zeroOnOneStrand, out merQual);
            //        if (alternateDepth >= maxDepth * 80 / 100)
            //            noAlternatives++;
            //    }
            //}

            //return noAlternatives > 1;
        }

        private static int IsPairPresent(Sequence healingRead, int m, int minDepth, bool merPairDirection)
        {
            ulong merPair;                      // a mer pair ending at 'm' (or possibly starting at 'm' if we're close to the start of the read)
            ulong anchorMer;                    // a full k-mer corresponding to the first fragment in the pair (or the last if we're going forwards)
            //textPair = "no pair";
            int merPairDepth = 0;
            bool forwardsOrBackwards = merPairBack;

            // try first to see if we can have a pair finishing with the end of this mer
            if ((m + merSize) >= merPairLength)
            {
                forwardsOrBackwards = merPairBack;
            }
            else
            // can we have a pair starting with the end of this mer (and are we allowed to go forwards anyway?)
            {
                if (merPairDirection == merPairBackOnly)
                    return pairNotPossible;
                if (m + merSize - 1 + merPairGap + merPairFragmentLength >= healingRead.Length)
                    return pairNotPossible;
                forwardsOrBackwards = merPairForward;
            }

            // generate the k-mer pair (separated by a 'gap' gap)
            if (!GenerateMerPairFromRead(healingRead, m, forwardsOrBackwards, out merPair, out anchorMer))
            {
                // found non-ACGT bases 
                //textPair = MerStrings.ExpandMer(merPair, 32);
                return pairNotPossible;
            }

            // say we can't say if the first fragment comes from a k-mer that appears itself to have an error
            int startingPlusCount = 0;
            int startingRCCount = 0;
            FindMerInUniqueMers(anchorMer, out startingPlusCount, out startingRCCount);
            if ((startingPlusCount + startingRCCount) < minReps)
                return pairNotPossible;

            //textPair = MerStrings.ExpandMer(merPair, 32);
            //string textPairRC = MerStrings.ExpandMer(MerStrings.ReverseComplementPacked(merPair, 32),32);

            // look up the pair (in canonical form) in the pairs table, and return not found if it isn't there or if there aren't enough repetitions
            ulong merPairRC = MerStrings.ReverseComplementPacked(merPair, 32);
            if (merPairRC < merPair)
                merPair = merPairRC;
            if (!merPairs.TryGetValue(merPair, out merPairDepth))
                return pairNotFound;
            if (merPairDepth < minDepth)
                return pairNotFound;

            return pairFound;
        }

        private static bool GenerateMerPairFromRead(Sequence read, int m, bool backwardsOrForwards, out ulong merPair, out ulong anchorMer)
        {
            // Generate a canonical packed k-mer pair from the read. 
            // The second 16-mer is normally the last part of the k-mer starting at 'm', the first 16-mer starts (16+gap) bases before this.
            // If 'm' is too close to the start of the read, we'll try for a pair starting with the k-mer starting at 'm'.
            // We also return an anchor k-mer - this is the k-mer that is not the one being checked/corrected by the caller 
            //             
            //           F                               S
            // ..........1234567890123456.......m........1234567890123456......... (for k=25 and gap=16) - normal 'backwards' case
            //                                  -------- current mer ----
            //           ==== anchor mer =========
            //
            // Early-in-read case - less reliable as it includes the more error-prone (and uncorrected) tail of the read. 
            //                   F                               S
            // .........m........1234567890123456................1234567890123456...... 
            //          -------- current mer ----.......==== anchor mer =========

            merPair = 0;
            anchorMer = 0;                          // a full k-mer based on the starting 16-mer - used to detect errors in first part of pair

            // accumulate lowest 16-mer into a ulong
            int startOfFirstFragment = 0;
            // usual case is to go backwards from m into the already corrected/checked start of the read. 
            if (backwardsOrForwards == merPairBack)
                startOfFirstFragment = m + merSize - merPairLength;                 // m+25-48
            else
                startOfFirstFragment = m + (merSize - merPairFragmentLength);       // m+(25-16)

            for (int i = 0; i < merPairFragmentLength; i++)
            {
                long packedBase = MerStrings.BaseCharToInt(read.Bases[startOfFirstFragment + i]);
                //ulong packedBase = (ulong)"ACGT".IndexOf(read[startOfFirstFragment + i]);
                if (packedBase < 0)
                {
                    return false;
                }
                merPair = (merPair << 2) | (ulong)packedBase;
            }

            // extend the anchor fragment to a full k-mer
            if (backwardsOrForwards == merPairBack)
            {
                anchorMer = merPair;
                for (int i = merPairFragmentLength; i < merSize; i++)
                {
                    long packedBase = MerStrings.BaseCharToInt(read.Bases[startOfFirstFragment + i]);
                    //ulong packedBase = (ulong)"ACGT".IndexOf(read[startOfFirstFragment + i]);
                    if (packedBase < 0)
                    {
                        return false;
                    }
                    anchorMer = (anchorMer << 2) | (ulong)packedBase;
                }
                anchorMer = anchorMer << (64 - merSize * 2);
            }

            // and shift in the second (later) 16-mer
            int startOfSecondFragment = 0;
            if (backwardsOrForwards == merPairBack)
                startOfSecondFragment = m + (merSize - merPairFragmentLength);  // m+(25-16)
            else
                startOfSecondFragment = m + merSize + merPairGap;               // m+25+16
            for (int i = 0; i < merPairFragmentLength; i++)
            {
                long packedBase = MerStrings.BaseCharToInt(read.Bases[startOfSecondFragment + i]);
                //ulong packedBase = (ulong)"ACGT".IndexOf(read[startOfSecondFragment + i]);
                if (packedBase < 0)
                {
                    return false;
                }
                merPair = (merPair << 2) | (ulong)packedBase;
            }

            if (backwardsOrForwards == merPairForward)
            {
                int startOfAnchor = startOfSecondFragment - (merSize - merPairFragmentLength);
                for (int i = 0; i < merSize; i++)
                {
                    long packedBase = MerStrings.BaseCharToInt(read.Bases[startOfAnchor + i]);
                    //ulong packedBase = (ulong)"ACGT".IndexOf(read[startOfAnchor + i]);
                    if (packedBase < 0)
                    {
                        return false;
                    }
                    anchorMer = (anchorMer << 2) | (ulong)packedBase;
                }
                anchorMer = anchorMer << (64 - merSize * 2);
            }

            return true;
        }


        private static bool MerIsHomopolymerEnd(ulong packedMer, char followingBaseChar)
        {
            // test whether we're at the end of a homopolymer (....XXX[y] or ....XXXy) (only used for 454 reads)

            // all possible homopolymer ends in packed form
            const ulong a3 = 0x0;               // AAA
            const ulong c3 = 0x15;              // CCC
            const ulong g3 = 0x2a;              // GGG
            const ulong t3 = 0x3f;              // TTT

            ulong shiftedMer = packedMer >> (32 - merSize) * 2;     // right shift mer
            ulong low3Bases = shiftedMer & 0x3f;                    // last 3 bases XXX

            // does current mer end in an HP (and followed by a non-HP base)
            if (low3Bases == a3 || low3Bases == c3 || low3Bases == g3 || low3Bases == t3)
            {
                // mer ends in an HP string
                ulong lastBase = low3Bases & 0x3;
                ulong followingBase = (ulong)MerStrings.BaseCharToInt(followingBaseChar);
                //ulong followingBase = (ulong)("ACGT".IndexOf(followingBaseChar));
                // so return true if the next base in the read isn't an extension of the HP string (ie we're at the end of the HP)
                return (lastBase != followingBase);
            }
            else
            {
                // mer doesn't end in an HP, are we just after the end of an HP (...XXXy)
                shiftedMer = shiftedMer >> 2;       // drop off the last base (we know it isn't the same as the bases before it)
                low3Bases = shiftedMer & 0x3f;      // get the XXX bases (from a possible initial XXXy mer)
                return (low3Bases == a3 || low3Bases == c3 || low3Bases == g3 || low3Bases == t3);
            }
        }

        private static bool MerIsHomopolymer(ulong packedMer)
        {
            if (errorModel == modelMostlySubs)
                return false;

            // test whether we're at the the end of a homopolymer (...XXX). 
            ulong shiftedMer = packedMer >> (32 - merSize) * 2;     // right shifted mer
            ulong low3Bases = shiftedMer & 0x3f;                    // bases XXX
            const ulong a3 = 0x0;
            const ulong c3 = 0x15;
            const ulong g3 = 0x2a;
            const ulong t3 = 0x3f;

            return (low3Bases == a3) | (low3Bases == c3) | (low3Bases == g3) | (low3Bases == t3);
        }

        private static bool MerIsLongHP(ulong packedMer)
        {
            // test whether we're at the the end of a long homopolymer (...XXXXXXXXXXXX). 
            ulong shiftedMer = packedMer >> (32 - merSize) * 2;     // right shifted mer
            ulong low12Bases = shiftedMer & 0xffffff;               // bases XXXXXXXXXXXX
            const ulong a12 = 0x000000;
            const ulong c12 = 0x555555;
            const ulong g12 = 0xaaaaaa;
            const ulong t12 = 0xffffff;

            return (low12Bases == a12) | (low12Bases == c12) | (low12Bases == g12) | (low12Bases == t12);
        }

        // Does this mer look like a low-rep error variant of a high-rep mer?
        private static bool MerIsLowRepVariant(ulong packedMer, int merDepth)
        {
            int merFill = 64 - merSize * 2;
            ulong rshiftedMer = packedMer >> merFill & 0xFFFFFFFFFFFFFFFC;
            int sumVariantCounts = 0;
            // generate all variants of this mer, look them up in the table and sum their counts
            for (ulong b = 0; b <= 3; b++)
            {
                ulong variantMer = (rshiftedMer | b) << merFill;
                int variantCountPlus = 0;
                int variantCountRC = 0;
                FindMerInUniqueMers(variantMer, out variantCountPlus, out variantCountRC);
                sumVariantCounts += variantCountPlus + variantCountRC;
            }
            // and reject the current mer if it's counts are tiny relative to the sum for the variants (3%)
            return (sumVariantCounts > 0) && (((merDepth * 1000) / sumVariantCounts) <= 20);
        }

        // Try healing the current mer (in packedMer). Returns true if the mer was changed. 
        // This function can be recursively called through CountFollowers, controlled by the repairsLeft parameter 
        static bool TryHealingMer(ulong packedMer, int merDepth, bool lowRepVariant, bool unbalanced, bool pairWasFound, CheckResult checkReason, int repairsLeft, int consecutivePoor,
                                  int minDepth, int OKDepth, ref bool subFixesOnly, Sequence healingContext, int readMerCount, int m, int startingM, int merCount,
                                  followerCache cachedFollowersFinal, followerCache cachedFollowersAbandoned,
                                  merProperties bestMerVariant, out bool cacheFinal, bool outer, int indentDepth, ref int thmCalls)
        {
            if (traceClump.Count > 1000)
                WriteTraceClump(traceClump);

            Stopwatch healMerTimer = null;
            if (tracing > traceOff)
            {
                healMerTimer = new Stopwatch();
                healMerTimer.Start();
            }

            thmCalls++;
            if (thmCalls > maxTHMAllowed)
                throw (new ApplicationException("thm=" + thmCalls));

            //List<merProperties> variantSet = new List<merProperties>(10);
            List<merProperties> variantSet = GetVariantSet();

            // the set of mer variants being chosen
            int viableVariantsInSet = 0;                                        // how many variants are viable 
            merProperties startingMer = null;                                   // remember starting values when the current mer is viable

            int bestChoiceIdx = 0;                                              // return this variant as the choice
            bool bestChoiceFound = false;                                       // set whenever we've decided on the best choice
            string bestChoiceReason = null;                                     // and why we chose it (for tracing)
            bool bestChoicePerfect = false;                                     // did we make a 'perfect' choice?

            string tracedStartingMer = null;
            if (tracing > traceOff)
            {
                bestChoiceReason = "";
                tracedStartingMer = MerStrings.ExpandMer(packedMer);
            }
            //if (startingMer == "GATTCCCATCGTCCACCCCATATTC" && tracing == traceFollowers)
            //    Debugger.Break();
            if (tracing == traceFollowers)
                TraceClumpAdd(Indent(indentDepth) + "heal " + tracedStartingMer + " @" + m);

            // always add the unchanged mer as the first variant if it's viable
            if (merDepth >= minDepth && !lowRepVariant)
            {
                if (tracing >= traceRead)
                    TraceFollowers(checkReason, FixTypes.fixNone, m, indentDepth, MerStrings.ExpandMer(packedMer), 0, merDepth, repairsLeft, consecutivePoor, healingContext.ToString(), merCount);

                startingMer = GetMerProperty();
                //startingMer = new merProperties();

                startingMer.variant = packedMer;
                startingMer.fixType = FixTypes.fixNone;
                startingMer.depth = merDepth;
                startingMer.unbalanced = unbalanced;
                startingMer.pairWasFound = pairWasFound;

                CountFollowers(startingMer, healingContext, readMerCount, m, startingM, merCount, repairsLeft, consecutivePoor, minDepth, OKDepth, ref subFixesOnly,
                               cachedFollowersFinal, cachedFollowersAbandoned, indentDepth, ref thmCalls);

                startingMer.sum += merDepth;
                variantSet.Add(startingMer);

                if (tracing >= traceRead)
                {
                    string cachedTag = startingMer.cached ? " (cached)" : "";
                    TraceClumpAdd(Indent(indentDepth) + "start @" + m + " " + MerStrings.ExpandMer(packedMer) +
                                    " sum=" + startingMer.sum + " fo=" + startingMer.goodFollowers + "/" + startingMer.allFollowers + "/" + startingMer.maxFollowers + "/" + startingMer.healedMerCount +
                                    " fx=" + startingMer.fixes +
                                    " with " + repairsLeft + " fixes left" + cachedTag);
                }
            }
            else
                if (tracing >= traceRead)
                {
                    TraceClumpAdd(Indent(indentDepth) + "start @" + m + " " + MerStrings.ExpandMer(packedMer) +
                                    " had no matches");
                }

            if (tracing >= traceRead)
                TraceClumpAdd(Indent(indentDepth) + "Finding variants for " + tracedStartingMer + "@" + m + " with " + repairsLeft + " fixes left");

            // add best (perhaps tied) variants for each repair type to the set of candidate mers/fix types
            //ulong[] merVariants = new ulong[merSize * 4];   // allocate here and pass in to save re-allocations
            ulong[] merVariants = GetMerVariants(merSize * 4);

            int subVars = FindPlausibleVariants(FixTypes.fixSub, varySingle, checkReason, minDepth, OKDepth, ref subFixesOnly, packedMer, healingContext, readMerCount, m, startingM, merCount,
                                                repairsLeft, consecutivePoor, cachedFollowersFinal, cachedFollowersAbandoned, variantSet, indentDepth, merVariants, ref thmCalls);
            // if we have a particularly bad start (no plausible variants or original mer) try a bit harder to catch multiple (Illumina) sub errors in first mer
            // (we could try to catch these on the reverse pass, but there's a chance with shorter Illumina data that we'll not find any good mer to start from...)
            if (subVars == 0 && m == 0)
                subVars = FindPlausibleVariants(FixTypes.fixSub, varyTwo, checkReason, minDepth, OKDepth, ref subFixesOnly, packedMer, healingContext, readMerCount, m, startingM, merCount,
                                                repairsLeft, consecutivePoor, cachedFollowersFinal, cachedFollowersAbandoned, variantSet, indentDepth, merVariants, ref thmCalls);

            // shift to sub-only correction if we know there is a very high depth mer in this read
            if (!subFixesOnly)
            {
                int delVars = FindPlausibleVariants(FixTypes.fixDel, varySingle, checkReason, minDepth, OKDepth, ref subFixesOnly, packedMer, healingContext, readMerCount, m, startingM, merCount,
                                                    repairsLeft, consecutivePoor, cachedFollowersFinal, cachedFollowersAbandoned, variantSet, indentDepth, merVariants, ref thmCalls);
                int insVars = FindPlausibleVariants(FixTypes.fixIns, varySingle, checkReason, minDepth, OKDepth, ref subFixesOnly, packedMer, healingContext, readMerCount, m, startingM, merCount,
                                                    repairsLeft, consecutivePoor, cachedFollowersFinal, cachedFollowersAbandoned, variantSet, indentDepth, merVariants, ref thmCalls);
            }
            // calculate highest possible follower count 
            int possibleFollowerCount = merCount - m - 1;
            // and how many fixes are viable for the first mer (stops us just rewriting the entire read when there are too many errors at the start) 
            int allowableFixes = possibleFollowerCount / 2;
            // used to give a preference to Dels when they are extending HP runs (454-like data only, -hp option set)
            bool HPExtendingDel = false;
            //  ???
            int maxPossibleFollower = readMerCount - (m + 1);

            viableVariantsInSet = variantSet.Count;

            if ((tracing >= traceChoices && outer) || tracing >= traceRead)
            {
                TraceVariants(m, indentDepth, variantSet, possibleFollowerCount, merCount, readMerCount);
            }

            // save the follower counts, see if the result will be cacheable and do anything else we can on this initial pass through the variants
            bool allCompleted = true;
            bool allIncomplete = true;
            bool allShortAbandoned = true;   // all variants result in abandoned paths, starting within a few bases of where we are now
            foreach (merProperties variant in variantSet)
            {
                variant.savedGoodFollowers = variant.goodFollowers;
                variant.savedAllFollowers = variant.allFollowers;
                variant.savedMaxFollowers = variant.maxFollowers;
                if (variant.nextFixType == FixTypes.fixAbandon)
                    allCompleted = false;
                else
                    allIncomplete = false;
                if (variant.fixType == FixTypes.fixDel && MerIsHomopolymer(variant.variant))
                    HPExtendingDel = true;
                allShortAbandoned = allShortAbandoned & (variant.nextFixType == FixTypes.fixAbandon && variant.mersToNextFix <= 2);
            }
            cacheFinal = allCompleted | allIncomplete;
            //cacheFinal = allCompleted;

            // disallow heroic fixes for the first mer - stops us effectively rewriting it (and everything that follows)
            if (m == 0)
            {
                foreach (merProperties variant in variantSet)
                {
                    if (variant.fixes >= allowableFixes && variant.mersToNextFix < 10)    // rewriting half the remaining string, and not just at end
                    {
                        variant.validVariant = false;                                     // mark this alternative to be ignored
                        viableVariantsInSet--;
                    }
                }
            }

            // don't allow del or ins variants if we appear to be in a rewriting patch for Illumina data
            if (allShortAbandoned && errorModel == modelMostlySubs)
            {
                foreach (merProperties variant in variantSet)
                {
                    if (variant.fixType == FixTypes.fixIns || variant.fixType == FixTypes.fixDel)           // drop ins/del variants if we're in a bad Illumina patch
                    {
                        variant.validVariant = false;                                     // mark this alternative to be ignored
                        viableVariantsInSet--;
                    }
                }
            }

            // only one variant looks possible 
            if (viableVariantsInSet == 1)
            {
                bestChoiceFound = true;
                bestChoiceIdx = 0;
                bestChoiceReason = "only";
            }

            // normalise the follower counts by seeing how much the context has changed
            foreach (merProperties variant in variantSet)
            {
                int adjustment = readMerCount - variant.healedMerCount;
                if (adjustment < 0)
                    adjustment = 0;
                if (adjustment > 0)
                {
                    if (variant.goodFollowers > 0)
                        variant.goodFollowers += adjustment;
                    if (variant.goodFollowers > maxPossibleFollower)
                        variant.goodFollowers = maxPossibleFollower;
                    if (variant.allFollowers > 0)
                        variant.allFollowers += adjustment;
                    if (variant.allFollowers > maxPossibleFollower)
                        variant.allFollowers = maxPossibleFollower;
                    if (variant.maxFollowers > 0)
                        variant.maxFollowers += adjustment;
                    if (variant.maxFollowers > maxPossibleFollower)
                        variant.maxFollowers = maxPossibleFollower;
                }
            }

            int highestGoodFollowerCount;
            int highestAllFollowerCount;
            GetHighestFollowers(variantSet, out highestGoodFollowerCount, out highestAllFollowerCount);
            int highestFollowerCount = highestGoodFollowerCount;
            if (highestAllFollowerCount > highestGoodFollowerCount)
                highestFollowerCount = highestAllFollowerCount;

            // adjust fixes count for HP-extending dels
            foreach (merProperties variant in variantSet)
            {
                if (!variant.validVariant)
                    continue;
                // treat HP extending dels that produce the best follower count specially... by reducing their fix count
                if (variant.fixType == FixTypes.fixDel && (variant.goodFollowers == highestFollowerCount || variant.allFollowers == highestFollowerCount) &&
                        MerIsHomopolymer(variant.variant))
                {
                    if (variant.fixes > 0)
                        variant.fixes--;
                }
            }

            // first try find the 'best' fix for this mer - or decide to leave it alone          
            if (viableVariantsInSet > 0)
                bestChoiceFound = ChooseBestVariant(variantSet, viableVariantsInSet, merDepth, minDepth, OKDepth, checkReason, highestFollowerCount, highestGoodFollowerCount,
                                                    possibleFollowerCount, HPExtendingDel, out bestChoiceIdx, out bestChoiceReason, out bestChoicePerfect);

            // no choices left... use the starting variant if it's in the list, otherwise create a starting variant and insert it at the start of the list
            if (!bestChoiceFound)
            {
                if (startingMer == null)
                {
                    startingMer = GetMerProperty();
                    startingMer.fixType = FixTypes.fixNone;
                    startingMer.variant = packedMer;
                    startingMer.sum = merDepth;
                    startingMer.goodFollowers = 0;
                    startingMer.allFollowers = 0;
                    startingMer.healedMerCount = merCount;
                    startingMer.fixes = 0;
                    startingMer.nextFixType = FixTypes.fixNone;
                    startingMer.perfectFix = false;
                    startingMer.mersToNextFix = 0;
                    variantSet.Insert(0, startingMer);
                }

                bestChoiceIdx = 0;
                bestChoiceReason = "orig";
            }

            merProperties chosenVariant = variantSet[bestChoiceIdx];
            bestMerVariant.fixType = chosenVariant.fixType;
            bestMerVariant.variant = chosenVariant.variant;
            bestMerVariant.depth = chosenVariant.depth;
            bestMerVariant.sum = chosenVariant.sum;
            bestMerVariant.goodFollowers = chosenVariant.savedGoodFollowers;
            bestMerVariant.allFollowers = chosenVariant.savedAllFollowers;
            bestMerVariant.maxFollowers = chosenVariant.savedMaxFollowers;
            bestMerVariant.healedMerCount = chosenVariant.healedMerCount;
            bestMerVariant.fixes = chosenVariant.fixes;
            bestMerVariant.nextFixType = chosenVariant.nextFixType;
            bestMerVariant.mersToNextFix = chosenVariant.mersToNextFix;
            bestMerVariant.perfectFix = chosenVariant.perfectFix;
            bestMerVariant.cached = chosenVariant.cached;
            bestMerVariant.pairWasFound = chosenVariant.pairWasFound;
            bestMerVariant.unbalanced = chosenVariant.unbalanced;

            // now return all the merProperty and other cached objects used here 
            foreach (merProperties mp in variantSet)
                ReturnMerProperty(mp);
            ReturnVariantSet(variantSet);
            ReturnMerVariants(merVariants);

            if ((tracing >= traceChoices && outer) || tracing >= traceRead)
            {
                healMerTimer.Stop();
                float timeTakenToHeal = (float)(healMerTimer.ElapsedMilliseconds) / 1000;
                int cacheSize = 0;
                if (cachedFollowersFinal != null)
                    cacheSize = cachedFollowersFinal.cachedKeys.Count;
                if (cachedFollowersAbandoned != null)
                    cacheSize += cachedFollowersAbandoned.cachedKeys.Count;

                TraceClumpAdd(Indent(indentDepth) +
                                "@" + m + " Choose " + bestChoiceReason + " " + MerStrings.ExpandMer(bestMerVariant.variant) + " " +
                                fixNames[(int)bestMerVariant.fixType] + " sum=" + bestMerVariant.sum + " fo=" + bestMerVariant.goodFollowers + "/" +
                                bestMerVariant.allFollowers + "/" + bestMerVariant.maxFollowers + "/" + bestMerVariant.healedMerCount + " " +
                                fixNames[(int)bestMerVariant.nextFixType] + "@" + bestMerVariant.mersToNextFix + " fx=" + bestMerVariant.fixes + (bestMerVariant.perfectFix ? " (perfect) " : "") +
                                " with " + repairsLeft + " fixes left. cached=" + cacheSize + " " + checkNames[(int)checkReason]);
            }

            return (bestMerVariant.fixType != FixTypes.fixNone);

        }

        private static bool ChooseBestVariant(List<merProperties> variantSet, int viableVariantsInSet, int merDepth, int minDepth, int OKDepth, CheckResult checkReason,
                                              int highestFollowerCount, int highestGoodFollowerCount, int possibleFollowerCount, bool HPExtendingDel,
                                              out int bestChoiceIdx, out string bestChoiceReason, out bool bestChoicePerfect)
        {
            bool bestChoiceFound = false;
            int variantsInSet = variantSet.Count;
            int chosenIdx = 0;
            bestChoicePerfect = false;

            // look for a 'perfect' fix - one that takes us to the max followers with no more changes
            int countPerfect = 0;
            int perfectMaxFollowers = -1;
            if (possibleFollowerCount > almostEnd)  // don't try for 'perfect if we're almost at the end
            {
                // if 'none' gives a perfect result, always choose it over any fixes
                bool noneIsPerfect = false;
                // first look for a 'good' path that will take us straight to the end
                for (int v = 0; v < variantsInSet; v++)
                {
                    merProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;

                    if (variant.mersToNextFix >= variant.maxFollowers && variant.nextFixType == FixTypes.fixNone &&
                        (variant.maxFollowers == variant.goodFollowers) && variant.maxFollowers > 0)
                    {
                        countPerfect++;
                        chosenIdx = v;
                        if (variant.maxFollowers > perfectMaxFollowers)      // remember the highest 'perfect' follower value
                            perfectMaxFollowers = variant.maxFollowers;
                        noneIsPerfect = noneIsPerfect | variant.fixType == FixTypes.fixNone;
                        //if (variant.fixType == FixTypes.fixNone && !HPExtendingDel) 
                        //    break;          // stop if 'no change' takes us to the end (and None will come before any changes, so count will be 1)
                    }
                }

                if (countPerfect == 0)
                {
                    // no 'good' perfect, try for a perfect 'all'
                    for (int v = 0; v < variantsInSet; v++)
                    {
                        merProperties variant = variantSet[v];
                        if (!variant.validVariant)
                            continue;

                        if (variant.mersToNextFix >= variant.maxFollowers && variant.nextFixType == FixTypes.fixNone &&
                            (variant.maxFollowers == variant.allFollowers && variant.maxFollowers > 0))
                        //(variant.maxFollowers == variant.allFollowers && variant.maxFollowers > 0 && !variant.zeroOnOneStrand))
                        {
                            countPerfect++;
                            chosenIdx = v;
                            if (variant.maxFollowers > perfectMaxFollowers)      // remember the highest 'perfect' follower value
                                perfectMaxFollowers = variant.maxFollowers;
                            noneIsPerfect = noneIsPerfect | variant.fixType == FixTypes.fixNone;
                            //if (variant.fixType == FixTypes.fixNone && !HPExtendingDel)     
                            //    break;
                        }
                    }
                }

                if (countPerfect > 1)
                {
                    // more than one 'perfect' - is there only one with the longest followers (handling the case where a shortened Ins is blocking something longer)
                    countPerfect = 0;
                    for (int v = 0; v < variantsInSet; v++)
                    {
                        merProperties variant = variantSet[v];
                        if (!variant.validVariant)
                            continue;

                        if (variant.mersToNextFix >= variant.maxFollowers && variant.nextFixType == FixTypes.fixNone &&
                            (variant.maxFollowers == perfectMaxFollowers || variant.maxFollowers == perfectMaxFollowers) && variant.maxFollowers > 0)
                        {
                            countPerfect++;
                            chosenIdx = v;
                        }
                    }
                }

                // exactly one 'perfect' result
                if (countPerfect == 1 | noneIsPerfect)
                {
                    if (noneIsPerfect)
                        bestChoiceIdx = 0;
                    bestChoiceFound = true;
                    bestChoiceIdx = chosenIdx;
                    bestChoiceReason = "perf";
                    bestChoicePerfect = true;
                    return bestChoiceFound;
                }

            }

            // get the best balanced depth from all of the variant variants (not fixNone)
            int bestVariantDepth = 0;
            foreach (merProperties variant in variantSet)
            {
                if (!variant.validVariant)
                    continue;
                if (variant.fixType == FixTypes.fixNone)
                    continue;

                if (variant.depth > bestVariantDepth && !variant.unbalanced)
                    bestVariantDepth = variant.depth;
            }

            // if there is at least one depth-viable balanced read, cull any unbalanced variants (including the starting k-mer)
            if (bestVariantDepth >= minDepth)
            {
                foreach (merProperties variant in variantSet)
                {
                    if (!variant.validVariant)
                        continue;

                    if (variant.unbalanced)
                    {
                        variant.validVariant = false;
                        viableVariantsInSet--;
                    }
                }
            }

            // are any of the variants 'perfect'? If so, delete all the non-perfect ones, and if there's only one left, choose it
            bool perfectFound = false;
            foreach (merProperties variant in variantSet)
            {
                if (!variant.validVariant)
                    continue;

                perfectFound = perfectFound | variant.perfectFix;
            }
            if (perfectFound)
            {
                countPerfect = 0;
                for (int v = 0; v < variantsInSet; v++)
                {
                    merProperties variant = variantSet[v];

                    if (!variant.validVariant)
                        continue;

                    if (!variant.perfectFix)
                    {
                        variant.validVariant = false;
                        viableVariantsInSet--;
                    }
                    else
                    {
                        countPerfect++;
                        chosenIdx = v;
                    }

                }

                if (countPerfect == 1)
                {
                    bestChoiceFound = true;
                    bestChoiceIdx = chosenIdx;
                    bestChoiceReason = "perfRight";
                    bestChoicePerfect = false;
                    return bestChoiceFound;
                }
            }

            // find the best gap and fix counts
            int highestGoodGapCount = 0;
            int highestAllGapCount = 0;
            int lowestGoodFixCount = 9999;
            int lowestAllFixCount = 9999;
            foreach (merProperties variant in variantSet)
            {
                if (!variant.validVariant)
                    continue;
                if (variant.mersToNextFix > highestGoodGapCount && variant.goodFollowers >= highestFollowerCount)
                    highestGoodGapCount = variant.mersToNextFix;
                if (variant.mersToNextFix > highestAllGapCount && variant.allFollowers >= highestFollowerCount)
                    highestAllGapCount = variant.mersToNextFix;
                if (variant.fixes < lowestGoodFixCount && variant.goodFollowers >= highestFollowerCount)
                    lowestGoodFixCount = variant.fixes;
                if (variant.fixes < lowestAllFixCount && variant.allFollowers >= highestFollowerCount)
                    lowestAllFixCount = variant.fixes;
            }

            // can we get this highest follower count @ lowest fixes from the 'good' followers list (and only once)?
            bool foundinGoodListTwice = false;
            if (possibleFollowerCount > almostEnd)
            {
                bool foundInGoodList = false;
                for (int v = 0; v < variantsInSet; v++)
                {
                    merProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;
                    if (variant.goodFollowers >= highestFollowerCount && highestFollowerCount > 0 && variant.fixes == lowestGoodFixCount)
                    {
                        if (foundInGoodList)
                            foundinGoodListTwice = true;
                        foundInGoodList = true;
                        chosenIdx = v;
                    }
                    else
                        variant.markedVariant = true;
                }

                if (foundInGoodList && !foundinGoodListTwice)
                {
                    bestChoiceFound = true;
                    bestChoiceIdx = chosenIdx;
                    bestChoiceReason = "good";
                    return bestChoiceFound;
                }

                // found multiple 'good' variants, so remove any others from the list
                foreach (merProperties variant in variantSet)
                    if (foundinGoodListTwice)
                    {
                        if (variant.markedVariant && variant.validVariant)
                        {
                            variant.validVariant = false;
                            viableVariantsInSet--;
                        }
                    }
                    else
                        variant.markedVariant = false;

            }

            // alternatively, can we get this highest follower count @ lowest fixes from the 'any' followers list (and only once)?
            if (!foundinGoodListTwice && possibleFollowerCount > almostEnd)
            {
                bool foundInAllList = false;
                bool foundinAllListTwice = false;
                for (int v = 0; v < variantsInSet; v++)
                {
                    merProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;
                    if (variant.allFollowers >= highestFollowerCount && highestFollowerCount > 0 && variant.fixes == lowestAllFixCount)
                    //if (variant.allFollowers >= highestFollower && variant.fixes == lowestAllFixCount && !variant.zeroOnOneStrand)
                    //if (allFollowers[v] >= highestFollower && mersToNextFix[v] == highestAllGapCount)
                    {
                        if (foundInAllList)
                            foundinAllListTwice = true;
                        foundInAllList = true;
                        chosenIdx = v;
                    }
                    else
                        variant.markedVariant = true;               // mark for possible culling
                }

                if (foundInAllList && !foundinAllListTwice)
                {
                    bestChoiceFound = true;
                    bestChoiceIdx = chosenIdx;
                    bestChoiceReason = "all";
                    return bestChoiceFound;
                }

                // found multiple 'all' variants, so remove any others from the list
                foreach (merProperties variant in variantSet)
                    if (foundinAllListTwice)
                    {
                        if (variant.markedVariant && variant.validVariant)
                        {
                            variant.validVariant = false;
                            viableVariantsInSet--;
                        }
                    }
                    // and reset mark if we're not culling
                    else
                        variant.markedVariant = false;
            }

            // no clear winner (with the 'best' followers at lowest fixes)

            // try first for the just highest (good/all) follower (disregarding the fixes need to get the followers)
            if (viableVariantsInSet > 0 && possibleFollowerCount > almostEnd)
            {
                int foundFollower = 0;
                for (int v = 0; v < variantsInSet; v++)
                {
                    merProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;

                    if ((variant.goodFollowers >= highestFollowerCount || variant.allFollowers >= highestFollowerCount) && highestFollowerCount > 0)
                    {
                        chosenIdx = v;
                        foundFollower++;
                    }
                    else
                    {
                        variant.markedVariant = true;
                    }

                }
                // and if there's only one variant with the best followers, choose that one
                if (foundFollower == 1)
                {
                    bestChoiceFound = true;
                    bestChoiceIdx = chosenIdx;
                    bestChoiceReason = "foll";
                    return bestChoiceFound;
                }
                else
                    foreach (merProperties variant in variantSet)
                        if (foundFollower > 1)
                        {
                            if (variant.markedVariant && variant.validVariant)
                            {
                                variant.validVariant = false;
                                viableVariantsInSet--;
                            }
                        }
                        else
                            variant.markedVariant = false;
            }

            // more than one with the highest follower count, try next for the highest 'good' follower out of the remaining variants
            if (viableVariantsInSet > 0 && possibleFollowerCount > almostEnd)
            {
                int foundTarget = 0;
                for (int v = 0; v < variantsInSet; v++)
                {
                    merProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;

                    if (variant.goodFollowers >= highestGoodFollowerCount && highestGoodFollowerCount > 0)
                    {
                        chosenIdx = v;
                        foundTarget++;
                    }
                    else
                    {
                        variant.markedVariant = true;
                    }
                }

                // and if there's only one variant with the best followers, choose that one
                if (foundTarget == 1)
                {
                    bestChoiceFound = true;
                    bestChoiceIdx = chosenIdx;
                    bestChoiceReason = "folg";
                    return bestChoiceFound;
                }
                else
                    foreach (merProperties variant in variantSet)
                        if (foundTarget > 1)
                        {
                            if (variant.markedVariant && variant.validVariant)
                            {
                                variant.validVariant = false;
                                viableVariantsInSet--;
                            }
                        }
                        else
                            variant.markedVariant = false;
            }

            // no best choice by looking at followers... try for the longest gap amongst the remaining candidates
            if (viableVariantsInSet > 0 && possibleFollowerCount > almostEnd)
            {
                int highestGaps = 0;
                if (possibleFollowerCount > 0)
                {
                    foreach (merProperties variant in variantSet)
                    {
                        if (!variant.validVariant)
                            continue;
                        if (variant.mersToNextFix > highestGaps)
                            highestGaps = variant.mersToNextFix;
                    }

                    // how many valid variants have this 'longest gap'?
                    int foundTarget = 0;
                    for (int v = 0; v < variantsInSet; v++)
                    {
                        merProperties variant = variantSet[v];
                        if (!variant.validVariant)
                            continue;

                        if (variant.mersToNextFix == highestGaps)
                        {
                            chosenIdx = v;
                            foundTarget++;
                        }
                        else
                        {
                            variant.markedVariant = true;
                        }
                    }

                    // and if there's only one variant with the best gap, choose that one
                    if (foundTarget == 1)
                    {
                        bestChoiceFound = true;
                        bestChoiceIdx = chosenIdx;
                        bestChoiceReason = "gaps";
                        return bestChoiceFound;
                    }
                    else
                        foreach (merProperties variant in variantSet)
                            if (foundTarget > 1)
                            {
                                if (variant.markedVariant && variant.validVariant)
                                {
                                    variant.validVariant = false;
                                    viableVariantsInSet--;
                                }
                            }
                            else
                                variant.markedVariant = false;
                }
            }

            // multiple choices with best followers, best gaps, choose the one with the lowest fixes
            if (!bestChoiceFound && viableVariantsInSet > 0 && possibleFollowerCount > almostEnd)
            {
                int lowestFixes = 9999;
                foreach (merProperties variant in variantSet)
                {
                    if (!variant.validVariant)
                        continue;
                    if (variant.fixes < lowestFixes)
                        lowestFixes = variant.fixes;
                }

                int foundTarget = 0;
                for (int v = 0; v < variantsInSet; v++)
                {
                    merProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;

                    if (variant.fixes == lowestFixes)
                    {
                        chosenIdx = v;
                        foundTarget++;
                    }
                    else
                    {
                        variant.markedVariant = true;
                    }
                }

                // and if there's only one variant with the best fixes, choose that one
                if (foundTarget == 1)
                {
                    bestChoiceFound = true;
                    bestChoiceIdx = chosenIdx;
                    bestChoiceReason = "fixs";
                    return bestChoiceFound;
                }
                else
                    foreach (merProperties variant in variantSet)
                        if (foundTarget > 1)
                        {
                            if (variant.markedVariant && variant.validVariant)
                            {
                                variant.validVariant = false;
                                viableVariantsInSet--;
                            }
                        }
                        else
                            variant.markedVariant = false;
            }

            // multiple choices with best followers, best gaps and best fixes
            // try for status quo if this is a possible choice
            if (variantSet[0].validVariant && variantSet[0].fixType == FixTypes.fixNone && !HPExtendingDel)
            {
                // unchanged mer will always be the first in the lists if it was viable 
                if (merDepth >= OKDepth || (merDepth > minDepth && variantSet[0].maxFollowers <= almostEnd))
                {
                    bestChoiceFound = true;
                    bestChoiceIdx = 0;
                    bestChoiceReason = "none";
                    return bestChoiceFound;
                }
            }

            // poor mer but multiple choices with best followers and best gaps and fewest fixes, try for best sum if one is clearly better than the others
            if (viableVariantsInSet > 0 && checkReason >= CheckResult.checkPoor)
            //if (viableVariantsInSet > 0 && checkReason >= checkPoor && possibleFollowerCount > almostEnd)
            {
                int highestSum = 0;
                foreach (merProperties variant in variantSet)
                {
                    if (!variant.validVariant)
                        continue;
                    if (variant.sum > highestSum)
                        highestSum = variant.sum;
                }

                int foundTarget = 0;
                for (int v = 0; v < variantsInSet; v++)
                {
                    merProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;

                    if (variant.sum >= highestSum * 9 / 10)
                    {
                        chosenIdx = v;
                        foundTarget++;
                    }
                    else
                    {
                        variant.markedVariant = true;
                    }
                }

                // and if there's only one variant with the high-enough sum, choose that one
                if (foundTarget == 1)
                {
                    bestChoiceFound = true;
                    bestChoiceIdx = chosenIdx;
                    bestChoiceReason = "sums";
                    return bestChoiceFound;
                }
                else
                    foreach (merProperties variant in variantSet)
                        if (foundTarget > 1)
                        {
                            if (variant.markedVariant && variant.validVariant)
                            {
                                variant.validVariant = false;
                                viableVariantsInSet--;
                            }
                        }
                        else
                            variant.markedVariant = false;
            }

            // make an almost arbitrary decision - preferring del/ins over sub and del over ins to avoid rewriting too much
            // but preferring sub if we're at the end of a read or we are correcting Illumina data
            if (viableVariantsInSet > 0 && (possibleFollowerCount <= almostEnd || errorModel == modelMostlySubs))
            {
                int bestSubSum = 0;
                bestChoiceIdx = -1;
                for (int v = 0; v < variantsInSet; v++)
                {
                    merProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;

                    // find the Sub variant with the highest sum (if there is one)
                    if (variant.fixType == FixTypes.fixSub && variant.sum > bestSubSum)
                    {
                        bestSubSum = variant.sum;
                        bestChoiceIdx = v;
                    }
                }
                if (bestChoiceIdx >= 0)
                {
                    bestChoiceFound = true;
                    bestChoiceReason = "pSub";
                    return bestChoiceFound;
                }
            }

            if (!bestChoiceFound && viableVariantsInSet > 0)
            {
                int bestDelSum = 0;
                bestChoiceIdx = -1;
                for (int v = 0; v < variantsInSet; v++)
                {
                    merProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;

                    if (variant.fixType == FixTypes.fixDel && variant.sum > bestDelSum)
                    {
                        bestDelSum = variant.sum;
                        bestChoiceIdx = v;
                    }
                }
                if (bestChoiceIdx >= 0)
                {
                    bestChoiceFound = true;
                    bestChoiceReason = "pDel";
                    return bestChoiceFound;
                }
            }

            if (!bestChoiceFound && viableVariantsInSet > 0)
            {
                int bestInsSum = 0;
                bestChoiceIdx = -1;
                for (int v = 0; v < variantsInSet; v++)
                {
                    merProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;

                    if (variant.fixType == FixTypes.fixIns && variant.sum > bestInsSum)
                    {
                        bestInsSum = variant.sum;
                        bestChoiceIdx = v;
                    }
                }
                if (bestChoiceIdx >= 0)
                {
                    bestChoiceFound = true;
                    bestChoiceReason = "pIns";
                    return bestChoiceFound;
                }
            }

            if (!bestChoiceFound && viableVariantsInSet > 0)
            {
                // no del or ins in list but list isn't empty 
                for (int v = 0; v < variantsInSet; v++)
                {
                    merProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;

                    if (variant.fixType == FixTypes.fixSub)
                    {
                        bestChoiceFound = true;
                        bestChoiceIdx = v;
                        bestChoiceReason = "fSub";
                        return bestChoiceFound;
                    }
                }
            }

            // couldn't find a suitable replacement
            bestChoiceIdx = -1;
            bestChoiceReason = "null";
            return false;
        }

        private static void GetHighestFollowers(List<merProperties> results, out int highestGoodFollower, out int highestAllFollower)
        {
            highestGoodFollower = 0;
            highestAllFollower = 0;
            for (int v = 0; v < results.Count; v++)
            {
                merProperties result = results[v];
                if (!result.validVariant)
                    continue;
                int follower = result.goodFollowers;
                if (follower > highestGoodFollower)
                    highestGoodFollower = follower;
                follower = result.allFollowers;
                if (follower > highestAllFollower)
                    highestAllFollower = follower;
            }
        }

        private static void TraceVariants(int m, int indentDepth, List<merProperties> variantSet,
                                          int maxPossible, int merCount, int readMerCount)
        {
            for (int v = 0; v < variantSet.Count; v++)
                TraceClumpAdd(Indent(indentDepth) + "@" + m + " " +
                                "pv " + v + " " + fixNames[(int)variantSet[v].fixType] + " s=" + variantSet[v].sum +
                               " fo=" + variantSet[v].goodFollowers + "/" + variantSet[v].allFollowers + "/" + variantSet[v].maxFollowers + "/" + variantSet[v].healedMerCount +
                               " fx=" + variantSet[v].fixes + " nf=" + fixNames[(int)variantSet[v].nextFixType] + "@" + variantSet[v].mersToNextFix +
                               " mxp=" + maxPossible + " mc=" + merCount + " rmc=" + readMerCount + " " + variantSet[v].perfectFix + " " +
                               MerStrings.ExpandMer(variantSet[v].variant) + (variantSet[v].cached ? " (cached)" : ""));
        }

        // Counts the followers for a mer. 
        // Possibly temporarily making a few repairs if needed to overcome nearly adjacent errors (recursively calls TryHealingMer).
        // These temporary fixes are cached to avoid repeated recalculations. 
        // 
        private static void CountFollowers(merProperties variantMer, Sequence healingContext, int readMerCount,
                                           int m, int startingM, int merCount, int repairsLeft, int consecutivePoor, int minDepth, int OKDepth, ref bool subFixesOnly,
                                           followerCache cachedFollowersFinal, followerCache cachedFollowersAbandoned, int indentDepth, ref int thmCalls)
        {
            Stopwatch countFollowersTimer = null;
            if (tracing > traceOff)
            {
                countFollowersTimer = new Stopwatch();
                countFollowersTimer.Start();
            }

            int remainingLength = healingContext.Length - m - merSize;
            int followerKeyLength = remainingLength;
            if (followerKeyLength > merSize)
                followerKeyLength = merSize;

            bool cacheFinal = true;                            // this result is 'final' and can be cached (not a mix of abandoned and completed paths)

            //ulong followerCacheKeyFinal = 0;
            string followerCacheKeyFinal = null;
            if (cachedFollowersFinal != null)
                //followerCacheKeyFinal = ((ulong)healingContext.ComputeHash() << 32) | ((ulong)m << 12) | ((ulong)variantMer.fixType << 6) + (ulong)repairsLeft;
                //followerCacheKeyFinal = ((ulong)healingContext.Bases.ToString().GetHashCode() << 32) | ((ulong)m << 12) | ((ulong)variantMer.fixType << 6) + (ulong)repairsLeft;
                followerCacheKeyFinal = (healingContext.Bases.ToString() + m + variantMer.fixType + repairsLeft);

            //if (cachedFollowersAbandoned != null)
            //    followerCacheKeyAbandoned = followerCacheKeyFinal + "+" + startingM;

            if (cachedFollowersFinal != null)
            {
                // look up the cache key in the relevant cache, remember which cache matched and the index inside the cache arrays
                int cacheIdx = -1;
                followerCache cachedFollowers = null;
                if (cachedFollowersFinal.cachedKeys.ContainsKey(followerCacheKeyFinal))
                {
                    cacheIdx = cachedFollowersFinal.cachedKeys[followerCacheKeyFinal];
                    cachedFollowers = cachedFollowersFinal;
                }
                //else
                //    if (cachedFollowersAbandoned != null)
                //        if (cachedFollowersAbandoned.cachedKeys.ContainsKey(followerCacheKeyAbandoned))
                //        {
                //            cacheIdx = cachedFollowersAbandoned.cachedKeys[followerCacheKeyAbandoned];
                //            cachedFollowers = cachedFollowersAbandoned;
                //        }

                // found this entry in one of the caches, so return the follower values immediately
                if (cacheIdx >= 0)
                {
                    variantMer.goodFollowers = cachedFollowers.cachedCountGoodFollowers[cacheIdx];
                    variantMer.allFollowers = cachedFollowers.cachedCountAllFollowers[cacheIdx];
                    variantMer.maxFollowers = cachedFollowers.cachedCountMaxFollowers[cacheIdx];
                    variantMer.healedMerCount = cachedFollowers.cachedHealedMerCount[cacheIdx];
                    variantMer.fixes = cachedFollowers.cachedCountFixes[cacheIdx];
                    variantMer.sum = cachedFollowers.cachedFixedSum[cacheIdx];
                    variantMer.mersToNextFix = cachedFollowers.cachedMersToFirstFix[cacheIdx];
                    variantMer.nextFixType = cachedFollowers.cachedFirstFixType[cacheIdx];
                    variantMer.perfectFix = cachedFollowers.cachedPerfectFix[cacheIdx];
                    variantMer.cached = true;
                    if (tracing >= traceRead)
                        TraceFollowerResults(indentDepth, m, variantMer.fixType, true, countFollowersTimer, cachedFollowersFinal, cachedFollowersAbandoned);
                    return;
                }
            }

            FixTypes locFirstFixType = FixTypes.fixNone;
            int locMersToFirstFix = 0;
            int locGoodFollowerCount = 0;
            int locAllFollowerCount = 0;
            int locMaxFollowerCount = 0;
            int locHealedMerCount = Math.Min(merCount, readMerCount);               // use this count if the follower loop terminates here
            int locFixesCount = 0;
            int locFollowerSum = 0;
            bool locPerfectFix = false;
            merProperties replacementMer = null;

            int previousDepth = variantMer.depth; ;

            ulong nextMer = variantMer.variant;
            int nextM = m + 1;
            char followingBase;

            while (GenerateNextMer(nextMer, healingContext, nextM, startingM, merCount, readMerCount, OKDepth, out nextMer, out followingBase))   // generate the next mer 
            {
                string watchMer = null;
                if (tracing >= traceRead)
                    watchMer = MerStrings.ExpandMer(nextMer);

                int nextMerDepth;
                bool nextLowRepVariant;
                bool nextMerUnbalanced;
                CheckResult checkReason = CheckResult.checkNone;
                bool nextMerHasPair;

                LookupMer(nextMer, subFixesOnly, OKDepth, out nextMerDepth, out nextLowRepVariant, out nextMerUnbalanced);

                // if the next mer (current mer + next base from read) looks OK, keep going 
                if (!MerIsHealingCandidate(healingContext, nextM, nextMer, followingBase, nextMerDepth, nextLowRepVariant, nextMerUnbalanced, minDepth, OKDepth, previousDepth,
                                           out checkReason, out nextMerHasPair))
                {
                    locGoodFollowerCount++;
                    locAllFollowerCount++;
                    locMaxFollowerCount++;
                    locFollowerSum += nextMerDepth;
                    locMersToFirstFix++;
                    nextM++;
                    previousDepth = nextMerDepth;
                    consecutivePoor = 0;
                    continue;
                }

                // reduce the number of repairs allowed if we're not just checking for an HP improvement
                if (checkReason == CheckResult.checkBad)
                    repairsLeft--;

                if (checkReason == CheckResult.checkPoor)
                {
                    consecutivePoor++;
                    if (consecutivePoor > maxConsecutivePoor)
                        repairsLeft--;
                }

                if (checkReason != CheckResult.checkBad && checkReason != CheckResult.checkPoor)
                    consecutivePoor = 0;

                // ran out of good followers (or we've hit an HP), perhaps the mer can/should be repaired and we can keep going.
                // Only a few repairs are allowed to stop us just rewriting the read from the consensus (and to limit the tree exploration depth)

                if (repairsLeft > 0)
                //if (repairsLeft > 0 && (nextM < (startingM + merSize)))
                {
                    // best replacement mer - returned from TryHealingMer
                    replacementMer = GetMerProperty();
                    //replacementMer = new merProperties();
                    bool merWasChanged = false;

                    merWasChanged = TryHealingMer(nextMer, nextMerDepth, nextLowRepVariant, nextMerUnbalanced, nextMerHasPair, checkReason, repairsLeft, consecutivePoor,
                                                  minDepth, OKDepth, ref subFixesOnly, healingContext, readMerCount, nextM, startingM, merCount,
                                                  cachedFollowersFinal, cachedFollowersAbandoned,
                                                  replacementMer, out cacheFinal, false, indentDepth + 1, ref thmCalls);

                    // did our attempt to repair a mer end up with the same mer as we started with?
                    if (replacementMer.variant == variantMer.variant && merWasChanged)
                    {
                        if (tracing >= traceRead)
                            TraceClumpAdd("Killing follower path at potential loop @" + nextM + " " +
                                            MerStrings.ExpandMer(nextMer));
                        nextMerDepth = 0;
                        replacementMer.allFollowers = 0;
                        replacementMer.goodFollowers = 0;
                        replacementMer.maxFollowers = 0;
                        replacementMer.fixes = 99;
                        replacementMer.mersToNextFix = 99;
                        replacementMer.perfectFix = false;
                        replacementMer.fixType = FixTypes.fixAbandon;
                    }

                    // if the corrected mer has a viable score and doesn't look to be an error, record *this* repair as producing a follower
                    if (replacementMer.depth >= minDepth)
                    {
                        if (replacementMer.depth >= OKDepth)
                            locGoodFollowerCount++;
                        locAllFollowerCount++;
                    }

                    // and add in the downstream followers to the follower counts we got from the loop prior to calling TryHealingMer
                    locGoodFollowerCount += replacementMer.goodFollowers;
                    locAllFollowerCount += replacementMer.allFollowers;
                    locMaxFollowerCount += replacementMer.maxFollowers + 1;
                    locFixesCount = replacementMer.fixes;
                    if (merWasChanged)
                    {
                        locFixesCount++;
                        locFirstFixType = replacementMer.fixType;
                    }
                    else
                    {
                        locMersToFirstFix += replacementMer.mersToNextFix + 1;
                        locFirstFixType = replacementMer.nextFixType;
                    }
                    locFollowerSum += replacementMer.sum;
                    locHealedMerCount = replacementMer.healedMerCount;
                    if (replacementMer.nextFixType == FixTypes.fixAbandon)
                        locFirstFixType = FixTypes.fixAbandon;

                    // no need to continue with this loop after trying a fix - recursion will have done it for us
                    break;
                }

                // not a good mer and no more fixes allowed
                locFirstFixType = FixTypes.fixAbandon;
                break;

            }

            if (replacementMer != null)
                ReturnMerProperty(replacementMer);

            // populate the 'downstream' field in the mer object passed in as a parameter
            variantMer.goodFollowers = locGoodFollowerCount;
            variantMer.allFollowers = locAllFollowerCount;
            variantMer.maxFollowers = locMaxFollowerCount;
            variantMer.healedMerCount = locHealedMerCount;
            variantMer.fixes = locFixesCount;
            variantMer.sum = locFollowerSum;
            variantMer.mersToNextFix = locMersToFirstFix;
            variantMer.nextFixType = locFirstFixType;
            variantMer.perfectFix = locPerfectFix;
            variantMer.cached = false;

            if (tracing >= traceRead)
                TraceFollowerResults(indentDepth, m, variantMer.fixType, false, countFollowersTimer, cachedFollowersFinal, cachedFollowersAbandoned);

            if (cachedFollowersFinal != null)
            {
                followerCache cachedFollowers = null;
                string followerCacheKey = null;
                //ulong followerCacheKey = 0;
                if (cacheFinal)
                {
                    cachedFollowers = cachedFollowersFinal;
                    followerCacheKey = followerCacheKeyFinal;
                }
                //else
                //if (cachedFollowersAbandoned != null)
                //{
                //    cachedFollowers = cachedFollowersAbandoned;
                //    followerCacheKey = followerCacheKeyAbandoned;
                //}

                if (cachedFollowers != null)
                {
                    int nextCacheIdx = cachedFollowers.nextCachedFollower;

                    if (nextCacheIdx == cachedFollowers.cachedFollowerMax)
                    {
                        if (nextCacheIdx == maxCacheAllowed)
                            throw (new ApplicationException("cache=" + nextCacheIdx));
                        int newCacheSize = maxCacheAllowed;

                        cachedFollowers.cachedFollowerMax = newCacheSize;
                        Array.Resize<int>(ref cachedFollowers.cachedCountGoodFollowers, newCacheSize);
                        Array.Resize<int>(ref cachedFollowers.cachedCountAllFollowers, newCacheSize);
                        Array.Resize<int>(ref cachedFollowers.cachedCountMaxFollowers, newCacheSize);
                        Array.Resize<int>(ref cachedFollowers.cachedHealedMerCount, newCacheSize);
                        Array.Resize<int>(ref cachedFollowers.cachedCountFixes, newCacheSize);
                        Array.Resize<int>(ref cachedFollowers.cachedFixedSum, newCacheSize);
                        Array.Resize<int>(ref cachedFollowers.cachedMersToFirstFix, newCacheSize);
                        Array.Resize<FixTypes>(ref cachedFollowers.cachedFirstFixType, newCacheSize);
                        Array.Resize<bool>(ref cachedFollowers.cachedPerfectFix, newCacheSize);
                        stats.cacheResizes++;
                    }

                    cachedFollowers.cachedKeys.Add(followerCacheKey, nextCacheIdx);
                    cachedFollowers.nextCachedFollower++;
                    cachedFollowers.cachedCountGoodFollowers[nextCacheIdx] = locGoodFollowerCount;
                    cachedFollowers.cachedCountAllFollowers[nextCacheIdx] = locAllFollowerCount;
                    cachedFollowers.cachedCountMaxFollowers[nextCacheIdx] = locMaxFollowerCount;
                    cachedFollowers.cachedHealedMerCount[nextCacheIdx] = locHealedMerCount;
                    cachedFollowers.cachedCountFixes[nextCacheIdx] = locFixesCount;
                    cachedFollowers.cachedFixedSum[nextCacheIdx] = locFollowerSum;
                    cachedFollowers.cachedMersToFirstFix[nextCacheIdx] = locMersToFirstFix;
                    cachedFollowers.cachedFirstFixType[nextCacheIdx] = locFirstFixType;
                    cachedFollowers.cachedPerfectFix[nextCacheIdx] = locPerfectFix;
                    if (tracing >= traceRead)
                    {
                        string whichCache = cacheFinal ? " (final) " : " (abandoned) ";
                        TraceClumpAdd(Indent(indentDepth) + "cached " + followerCacheKey + " g=" + locGoodFollowerCount + " a=" + locAllFollowerCount + " m=" + locMaxFollowerCount + whichCache + "(" + cachedFollowers.cachedKeys.Count + ") ");
                    }
                }
            }
        }

        private static void TraceFollowerResults(int indentDepth, int m, FixTypes variantType, bool cachedResult,
                                                 Stopwatch countFollowersTimer, followerCache followerCacheFinal, followerCache followerCacheAbandoned)
        {
            countFollowersTimer.Stop();
            float timeTakenToCount = (float)(countFollowersTimer.ElapsedMilliseconds) / 1000;

            string cacheSizes = "";
            if (followerCacheFinal != null)
                cacheSizes = followerCacheFinal.cachedKeys.Count.ToString(); ;
            if (followerCacheAbandoned != null)
                cacheSizes += "+" + followerCacheAbandoned.cachedKeys.Count;
            string cached = "";
            if (cachedResult)
                cached = " (cached) ";
            TraceClumpAdd(Indent(indentDepth) + "counted followers for " + fixNames[(int)variantType] + " @" + m + cached + " in " + timeTakenToCount.ToString("F3") + " cache=" + cacheSizes);
        }

        private static void UpdateHealingRead(FixTypes repType, Sequence healingRead, int m, ulong replacementMer, int readMerCount,
                                              ref int merCount, ref char repBase)
        {
            if (repType == FixTypes.fixDel)
            {
                // the replacement 'del' mer has had a base inserted somewhere. This means that we have to insert 
                // a placeholder base in the read string, prior to replacing the healed mer  

                healingRead.Insert(m + merSize - 1, 'X');
                merCount++;
                healingRead.Replace(m, replacementMer, merSize);
                repBase = healingRead.Bases[m + merSize - 1];
            }

            if (repType == FixTypes.fixSub)
            {
                // for substitutions, the structure of the healed read doesn't change, 
                // so just overwrite the part of the read corresponding to the mer being healed
                healingRead.Replace(m, replacementMer, merSize);
                repBase = healingRead.Bases[m + merSize - 1];
            }

            if (repType == FixTypes.fixIns)
            {
                // the replacement 'ins' mer has had a base deleted somewhere. This means that we have to delete 
                // a base from the read string, prior to replacing the healed mer (which includes the following base from the read)
                merCount--;
                healingRead.Remove(m);
                healingRead.Replace(m, replacementMer, merSize);
                repBase = healingRead.Bases[m + merSize - 1];
            }
        }

        private static int FindPlausibleVariants(FixTypes variantType, int variantsWanted, CheckResult checkReason, int minDepth, int OKDepth, ref bool subFixesOnly,
                                                 ulong startingMer, Sequence healingContext, int readMerCount, int m, int startingM, int merCount,
                                                 int repairsLeft, int consecutivePoor, followerCache cachedFollowersFinal, followerCache cachedFollowersAbandoned,
                                                 List<merProperties> variantMers, int indentDepth, ulong[] merVariants, ref int thmCalls)
        {
            int variantsAdded = 0;                  // count of variants added here (return value)

            if (m > 0)                              // only vary last base once past the first mer
                variantsWanted = varyLast;

            int noMerVariants = 0;
            //ulong[] merVariants = new ulong[merSize * 4];             // passed in to save on allocation overhead
            ulong[] generatedMerVariants = merVariants;                 // local pointer to mer variant array to allow use of locally allocated array as well
            string[] merVariantsTrace = null;

            // generate all the variants of the current mer (including current mer itself)
            if (variantType == FixTypes.fixSub)
            {
                noMerVariants = GenerateMerSubVariants(startingMer, merVariants, 0, variantsWanted);
                if (variantsWanted == varyTwo)
                {
                    // collected variants of variants (unsorted, with duplicates) 
                    ulong[] allDoubleMerVariants = new ulong[noMerVariants * merSize * 4];
                    int allIdx = 0;
                    int totalDoubleVariants = 0;

                    for (int v = 0; v < noMerVariants; v++)
                    {
                        // generate variants of next variant
                        int noNewSubVariants = GenerateMerSubVariants(merVariants[v], allDoubleMerVariants, allIdx, variantsWanted);
                        totalDoubleVariants += noNewSubVariants;
                        allIdx += noNewSubVariants;
                    }

                    noMerVariants = totalDoubleVariants;
                    generatedMerVariants = allDoubleMerVariants;
                }
            }

            if (variantType == FixTypes.fixDel)
            {
                noMerVariants = GenerateMerDelVariants(startingMer, merVariants, variantsWanted);
            }

            if (variantType == FixTypes.fixIns)
            {
                char nextBase = 'N';
                if (m + merSize < healingContext.Length)
                    nextBase = healingContext.Bases[m + merSize];  // if we're not right at the end of the read, find the base after this mer
                if (nextBase != 'N')
                    noMerVariants = GenerateMerInsVariants(startingMer, nextBase, merVariants, variantsWanted);
            }

            if (tracing >= traceChoices)
            {
                merVariantsTrace = new string[noMerVariants];
                for (int v = 0; v < noMerVariants; v++)
                    merVariantsTrace[v] = MerStrings.ExpandMer(generatedMerVariants[v]);
            }

            // look up the scores (and followers) of each variant and add them to the return list if they meet the requirements
            //Depth[] merVariantDepths = new Depth[Math.Max(merSize * 4, noMerVariants)];
            Depth[] merVariantDepths = GetDepths(noMerVariants);
            int sumVariantDepths = 0;

            // get the depths of each variant (and the sum of all depths)
            for (int v = 0; v < noMerVariants; v++)
            {
                int variantDepthPlus = 0;
                int variantDepthRC = 0;
                ulong variantMer = generatedMerVariants[v];
                FindMerInUniqueMers(variantMer, out variantDepthPlus, out variantDepthRC);
                int sumDepth = variantDepthPlus + variantDepthRC;
                merVariantDepths[v].depth = variantDepthPlus + variantDepthRC;
                if (sumDepth > 0)
                    merVariantDepths[v].unbalanced = (variantDepthPlus * 1000 / sumDepth > 950) | (variantDepthRC * 1000 / sumDepth > 950);
                else
                    merVariantDepths[v].unbalanced = false;
                sumVariantDepths += variantDepthPlus + variantDepthRC;
            }

            // condense the list of variants by copying any likely-looking alternatives towards the start of the merVariants array,
            // discarding 'poor' alternatives (and any ...XXXXXXXXXXXX (long HP) alternatives as these attractors cause performance problems and incorrect healing)
            // also discard the startingMer as this will always be in the list somewhere
            int viableMerVariants = 0;
            for (int v = 0; v < noMerVariants; v++)
            {
                ulong variantMer = generatedMerVariants[v];
                Depth variantDepth = merVariantDepths[v];
                bool lowRepVariant = false;
                if (sumVariantDepths > 0)
                    lowRepVariant = ((variantDepth.depth * 1000) / sumVariantDepths) <= 50;
                if (variantMer != startingMer && variantDepth.depth >= minDepth && !lowRepVariant && !MerIsLongHP(variantMer))
                {
                    generatedMerVariants[viableMerVariants] = variantMer;
                    merVariantDepths[viableMerVariants] = variantDepth;
                    viableMerVariants++;
                }
            }

            if ((sumVariantDepths >= averageDepth * highDepthFactor) && errorModel == modelMostlySubs)
                subFixesOnly = true;

            // multiple viable variants, so sort them so that the next loop can ignore duplicates (can't get duplicates if we're only varying the last base)
            if (viableMerVariants > 1 && variantsWanted != varyLast)
                Array.Sort<ulong, Depth>(generatedMerVariants, merVariantDepths, 0, viableMerVariants);

            ulong previousMer = 0xffffffffffffffff;
            for (int v = 0; v < viableMerVariants; v++)
            {
                ulong currentVariant = generatedMerVariants[v];
                int variantRepairsLeft = repairsLeft;
                int variantMerCount = merCount;

                // it's possible to have duplicates in this list, so just ignore them in this pass
                if (previousMer == currentVariant)
                    continue;
                previousMer = currentVariant;

                // do the next variants look likely (have suitable scores)
                Depth variantDepth = merVariantDepths[v];
                if (variantDepth.depth >= minDepth)
                {
                    Sequence variantReadContext = GetSequence();
                    healingContext.CopyTo(variantReadContext);
                    char variantBase = ' ';

                    UpdateHealingRead(variantType, variantReadContext, m, currentVariant, readMerCount, ref variantMerCount, ref variantBase);

                    // don't continue with a variant if we know it does not fit into a longer context, and only if the read is stable
                    bool forwardsOrBackwards = (errorModel == modelMostlySubs && variantType == FixTypes.fixSub) ? merPairEither : merPairBackOnly;
                    int pairPresent = pairNotPossible;
                    if (merPairs != null)
                        pairPresent = IsPairPresent(variantReadContext, m, minDepth, forwardsOrBackwards);
                    if (pairPresent == pairNotFound)
                        continue;

                    merProperties variantMer = GetMerProperty();
                    //merProperties variantMer = new merProperties();
                    variantMer.fixType = variantType;
                    variantMer.variant = currentVariant;
                    variantMer.depth = variantDepth.depth;
                    variantMer.unbalanced = variantDepth.unbalanced;
                    variantMer.pairWasFound = (pairPresent == pairFound);

                    //if (tracing == traceFollowers)
                    //    TraceClumpAdd(Indent(indentDepth) +
                    //                     "var " + v + " " + fixTypes[variantType] + " " + tempHealingRead + " @" + m);

                    if (tracing > traceChoices)
                        TraceFollowers(checkReason, variantType, m, indentDepth, MerStrings.ExpandMer(currentVariant), v, variantDepth.depth, variantRepairsLeft, consecutivePoor, variantReadContext.ToString(), variantMerCount);

                    CountFollowers(variantMer, variantReadContext, readMerCount, m, startingM, variantMerCount, variantRepairsLeft, consecutivePoor,
                                   minDepth, OKDepth, ref subFixesOnly, cachedFollowersFinal, cachedFollowersAbandoned, indentDepth, ref thmCalls);

                    variantMer.sum += variantDepth.depth;
                    variantMers.Add(variantMer);
                    variantsAdded++;

                    if (tracing == traceFollowers)
                    {
                        string cachedTag = variantMer.cached ? " (cached)" : "";
                        TraceClumpAdd(Indent(indentDepth) +
                                          "var:" + fixNames[(int)variantType] + "@" + m + " " + MerStrings.ExpandMer(currentVariant) +
                                          " sum=" + variantMer.sum + " fx " + variantMer.fixes +
                                          " fo=" + variantMer.goodFollowers + "/" + variantMer.allFollowers + "/" + variantMer.maxFollowers + "/" + variantMerCount + " gap " + variantMer.mersToNextFix + cachedTag);
                    }

                    ReturnSequence(variantReadContext);
                }
                else
                {
                    //if (tracing == traceFollowers)
                    //    TraceClumpAdd(Indent(indentDepth) + "skipped variant " + fixTypes[variantType] + " " + v + "  " + MerStrings.ExpandMer(merVariant));
                }

            }

            ReturnDepths(merVariantDepths);

            return variantsAdded;
        }

        private static void ExtendShortFixedRead(Sequence read, int wantedLength, Sequence quals)
        {
            int readLength = read.Length;

            // try adding non-ambiguous bases to required length
            while (readLength < wantedLength)
            {
                bool extended = false;
                extended = ExtendReadByOneBase(read, quals);
                if (!extended)
                    break;
                readLength++;
            }

            // if the read is still short, just pad read out with Ns
            while (readLength < wantedLength)
            {
                read.Append('N');
                if (quals.Length != 0)
                    quals.Append((char)1);
                readLength++;
            }
        }

        private static bool ExtendReadByOneBase(Sequence read, Sequence quals)
        {
            ulong packedLastMer;

            // don't extend if the trailing k-mer contains N
            if (!MerStrings.CondenseMer(read, (read.Length - merSize), merSize, out packedLastMer))
                return false;

            // pack the mer at the end of the read
            int depthPlus;
            int depthRC;
            FindMerInUniqueMers(packedLastMer, out depthPlus, out depthRC);

            // only extend from a good mer
            if ((depthPlus + depthRC) < minReps)
                return false;

            // generate a packed right-shifted k-mer with a hole at the end for a base to be inserted
            int merFill = 64 - merSize * 2;
            ulong rshiftedMer = (packedLastMer >> merFill) << 2;
            int depthSum = 0;
            int totalDepth = 0;
            int deepestDepth = 0;
            int deepestBase = -1;       // (A=0, ...)

            // find the total depth and the deepest individual depth
            for (ulong b = 0; b < 4; b++)
            {
                ulong possibleMer = (rshiftedMer | b) << merFill;
                FindMerInUniqueMers(possibleMer, out depthPlus, out depthRC);
                depthSum = depthPlus + depthRC;
                totalDepth += depthSum;
                if (depthSum > deepestDepth)
                {
                    deepestDepth = depthSum;
                    deepestBase = (int)b;
                }
            }

            // and choose the highest base if it looks to be the only sensible candidate (if it even exists)
            if (deepestBase < 0)
                return false;                       // none of the alternatives were viable (will add an N in the caller)

            if ((deepestDepth * 100 / totalDepth) < 90)
                return false;                       // ambiguity - want the candidate to contribute 90% to the total

            // found a likely extending base so concatenate it to the starting read
            char appendingBase = "ACGT"[deepestBase];
            read.Append(appendingBase);
            if (quals.Length != 0)
                quals.Append((char)replacementQual);

            return true;
        }

        private static void TraceFollowers(CheckResult checkReason, FixTypes variantType, int m, int indentDepth, string merVariant, int v, int depth, int repairsLeft,
                                           int consecutivePoor, String readContext, int merCount)
        {
            int startOfRest = m + merSize;
            int contextLength = readContext.Length - startOfRest;
            if (contextLength > merSize)
                contextLength = merSize;
            string context = "(end)";
            if (contextLength > 0)
                context = readContext.Substring(startOfRest, contextLength);

            TraceClumpAdd(Indent(indentDepth) + "@" + m + " counting followers for " +
                            "[" + v + "] " + fixNames[(int)variantType] + " " + merVariant + "+" + context +
                            " with " + repairsLeft + " fixes left; sum=" + depth + "; mc=" + merCount + "; cp=" + consecutivePoor +
                            " " + checkNames[(int)checkReason]);
        }

        // Generates a 'following' mer by taking the base following the end of a mer variant and 
        // concatenating it to the end of the variant. Returns true if such a mer could be generated.
        // Only used when counting followers. 
        private static bool GenerateNextMer(ulong packedMer, Sequence healingRead, int m, int startingM,
                                            int merCount, int readMerCount, int OKDepth, out ulong nextPackedMer, out char followingBase)
        {
            //if ((m >= merCount) || (m > startingM + merSize))
            if (m >= Math.Min(merCount, readMerCount))      // go as far as we can to maximise followers 
            {
                nextPackedMer = 0;
                followingBase = '\0';
                return false;
            }
            else
            {
                char nextBase = healingRead.Bases[m + merSize - 1];
                if (m + merSize < healingRead.Length)
                    followingBase = healingRead.Bases[m + merSize];
                else
                    followingBase = '\0';   // avoid triggering end-of-homopolymer at end of read
                long nextBasePacked = MerStrings.BaseCharToInt(nextBase);
                //int nextBasePacked = "ACGT".IndexOf(nextBase);
                if (nextBasePacked < 0)
                {
                    return FindBestReplacementForNs(healingRead, m, OKDepth, out nextPackedMer);
                }
                else
                    nextPackedMer = packedMer << 2 | (ulong)nextBasePacked << (64 - merSize * 2);
                //string merIn = MerStrings.ExpandMer(packedMer);
                //string merOut = MerStrings.ExpandMer(nextPackedMer);
                return true;
            }
        }

        private static int GenerateMerSubVariants(ulong mer, ulong[] merVariants, int startingIdx, int variantsWanted)
        {
            int start = 0;
            if (variantsWanted == varyLast)
                start = merSize - 1;
            int variantsAdded = 0;

            ulong baseMask = 0xc000000000000000;

            for (int m = start; m < merSize; m++)
            {
                ulong merWithHole = mer & ~(baseMask >> (m * 2));
                for (ulong b = 0; b <= 3; b++)
                {
                    ulong newBase = b << (64 - (m + 1) * 2);
                    ulong merVariant = merWithHole | newBase;
                    merVariants[startingIdx + variantsAdded] = merVariant;
                    variantsAdded++;
                }
            }
            return variantsAdded;
        }

        private static int GenerateMerDelVariants(ulong mer, ulong[] merVariants, int variantsWanted)
        {
            int start = 1;                                          // never try inserting at mer[0] as this is effectively shifting the read 1 base to the left
            if (variantsWanted == varyLast)
                start = merSize - 1;
            int variantsAdded = 0;

            ulong maskMer = 0xffffffffffffffff << (64 - merSize * 2);   // just the RHS merSize bits

            for (int m = start; m < merSize; m++)
            {
                ulong maskLHS = 0xffffffffffffffff << (64 - (m * 2));   // the retained LHS part of the mer
                ulong merLHS = mer & maskLHS;
                ulong maskRHS = ~maskLHS;                               // the retained RHS part of the mer
                ulong merRHS = ((mer & maskRHS) >> 2) & maskMer;
                for (ulong b = 0; b <= 3; b++)
                {
                    ulong newBase = b << (64 - (m + 1) * 2);
                    ulong merVariant = merLHS | newBase | merRHS;
                    if (variantsWanted == varySingle)                           // being careful about dels that effectively shift the read left one base in the genome
                    {
                        ulong topBase = merVariant & 0xc000000000000000;        // just the first base of the variant
                        ulong shiftedMer = ((mer >> 2) | topBase) & maskMer;    // starting mer with first base from variant inserted at start

                        if (merVariant == shiftedMer)                           // this variant didn't really change anything
                            continue;
                    }
                    merVariants[variantsAdded] = merVariant;
                    variantsAdded++;
                }
            }
            return variantsAdded;
        }

        private static int GenerateMerInsVariants(ulong mer, char followingBase, ulong[] merVariants, int variantsWanted)
        {
            int start = 1;                                          // never try deleting the base at mer[0] as this is effectively shifting the read 1 base to the right
            if (variantsWanted == varyLast)
                start = merSize - 1;
            int variantsAdded = 0;
            ulong newBase = (ulong)MerStrings.BaseCharToInt(followingBase) << (64 - merSize * 2);    

            for (int m = start; m < merSize; m++)
            {
                ulong maskLHS = 0xffffffffffffffff << (64 - (m * 2));   // the retained LHS part of the mer
                ulong merLHS = mer & maskLHS;
                ulong maskRHS = ~maskLHS;                               // the RHS part of the mer (including the base to be excised)
                ulong merRHS = ((mer << 2) & maskRHS);                  // and then without the right-most base
                ulong merVariant = merLHS | newBase | merRHS;
                if (merVariant == mer)                                  // don't include starting mer as one of its own variants
                    continue;
                merVariants[variantsAdded] = merVariant;
                variantsAdded++;
            }
            return variantsAdded;
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

                    if (inLoadingPhase)
                        Console.WriteLine("loaded: " + MerStrings.progressUniqueMers + "(" + MerStrings.progressLoadedUniqueMers + ")/" +
                                           MerStrings.progressTotalMers + "(" + MerStrings.progressLoadedTotalMers + ") mers");
                    if (inHealingPhase)
                        Console.WriteLine("healed " + progressHealedReads + " from " + currentReadsCount + " reads at " + readsRate + " rps");
                }
            }
        }
    }

    public class healingThreadParams
    {
        public int threadNumber;
        public BufferedReader[] bufferedReadsFiles;
        public StreamWriter[] healedReads;
        public StreamWriter[] healedQuals;
        public StreamWriter[] problemReads;
        public StreamWriter[] problemQuals;
        public StreamWriter[] singleReads;
        public StreamWriter[] singleQuals;
    }

    public class followerCache
    {
        public int nextCachedFollower = 0;
        public int cachedFollowerMax = 1000;
        public Dictionary<string, int> cachedKeys;                   // cached result key --> array index
        //public Dictionary<ulong, int> cachedKeys;                   // cached result key --> array index
        public FixTypes[] cachedFirstFixType;
        public int[] cachedMersToFirstFix;                          // should rewrite using struct
        public int[] cachedCountGoodFollowers;
        public int[] cachedCountAllFollowers;
        public int[] cachedCountMaxFollowers;
        public int[] cachedHealedMerCount;
        public int[] cachedCountFixes;
        public int[] cachedFixedSum;
        public bool[] cachedPerfectFix;

        public followerCache()
        {
            cachedKeys = new Dictionary<string, int>(cachedFollowerMax);
            //cachedKeys = new Dictionary<ulong, int>(cachedFollowerMax);
            cachedFirstFixType = new FixTypes[cachedFollowerMax];
            cachedMersToFirstFix = new int[cachedFollowerMax];
            cachedCountGoodFollowers = new int[cachedFollowerMax];
            cachedCountAllFollowers = new int[cachedFollowerMax];
            cachedCountMaxFollowers = new int[cachedFollowerMax];
            cachedHealedMerCount = new int[cachedFollowerMax];
            cachedCountFixes = new int[cachedFollowerMax];
            cachedFixedSum = new int[cachedFollowerMax];
            cachedPerfectFix = new bool[cachedFollowerMax];
        }
    }

    public class merProperties
    {
        public bool validVariant = true;            // this variant can be considered (or has been deleted)
        public bool markedVariant = false;          // this variant marked for possible deletion
        public FixTypes fixType = FixTypes.fixNone;
        public FixTypes nextFixType = FixTypes.fixNone;
        public int mersToNextFix = 0;
        public bool perfectFix = false;
        public ulong variant = 0;
        public int depth = 0;
        public bool unbalanced = false;
        public bool pairWasFound = false;
        public int sum = 0;
        public int goodFollowers = 0;
        public int allFollowers = 0;
        public int maxFollowers = 0;
        public int healedMerCount = 0;
        public int fixes = 0;
        public bool cached = false;
        public int savedGoodFollowers = 0;
        public int savedAllFollowers = 0;
        public int savedMaxFollowers = 0;
    }

    public struct Depth
    {
        public int depth;
        public bool unbalanced;
    }

    public class healingStats
    {
        const int maxFixTypes = 6;                                   // kludge - defined twice
        public long reads = 0;                                       // how many reads have been processed so far
        public long shortReads = 0;                                  // how many short reads were passed over
        public long longReads = 0;                                   // how many overly-long reads were passed over
        public long OKReads = 0;                                     // how many reads did not need healing (on initial look)
        public long OKReadsChecked = 0;                              // and how many reads looked like they might have a problem but proved to be OK
        public long healedReads = 0;                                 // how many reads were healed
        public long readsNotHealedFully = 0;                         // how many reads needed healing but were not completely healed
        public long readsNotHealedAtAll = 0;                         // how many reads needed healing but could not be healed at all 
        public long singletonReads = 0;                              // how many reads with no shared mers
        public long unbalancedReads = 0;                             // how many reads with one strand very high and the other effectively zero
        public long lowComplexityReads = 0;                          // how many reads deemed to be low complexity (and not corrected)
        public long tooHighDepthReads = 0;                           // how many reads were not corrected because they seemed to be too deep
        public long abandonedReads = 0;                              // how many reads were abandoned (poor quality, rewriting, explosion or too many Ns)
        public long abandonedRewriting = 0;                          // ... decided we were just rewriting the read
        public long abandonedTree = 0;                               // ... tree got too big (cache or thm calls)
        public long abandonedNs = 0;                                 // ... read had too many Ns
        public long cacheResizes = 0;                                // how many times we resized the cache arrays
        public float abandonedTime = 0.0f;                           // how much time was spent on reads that ended up being abandoned (cache or thm calls)
        public long RCHealedReads = 0;                               // reads that were healed by trying again in the reverse direction
        public long fwdHealedReads = 0;                              // reads that were changed by trying the third (forward) pass - may be none...
        public long mers = 0;                                        // how many mers were in these reads
        public long replacedMers = 0;                                // how many of these mers were healed
        public long slowHealings = 0;                                // how many slow healing reads were encountered
        public long discarded = 0;                                   // how many 'poor' reads were discarded
        public long extended = 0;                                    // how many shortened Illumina reads were restored their fixed length
        public long[] fixesByType = new long[maxFixTypes];           // recording types of fixes used
    }

}


