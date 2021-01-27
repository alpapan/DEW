using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading;
using System.Diagnostics;

namespace CommonCode
{
    public static class MerStrings
    {
        static int merSize = 0;
        public const int formatNone = 0;            // no format specified
        public const int formatFNA = 1;             // fa with multiple data lines (.fna)
        public const int formatFASTQ = 2;           // fastq (with multiple data lines per read)

        public const int tfaType = 0;
        public const int bfaType = 1;
        public const int bfqType = 2;

        public static long progressUniqueMers = 0;
        public static long progressTotalMers = 0;
        public static long progressLoadedUniqueMers = 0;
        public static long progressLoadedTotalMers = 0;

        static char[] spaceDelimiter = new char[] { ' ' };

        const int bytesPerCBTMer = 8 * 2;

        // Mapping low 5 bits of ASCII chars to 2-bit base encoding (A=00, C=01, G=10, T=11). Same mapping works for lower case bases.
        // Used by BaseCharToInt(char baseChar) to convert char bases to 2-bit longs
        private static long[] baseToInt = new long[] { -1, 0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
        //                                              @, A,  B, C,  D,  E,  F, G,  H,  I,  J,  K,  L,  M,  N,  O,  P,  Q,  R,  S, T,  U,  V,  W,  X,  Y,  Z,  [, \\,  ],  ^,  _
        // mapping 2-bit bases back to char form
        public static char[] baseToChar = new char[] { 'A', 'C', 'G', 'T' };

        public static void Initialise(int merSizeP)
        {
            merSize = merSizeP;
        }

        public static string ExpandMer(ulong packedMer)
        {
            if (merSize <= 0)
                throw new ArgumentException("MerStrings: merSize not set");

            return ExpandMer(packedMer, merSize);
        }

        public static string ExpandMer(ulong packedMer, int merSize)
        {
            char[] charMer = new char[merSize];             // char mer string under construction
            ulong tempMer = packedMer;                      // copy of packedMer for shifting
            ulong intBase = 0;                              // current base in binary form

            for (int i = 0; i < merSize; i++)
            {
                intBase = (ulong)tempMer >> 62;
                tempMer = (ulong)(tempMer << 2);

                charMer[i] = baseToChar[intBase];
            }

            return new string(charMer, 0, merSize);
        }

        // delete after MapDraftGenomeOntoDraftGenome changed to use Sequence
        public static void ReplaceMer(StringBuilder read, int m, ulong packedMer, int merSize)
        {
            ulong intBase = 0;                              // current base in binary form
            ulong tempMer = packedMer;                      // progressively shifted left
            char ACGT;                                      // char form of base

            for (int i = 0; i < merSize; i++)
            {
                intBase = (ulong)tempMer >> 62;
                tempMer = (ulong)(tempMer << 2);

                ACGT = baseToChar[intBase];
                read[m + i] = ACGT;
            }
        }

        public static ulong CondenseMer(string mer)
        {
            ulong packedMer;
            CondenseMer(mer, out packedMer);
            return packedMer;
        }

        // Convert a text-form k-mer into its 2-bit binary form. 
        // The function will return false if any non-ACGT bases are found. Lower case acgt are also accepted.
        // 
        public static bool CondenseMer(string seq, out ulong packedMer)
        {
            packedMer = 0;

            if (seq.Length > 32)
                return false;

            for (int m = 0; m < seq.Length; m++)
            {
                char nextBase = seq[m];
                long packedBase = BaseCharToInt(nextBase);

                if (packedBase < 0)
                {
                    packedMer = 0;
                    return false;
                }

                packedMer = (packedMer << 2) | (ulong)packedBase;
            }
            packedMer = packedMer << (64 - seq.Length * 2);
            return true;
        }

        public static bool CondenseMer(string seq, int start, int merSize, out ulong packedMer)
        {
            packedMer = 0;

            for (int m = 0; m < merSize; m++)
            {
                char nextBase = seq[m + start];
                long packedBase = BaseCharToInt(nextBase);

                if (packedBase < 0)
                {
                    packedMer = 0;
                    return false;
                }

                packedMer = (packedMer << 2) | (ulong)packedBase;
            }
            packedMer = packedMer << (64 - merSize * 2);
            return true;
        }

        public static bool CondenseMer(StringBuilder seq, int start, int merSize, out ulong packedMer)
        {
            packedMer = 0;

            for (int m = 0; m < merSize; m++)
            {
                char nextBase = seq[m + start];
                long packedBase = BaseCharToInt(nextBase);

                if (packedBase < 0)
                {
                    packedMer = 0;
                    return false;
                }

                packedMer = (packedMer << 2) | (ulong)packedBase;
            }
            packedMer = packedMer << (64 - merSize * 2);
            return true;
        }

        public static bool CondenseMer(Sequence seq, int start, int merSize, out ulong packedMer)
        {
            packedMer = 0;

            for (int m = 0; m < merSize; m++)
            {
                char nextBase = seq.Bases[m + start];
                long packedBase = BaseCharToInt(nextBase);

                if (packedBase < 0)
                {
                    packedMer = 0;
                    return false;
                }

                packedMer = (packedMer << 2) | (ulong)packedBase;
            }
            packedMer = packedMer << (64 - merSize * 2);
            return true;
        }

        public static bool CondenseMerIncremental(int merSize, ulong previousMer, string line, int m, out ulong nextMer)
        {
            char nextBase = line[m + merSize - 1];
            long packedBase = BaseCharToInt(nextBase);
            nextMer = (previousMer << 2) | (ulong)packedBase << (64 - merSize * 2);
            if (packedBase < 0)
            {
                nextMer = 0;
                return false;
            }
            return true;
        }

        public static bool CondenseMerIncremental(int merSize, ulong previousMer, Sequence read, int m, out ulong nextMer)
        {
            char nextBase = read.Bases[m + merSize - 1];
            long packedBase = BaseCharToInt(nextBase);
            nextMer = (previousMer << 2) | (ulong)packedBase << (64 - merSize * 2);
            if (packedBase < 0)
            {
                nextMer = 0;
                return false;
            }
            return true;
        }

        public static long BaseCharToInt(char baseChar)
        {
            return baseToInt[baseChar & 0x1f];
        }

        public static string ReverseComplement(string s)
        {
            int sl = s.Length;
            char[] rcs = new char[sl];

            for (int si = 0; si < sl; si++)
            {
                long b = BaseCharToInt(s[si]);
                switch (b)
                {
                    case 0:
                        rcs[sl - si - 1] = 'T';
                        break;
                    case 1:
                        rcs[sl - si - 1] = 'G';
                        break;
                    case 2:
                        rcs[sl - si - 1] = 'C';
                        break;
                    case 3:
                        rcs[sl - si - 1] = 'A';
                        break;
                    case -1:
                        rcs[sl - si - 1] = 'N';
                        break;
                }
            }
            return new string(rcs);
        }

        public static StringBuilder ReverseComplement(StringBuilder s)
        {
            int sl = s.Length;
            StringBuilder rcs = new StringBuilder(sl);
            rcs.Length = sl;

            for (int si = 0; si < sl; si++)
            {
                long b = BaseCharToInt(s[si]);
                switch (b)
                {
                    case 0:
                        rcs[sl - si - 1] = 'T';
                        break;
                    case 1:
                        rcs[sl - si - 1] = 'G';
                        break;
                    case 2:
                        rcs[sl - si - 1] = 'C';
                        break;
                    case 3:
                        rcs[sl - si - 1] = 'A';
                        break;
                    case -1:
                        rcs[sl - si - 1] = 'N';
                        break;
                }
            }
            return rcs;
        }

        public static void ReverseComplement(Sequence s)
        {
            int sl = s.Length;
            char[] rcs = new char[sl];

            for (int si = 0; si < sl; si++)
            {
                long b = BaseCharToInt(s.Bases[si]);
                switch (b)
                {
                    case 0:
                        rcs[sl - si - 1] = 'T';
                        break;
                    case 1:
                        rcs[sl - si - 1] = 'G';
                        break;
                    case 2:
                        rcs[sl - si - 1] = 'C';
                        break;
                    case 3:
                        rcs[sl - si - 1] = 'A';
                        break;
                    case -1:
                        rcs[sl - si - 1] = 'N';
                        break;
                }
            }
            rcs.CopyTo(s.Bases, 0);
        }

        public static ulong ReverseComplementPacked(ulong mer)
        {
            if (merSize <= 0)
                throw new ArgumentException("MerStrings: merSize not set");

            return ReverseComplementPacked(mer, merSize);
        }

        public static ulong ReverseComplement(ulong mer)
        {
            if (merSize <= 0)
                throw new ArgumentException("MerStrings: merSize not set");

            return ReverseComplementPacked(mer, merSize);
        }

        // fast reverse complement for a packed k-mer (delete after changing all callers)
        public static ulong ReverseComplementPacked(ulong mer, int merSize)
        {
            mer = ((mer >> 2) & 0x3333333333333333) | ((mer & 0x3333333333333333) << 2);
            mer = ((mer >> 4) & 0x0F0F0F0F0F0F0F0F) | ((mer & 0x0F0F0F0F0F0F0F0F) << 4);
            mer = ((mer >> 8) & 0x00FF00FF00FF00FF) | ((mer & 0x00FF00FF00FF00FF) << 8);
            mer = ((mer >> 16) & 0x0000FFFF0000FFFF) | ((mer & 0x0000FFFF0000FFFF) << 16);
            mer = (mer >> 32) | (mer << 32);  // reversed and right-shifted after this line. Need to complement then left-shift.

            if (merSize < 32)
                return ~mer << (64 - merSize * 2);
            else
                return ~mer;
        }

        // fast reverse complement for a packed k-mer 
        public static ulong ReverseComplement(ulong mer, int merSize)
        {
            mer = ((mer >> 2) & 0x3333333333333333) | ((mer & 0x3333333333333333) << 2);
            mer = ((mer >> 4) & 0x0F0F0F0F0F0F0F0F) | ((mer & 0x0F0F0F0F0F0F0F0F) << 4);
            mer = ((mer >> 8) & 0x00FF00FF00FF00FF) | ((mer & 0x00FF00FF00FF00FF) << 8);
            mer = ((mer >> 16) & 0x0000FFFF0000FFFF) | ((mer & 0x0000FFFF0000FFFF) << 16);
            mer = (mer >> 32) | (mer << 32);  // reversed and right-shifted after this line. Need to complement then left-shift.

            if (merSize < 32)
                return ~mer << (64 - merSize * 2);
            else
                return ~mer;
        }

        public static int GenerateMerSubVariants(ulong mer, List<ulong> merVariants, int merSize)
        {
            int start = 0;
            int variantsAdded = 0;

            ulong baseMask = 0xc000000000000000;

            for (int m = start; m < merSize; m++)
            {
                ulong merWithHole = mer & ~(baseMask >> (m * 2));
                for (ulong b = 0; b <= 3; b++)
                {
                    ulong newBase = b << (64 - (m + 1) * 2);
                    ulong merVariant = merWithHole | newBase;
                    if (merVariant == mer)
                        continue;
                    merVariants.Add(merVariant);
                    variantsAdded++;
                }
            }
            return variantsAdded;
        }

        // The qual scores used with Fastq reads are ambiguous. Illumina data may use a base of 64 ('@') while Sanger formatted data
        // uses a base of 33 ('!'). This methods reads quals from a reads file until it finds a value that can only be either Sanger or Illumina
        // and then sets the static qual parameters appropriately.
        public static int ResolveFastqQualAmbiguity(string fastqReadsFN, out bool fullQualHeader)
        {
            fullQualHeader = true;
            StreamReader fastqReads = new StreamReader(fastqReadsFN);
            int qualBase = 64;                          // assume Illumina format until we find a Sanger-only value - this forces Illumina to be the default
            const char minIlluminaQual = '@';           // lowest possible Illumina qual
            const int readsToExamine = 100;             // only look at the first 100 reads
            int readsExamined = 0;
            bool foundIndicativeBase = false;           // found a base that can only appear in Sanger

            while (readsExamined < readsToExamine)
            {
                string readHeader = fastqReads.ReadLine();
                if (readHeader == null)
                {
                    Console.WriteLine("failed to resolve fastq qual ambiguity before reaching EOF on " + fastqReadsFN);
                    break;
                }
                readsExamined++;
                string read = fastqReads.ReadLine();
                string qualHeader = fastqReads.ReadLine();
                if (qualHeader == "+")
                    fullQualHeader = false;
                string quals = fastqReads.ReadLine();

                for (int i = 0; i < quals.Length; i++)
                {
                    if (quals[i] < minIlluminaQual)         // lower than any possible Illumina qual --> Sanger formatted data
                    {
                        qualBase = 33;
                        // Console.WriteLine("Sanger");
                        break;
                    }
                }

                if (foundIndicativeBase)
                    break;
            }
            fastqReads.Close();
            Console.WriteLine("FASTQ file found to be have " + qualBase + " base");
            return qualBase;
        }

        public static string LFConvention(string readsFN)
        {
            StreamReader reads = new StreamReader(readsFN);
            bool foundCR = false;
            bool foundLF = false;

            while (!foundLF)
            {
                int c = reads.Read();
                if (c == '\r')
                    foundCR = true;
                if (c == '\n')
                    break;
            }

            reads.Close();
            return foundCR ? "\r\n" : "\n";
        }

        public static int DetermineFileFormat(string fn)
        {
            int fileFormat = formatNone;
            StreamReader file = new StreamReader(fn);
            string firstLine = file.ReadLine();
            if (firstLine[0] == '@')
                fileFormat = formatFASTQ;
            if (firstLine[0] == '>')
                fileFormat = formatFNA;
            return fileFormat;
        }

        public static int GenerateMersFromRead(string line, int merSize, string[] merSet)
        {
            int lineLength = line.Length;
            int mersInRead = lineLength - merSize + 1;
            int mersReturned = 0;

            for (int i = 0; i < mersInRead; i++)
            {
                merSet[i] = line.Substring(i, merSize);
                mersReturned++;
            }

            return mersReturned;
        }

        public static int GenerateMersFromReadPacked(string line, int merSize, ref ulong[] merSet, ref bool[] merValid)
        {
            int lineLength = line.Length;
            int mersInRead = lineLength - merSize + 1;
            bool merIsValid;

            if (mersInRead < 1)
                return 0;

            if (merSet.Length < mersInRead)
            {
                Array.Resize<ulong>(ref merSet, mersInRead + 100);
                Array.Resize<bool>(ref merValid, mersInRead + 100);
            }

            // try the incremental (and faster) approach first

            // condense the first k-mer 
            int m = 0;
            merIsValid = CondenseMer(line.Substring(0, merSize), out merSet[0]);
            m++;
            // condense the rest of the mers incrementally if the first one was OK
            if (merIsValid)
            {
                merValid[0] = true;
                for (; m < mersInRead; m++)
                {
                    merIsValid = MerStrings.CondenseMerIncremental(merSize, merSet[m - 1], line, m, out merSet[m]);
                    if (merIsValid)
                        merValid[m] = true;
                    else
                        break;                  // found the first N or non-ACGT so restart the tiling (below)
                }
            }

            // if the line contains any non-ACGT bases, scan and generate the mers non-incrementally
            if (!merIsValid)
            {
                for (int i = 0; i < mersInRead; i++)
                {
                    string nextMer = line.Substring(i, merSize);
                    merValid[i] = MerStrings.CondenseMer(nextMer, out merSet[i]);
                }
            }

            return mersInRead;
        }

        public static int GenerateMersFromRead(string read, int merSize, ref ulong[] merSet, ref bool[] merValid)
        {
            int readLength = read.Length;
            int mersInRead = readLength - merSize + 1;
            bool merIsValid = false;
            ulong lastMer = 0;

            if (mersInRead < 1)
                return 0;

            if (merSet.Length < mersInRead)
            {
                Array.Resize<ulong>(ref merSet, mersInRead + 100);
                Array.Resize<bool>(ref merValid, mersInRead + 100);
            }

            for (int i = 0; i < mersInRead; i++)
            {
                if (merIsValid)
                {
                    merIsValid = CondenseMerIncremental(merSize, lastMer, read, i, out lastMer);
                    merValid[i] = merIsValid;
                    merSet[i] = lastMer;
                }
                else
                {
                    merIsValid = CondenseMer(read.Substring(i, merSize), out lastMer);
                    merValid[i] = merIsValid;
                    merSet[i] = lastMer;
                }
            }

            return mersInRead;
        }

        public static int GenerateMersFromRead(Sequence read, int merSize, ref ulong[] merSet, ref bool[] merValid)
        {
            int readLength = read.Length;
            int mersInRead = readLength - merSize + 1;
            bool merIsValid = false;
            ulong lastMer = 0;

            if (mersInRead < 1)
                return 0;

            if (merSet.Length < mersInRead)
            {
                Array.Resize<ulong>(ref merSet, mersInRead + 100);
                Array.Resize<bool>(ref merValid, mersInRead + 100);
            }

            for (int i = 0; i < mersInRead; i++)
            {
                if (merIsValid)
                {
                    merIsValid = CondenseMerIncremental(merSize, lastMer, read, i, out lastMer);
                    merValid[i] = merIsValid;
                    merSet[i] = lastMer;
                }
                else
                {
                    merIsValid = CondenseMer(read, i, merSize, out lastMer);
                    merValid[i] = merIsValid;
                    merSet[i] = lastMer;
                }
            }

            return mersInRead;
        }
 
        public static string ReadRead(StreamReader reads, int readFormat)
        {
            string read = null;
            string readHeader;

            if (readFormat == formatFASTQ)
            {
                readHeader = reads.ReadLine();
                if (readHeader == null)
                    return null;
                read = reads.ReadLine();
                reads.ReadLine();
                reads.ReadLine();
            }

            if (readFormat == formatFNA)
            {
                read = ReadFASTA(reads, false, out readHeader);
            }

            return read;
        }

        public static string ReadRead(StreamReader reads, int readFormat, out string header)
        {
            string read = null;
            header = "";

            if (readFormat == formatFASTQ)
            {
                header = reads.ReadLine();
                if (header == null)
                    return null;
                read = reads.ReadLine();
                reads.ReadLine();
                reads.ReadLine();
            }
            if (readFormat == formatFNA)
            {
                read = ReadFASTA(reads, false, out header);
            }

            return read;
        }

        public static string ReadRead(StreamReader reads, int readFormat, out string readHeader, out string qualHeader, out string quals)
        {
            string read = null;
            readHeader = null;
            qualHeader = null;
            quals = null;

            if (readFormat == formatFASTQ)
            {
                readHeader = reads.ReadLine();
                if (readHeader == null)
                    return null;
                read = reads.ReadLine();
                qualHeader = reads.ReadLine();
                quals = reads.ReadLine();
            }

            if (readFormat == formatFNA)
            {
                read = ReadFASTA(reads, false, out readHeader);
            }

            return read;
        }

        public static string ReadRead(StreamReader readsFile, StreamReader qualsFile, int readFormat, out string readHeader, out string qualHeader, out string quals)
        {
            string read = null;
            readHeader = null;
            qualHeader = null;
            quals = null;

            if (readFormat == formatFASTQ)
            {
                readHeader = readsFile.ReadLine();
                if (readHeader == null)
                    return null;
                read = readsFile.ReadLine();
                qualHeader = readsFile.ReadLine();
                quals = readsFile.ReadLine();
            }

            if (readFormat == formatFNA)
            {
                read = ReadFASTA(readsFile, false, out readHeader);
                if (qualsFile != null)
                {
                    quals = ReadFASTA(qualsFile, true, out qualHeader);
                    if (quals != null)
                    {
                        Sequence qualsSeq = new Sequence(quals);
                        ConvertQualsToCanonicalForm(qualsSeq, readFormat, 0);
                        quals = qualsSeq.ToString();
                    }
                }
            }

            return read;
        }

        public static bool ReadRead(StreamReader reads, StreamReader quals, int readFormat, Sequence readHeader, Sequence read, Sequence qualHeader, Sequence qual)
        {
            string readString = null;
            string readHeaderString = null;
            string qualHeaderString = null;        
            string qualString = null;

            readString = ReadRead(reads, quals, readFormat, out readHeaderString, out qualHeaderString, out qualString);

            if (readString != null)
            {
                readHeader.CopyFrom(readString);
                read.CopyFrom(readHeaderString);
                qualHeader.CopyFrom(qualHeaderString);
                qual.CopyFrom(qualString);
            }
            else
            {
                readHeader.Length = 0;
                read.Length = 0;
                qualHeader.Length = 0;
                qual.Length = 0;
            }

            return read.Length > 0;
        }

        public static int ReadReads(int readsWanted, StreamReader readsFile, StreamReader qualsFile, int readsFormat, Sequence[] readHeaders, Sequence[] readSeqs, Sequence[] qualHeaders, Sequence[] quals)
        {
            int readsRead = 0;

            for (int i = 0; i < readsWanted; i++)
            {
                string readHeader;
                string readSeq;
                string qualHeader;
                string qual;

                readSeq = ReadRead(readsFile, qualsFile, readsFormat, out readHeader, out qualHeader, out qual);

                if (readSeq == null)
                    break;

                readHeaders[i].CopyFrom(readHeader);
                readSeqs[i].CopyFrom(readSeq);
                if (qualHeaders != null)
                {
                    qualHeaders[i].CopyFrom(qualHeader);
                    quals[i].CopyFrom(qual);
                }

                readsRead++;
            }

            return readsRead;
        }

        // Read a FASTA read. The 'addSpace' parameter controls whether the lines of data are concatenated with a gap. For base data, no gap is wanted
        // but the same function is called for quals (for the 454 separate-qual format) and a space is needed to separate the two adjacent quals.
        private static string ReadFASTA(StreamReader reads, bool addSpace, out string readHeader)
        {
            readHeader = "";
            string currentLine = reads.ReadLine();          // read the next line from the reads file
            if (currentLine == null)                        // and return if we're at EOF
                return null;

            if (currentLine[currentLine.Length - 1] == ' ')
                currentLine = currentLine.TrimEnd();

            readHeader = currentLine;

            // multiple file lines per read
            StringBuilder currentRead = new StringBuilder(1000);
            bool keepReading = true;
            int peekChar = 0;
            currentLine = reads.ReadLine();

            while (keepReading)
            { 
                if (currentLine.Length != 0)
                {
                    currentRead.Append(currentLine);        // add the just read line to the read under construction
                    if (addSpace)
                        currentRead.Append(' ');            // leave space between lines if needed
                }

                peekChar = reads.Peek();                    // and peek at the next line to see if we're at the last line for this current read
                if (peekChar < 0 || peekChar == '>')
                    keepReading = false;                    // yes... either EOF or the header for the next read
                else
                {
                    currentLine = reads.ReadLine();         // no... next line is a read line
                    if (currentLine.Length != 0 && currentLine[currentLine.Length - 1] == ' ')
                        currentLine = currentLine.TrimEnd();
                }
            }

            return currentRead.ToString();
        }

 
        // Reads a read (headers+sequence+quals) from either a Fasta or Fastq source. Quals (if any) are always returned in one-qual-per-char form.
        //
        // Fasta files do not include qual scores, but these can come in a separate file in space-separated numeric format (old 454 data).
        // The QualFile pointer is only used for Fasta data, and will be null if there is no separate qual file associated with the reads.
        // Fastq files always come with qual data, the 3rd and 4th lines in each read.
        //
        // Fastq quals are returned as-read - in either Sanger or Solexa format.
        // Fasta quals are converted from nnn form to 0-based character form (canonical)
        //
        public static bool NextReadPlusQual(StreamReader readsFile, StreamReader qualFile, int readFormat, Sequence readHeader, Sequence read, Sequence qualHeader, Sequence quals)
        {
            string readHeaderString = null;
            string readString = null;
            string qualHeaderString = null;
            string qualsString = null;

            readString = ReadRead(readsFile, qualFile, readFormat, out readHeaderString, out qualHeaderString, out qualsString);

            if (readString == null)
            {
                readHeader.Length = 0;
                read.Length = 0;
                qualHeader.Length = 0;
                quals.Length = 0;
                return false;
            }

            readHeader.CopyFrom(readHeaderString);                
            read.CopyFrom(readString);
            qualHeader.CopyFrom(qualHeaderString);
            quals.CopyFrom(qualsString);
            if (readFormat == formatFNA && quals.Length != 0)
                ConvertQualsToCanonicalForm(quals, readFormat, 0);

            return true;
        }

        public static string NextReadPlusQual(StreamReader readsFile, StreamReader qualFile, int readFormat, out string readHeader, out string qualHeader, out string quals)
        {
            string read = null;

            read = NextReadPlusQual(readsFile, qualFile, readFormat, out readHeader, out qualHeader, out quals);

            return read;
        }

        public static void ConvertQualsToCanonicalForm(Sequence quals, int readsFormat, int qualBase)
        {
            if (readsFormat == formatFASTQ)
            {
                for (int i = 0; i < quals.Length; i++)
                    quals.Bases[i] -= (char)qualBase;
                return;
            }

            if (readsFormat == formatFNA && quals.Length != 0)
            {
                string qualString = quals.ToString();
                string[] qualStrings = qualString.Split(spaceDelimiter, StringSplitOptions.RemoveEmptyEntries);
                for (int i = 0; i < qualStrings.Length; i++)
                    quals.Bases[i] = (char)(Convert.ToInt16(qualStrings[i]));
                quals.Length = qualStrings.Length;
            }
        }

        // deprecated. Should be replaced with calls on following WriteRead (with quals parameters)
        public static void WriteReadAndQual(StreamWriter readsFile, string header, string read, string qualHeader, string quals, int readsFormat)
        {
            WriteRead(readsFile, null, header, read, qualHeader, quals, readsFormat);
        }

        public static void WriteRead(StreamWriter readsFile, StreamWriter qualsFile, string header, string read, string qualHeader, string quals, int readsFormat)
        {
            if (readsFormat == formatFASTQ)
            {
                readsFile.WriteLine(header);
                readsFile.WriteLine(read);
                readsFile.WriteLine(qualHeader);
                readsFile.WriteLine(quals);
            }
            else
            {
                WriteRead(readsFile, header, read, readsFormat);
                WriteQual(qualsFile, qualHeader, quals, readsFormat);
            }
        }

        public static void WriteRead(StreamWriter readsFile, string header, string read, string qualHeader, string quals, int readsFormat)
        {
            if (readsFormat == formatFASTQ)
            {
                readsFile.WriteLine(header);
                readsFile.WriteLine(read);
                readsFile.WriteLine(qualHeader);
                readsFile.WriteLine(quals);
            }
            else
            {
                WriteRead(readsFile, header, read, readsFormat);
            }
        }

        public static void WriteRead(StreamWriter readsFile, Sequence header, Sequence read, Sequence qualHeader, Sequence quals, int readsFormat)
        {
            if (readsFormat == formatFASTQ)
            {
                readsFile.WriteLine(header.Bases, 0, header.Length);
                readsFile.WriteLine(read.Bases, 0, read.Length);
                readsFile.WriteLine(qualHeader.Bases, 0, qualHeader.Length);
                readsFile.WriteLine(quals.Bases, 0, quals.Length);
            }
            else
            {
                WriteRead(readsFile, header.ToString(), read.ToString(), readsFormat);
            }
        }

        public static void WriteRead(StreamWriter readsFile, string readHeader, string read, int readFormat)
        {
            // @1:1:0:686#0/1
            // NTGGAGAATTCTGGATCCTCGGACTAAAACAATAGCAGTTGATTCGCTCACAGTTCTGGAGGCTAGAGGTATGAAA
            // +1:1:0:686#0/1
            // @hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh\hhhhhhhhhhhhhhhhVg[hhhU^hhWhfgVhc^Dhh`_V
            //
            // >FUOCO2J02IPX2K rank=0000260 x=3458.0 y=3626.0 length=70
            // GAGCAGCTCCTCAAGCAACTGCAACTGGATTGAGCAAGGAATGTTCGCAGCTACCCGACT
            // GACCCGTCTT

            if (readsFile == null)
                return;

            // adjust length for 454 reads
            if (readFormat == formatFNA && readHeader != null && readHeader != "")
            {
                int lengthIdx = readHeader.IndexOf("length=");
                if (lengthIdx > 0)
                {
                    int numIdx = lengthIdx + "length=".Length;
                    readHeader = readHeader.Substring(0, numIdx) + read.Length;
                }
            }

            // write out the header if it exists
            if (readHeader != null)
                readsFile.WriteLine(readHeader);

            if (readFormat == formatFNA)
            {
                // if we concatenated shorter source lines (e.g. from 454 reads) when the data was read,
                // we'll now split it again as it's being written
                int m = 0;
                int hrLen = read.Length;
                while (m < hrLen)
                {
                    int wLen = Math.Min(60, hrLen - m);
                    readsFile.WriteLine(read.Substring(m, wLen));
                    m += 60;
                }
            }
            else
                readsFile.WriteLine(read);
        }

        public static void WriteQual(StreamWriter qualFile, string qualHeader, string quals, int readFormat)
        {
            if (qualFile == null)
                return;

            int qualsLength = quals.Length;

            if (readFormat == formatFNA && qualHeader != null && qualHeader != "")
            {
                int lengthIdx = qualHeader.IndexOf("length=");
                if (lengthIdx > 0)
                {
                    int numIdx = lengthIdx + "length=".Length;
                    qualHeader = qualHeader.Substring(0, numIdx) + qualsLength;
                }
            }

            if (qualHeader != null)
                qualFile.WriteLine(qualHeader);

            if (readFormat == formatFNA)
            {
                // if we concatenated shorter qual lines (e.g. from 454 reads) when they were being read,
                // we'll now reformat and split them again as they're being written
                int m = 0;
                while (m < quals.Length)
                {
                    int wLen = Math.Min(20, qualsLength - m);
                    for (int i = m; i < m + wLen; i++)
                    {
                        int qualInt = (int)quals[i];
                        qualFile.Write(qualInt);
                        qualFile.Write(' ');
                    }
                    qualFile.WriteLine();
                    m += 20;
                }
            }

            if (readFormat == formatFASTQ)
            {
                qualFile.WriteLine(quals);
            }
        }

        // load a single tile file into a sorted MerTable
        public static long LoadTileFile(string tileFN, int minSum, int minAsRead, int minRC, int minReps,
                                        out MerTable<ulong> merArray, out int tableMerSize, out int averageDepth, out bool canonical,
                                        out long loadedUniqueMers, out long loadedTotalMers)
        {
            Console.WriteLine("Loading saved tiles from " + tileFN + " (min=" + minSum + ")");
            DateTime loadStart = DateTime.Now;

            FileInfo tileFI = new FileInfo(tileFN);
            long merArraySize = 0;
            long tileFileLength = tileFI.Length;
            merArraySize = tileFileLength / bytesPerCBTMer;     // file size in mers
            merArraySize = merArraySize / 3;                    // allow for all the low-rep mers being discarded on load

            // allocate the mer arrays
            merArray = new MerTable<ulong>(merArraySize);

            loadedUniqueMers = 0;
            loadedTotalMers = 0;
            long passableMers = 0;
            long passableTotalMers = 0;

            int tileFileType = -1;
            canonical = false;

            if (tileFN.EndsWith(".bfa"))
                tileFileType = bfaType;
            if (tileFN.EndsWith(".cbt"))
            {
                tileFileType = bfaType;
                canonical = true;
            }
            if (tileFN.EndsWith(".bfq"))
                tileFileType = bfqType;

            BinaryReader bfaFile = null;
            if (tileFileType == bfaType || tileFileType == bfqType)
                bfaFile = new BinaryReader(File.Open(tileFN, FileMode.Open, FileAccess.Read, FileShare.Read));

            // handle the header on .cbt files - source of mers and merSize
            if (canonical)
            {
                //string sourceOfMers = bfaFile.ReadString();
                merSize = bfaFile.ReadInt32();
                if (merSize <= 0 || merSize > 32)
                {
                    tableMerSize = 0;
                    averageDepth = 0;
                    return 0;
                }
            }
            else
            {
                if (merSize == 0)
                    merSize = 25;
            }

            LoadTileFile(tileFileType, bfaFile, merArray, minSum, minAsRead, minRC, minReps,
                         out loadedUniqueMers, out loadedTotalMers, out passableMers, out passableTotalMers);

            bfaFile.Close();

            averageDepth = (int)(passableTotalMers / passableMers);

            Console.WriteLine("finished reload of unique mers dictionary. " +
                                  loadedUniqueMers + "/" + loadedTotalMers + " mers, average depth " + averageDepth +
                                  " in " + (DateTime.Now - loadStart).TotalSeconds.ToString("#.0") + "s");

            tableMerSize = merSize;
            return loadedUniqueMers;
        }

        // load a single canonical Binary tiles file into a sorted MerTable
        public static long LoadCBTFile(string tileFN, int minLoadReps, int minAsRead, int minRC, int minReps,
                                        out MerTable<ulong> merArray, out int tableMerSize, out int averageDepth,
                                        out long loadedUniqueMers, out long loadedTotalMers)
        {
            Console.WriteLine("Loading saved tiles from " + tileFN + " (min=" + minLoadReps + ")");
            DateTime loadStart = DateTime.Now;

            FileInfo tileFI = new FileInfo(tileFN);
            long merArraySize = 0;
            long tileFileLength = tileFI.Length;
            merArraySize = tileFileLength / bytesPerCBTMer;     // file size in mers
            //merArraySize = merArraySize / 2;                    // allow for all the low-rep mers being discarded on load

            // allocate the mer arrays
            merArray = new MerTable<ulong>(merArraySize);

            loadedUniqueMers = 0;
            loadedTotalMers = 0;
            long passableMers = 0;
            long passableTotalMers = 0;

            BinaryReader cbtFile = null;
            cbtFile = new BinaryReader(File.Open(tileFN, FileMode.Open, FileAccess.Read, FileShare.Read));

            merSize = cbtFile.ReadInt32();
            if (merSize <= 0 || merSize > 32)
            {
                tableMerSize = 0;
                averageDepth = 0;
                return 0;
            }

            LoadCBTTiles(cbtFile, merArray, minLoadReps, minAsRead, minRC, minReps,
                         out loadedUniqueMers, out loadedTotalMers, out passableMers, out passableTotalMers);

            cbtFile.Close();

            averageDepth = (int)(passableTotalMers / passableMers);

            Console.WriteLine("Finished reload of unique mers dictionary. " +
                                  loadedUniqueMers + "/" + loadedTotalMers + " mers, average depth " + averageDepth +
                                  " in " + (DateTime.Now - loadStart).TotalSeconds.ToString("#.0") + "s");

            tableMerSize = merSize;
            return loadedUniqueMers;
        }

        private static void LoadCBTTiles(BinaryReader cbtFile, MerTable<ulong> mers,
                                         int minSum, int minAsRead, int minRC, int minReps,
                                         out long uniqueMersRead, out long totalMersRead,
                                         out long passableUniqueMersRead, out long passableTotalMersRead)
        {
            ulong packedMer = 0;                                    // packed form of next read mer
            int merAsReadCount = 0;                                 // as-read count for the next read mer
            int merRCCount = 0;                                     // RC count for the next read mer
            ulong countPair = 0;                                    // a packed count pair
            int merSumCount;                                        // and the summed count

            bool EOF = false;

            long totalMerCount = 0;
            long uniqueMerCount = 0;
            long nonSingleMerCount = 0;
            long nonSingleTotalCount = 0;

            while (!EOF)
            {
                try
                {
                    packedMer = cbtFile.ReadUInt64();
                    merAsReadCount = cbtFile.ReadInt32();
                    merRCCount = cbtFile.ReadInt32();
                    countPair = ((ulong)merAsReadCount << 32) | (ulong)merRCCount;
                }
                catch
                {
                    EOF = true;
                }
                if (EOF)
                    break;

                progressUniqueMers++;
                progressTotalMers += merAsReadCount;
                merSumCount = merAsReadCount + merRCCount;

                // add all mers tiled from the reads into our hash table
                if (merSumCount >= minSum && merAsReadCount >= minAsRead && merRCCount >= minRC)
                {
                    mers.Add(packedMer, countPair);
                    uniqueMerCount++;
                    progressLoadedUniqueMers++;
                    totalMerCount += merAsReadCount + merRCCount;
                    progressLoadedTotalMers += merAsReadCount;
                    if (merSumCount >= minReps)
                    {
                        nonSingleMerCount++;
                        nonSingleTotalCount += merSumCount;
                    }
                }

            } // all cbt tuples in the file

            mers.LoadFinished();
            uniqueMersRead = uniqueMerCount;
            totalMersRead = totalMerCount;
            passableUniqueMersRead = nonSingleMerCount;
            passableTotalMersRead = nonSingleTotalCount;
        }

        private static void LoadTileFile(int tileFileType, BinaryReader bfaFile, MerTable<ulong> mers,
                                         int minSum, int minAsRead, int minRC, int minReps,
                                         out long uniqueMersRead, out long totalMersRead,
                                         out long passableUniqueMersRead, out long passableTotalMersRead)
        {
            ulong packedMer = 0;                                    // packed form of next read mer
            int merAsReadCount = 0;                                 // as-read count for the next read mer
            int merRCCount = 0;                                     // RC count for the next read mer
            int qualScore = 0;                                      // the qual score for the mer
            ulong countPair = 0;                                    // a packed count pair
            int merSumCount;                                        // and the summed count

            bool EOF = false;

            long totalMerCount = 0;
            long uniqueMerCount = 0;
            long nonSingleMerCount = 0;
            long nonSingleTotalCount = 0;

            while (!EOF)
            {
                // string watchMer = "";

                try
                {
                    packedMer = bfaFile.ReadUInt64();
                    merAsReadCount = bfaFile.ReadInt32();
                    merRCCount = bfaFile.ReadInt32();
                    if (tileFileType == bfqType)
                        qualScore = bfaFile.ReadInt32();
                    countPair = ((ulong)merAsReadCount << 36) | ((ulong)merRCCount << 8) | (ulong)qualScore;
                    //watchMer = MerStrings.ExpandMer(packedMer);
                }
                catch
                {
                    EOF = true;
                }
                if (EOF)
                    break;

                progressUniqueMers++;
                progressTotalMers += merAsReadCount;
                merSumCount = merAsReadCount + merRCCount;

                // add all mers tiled from the reads into our hash table
                if (merSumCount >= minSum && merAsReadCount >= minAsRead && merRCCount >= minRC)
                {
                    mers.Add(packedMer, countPair);
                    uniqueMerCount++;
                    progressLoadedUniqueMers++;
                    totalMerCount += merAsReadCount + merRCCount;
                    progressLoadedTotalMers += merAsReadCount;
                    if (merSumCount >= minReps)
                    {
                        nonSingleMerCount++;
                        nonSingleTotalCount += merSumCount;
                    }
                }

            } // all tfa/bfa tuples in the file

            mers.LoadFinished();
            uniqueMersRead = uniqueMerCount;
            totalMersRead = totalMerCount;
            passableUniqueMersRead = nonSingleMerCount;
            passableTotalMersRead = nonSingleTotalCount;
        }


    }

    public class loaderThreadParams
    {
        public int tileFileType;
        public StreamReader tfaFile;
        public BinaryReader bfaFile;
        public MerDictionary<ulong> merPartition;
        public int merPartitionSize;
        public int keySize;
        public int minSum;
        public int minAsRead;
        public int minRC;
        public long partitionUniqueMers;
        public long partitionTotalMers;
    }

    // MerDictionary: a partitioned dictionary class. 
    // ----------------------------------------------
    //
    // Creates an array of dictionaries, each of which should be smaller than the maximum safe size for a Dictionary (~50M)
    //
    public class MerDictionary<TV>
    {
        // number of internal Dictionary partitions is 4 ** keySize
        int keySize = 0;
        public int noOfPartitions = 0;
        Dictionary<ulong, TV>[] dictionaryPartitions = null;
        const int maxTable = 50000000;

        public MerDictionary(long dictionarySize)
        {
            long dictionarySizePlusSome = dictionarySize + dictionarySize / 5;              // allow for some extra space in the hash table
            keySize = (int)Math.Ceiling(Math.Log(dictionarySizePlusSome / maxTable, 4));    // find minimum key size
            if (keySize < 1)
                keySize = 1;                                                                // must partition on at least one base
            noOfPartitions = (int)Math.Pow(4, keySize);                                     // giving this many partitions
            int partitionLength = (int)(dictionarySizePlusSome / noOfPartitions);           // with this average length (but scaled to reflect canonical distributions)

            dictionaryPartitions = new Dictionary<ulong, TV>[noOfPartitions];
            for (int i = 0; i < noOfPartitions; i++)
            {
                int scaledPartitionLength = 2 * partitionLength * (noOfPartitions - i) / noOfPartitions; // =(16-I2)/16*2*(F19/16)
                dictionaryPartitions[i] = new Dictionary<ulong, TV>(scaledPartitionLength);
            }
        }

        public void Add(ulong key, TV value)
        {
            int partition = 0;
            if (keySize > 0)
                partition = (int)(key >> (64 - (keySize * 2)));
            dictionaryPartitions[partition].Add(key, value);
        }

        public void Remove(ulong key)
        {
            int partition = 0;
            if (keySize > 0)
                partition = (int)(key >> (64 - (keySize * 2)));
            dictionaryPartitions[partition].Remove(key);
        }

        public bool ContainsKey(ulong key)
        {
            int partition = (int)(key >> (64 - (keySize * 2)));
            return dictionaryPartitions[partition].ContainsKey(key);
        }

        public int Count
        {
            get
            {
                int totalCount = 0;
                for (int i = 0; i < noOfPartitions; i++)
                    totalCount += dictionaryPartitions[i].Count;
                return totalCount;
            }
        }

        public void Clear()
        {
            for (int i = 0; i < noOfPartitions; i++)
                dictionaryPartitions[i].Clear();
        }

        public bool TryGetValue(ulong key, out TV value)
        {
            int partition = (int)(key >> (64 - (keySize * 2)));
            return dictionaryPartitions[partition].TryGetValue(key, out value);
        }

        public TV this[ulong key]
        {
            get
            {
                int partition = (int)(key >> (64 - (keySize * 2)));
                return dictionaryPartitions[partition][key];
            }
            set
            {
                int partition = (int)(key >> (64 - (keySize * 2)));
                dictionaryPartitions[partition][key] = value;
            }
        }

        public IEnumerator<KeyValuePair<ulong, TV>> GetEnumerator()
        {
            for (int i = 0; i < noOfPartitions; i++)
            {
                foreach (KeyValuePair<ulong, TV> kvp in dictionaryPartitions[i])
                {
                    yield return kvp;
                }
            }
        }

        public int GetPartition(int partitionNo, out ulong[] mers, out TV[] values)
        {
            mers = new ulong[dictionaryPartitions[partitionNo].Count];
            values = new TV[dictionaryPartitions[partitionNo].Count];
            int idx = 0;

            foreach (KeyValuePair<ulong, TV> kvp in dictionaryPartitions[partitionNo])
            {
                //if (kvp.Key == 0x02cb7d7c0233df70)
                //    Debugger.Break();

                mers[idx] = kvp.Key;
                values[idx] = kvp.Value;
                idx++;
            }

            Array.Sort<ulong, TV>(mers, values);

            return idx;
        }
    }

    public class MerHashSet
    {
        // number of internal HashSet partitions is 4 ** keySize
        int keySize = 3;
        int noOfPartitions = 0;
        HashSet<ulong>[] hashSetPartitions = null;

        public MerHashSet(int keySize)
        {
            this.keySize = keySize;
            this.noOfPartitions = (int)Math.Pow(4, keySize);
            hashSetPartitions = new HashSet<ulong>[noOfPartitions];
            for (int i = 0; i < noOfPartitions; i++)
                hashSetPartitions[i] = new HashSet<ulong>();
        }

        public void Add(ulong key)
        {
            int partition = 0;
            if (keySize > 0)
                partition = (int)(key >> (64 - (keySize * 2)));
            hashSetPartitions[partition].Add(key);
        }

        public void Remove(ulong key)
        {
            int partition = 0;
            if (keySize > 0)
                partition = (int)(key >> (64 - (keySize * 2)));
            hashSetPartitions[partition].Remove(key);
        }

        public bool Contains(ulong key)
        {
            int partition = (int)(key >> (64 - (keySize * 2)));
            return hashSetPartitions[partition].Contains(key);
        }

        public int Count
        {
            get
            {
                int totalCount = 0;
                for (int i = 0; i < noOfPartitions; i++)
                    totalCount += hashSetPartitions[i].Count;
                return totalCount;
            }
        }

        public void Clear()
        {
            for (int i = 0; i < noOfPartitions; i++)
                hashSetPartitions[i].Clear();
        }

        public IEnumerator<ulong> GetEnumerator()
        {
            for (int i = 0; i < noOfPartitions; i++)
            {
                foreach (ulong v in hashSetPartitions[i])
                {
                    yield return v;
                }
            }
        }
    }

    // MerTable: a potentially very large MerDictionary (read-only once loaded) 
    // -------------------------------------------------
    // Implemented using a single Dictionary for small tables (faster) or an array of length-optimised Dictionaries for larger tables (slightly slower to load).
    // These tables are assumed to be read-only once loaded (indicated by LoadFinished) being called. The mers to be loaded are assumed to be ordered. 
    // 
    public class MerTable<TV>
    {
        bool partitioned = false;           // single Dictionary or multiple
        const int maxTable = 25000000;      // safe size for a single Dictionary (can be doubled safely in scaling)
        bool dataNeedsSorting = false;      // set if out-of-order load detected
        ulong previousMer = 0;              // loads should be monotonically increasing

        Dictionary<ulong, TV> mers = null;      // small sets of mers in a single dictionary (a bit faster to load)
        Dictionary<ulong, TV>[] pmers = null;   // large sets of mers in length-optimised Dictionaries
        ulong[][] orderedMers = null;         // arrays used (and re-used) during partition loading
        TV[][] orderedValues = null;          // values are loaded into these arrays prior to the allocation of the corresponding dictionary partition
        ulong[] partitionBoundaries;        // upper possible k-mer in each partition
        int currentPartition = 0;           // current partition being loaded
        int cpi = 0;                        // index into partition buffer currently being filled (via Add calls)
        int currentBuffer = 0;
        IAsyncResult iarCopyToTable;
        CopyToTableDelegate ctd;

        int noParts = 0;                    // how many partitions
        int keyBaseShift = 0;               // # of bits to shift to get partition-key bases (bits)

        public MerTable(long tableSize)
        {
            if (tableSize > maxTable)
            {
                partitioned = true;

                int keySizeBases = (int)Math.Ceiling(Math.Log(tableSize / maxTable, 4));    // find minimum key size
                if (keySizeBases < 1)
                    keySizeBases = 1;                                                       // must partition on at least one base
                keyBaseShift = 64 - keySizeBases * 2;
                noParts = (int)Math.Pow(4, keySizeBases);                                   // giving this many partitions
                int partitionLength = (int)(tableSize / noParts);                           // and each partition starts off being this big (and then scaled to reflect the canonical
                int scaledPartitionLength = 2 * partitionLength;                            // first partition will be bigger than all the others (canonical)

                orderedMers = new ulong[2][];
                orderedValues = new TV[2][];

                for (int b = 0; b < 2; b++)
                {
                    orderedMers[b] = new ulong[scaledPartitionLength];
                    orderedValues[b] = new TV[scaledPartitionLength];
                }
                ctd = new CopyToTableDelegate(CopyToTable);

                pmers = new Dictionary<ulong, TV>[noParts];

                ulong fillBases = 0xffffffffffffffff >> (keySizeBases * 2);
                partitionBoundaries = new ulong[noParts];
                for (ulong k = 0; k < (ulong)noParts; k++)
                    partitionBoundaries[k] = k << keyBaseShift | fillBases;
            }
            else
            {
                mers = new Dictionary<ulong, TV>((int)tableSize, new merEqualityComparer());
            }
        }

        public void Add(ulong key, TV value)
        {
            if (partitioned)
            {
                // is this the first kmer to go into the next partition?
                if (key > partitionBoundaries[currentPartition])
                {
                    //Console.WriteLine(cpi + " mers added to partition " + currentPartition);
                    // wait for previous buffer copy to complete
                    if (iarCopyToTable != null && !iarCopyToTable.IsCompleted)
                    {
                        //Console.WriteLine("waiting for copy to finish");
                        iarCopyToTable.AsyncWaitHandle.WaitOne();
                    }
                    //Console.WriteLine("calling copy for partition " + currentPartition + " for buffer " + currentBuffer);
                    iarCopyToTable = ctd.BeginInvoke(currentPartition, orderedMers[currentBuffer], orderedValues[currentBuffer], cpi, null, null);
                    //CopyToTable(currentPartition, orderedMers[currentBuffer], orderedValues[currentBuffer], cpi);

                    if (currentBuffer == 0)
                        currentBuffer = 1;
                    else
                        currentBuffer = 0;

                    currentPartition++;
                    cpi = 0;
                }

                orderedMers[currentBuffer][cpi] = key;
                orderedValues[currentBuffer][cpi] = value;
                cpi++;

                if (cpi == orderedMers[currentBuffer].Length)
                {
                    int newLength = cpi + cpi / 10;
                    Array.Resize<ulong>(ref orderedMers[currentBuffer], newLength);
                    Array.Resize<TV>(ref orderedValues[currentBuffer], newLength);
                }
                if (key < previousMer)
                    dataNeedsSorting = true;
                previousMer = key;
            }
            else
            {
                mers.Add(key, value);
            }
        }

        private delegate void CopyToTableDelegate(int partitionNo, ulong[] orderedMers, TV[] orderedValues, int merCount);

        private void CopyToTable(int partitionNo, ulong[] orderedMers, TV[] orderedValues, int merCount)
        {
            //Console.WriteLine("starting copy for partition " + partitionNo);
            Dictionary<ulong, TV> pmp = new Dictionary<ulong, TV>(merCount, new merEqualityComparer());
            pmers[partitionNo] = pmp;
            for (int i = 0; i < merCount; i++)
            {
                pmp.Add(orderedMers[i], orderedValues[i]);
                //Interlocked.Increment(ref mersCopied);
            }
            //Console.WriteLine(merCount + " mers copied for partition " + partitionNo);
        }

        public bool LoadFinished()
        {
            if (dataNeedsSorting)
            {
                Console.WriteLine("kmers being loaded into table were not in sorted order");
                return false;
            }

            if (partitioned)
            {
                // wait for previous buffer copy to complete
                if (iarCopyToTable != null && !iarCopyToTable.IsCompleted)
                {
                    //Console.WriteLine("waiting for copy to finish before final copy");
                    iarCopyToTable.AsyncWaitHandle.WaitOne();
                }
                //Console.WriteLine("calling copy from LoadFinished for partition " + currentPartition + " for buffer " + currentBuffer);
                CopyToTable(currentPartition, orderedMers[currentBuffer], orderedValues[currentBuffer], cpi);
            }

            //Console.WriteLine(mersAdded + " mers added, " + mersCopied + " mers copied");
            return true;
        }

        public bool ContainsKey(ulong key)
        {
            if (partitioned)
            {
                int partIdx = (int)(key >> keyBaseShift);
                return pmers[partIdx].ContainsKey(key);
            }
            else
                return mers.ContainsKey(key);
        }

        public int Count
        {
            get
            {
                if (partitioned)
                {
                    int sum = 0;
                    for (int partIdx = 0; partIdx < noParts; partIdx++)
                        sum += pmers[partIdx].Count;
                    return sum;
                }
                else
                    return mers.Count;
            }
        }

        public void Clear()
        {
            if (partitioned)
            {
                for (int p = 0; p < noParts; p++)
                    pmers[p].Clear();
            }
            else
                mers.Clear();
        }

        public bool TryGetValue(ulong key, out TV value)
        {
            if (partitioned)
            {
                int partIdx = (int)(key >> keyBaseShift);
                return pmers[partIdx].TryGetValue(key, out value);
            }
            else
                return mers.TryGetValue(key, out value);
        }

        public TV this[ulong key]
        {
            get
            {
                if (partitioned)
                {
                    int partIdx = (int)(key >> keyBaseShift);
                    return pmers[partIdx][key];
                }
                else
                    return mers[key];
            }
            set
            {
                if (partitioned)
                {
                    int partIdx = (int)(key >> keyBaseShift);
                    pmers[partIdx][key] = value;
                }
                else
                    mers[key] = value;
            }
        }

        public IEnumerator<KeyValuePair<ulong, TV>> GetEnumerator()
        {
            if (partitioned)
            {
                for (int p = 0; p < noParts; p++)
                    foreach (KeyValuePair<ulong, TV> kvp in pmers[p])
                        yield return kvp;
            }
            else
            {
                foreach (KeyValuePair<ulong, TV> kvp in mers)
                    yield return kvp;
            }
        }
    }

    public class merEqualityComparer : IEqualityComparer<ulong>
    {
        public bool Equals(ulong x, ulong y)
        {
            return x == y;
        }

        public int GetHashCode(ulong mer)
        {
            return (mer / 7).GetHashCode();           // divide packed mer by 7 to shuffle the bits
        }
    }

    public class Sequence
    {
        public int Length;
        public int Capacity;
        public char[] Bases;

        public Sequence(int capacity)
        {
            Length = 0;
            Capacity = capacity;
            Bases = new char[capacity];
        }

        public Sequence(string s)
        {
            Length = s.Length;
            Capacity = s.Length;
            Bases = s.ToCharArray();
        }

        public void Resize(int newCapacity)
        {
            Array.Resize<char>(ref this.Bases, newCapacity);
            this.Capacity = newCapacity;
        }

        public bool Matches(string s)
        {
            if (s.Length != Length)
                return false;

            for (int i = 0; i < Length; i++)
            {
                if (Bases[i] != s[i])
                    return false;
            }

            return true;
        }

        public new string ToString()
        {
            return new string(Bases, 0, Length);
        }

        public string ToString(int start, int length)
        {
            return new string(Bases, start, length);
        }

        public void Reverse()
        {
            int halfWay = this.Length / 2;
            for (int i = 0; i < halfWay; i++)
            {
                int rhsIdx = this.Length - 1 - i;
                char savedChar = this.Bases[rhsIdx];
                this.Bases[rhsIdx] = this.Bases[i];
                this.Bases[i] = savedChar;
            }
        }

        public void Append(char b)
        {
            if (Length == Capacity)
            {
                Capacity += 50;
                Array.Resize(ref Bases, Capacity);
            }

            Bases[Length] = b;
            Length++;
        }

        public void CopyTo(Sequence copy)
        {
            if (this.Capacity > copy.Capacity)
            {
                Array.Resize<char>(ref copy.Bases, this.Capacity);
                copy.Capacity = this.Capacity;
            }
            copy.Length = this.Length;
            Array.Copy(this.Bases, 0, copy.Bases, 0, this.Length);
        }

        public void CopyFrom(string source)
        {
            if (source == null)
            {
                this.Length = 0;
                return;
            }

            if (source.Length > this.Capacity)
            {
                Array.Resize<char>(ref this.Bases, source.Length);
                this.Capacity = source.Length;
            }

            this.Length = source.Length;
            source.CopyTo(0, this.Bases, 0, source.Length);
        }

        public void Insert(int m, char c)
        {
            if (this.Length == this.Capacity)
            {
                Array.Resize(ref this.Bases, this.Capacity + 20);
                this.Capacity += 20;
            }

            for (int i = this.Length; i > m; i--)
                this.Bases[i] = this.Bases[i - 1];

            this.Bases[m] = c;
            this.Length++;
        }

        public void Replace(int m, ulong packedMer, int merSize)
        {
            ulong intBase = 0;                              // current base in binary form
            ulong tempMer = packedMer;                      // progressively shifted left
            char ACGT;                                      // char form of base

            for (int i = 0; i < merSize; i++)
            {
                intBase = (ulong)tempMer >> 62;
                tempMer = (ulong)(tempMer << 2);

                ACGT = MerStrings.baseToChar[intBase];
                this.Bases[m + i] = ACGT;
            }
        }

        public void Remove(int m)
        {
            for (int i = m; i < this.Length - 1; i++)
                this.Bases[i] = this.Bases[i + 1];
            this.Length--;
        }

        public int ComputeHash()
        {
            // not used for now. Could possibly do with a better hash algorithm as well. 
            unchecked
            {
                const int p = 16777619;
                int hash = (int)2166136261;

                for (int i = 0; i < this.Bases.Length; i++)
                    hash = (hash ^ this.Bases[i]) * p;

                hash += hash << 13;
                hash ^= hash >> 7;
                hash += hash << 3;
                hash ^= hash >> 17;
                hash += hash << 5;
                return hash;
            }
        }
    }

    public class BufferedReader
    {
        int fileFormat = MerStrings.formatNone;
        StreamReader readsFile = null;
        StreamReader qualsFile = null;

        FillBufferDelegate fbd = null;
        IAsyncResult iarFillBuffer;
        char[][] buffers = new char[2][];
        const int bufferLength = 1000000;
        int[] endOfLastRead = new int[2];
        int currentBufferIdx;
        int previousBufferIdx;
        bool EOF;
        int nextReadIdx = 0;

        public BufferedReader(int fileFormat, StreamReader readsFile, StreamReader qualsFile)
        {
            this.fileFormat = fileFormat;
            this.readsFile = readsFile;
            this.qualsFile = qualsFile;

            if (fileFormat == MerStrings.formatFASTQ)
            {
                buffers[0] = new char[bufferLength];
                buffers[1] = new char[bufferLength];
                endOfLastRead[0] = -1;
                endOfLastRead[1] = -1;
                currentBufferIdx = 0;
                previousBufferIdx = 1;
                nextReadIdx = 0;
                EOF = false;
                // start filling the initial buffer
                fbd = new FillBufferDelegate(FillBuffer);
                iarFillBuffer = fbd.BeginInvoke(readsFile, previousBufferIdx, buffers[previousBufferIdx], buffers[currentBufferIdx], 0, out endOfLastRead[previousBufferIdx], out EOF, null, null);
            }
         }

        public void Close()
        {
            if (iarFillBuffer != null && !EOF)
            {
                // wait for any in-process buffer fill to complete
                if (!iarFillBuffer.IsCompleted)
                    iarFillBuffer.AsyncWaitHandle.WaitOne();
                // get out parameters from the async call
                fbd.EndInvoke(out endOfLastRead[previousBufferIdx], out EOF, iarFillBuffer);
                // and reset the wait handle
                iarFillBuffer.AsyncWaitHandle.Close();
            }
            readsFile.Close();
            if (qualsFile != null)
                qualsFile.Close();
        }

        public int ReadReads(int readsWanted, Sequence[] readHeaders, Sequence[] readSeqs, Sequence[] qualHeaders, Sequence[] quals)
        {
            int readsReturned = 0;

            if (fileFormat == MerStrings.formatFASTQ)
            {
                while (true)
                {
                    if (nextReadIdx >= endOfLastRead[currentBufferIdx])
                    {
                        if (EOF)
                            break;

                        // wait for previous buffer fill to complete
                        if (!iarFillBuffer.IsCompleted)
                        {
                            iarFillBuffer.AsyncWaitHandle.WaitOne();
                        }

                        // get out parameters from the async call
                        fbd.EndInvoke(out endOfLastRead[previousBufferIdx], out EOF, iarFillBuffer);
                        // and reset the wait handle
                        iarFillBuffer.AsyncWaitHandle.Close();

                        if (!EOF)
                        {
                            iarFillBuffer = fbd.BeginInvoke(readsFile, currentBufferIdx, buffers[currentBufferIdx], buffers[previousBufferIdx], endOfLastRead[previousBufferIdx], out endOfLastRead[currentBufferIdx], out EOF, null, null);
                        }

                        previousBufferIdx = currentBufferIdx;
                        if (currentBufferIdx == 0)
                            currentBufferIdx = 1;
                        else
                            currentBufferIdx = 0;

                        nextReadIdx = 0;
                    }

                    int readHeaderIdx;
                    int readHeaderLen;
                    int readIdx;
                    int readLen;
                    int qualHeaderIdx;
                    int qualHeaderLen;
                    int qualIdx;
                    int qualLen;

                    char[] currentBuffer = buffers[currentBufferIdx];

                    nextReadIdx = GetNextRead(currentBuffer, nextReadIdx, out readHeaderIdx, out readHeaderLen, out readIdx, out readLen,
                                                                          out qualHeaderIdx, out qualHeaderLen, out qualIdx, out qualLen);

                    if (readHeaderLen > readHeaders[readsReturned].Capacity)
                        readHeaders[readsReturned].Resize(readHeaderLen + 50);
                    Array.Copy(currentBuffer, readHeaderIdx, readHeaders[readsReturned].Bases, 0, readHeaderLen);
                    readHeaders[readsReturned].Length = readHeaderLen;

                    if (readLen > readSeqs[readsReturned].Capacity)
                        readSeqs[readsReturned].Resize(readLen + 50);
                    Array.Copy(currentBuffer, readIdx, readSeqs[readsReturned].Bases, 0, readLen);
                    readSeqs[readsReturned].Length = readLen;

                    if (qualHeaders != null)
                    {
                        if (qualHeaderLen > qualHeaders[readsReturned].Capacity)
                            qualHeaders[readsReturned].Resize(qualHeaderLen + 50);
                        Array.Copy(currentBuffer, qualHeaderIdx, qualHeaders[readsReturned].Bases, 0, qualHeaderLen);
                        qualHeaders[readsReturned].Length = qualHeaderLen;

                        if (qualLen > quals[readsReturned].Capacity)
                            quals[readsReturned].Resize(qualLen + 50);
                        Array.Copy(currentBuffer, qualIdx, quals[readsReturned].Bases, 0, qualLen);
                        quals[readsReturned].Length = qualLen;
                    }

                    readsReturned++;
                    if (readsReturned == readsWanted)
                        break;
                }
            }
            else
                readsReturned = MerStrings.ReadReads(readsWanted, readsFile, qualsFile, fileFormat, readHeaders, readSeqs, qualHeaders, quals);

            return readsReturned;
        }

        private int GetNextRead(char[] buffer, int bi, out int readHeaderIdx, out int readHeaderLen, out int readIdx, out int readLen,
                                                              out int qualHeaderIdx, out int qualHeaderLen, out int qualIdx, out int qualLen)
        {
            int nextLF = -1;
            //starting at the initial @ of a read
            readHeaderIdx = bi;
            nextLF = Array.IndexOf<char>(buffer, '\n', readHeaderIdx);
            readHeaderLen = nextLF - readHeaderIdx;
            if (buffer[nextLF - 1] == '\r')
                readHeaderLen--;

            readIdx = nextLF + 1;
            nextLF = Array.IndexOf<char>(buffer, '\n', readIdx);
            readLen = nextLF - readIdx;
            if (buffer[nextLF - 1] == '\r')
                readLen--;

            qualHeaderIdx = nextLF + 1;
            nextLF = Array.IndexOf<char>(buffer, '\n', qualHeaderIdx);
            qualHeaderLen = nextLF - qualHeaderIdx;
            if (buffer[nextLF - 1] == '\r')
                qualHeaderLen--;

            qualIdx = nextLF + 1;
            nextLF = Array.IndexOf<char>(buffer, '\n', qualIdx);
            qualLen = nextLF - qualIdx;
            if (buffer[nextLF - 1] == '\r')
                qualLen--;

            return nextLF + 1;
        }

        private delegate void FillBufferDelegate(StreamReader reads, int bufferIdx, char[] currentBuffer, char[] previousBuffer, int leftoverIdx, out int endOfLastRead, out bool EOF);

        private void FillBuffer(StreamReader reads, int bufferIdx, char[] currentBuffer, char[] previousBuffer, int leftoverIdx, out int endOfLastRead, out bool EOF)
        {
            int bi = 0;

            if (leftoverIdx > 0)
            {
                while (leftoverIdx < bufferLength)
                {
                    currentBuffer[bi] = previousBuffer[leftoverIdx];
                    bi++;
                    leftoverIdx++;
                }
            }

            int charsToRead = bufferLength - bi;
            int charsRead = reads.ReadBlock(currentBuffer, bi, charsToRead);
            EOF = charsRead != charsToRead;
            int charsToProcess = bi + charsRead;

            bi = 0;
            endOfLastRead = -1;

            while (bi < charsToProcess)
                bi = SkipOverNextRead(currentBuffer, charsToProcess, bi, ref endOfLastRead);

            //if (currentBuffer[endOfLastRead] != '@')
            //    Debugger.Break();
        }

        private int SkipOverNextRead(char[] buffer, int bufferLength, int bi, ref int endOfLastRead)
        {
            // starts with bi pointing at @ and ends with it pointing to the @ at the start of the next read (if there was a  whole read to skip)
            int nextLF = Array.IndexOf<char>(buffer, '\n', bi);
            if (nextLF == -1)
                return bufferLength;
            bi = nextLF + 1;
            if (bi == bufferLength)
                return bufferLength;

            nextLF = Array.IndexOf<char>(buffer, '\n', bi);
            if (nextLF == -1)
                return bufferLength;
            bi = nextLF + 1;
            if (bi == bufferLength)
                return bufferLength;

            nextLF = Array.IndexOf<char>(buffer, '\n', bi);
            if (nextLF == -1)
                return bufferLength;
            bi = nextLF + 1;
            if (bi == bufferLength)
                return bufferLength;

            nextLF = Array.IndexOf<char>(buffer, '\n', bi);
            if (nextLF == -1)
                return bufferLength;
            bi = nextLF + 1;

            endOfLastRead = bi;
            return bi;
        }
    }
}
