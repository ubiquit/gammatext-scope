using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GammaTextScope
{
    public static class SequenceHelper
    {
        /// <summary>
        /// Creates lexicon from the files in the given path. The words need to be separated by the given separator
        /// </summary>
        /// <param name="catName"></param>
        public static HashSet<string> ReadLexiconFromStream(string catName, string separator)
        {
            string pathLocal = PipelineConfig.Default.StreamDir(catName);
            string[] paths = System.IO.Directory.GetFiles(pathLocal);
            HashSet<string> lexicon = new HashSet<string>();

            foreach (string path in paths)
            {
                System.IO.FileInfo file = new System.IO.FileInfo(path);
                using (var reader = new StreamReader(path))
                {
                    string geneSeq = "";
                    while (!reader.EndOfStream)
                    {
                        string line = reader.ReadLine();

                        if (line.StartsWith(separator))
                        {
                            lexicon.Add(geneSeq);
                            geneSeq = "";

                            continue;
                        }

                        geneSeq += line;
                    }
                }
            }

            return lexicon;
        }

        /// <summary>
        /// Creates gene index from GFF3 file. Downloaded from: 
        /// ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/
        /// </summary>
        /// <param name="catName"></param>
        public static Dictionary<string, List<EnsembleGeneEntry>> ReadGeneIndexFromStream(string catName)
        {
            string pathLocal = PipelineConfig.Default.StreamDir(catName);
            string[] paths = System.IO.Directory.GetFiles(pathLocal);
            Dictionary<string, List<EnsembleGeneEntry>> geneIndex = new Dictionary<string, List<EnsembleGeneEntry>>();

            foreach (string path in paths)
            {
                System.IO.FileInfo file = new System.IO.FileInfo(path);
                using (var reader = new StreamReader(path))
                {
                    while (!reader.EndOfStream)
                    {
                        string line = reader.ReadLine();
                        var gene = CreateEnsembleGeneFromLine(line);

                        if (gene != null)
                        {
                            UpdateGeneIndex(gene, ref geneIndex);
                        }
                    }
                }
            }

            return geneIndex;
        }

        /// <summary>
        /// Creates gene index from GFF3 file. Downloaded from: 
        /// ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/
        /// </summary>
        /// <param name="catName"></param>
        public static Dictionary<string, List<EnsembleGeneEntry>> ReadCanisLupusFamiliarisGeneIndexFromStream(string catName)
        {
            string pathLocal = PipelineConfig.Default.StreamDir(catName);
            string[] paths = System.IO.Directory.GetFiles(pathLocal);
            Dictionary<string, List<EnsembleGeneEntry>> geneIndex = new Dictionary<string, List<EnsembleGeneEntry>>();
            int counter = 0;

            foreach (string path in paths)
            {
                System.IO.FileInfo file = new System.IO.FileInfo(path);
                using (var reader = new StreamReader(path))
                {
                    while (!reader.EndOfStream)
                    {
                        string line = reader.ReadLine();
                        var gene = CreateGnomonGeneFromLine(line);

                        if (gene != null)
                        {
                            counter++;
                            UpdateGeneIndex(gene, ref geneIndex);
                        }
                    }
                }
            }

            return geneIndex;
        }

        /// <summary>
        /// Creates gene index from GFF3 file. Downloaded from: 
        /// ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/
        /// </summary>
        /// <param name="catName"></param>
        public static Dictionary<string, List<EnsembleGeneEntry>> ReadGeneIndexFromStream(string catName, FeatureType featureType)
        {
            string pathLocal = catName;
            string[] paths = System.IO.Directory.GetFiles(pathLocal);
            Dictionary<string, List<EnsembleGeneEntry>> geneIndex = new Dictionary<string, List<EnsembleGeneEntry>>();
            int counter = 0;

            if (featureType.Equals(FeatureType.WholeGenome) && paths.Length < 1)
                return geneIndex;

            //Created to idenify correct exons
            Dictionary<string, List<EnsembleGeneEntry>> geneIndexForExons = featureType == FeatureType.Exon ? ReadGeneIndexFromStream(catName, FeatureType.Gene) : null;

            foreach (string path in paths)
            {
                System.IO.FileInfo file = new System.IO.FileInfo(path);
                using (var reader = new StreamReader(path))
                {
                    while (!reader.EndOfStream)
                    {
                        string line = reader.ReadLine();
                        var gene = ReadGeneFeatureLine(line, featureType, geneIndexForExons);

                        if (gene != null)
                        {
                            counter++;
                            UpdateGeneIndex(gene, ref geneIndex);
                        }
                    }
                }
            }

            return geneIndex;
        }

        public static Dictionary<string, List<double>> MakeDistances(Dictionary<string, List<EnsembleGeneEntry>> index)
        {
            var iiter = index.GetEnumerator();
            Dictionary<string, List<double>> geneDists = new Dictionary<string, List<double>>();

            while (iiter.MoveNext())
            {
                if (iiter.Current.Value.Count <= 1)
                    continue;

                if (!geneDists.ContainsKey(iiter.Current.Key))
                    geneDists.Add(iiter.Current.Key, new List<double>());

                for (int i = 1; i < iiter.Current.Value.Count; i++) //one indexed to not retrieve index -1
                {
                    if (iiter.Current.Value[i].StartPosition > iiter.Current.Value[i - 1].EndPosition)
                    {
                        geneDists[iiter.Current.Key].Add(iiter.Current.Value[i].StartPosition - iiter.Current.Value[i - 1].EndPosition);
                    }
                }
            }

            return geneDists;
        }

        public static bool TestCreateGeneticDistances()
        {
            var gene1 = new EnsembleGeneEntry { StartPosition = 20, EndPosition = 50, Chromosome = "chr1" };
            var gene2 = new EnsembleGeneEntry { StartPosition = 30, EndPosition = 50, Chromosome = "chr1" };
            var gene3 = new EnsembleGeneEntry { StartPosition = 30, EndPosition = 100, Chromosome = "chr1" };
            var gene4 = new EnsembleGeneEntry { StartPosition = 120, EndPosition = 150, Chromosome = "chr1" };
            var gene5 = new EnsembleGeneEntry { StartPosition = 130, EndPosition = 170, Chromosome = "chr1" };
            var gene6 = new EnsembleGeneEntry { StartPosition = 171, EndPosition = 180, Chromosome = "chr1" };
            var gene7 = new EnsembleGeneEntry { StartPosition = 190, EndPosition = 200, Chromosome = "chr1" };
            var gene8 = new EnsembleGeneEntry { StartPosition = 230, EndPosition = 245, Chromosome = "chr1" };
            var gene9 = new EnsembleGeneEntry { StartPosition = 247, EndPosition = 249, Chromosome = "chr1" };
            var gene10 = new EnsembleGeneEntry { StartPosition = 250, EndPosition = 300, Chromosome = "chr1" };
            var gene11 = new EnsembleGeneEntry { StartPosition = 120, EndPosition = 150, Chromosome = "chr2" };
            var gene12 = new EnsembleGeneEntry { StartPosition = 100, EndPosition = 300, Chromosome = "chr2" };

            Dictionary<string, List<EnsembleGeneEntry>> index = new Dictionary<string, List<EnsembleGeneEntry>>();
            index.Add("gene1", new List<EnsembleGeneEntry> { gene1 });
            index.Add("gene2", new List<EnsembleGeneEntry> { gene2 });
            index.Add("gene3", new List<EnsembleGeneEntry> { gene3 });
            index.Add("gene4", new List<EnsembleGeneEntry> { gene4 });
            index.Add("gene5", new List<EnsembleGeneEntry> { gene5 });
            index.Add("gene6", new List<EnsembleGeneEntry> { gene6 });
            index.Add("gene7", new List<EnsembleGeneEntry> { gene7 });
            index.Add("gene8", new List<EnsembleGeneEntry> { gene8 });
            index.Add("gene9", new List<EnsembleGeneEntry> { gene9 });
            index.Add("gene10", new List<EnsembleGeneEntry> { gene10 });
            index.Add("gene11", new List<EnsembleGeneEntry> { gene11 });

            List<Tuple<int, int, string>> interDists = CreateGeneticDistances(index, "chr1", 0, 240);

            if (interDists.Count != 5)
                return false;

            interDists = CreateGeneticDistances(index, "chr1", 241, 290);

            if (interDists.Count != 3)
                return false;

            return true;
        }

        public static bool TestCreateGeneticText()
        {
            var gene1 = new EnsembleGeneEntry { StartPosition = 6, EndPosition = 10, Chromosome = "chr1" };
            var gene2 = new EnsembleGeneEntry { StartPosition = 8, EndPosition = 15, Chromosome = "chr1" };
            var gene3 = new EnsembleGeneEntry { StartPosition = 26, EndPosition = 30, Chromosome = "chr1" };
            var gene4 = new EnsembleGeneEntry { StartPosition = 38, EndPosition = 48, Chromosome = "chr1" };
            var gene5 = new EnsembleGeneEntry { StartPosition = 135, EndPosition = 150, Chromosome = "chr1" };
            var gene6 = new EnsembleGeneEntry { StartPosition = 290, EndPosition = 305, Chromosome = "chr1" };
            var gene7 = new EnsembleGeneEntry { StartPosition = 311, EndPosition = 315, Chromosome = "chr1" };
            var gene8 = new EnsembleGeneEntry { StartPosition = 321, EndPosition = 325, Chromosome = "chr1" };
            var gene9 = new EnsembleGeneEntry { StartPosition = 331, EndPosition = 340, Chromosome = "chr1" };

            var gene10 = new EnsembleGeneEntry { StartPosition = 15, EndPosition = 30, Chromosome = "chr2" };


            Dictionary<string, List<EnsembleGeneEntry>> index = new Dictionary<string, List<EnsembleGeneEntry>>();
            index.Add("gene1", new List<EnsembleGeneEntry> { gene1 });
            index.Add("gene2", new List<EnsembleGeneEntry> { gene2 });
            index.Add("gene3", new List<EnsembleGeneEntry> { gene3 });
            index.Add("gene4", new List<EnsembleGeneEntry> { gene4 });
            index.Add("gene5", new List<EnsembleGeneEntry> { gene5 });
            index.Add("gene6", new List<EnsembleGeneEntry> { gene6 });
            index.Add("gene7", new List<EnsembleGeneEntry> { gene7 });
            index.Add("gene8", new List<EnsembleGeneEntry> { gene8 });
            index.Add("gene9", new List<EnsembleGeneEntry> { gene9 });
            index.Add("gene10", new List<EnsembleGeneEntry> { gene10 });

            string windowSeq = "AAAAATTTTTTTTTTAAAAAAAAAATTTTTAAAAA";

            var strings = CreateIntergenicText(index, windowSeq, 1, 36, "chr1");

            if (!strings[0].Equals("AAAAA"))
                return false;

            if (!strings[1].Equals("AAAAAAAAAA"))
                return false;

            if (!strings[2].Equals("AAAAA"))
                return false;

            windowSeq = "AAATTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAA";
            strings = CreateIntergenicText(index, windowSeq, 35, 70, "chr1");

            if (!strings[0].Equals("AAA"))
                return false;

            if (!strings[1].Equals("AAAAAAAAAAAAAAAAAAAAAA"))
                return false;

            windowSeq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTT";
            strings = CreateIntergenicText(index, windowSeq, 105, 140, "chr1");

            if (!strings[0].Equals("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"))
                return false;

            windowSeq = "TTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAA";
            strings = CreateIntergenicText(index, windowSeq, 140, 175, "chr1");

            if (!strings[0].Equals("AAAAAAAAAAAAAAAAAAAAAAAAA"))
                return false;

            windowSeq = "TTTTTAAAAATTTTTAAAAATTTTTAAAAATTTTT";
            strings = CreateIntergenicText(index, windowSeq, 301, 336, "chr1");

            if (!strings[0].Equals("AAAAA"))
                return false;

            if (!strings[1].Equals("AAAAA"))
                return false;

            if (!strings[2].Equals("AAAAA"))
                return false;

            return true;
        }

        public static bool TestCreateGeneticTextNew()
        {
            var gene1 = new EnsembleGeneEntry { StartPosition = 5, EndPosition = 11, Chromosome = "chr1" }; //test 1
            var gene2 = new EnsembleGeneEntry { StartPosition = 7, EndPosition = 15, Chromosome = "chr1" };
            var gene3 = new EnsembleGeneEntry { StartPosition = 25, EndPosition = 30, Chromosome = "chr1" };
            var gene4 = new EnsembleGeneEntry { StartPosition = 38, EndPosition = 48, Chromosome = "chr1" };  //test 2
            var gene5 = new EnsembleGeneEntry { StartPosition = 135, EndPosition = 150, Chromosome = "chr1" }; //test 3
            var gene6 = new EnsembleGeneEntry { StartPosition = 290, EndPosition = 305, Chromosome = "chr1" };
            var gene7 = new EnsembleGeneEntry { StartPosition = 310, EndPosition = 315, Chromosome = "chr1" };
            var gene8 = new EnsembleGeneEntry { StartPosition = 320, EndPosition = 325, Chromosome = "chr1" };
            var gene9 = new EnsembleGeneEntry { StartPosition = 330, EndPosition = 340, Chromosome = "chr1" };

            var gene10 = new EnsembleGeneEntry { StartPosition = 15, EndPosition = 30, Chromosome = "chr2" };


            Dictionary<string, List<EnsembleGeneEntry>> index = new Dictionary<string, List<EnsembleGeneEntry>>();
            index.Add("gene1", new List<EnsembleGeneEntry> { gene1 });
            index.Add("gene2", new List<EnsembleGeneEntry> { gene2 });
            index.Add("gene3", new List<EnsembleGeneEntry> { gene3 });
            index.Add("gene4", new List<EnsembleGeneEntry> { gene4 });
            index.Add("gene5", new List<EnsembleGeneEntry> { gene5 });
            index.Add("gene6", new List<EnsembleGeneEntry> { gene6 });
            index.Add("gene7", new List<EnsembleGeneEntry> { gene7 });
            index.Add("gene8", new List<EnsembleGeneEntry> { gene8 });
            index.Add("gene9", new List<EnsembleGeneEntry> { gene9 });
            index.Add("gene10", new List<EnsembleGeneEntry> { gene10 });

            string windowSeq = "AAAAATTTTTTTTTTAAAAAAAAAATTTTTAAAAA";
            int windowstart = 0;
            var strings = SequenceInverter(index, windowSeq, windowstart, windowstart + windowSeq.Length, "chr1");

            if (!strings[0].Equals("AAAAA"))
                return false;

            if (!strings[1].Equals("AAAAAAAAAA"))
                return false;

            if (!strings[2].Equals("AAAAA"))
                return false;

            windowSeq = "AAATTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAA";
            windowstart = 35;
            strings = SequenceInverter(index, windowSeq, windowstart, windowstart + windowSeq.Length, "chr1");

            if (!strings[0].Equals("AAA"))
                return false;

            if (!strings[1].Equals("AAAAAAAAAAAAAAAAAAAAAA"))
                return false;

            windowSeq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTT";
            windowstart = 105;
            strings = SequenceInverter(index, windowSeq, windowstart, windowstart + windowSeq.Length, "chr1");

            if (!strings[0].Equals("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"))
                return false;

            windowSeq = "TTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAA";
            windowstart = 140;
            strings = SequenceInverter(index, windowSeq, windowstart, windowstart + windowSeq.Length, "chr1");

            if (!strings[0].Equals("AAAAAAAAAAAAAAAAAAAAAAAAA"))
                return false;

            windowSeq = "TTTTTAAAAATTTTTAAAAATTTTTAAAAATTTTT";
            windowstart = 300;
            strings = SequenceInverter(index, windowSeq, windowstart, windowstart + windowSeq.Length, "chr1");

            if (!strings[0].Equals("AAAAA"))
                return false;

            if (!strings[1].Equals("AAAAA"))
                return false;

            if (!strings[2].Equals("AAAAA"))
                return false;

            return true;
        }

        /// <summary>
        /// Creates the intergenig sequence acording to the genes passed in on the index. Genes inside this method is 1 indexed, according to Gencode coordinates 
        /// </summary>
        /// <param name="index"></param>
        /// <param name="windowSeq"></param>
        /// <param name="windowStart"></param>
        /// <param name="windowEnd"></param>
        /// <param name="currentChromosome"></param>
        /// <returns></returns>
        public static List<string> CreateIntergenicText(Dictionary<string, List<EnsembleGeneEntry>> index, string windowSeq, int windowStart, int windowEnd, string currentChromosome)
        {
            List<string> interTexts = new List<string>();
            List<Tuple<int, int, string>> interDists = CreateGeneticDistances(index, currentChromosome, windowStart, windowEnd);

            //defent null or empty collection
            if (interDists == null || interDists.Count < 1)
                return interTexts;

            //this is to handle the case of only one distance in interdists
            if (interDists.Count < 2)
            {
                var t = interDists.First();
                if (t.Item1 > windowStart)
                    interTexts.Add(new string(windowSeq.ToCharArray().SubArray(0, t.Item1 - windowStart)));

                if (t.Item2 < windowEnd)
                    interTexts.Add(new string(windowSeq.ToCharArray().SubArray(t.Item2 - windowStart, windowEnd - t.Item2)));

                return interTexts;
            }

            for (int i = 1; i < interDists.Count; i++)
            {
                var t2 = interDists[i];
                var t1 = interDists[i - 1];

                //append first section on first iteration
                if (i == 1 && t1.Item1 > windowStart)
                    interTexts.Add(new string(windowSeq.ToCharArray().SubArray(0, t1.Item1 - windowStart)));

                if (t2.Item1 > t1.Item2)
                {
                    //subtracting window start here to compensate for sequence index of the intergenig positions
                    string s = new string(windowSeq.ToCharArray().SubArray((t1.Item2 - windowStart) + 1, t2.Item1 - t1.Item2 - 1));
                    interTexts.Add(s);
                }
            }

            //+1 to compensate for the zero indexed window positions when using it in the number ogf elements parameter
            if (interDists.Last().Item2 < windowEnd)
                interTexts.Add(new string(windowSeq.ToCharArray().SubArray((interDists.Last().Item2 - windowStart) + 1, (windowEnd - interDists.Last().Item2) - 1)));

            return interTexts;
        }

        /// <summary>
        /// Trying new sequence inverting method
        /// </summary>
        /// <param name="text"></param>
        /// <param name="dictionary"></param>
        /// <param name="lengthOfSeq"></param>
        /// <returns></returns>
        public static List<string> SequenceInverter(Dictionary<string, List<EnsembleGeneEntry>> index, string windowSeq, int windowStart, int windowEnd, string currentChromosome)
        {
            //windowEnd -= 1;

            // start, end, chromosome name
            List<Tuple<int, int, string>> interDists = CreateGeneticDistances(index, currentChromosome, windowStart, windowEnd);
            // base symbol, sequence position, is intragenic
            List<IntragenicSequenceMarke> sequenceIndex = new List<IntragenicSequenceMarke>();

            //in case there are no genes in the window sequence
            if (interDists == null)
                return new List<string> { windowSeq };

            //create array of entries to keep window position and chromosome position in sync
            char[] baseArray = windowSeq.ToCharArray();
            int windowStartCounter = windowStart;
            foreach (char c in baseArray)
            {
                sequenceIndex.Add(new IntragenicSequenceMarke
                {
                    BaseSymbol = c.ToString(),
                    SequencePosition = windowStartCounter,
                    IsIntragenic = true
                });

                windowStartCounter++;
            }

            //mark genetic base symbols
            foreach (Tuple<int, int, string> tuple in interDists)
            {
                int start = tuple.Item1;
                int distanceDiff = (tuple.Item2 - tuple.Item1);
                var startingMarker = sequenceIndex.Find(x => x.SequencePosition == start);
                int indexOfStartingPosition = sequenceIndex.IndexOf(startingMarker);
                int end = (indexOfStartingPosition + distanceDiff);

                for (; indexOfStartingPosition < end; indexOfStartingPosition++)
                {
                    //REMOVE after debug
                    //if (indexOfStartingPosition == sequenceIndex.Count)
                    //    continue;

                    sequenceIndex[indexOfStartingPosition].IsIntragenic = false;
                }
            }

            StringBuilder sb = new StringBuilder();
            ////TEST
            //foreach (var sequenceMarker in sequenceIndex)
            //{
            //    if (!sequenceMarker.IsIntragenic)
            //    {
            //        continue;
            //    }

            //    sb.Append(sequenceMarker.BaseSymbol);
            //}

            //string lal = sb.ToString();
            //sb = new StringBuilder();

            //build intragenic strings for return
            List<string> intragenicStrings = new List<string>();
            int cou = 0;
            foreach (var sequenceMarker in sequenceIndex)
            {
                if (!sequenceMarker.IsIntragenic)
                {
                    if (sb.Length > 0)
                    {
                        intragenicStrings.Add(sb.ToString());
                        sb = new StringBuilder();
                    }
                    else
                    {
                        cou++;
                    }

                    continue;
                }

                sb.Append(sequenceMarker.BaseSymbol);
            }

            //add last intragenic string
            if (sb.Length > 0)
                intragenicStrings.Add(sb.ToString());

            return intragenicStrings;
        }

        public static Dictionary<string, List<int>> ReadTimeSeriesFromString(string text, HashSet<string> dictionary, int lengthOfSeq)
        {
            char[] textChars = text.ToCharArray();
            int lastEventDist = 0;
            Dictionary<string, List<int>> seqDistances = new Dictionary<string, List<int>>();

            for (int i = 0; i < textChars.Length; i++)
            {
                if (textChars.Length < (i + lengthOfSeq + 1))
                    continue;

                var seqArr = textChars.SubArray(i, lengthOfSeq);
                string seq = string.Join("", seqArr);

                foreach (string seqWord in dictionary)
                {
                    if (seqWord.Equals(seq))
                    {
                        //verify that there should be 1-indexing here
                        if (seqDistances.ContainsKey(seqWord))
                        {
                            int sinceLastOccurrence = (i + 1) - seqDistances[seqWord].Sum();
                            seqDistances[seqWord].Add(sinceLastOccurrence);
                        }
                        else
                        {
                            seqDistances.Add(seqWord, new List<int> { i + 1 });
                        }
                    }
                }
            }

            return seqDistances;
        }

        public static List<string> FindNextSlidingWindow(string sequence, int lastStartPosition, int retryCounter, ref int lastFinishIndex, out bool isEnd)
        {
            string endString = sequence.Substring(lastStartPosition);
            var tokenizedProcessingString = Tokenizer.TokenizeOnStopCodons(endString, GenomeConstants.stopWords, false);
            isEnd = false;

            Dictionary<string, int> WordOccurrenceStore = new Dictionary<string, int>();

            for (int i = 0; i < tokenizedProcessingString.Count; i++)
            {
                string word = tokenizedProcessingString[i];

                if (i < lastFinishIndex)
                    continue;

                if (word.Equals(string.Empty))
                    continue;

                lastFinishIndex = i;

                if (WordOccurrenceStore.ContainsKey(word))
                    WordOccurrenceStore[word]++;
                else
                    WordOccurrenceStore.Add(word, 1);

                if (IsNumberOfKWordsReached(ref WordOccurrenceStore))
                {
                    var tokenizedList = tokenizedProcessingString.GetRange(0, i + 1);

                    return tokenizedList;
                }
            }

            //for unqualified end string
            isEnd = true;
            return tokenizedProcessingString;
        }

        public static List<string> FindNextDictionaryWindow(string sequence, int lastStartPosition, ref HashSet<string> dictionary, ref int lastDictionaryWordCount, ref int lastFinishIndex, out bool isEnd, ref int retry)
        {
            string endString = sequence.Substring(lastStartPosition);
            var tokenizedProcessingString = Tokenizer.TokenizeOnStopCodons(endString, GenomeConstants.stopWords, false);
            isEnd = false;
            int dictionaryWordCount = lastDictionaryWordCount;

            Dictionary<string, int> WordOccurrenceStore = new Dictionary<string, int>();
            for (int i = 0; i < tokenizedProcessingString.Count; i++)
            {
                string word = tokenizedProcessingString[i];

                if (i < dictionary.Count * retry)
                    continue;

                if (dictionary.Count * retry > tokenizedProcessingString.Count)
                    break;

                if (i <= lastFinishIndex)
                    continue;

                if (word.Equals(string.Empty))
                    continue;

                if (dictionary.Contains(word))
                    dictionaryWordCount++;

                lastFinishIndex = i;

                if (IsDictionaryWordRatioReached(i + 1, dictionaryWordCount))
                {
                    var tokenizedList = tokenizedProcessingString.GetRange(0, i + 1);
                    lastDictionaryWordCount = dictionaryWordCount;
                    lastFinishIndex += dictionary.Count; //TODO iffy
                    return tokenizedList;
                }
            }

            //for unqualified end string
            isEnd = true;
            lastDictionaryWordCount = dictionaryWordCount;
            return tokenizedProcessingString;
        }

        public static List<string> FindNextLengthWindow(string sequence, int lastStartPosition, ref int lastFinishIndex, int length, out bool isEnd, ref int retry, int windowWidth)
        {
            string endString = sequence.Substring(lastStartPosition);
            var tokenizedProcessingString = Tokenizer.TokenizeOnStopCodons(endString, GenomeConstants.stopWords, false);
            isEnd = false;
            //int buffer = (int)(length * windowWidth * 0.2); //delete buffer reference after extension experiment

            for (int i = 0; i < tokenizedProcessingString.Count; i++)
            {
                if (i < length * retry)
                    continue;

                if (length * retry > tokenizedProcessingString.Count)
                    break;

                if (i <= lastFinishIndex)
                    continue;

                string word = tokenizedProcessingString[i];
                if (word.Equals(string.Empty))
                    continue;

                var tokenizedList = tokenizedProcessingString.GetRange(0, i + 1);
                lastFinishIndex += length; //TODO iffy

                return tokenizedList;
            }

            //for unqualified end string
            isEnd = true;
            return tokenizedProcessingString;
        }

        public static int GetTokenizedSequenceLength(List<string> tokenized, bool isEnd)
        {
            int length = 0;
            int emptyCounter = 1; // emptyCounter is to compensate for areas where there are several stop codons in a row
            for (int i = 0; i < tokenized.Count; i++)
            {
                string w = tokenized[i];

                if (w == string.Empty)
                {
                    emptyCounter++;

                    if (i == tokenized.Count - 1)
                        length += (emptyCounter * 3); // addition here to add possible empty window lengths
                }
                else
                {
                    length += (w.Length + (emptyCounter * 3)); // + 3 to compensate for lost stop codon length
                    emptyCounter = 1;
                }
            }

            if (isEnd)
                length -= 3; //subtract here if this si the last window. there is no end element

            return length;
        }

        public static List<string> ConvertTokenizedWithEptyPositions(List<string> tokenizedwithEmptyPositions)
        {
            List<string> tokenized = new List<string>();

            foreach (var w in tokenizedwithEmptyPositions)
                if (w != string.Empty)
                    tokenized.Add(w);

            return tokenized;
        }

        public static Dictionary<string, int> CreateWordCounts(List<string> tokenizedSequence)
        {
            Dictionary<string, int> countDict = new Dictionary<string, int>();

            foreach (string word in tokenizedSequence)
                if (countDict.ContainsKey(word))
                    countDict[word]++;
                else
                    countDict.Add(word, 1);

            return countDict;
        }

        private static bool IsNumberOfKWordsReached(ref Dictionary<string, int> WordOccurrenceStore)
        {
            int thresholdCounter = 0;
            foreach (var kvp in WordOccurrenceStore)
                if (kvp.Value >= 30)
                    thresholdCounter++;

            if (thresholdCounter >= GenomeConstants.SlidingWindowKWordThreshold)
                return true;

            return false;
        }

        private static bool IsDictionaryWordRatioReached(int numberOfWords, int numberOfDictionaryWords)
        {
            if ((double)numberOfDictionaryWords / numberOfWords > 0.3)
                return true;

            return false;
        }

        /// <summary>
        /// created the genetic distances of transcribing sequence.
        /// </summary>
        /// <param name="index"></param>
        /// <param name="currentChromosome"></param>
        /// <param name="windowStart"></param>
        /// <param name="windowEnd"></param>
        /// <returns></returns>
        private static List<Tuple<int, int, string>> CreateGeneticDistances(Dictionary<string, List<EnsembleGeneEntry>> index, string currentChromosome, int windowStart, int windowEnd)
        {
            List<Tuple<int, int, string>> distsList = new List<Tuple<int, int, string>>();

            int test = 0;
            int totGeneticSeq = 0;
            foreach (var entry in index)
            {
                foreach (var gene in entry.Value)
                {
                    if (gene.Chromosome == "chr1")
                    {
                        totGeneticSeq += (gene.EndPosition - gene.StartPosition);
                        test++;
                    }
                }
            }

            var iiter = index.GetEnumerator();
            while (iiter.MoveNext())
            {
                foreach (var geneEntry in iiter.Current.Value)
                {
                    //outside of window
                    if (geneEntry.EndPosition < windowStart || geneEntry.StartPosition > windowEnd || !geneEntry.Chromosome.Equals(currentChromosome, StringComparison.InvariantCultureIgnoreCase))
                        continue;

                    distsList.Add(new Tuple<int, int, string>(geneEntry.StartPosition, geneEntry.EndPosition, geneEntry.Chromosome));
                }
            }

            distsList.OrderBy(x => x.Item3).ThenBy(x => x.Item1);
            List<Tuple<int, int, string>> retList = new List<Tuple<int, int, string>>();

            if (distsList.Count < 1)
                return null;

            int start = windowStart;
            //runs through the genetic distances, and attempts to merge the distances which are overlapping
            for (int i = 0; i < distsList.Count;)
            {
                var tuple = distsList[i];

                Tuple<int, int, int> intervalProbe = FindEndOfInterval(distsList.GetRange(i, (distsList.Count - i)));

                int end = intervalProbe.Item2 > windowEnd ? windowEnd : intervalProbe.Item2; //to not exceed window
                int beginning = intervalProbe.Item1 < windowStart ? windowStart : intervalProbe.Item1; //to not exceed window

                retList.Add(new Tuple<int, int, string>(beginning, end, distsList[i].Item3));

                i += intervalProbe.Item3;
            }

            return retList;
        }

        /// <summary>
        /// To probe the length of a interval. Selects the longest when given a choice of two ends
        /// </summary>
        /// <param name="distsList"></param>
        /// <returns></returns>
        private static Tuple<int, int, int> FindEndOfInterval(List<Tuple<int, int, string>> distsList)
        {
            if (distsList.Count < 1)
                return null;

            int index = 1;
            int start = distsList[0].Item1;
            int end = distsList[0].Item2;
            for (int i = 1; i < distsList.Count; i++)
            {
                if (end < distsList[i].Item1)
                    return new Tuple<int, int, int>(start, end, index);

                index++;
                end = distsList[i].Item2 > end ? distsList[i].Item2 : end;
            }

            //final distance
            return new Tuple<int, int, int>(start, end, index);
        }

        private static EnsembleGeneEntry CreateEnsembleGeneFromLine(string lineIn)
        {
            string[] line = lineIn.Split('\t');

            if (line.Length <= 1)
                return null;

            int startPos = -1;
            int endPos = -1;

            if (line.Length > 2 && line[2].Equals("gene") && int.TryParse(line[3], out startPos) && int.TryParse(line[4], out endPos))
            {
                var geneData = line[8].Split(';');
                var geneType = geneData[2].Split('=')[1];
                var level = int.Parse(geneData[4].Split('=')[1]);
                var transcriptType = line[1];
                if (geneType.Equals("protein_coding") && level == 1)
                {
                    return new EnsembleGeneEntry
                    {
                        ID = geneData[0].Split('=')[1],
                        GeneName = geneData[3].Split('=')[1],
                        Level = level,
                        GeneType = geneType,
                        TranscriptType = transcriptType,
                        StartPosition = startPos - 1,
                        EndPosition = endPos - 1,
                        Chromosome = line[0]
                    };
                }
            }

            return null;
        }

        private static EnsembleGeneEntry CreateGnomonGeneFromLine(string lineIn)
        {
            string[] line = lineIn.Split('\t');

            if (line.Length <= 1)
                return null;

            int startPos = -1;
            int endPos = -1;

            if (line.Length > 2 && line[2].Equals("gene") && int.TryParse(line[3], out startPos) && int.TryParse(line[4], out endPos))
            {
                var geneData = line[8].Split(';');
                var geneType = geneData.FirstOrDefault(x => x.StartsWith("gene_biotype"))?.Split("=")[1];
                if (geneType == null) geneType = geneData.FirstOrDefault(x => x.StartsWith("gene_type"))?.Split("=")[1];

                var transcriptType = line[1];

                if (lineIn.Contains("protein_coding"))
                {
                    int i = 0;
                }

                if (geneType != null && geneType.Equals("protein_coding"))
                {
                    return new EnsembleGeneEntry
                    {
                        ID = geneData[0].Split('=')[1],
                        GeneName = geneData[0].Split('=')[1],
                        GeneType = geneType,
                        TranscriptType = transcriptType,
                        StartPosition = startPos - 1,
                        EndPosition = endPos - 1,
                        Chromosome = line[0]
                    };
                }
            }

            return null;
        }

        private static EnsembleGeneEntry ReadGeneFeatureLine(string lineIn, FeatureType featureType, Dictionary<string, List<EnsembleGeneEntry>> geneIndex = null)
        {
            string[] line = lineIn.Split('\t');

            if (line.Length <= 1)
                return null;

            EnsembleGeneEntry entry;

            switch (featureType)
            {
                case FeatureType.Gene:
                    entry = CreateGeneFromLine(line);
                    break;
                case FeatureType.Exon:
                    entry = CreateExonFromLine(line, geneIndex);
                    break;
                case FeatureType.Intergenic:
                    entry = CreateGeneFromLine(line);
                    break;
                case FeatureType.WholeGenome:
                    entry = CreateGeneFromLine(line);
                    break;
                default:
                    entry = CreateGeneFromLine(line);
                    break;
            }

            return entry;
        }

        private static EnsembleGeneEntry CreateGeneFromLine(string[] line)
        {
            int startPos = -1;
            int endPos = -1;

            if (line.Length > 2 && line[2].Equals("gene") && int.TryParse(line[3], out startPos) && int.TryParse(line[4], out endPos))
            {
                var geneData = line[8].Split(';');
                var geneType = geneData.FirstOrDefault(x => x.StartsWith("gene_biotype"))?.Split("=")[1];
                if(geneType == null) geneType = geneData.FirstOrDefault(x => x.StartsWith("gene_type"))?.Split("=")[1];
                if(geneType == null) geneType = geneData.FirstOrDefault(x => x.StartsWith("biotype"))?.Split("=")[1];
                var transcriptType = line[1];

                if (geneType != null && geneType.Equals("protein_coding"))
                {
                    return new EnsembleGeneEntry
                    {
                        ID = geneData[0].Split('=')[1],
                        GeneName = geneData[0].Split('=')[1],
                        GeneType = geneType,
                        TranscriptType = transcriptType,
                        StartPosition = startPos - 1,
                        EndPosition = endPos - 1,
                        Chromosome = line[0]
                    };
                }
            }

            return null;
        }

        private static EnsembleGeneEntry CreateExonFromLine(string[] line, Dictionary<string, List<EnsembleGeneEntry>> workingGeneIndex)
        {
            int startPos = -1;
            int endPos = -1;

            if (line.Length > 2 && line[2].Equals("exon") && int.TryParse(line[3], out startPos) && int.TryParse(line[4], out endPos))
            {
                var geneData = line[8].Split(';');
                var geneDataField = geneData.FirstOrDefault(x => x.StartsWith("gene="));
                string parentName = null;
                if (geneDataField != null)
                {
                    parentName = geneDataField.Split("=")[1];
                }


                //declaring thios to handle exons where locus tag is used
                string locus_tag = null;
                if (geneData.Any(x => x.StartsWith("locus_tag")))
                    locus_tag = geneData.FirstOrDefault(x => x.StartsWith("locus_tag")).Split("=")[1];

                var transcriptType = line[1];

                bool hasEnsembleGene = locus_tag == null ? false : workingGeneIndex.Any(x => x.Value.Any(y => y.GeneName.Equals("gene-" + locus_tag)));

                if (parentName != null || (workingGeneIndex.ContainsKey("gene-" + parentName) || hasEnsembleGene))
                {
                    try
                    {
                        var gene = hasEnsembleGene ? workingGeneIndex["gene-" + locus_tag].First() : workingGeneIndex["gene-" + parentName].First();
                        return new EnsembleGeneEntry
                        {
                            ID = geneData[0].Split('=')[1],
                            GeneName = geneData[0].Split('=')[1],
                            GeneType = gene.GeneType,
                            TranscriptType = transcriptType,
                            StartPosition = startPos - 1,
                            EndPosition = endPos - 1,
                            Chromosome = line[0]
                        };
                    }
                    catch (Exception e)
                    {
                        return null;
                    }
                }

                // Ensembl format: exons have Parent=transcript:... instead of gene=...
                // Match exon to a gene by checking if it falls within a gene's position range on the same chromosome
                if (parentName == null && !hasEnsembleGene)
                {
                    string chromosome = line[0];
                    int exonStart = startPos - 1;
                    int exonEnd = endPos - 1;

                    foreach (var kvp in workingGeneIndex)
                    {
                        var gene = kvp.Value.FirstOrDefault(g => g.Chromosome == chromosome && exonStart >= g.StartPosition && exonEnd <= g.EndPosition);
                        if (gene != null)
                        {
                            return new EnsembleGeneEntry
                            {
                                ID = geneData[0].Split('=')[1],
                                GeneName = geneData[0].Split('=')[1],
                                GeneType = gene.GeneType,
                                TranscriptType = transcriptType,
                                StartPosition = exonStart,
                                EndPosition = exonEnd,
                                Chromosome = chromosome
                            };
                        }
                    }
                }
            }

            return null;
        }

        private static void UpdateGeneIndex(EnsembleGeneEntry entry, ref Dictionary<string, List<EnsembleGeneEntry>> index)
        {
            if (index.ContainsKey(entry.GeneName))
            {
                index[entry.GeneName].Add(entry);

                return;
            }

            index.Add(entry.GeneName, new List<EnsembleGeneEntry> { entry });

            return;
        }


    }
}

public class IntragenicSequenceMarke
{
    public string BaseSymbol { get; set; }
    public int SequencePosition { get; set; }
    public bool IsIntragenic { get; set; }
}

public enum NucleobaseNames { A, C, G, T, U };

public class EnsembleGeneEntry
{
    public string ID { get; set; }
    public string GeneName { get; set; }
    public int Level { get; set; }
    public string Chromosome { get; set; }
    public string GeneType { get; set; }
    public string TranscriptType { get; set; }
    public int StartPosition { get; set; }
    public int EndPosition { get; set; }
}

public static class extensions
{
    public static T[] SubArray<T>(this T[] data, int index, int length)
    {
        T[] result = new T[length];
        Array.Copy(data, index, result, 0, length);
        return result;
    }
}


