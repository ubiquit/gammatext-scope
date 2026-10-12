using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;

namespace GammaTextScope
{
    public class WindowClassifier
    {
        public void RunSimplified(string species, int windowWidth, int commonWordsCount, FeatureType featureType, string printPath)
        {
            string printpath = Path.Combine(printPath, "WindowWidth" + windowWidth);
            Directory.CreateDirectory(printpath);

            bool skipper = true;

            string pathLocal = PipelineConfig.Default.StreamDir(species);
            string[] paths = System.IO.Directory.GetFiles(pathLocal);
            string chromosomeNameString = "";

            //List<double> timeStampsRegistry = new List<double>();

            var geneIndex = SequenceHelper.ReadGeneIndexFromStream(PipelineConfig.Default.StreamDir(species + "_GeneInfo"), featureType);

            if (!SequenceHelper.TestCreateGeneticDistances() || !SequenceHelper.TestCreateGeneticTextNew())
                throw new Exception("genetic distance algorithm faliure");

            foreach (string path in paths)
            {
                Dictionary<string, List<Tuple<string, int>>> WordMapping = new Dictionary<string, List<Tuple<string, int>>>();
                System.IO.FileInfo file = new System.IO.FileInfo(path);
                string chromosome = file.Name.Split('_').Last().Split('.').First();
                if (!PipelineConfig.Default.IncludesChromosome(chromosome)) continue;

                //skipping
                //if (chromosome.Equals("chr22") && windowWidth == 26 && skipper)
                //{
                //    skipper = false;
                //}
                //else if (windowWidth != 26)
                //{
                //}
                //else if (skipper)
                //{
                //    continue;
                //}

                //skipping chromosomes
                //if (chromosome.Equals("chr2"))
                //    continue;

                string positionString = "";

                var sequence = WindowHelper.ReadSequenceFromFile(path, file.Name, ref geneIndex, featureType);
                bool longGeneticSequence = true;

                //todo, try this to not have repetitive WWs
                if (sequence.Length <= (10 * commonWordsCount * windowWidth))
                {
                    var tokenizedWithEmptyPos = Tokenizer.TokenizeOnStopCodons(sequence, GenomeConstants.stopWords, false);

                    if (tokenizedWithEmptyPos.Count <= (commonWordsCount * windowWidth))
                    {
                        Console.WriteLine("Sequence too short");
                        longGeneticSequence = false;
                    }
                }

                int lastStartPosition = 0;
                int retryCounter = 1;
                int lastFinishIndex = 0;
                int windowCounter = 1; //just for file name
                bool isEnd = false;
                bool isNullableWindow = false;

                //if(chromosome != "chr37")
                //{
                //    int i = 0;
                //    continue;
                //}

                Stopwatch stopWatch = new Stopwatch();
                stopWatch.Start();

                while (lastStartPosition < sequence.Length)
                {
                    var tokenizedWithEmptyPos = SequenceHelper.FindNextLengthWindow(sequence, lastStartPosition, ref lastFinishIndex, commonWordsCount * windowWidth, out isEnd, ref retryCounter, windowWidth);
                    var tokenized = SequenceHelper.ConvertTokenizedWithEptyPositions(tokenizedWithEmptyPos);

                    bool isPrintableWindow = false;

                    var countDict = SequenceHelper.CreateWordCounts(tokenized);
                    var isZipfTuple1 = ZipfHelper.ClassifyDistributionDeltaCenter(countDict);
                    //var isZipfTupl2 = ZipfHelper.ClassifyDistributionSlope(countDict);
                    string windowName = chromosome + "_WW" + windowWidth + "_N" + windowCounter;
                    Tuple<bool, Dictionary<string, double>> isQualifiedDoc = null;

                    if (isZipfTuple1.Item1 && tokenizedWithEmptyPos.Count >= (commonWordsCount * windowWidth))
                    {
                        isQualifiedDoc = IsQualifiedAlphaDistribution(ref tokenizedWithEmptyPos, chromosomeNameString);
                        var lexicon = GetLexicon(tokenized).Count;
                        windowName = "Q_" + windowName;
                        int length = SequenceHelper.GetTokenizedSequenceLength(tokenizedWithEmptyPos, isEnd);
                        positionString += windowName + "\t" + lastStartPosition + "\t" + length + "\tQUALIFIED\t" + isQualifiedDoc.Item1 + "\t" + lexicon + "\t" + stopWatch.Elapsed.TotalMilliseconds + Environment.NewLine;
                        isPrintableWindow = true;
                        isNullableWindow = true;
                    }
                    else
                    {
                        windowName = "U_" + windowName;
                        int length = SequenceHelper.GetTokenizedSequenceLength(tokenizedWithEmptyPos, isEnd); // using previous to retain last window
                        var lexicon = GetLexicon(tokenized).Count;
                        positionString += windowName + "\t" + lastStartPosition + "\t" + length + "\tUNQUALIFIED\t\t" + lexicon + "\t" + stopWatch.Elapsed.TotalMilliseconds + Environment.NewLine;
                        lastStartPosition += (length);
                        isPrintableWindow = false;
                        isNullableWindow = true;
                    }

                    if (isPrintableWindow)
                    {
                        int length = SequenceHelper.GetTokenizedSequenceLength(tokenizedWithEmptyPos, isEnd);
                        lastStartPosition += length;


                        isPrintableWindow = false;
                        isNullableWindow = true;
                    }

                    if (isNullableWindow)
                    {
                        //timeStampsRegistry.Add(stopWatch.Elapsed.TotalMilliseconds);
                        stopWatch.Restart();
                        retryCounter = 0;
                        lastFinishIndex = 0;
                        windowCounter++;
                        isEnd = false;
                        isNullableWindow = false;
                        isQualifiedDoc = null;
                    }

                    retryCounter++;
                }

                //writing indow positions for chromosome
                using System.IO.StreamWriter sw = new System.IO.StreamWriter(Path.Combine(printpath, chromosome + "WW" + windowWidth + "_window_positions_indexLength.txt"));
                sw.WriteLine(positionString);
            }
        }

        private HashSet<string> GetLexicon(List<string> tokenized)
        {
            HashSet<string> lexicon = new HashSet<string>();

            foreach (var w in tokenized)
            {
                lexicon.Add(w);
            }

            return lexicon;
        }

        private static void EnsureDirectoryExists(string filePath)
        {
            FileInfo fi = new FileInfo(filePath);
            if (!fi.Directory.Exists)
            {
                System.IO.Directory.CreateDirectory(fi.DirectoryName);
            }
        }

        private Tuple<bool, Dictionary<string, double>> IsQualifiedAlphaDistribution(ref List<string> tokenizedWithEmptyPos, string chromosome)
        {
            var tokenized = SequenceHelper.ConvertTokenizedWithEptyPositions(tokenizedWithEmptyPos);
            var wordFluctuations = WordFluctuations.Compute(tokenized.ToArray(), GenomeConstants.AlphaWordThreshold);
            var groupings = DistributionHelper.CountNumberGroupings(wordFluctuations.Values.OrderBy(x => x).ToList(), 0.1);
            var tuple = ClassifyDistribution(groupings);
            return new Tuple<bool, Dictionary<string, double>>(tuple.Item1, wordFluctuations);
        }

        private static Tuple<bool, double> ClassifyDistribution(Dictionary<double, double> distribution)
        {
            double totalCounts = 0;
            var largestPeak = new KeyValuePair<double, double>(0.0, 0);

            foreach (KeyValuePair<double, double> kvp in distribution)
            {
                totalCounts += kvp.Value;

                if (kvp.Value > largestPeak.Value)
                {
                    largestPeak = kvp;
                }
            }

            if (largestPeak.Value == 0.0)
                return new Tuple<bool, double>(false, double.NaN);

            if (Math.Round(largestPeak.Key, 2) != 0.6 && Math.Round(largestPeak.Key, 2) != 0.7)
                return new Tuple<bool, double>(false, Math.Round(largestPeak.Key, 2));

            if (Math.Round(largestPeak.Key, 2) == 0.7)
            {
                int i = 0;
            }

            if ((distribution[1.0] / largestPeak.Value) < 0.07)
                return new Tuple<bool, double>(false, Math.Round(largestPeak.Key, 2));

            return new Tuple<bool, double>(true, Math.Round(largestPeak.Key, 2));
        }

    }
}
