using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;

namespace GammaTextScope
{
    /// <summary>
    /// The classify step. Each chromosome of a sample is cut into consecutive windows of
    /// commonWordsCount γwords (stop-codon delimited words) and one index line is written per
    /// window: name, start, length, QUALIFIED/UNQUALIFIED (Zipf test, <see cref="ZipfHelper"/>),
    /// isAlpha (wizard-hat test of the per-word DFA exponents), lexicon size, elapsed ms.
    /// The experiments read columns 1–5.
    /// </summary>
    public class WindowClassifier
    {
        public void RunSimplified(string species, int windowWidth, int commonWordsCount, FeatureType featureType, string printPath)
        {
            string printpath = Path.Combine(printPath, "WindowWidth" + windowWidth);
            Directory.CreateDirectory(printpath);

            string[] paths = Directory.GetFiles(PipelineConfig.Default.StreamDir(species));
            var geneIndex = SequenceHelper.ReadGeneIndexFromStream(PipelineConfig.Default.StreamDir(species + "_GeneInfo"), featureType);

            foreach (string path in paths)
            {
                var file = new FileInfo(path);
                string chromosome = file.Name.Split('_').Last().Split('.').First();
                if (!PipelineConfig.Default.IncludesChromosome(chromosome)) continue;

                var sequence = WindowHelper.ReadSequenceFromFile(path, file.Name, geneIndex, featureType);

                if (sequence.Length <= 10 * commonWordsCount * windowWidth
                    && Tokenizer.TokenizeOnStopCodons(sequence, GenomeConstants.stopWords, false).Count <= commonWordsCount * windowWidth)
                    Console.WriteLine("Sequence too short");

                var positions = new StringBuilder();
                int lastStartPosition = 0;
                int lastFinishIndex = 0;
                int retryCounter = 1;
                int windowCounter = 1; // only for the window name
                var stopWatch = Stopwatch.StartNew();

                while (lastStartPosition < sequence.Length)
                {
                    var tokenizedWithEmptyPos = SequenceHelper.FindNextLengthWindow(sequence, lastStartPosition, ref lastFinishIndex, commonWordsCount * windowWidth, out bool isEnd, ref retryCounter, windowWidth);
                    var tokenized = SequenceHelper.ConvertTokenizedWithEptyPositions(tokenizedWithEmptyPos);
                    var countDict = SequenceHelper.CreateWordCounts(tokenized);
                    bool isZipf = ZipfHelper.ClassifyDistributionDeltaCenter(countDict).Item1;
                    int length = SequenceHelper.GetTokenizedSequenceLength(tokenizedWithEmptyPos, isEnd);
                    string windowName = chromosome + "_WW" + windowWidth + "_N" + windowCounter;

                    if (isZipf && tokenizedWithEmptyPos.Count >= commonWordsCount * windowWidth)
                    {
                        bool isAlpha = IsQualifiedAlphaDistribution(tokenized);
                        positions.Append("Q_" + windowName + "\t" + lastStartPosition + "\t" + length + "\tQUALIFIED\t" + isAlpha + "\t" + countDict.Count + "\t" + stopWatch.Elapsed.TotalMilliseconds + GenomeConstants.OutputNewLine);
                    }
                    else
                    {
                        positions.Append("U_" + windowName + "\t" + lastStartPosition + "\t" + length + "\tUNQUALIFIED\t\t" + countDict.Count + "\t" + stopWatch.Elapsed.TotalMilliseconds + GenomeConstants.OutputNewLine);
                    }

                    lastStartPosition += length;
                    stopWatch.Restart();
                    retryCounter = 1;
                    lastFinishIndex = 0;
                    windowCounter++;
                }

                using var sw = new StreamWriter(Path.Combine(printpath, chromosome + "WW" + windowWidth + "_window_positions_indexLength.txt")) { NewLine = GenomeConstants.OutputNewLine };
                sw.WriteLine(positions.ToString());
            }
        }

        /// <summary>
        /// Wizard-hat ("alfa") test: the histogram (bin width 0.1) of the per-word DFA exponents
        /// must peak in the bin centred on 0.6 or 0.7, and the bin centred on 1.0 must hold at
        /// least 7 % of the peak count.
        /// </summary>
        private static bool IsQualifiedAlphaDistribution(List<string> tokenized)
        {
            var wordFluctuations = WordFluctuations.Compute(tokenized.ToArray(), GenomeConstants.AlphaWordThreshold);
            var groupings = DistributionHelper.CountNumberGroupings(wordFluctuations.Values.OrderBy(x => x).ToList(), 0.1);
            return ClassifyDistribution(groupings);
        }

        private static bool ClassifyDistribution(Dictionary<double, double> distribution)
        {
            var largestPeak = new KeyValuePair<double, double>(0.0, 0);
            foreach (var kvp in distribution)
                if (kvp.Value > largestPeak.Value)
                    largestPeak = kvp;

            if (largestPeak.Value == 0.0)
                return false;

            double peakBin = Math.Round(largestPeak.Key, 2);
            if (peakBin != 0.6 && peakBin != 0.7)
                return false;

            return distribution[1.0] / largestPeak.Value >= 0.07;
        }
    }
}
