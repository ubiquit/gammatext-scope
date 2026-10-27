using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace GammaTextScope
{
    /// <summary>
    /// Analyses the spacing between alfa-decay markers from StopCodonAlphaComparison.
    ///
    /// For the genome-wide window sequence, this experiment replicates the decay detection
    /// logic (sliding window of 10, threshold of 3) to find where decay dots would be placed.
    /// It then measures the number of windows between consecutive dots, producing a time series
    /// of inter-dot gaps. This gap series is submitted to Detrended Fluctuation Analysis (DFA)
    /// to determine whether the spacing of decay events exhibits fractal/long-range correlation.
    ///
    /// Separate analyses are run for:
    ///   - Normal-side decay dots (tumor has lost qualification relative to normal)
    ///   - Tumor-side decay dots (normal has lost qualification relative to tumor)
    ///   - Each further split into alfa-dominant and zipf-dominant dots
    /// </summary>
    public class DecayGapDfa
    {
        private const int SlidingSize = 10;
        private const int MinDelta = 3;
        private const int DfaThreshold = 30;

        public void Run(string pairName, string featureLabel)
        {
            string scaDir = PipelineConfig.Default.TumorOutDir(pairName, "StopCodonAlphaComparison", featureLabel) + System.IO.Path.DirectorySeparatorChar;
            string tsvPath = Path.Combine(scaDir, $"StopCodonAlpha_{pairName}_AllChromosomes.tsv");

            if (!File.Exists(tsvPath))
            {
                Console.WriteLine($"  DecayGapDFA: TSV not found: {tsvPath}");
                return;
            }

            string outDir = PipelineConfig.Default.TumorOutDir(pairName, "DecayGapDFA", featureLabel) + System.IO.Path.DirectorySeparatorChar;
            if (!Directory.Exists(outDir)) Directory.CreateDirectory(outDir);

            Console.WriteLine($"  DecayGapDFA ({featureLabel}): reading {Path.GetFileName(tsvPath)}...");

            // Parse the genome-wide TSV
            var (normalClasses, tumorClasses) = ParseGenomeTsv(tsvPath);
            int n = Math.Min(normalClasses.Count, tumorClasses.Count);

            Console.WriteLine($"  DecayGapDFA: {n} windows loaded");

            if (n < SlidingSize)
            {
                Console.WriteLine("  DecayGapDFA: too few windows, skipping.");
                return;
            }

            // Detect decay dot positions (replicating ChartPrinter.AddDecayDeltaMarkers logic)
            var normalSideAlfa = new List<int>();  // positions where Normal has more qualified (alfa-dominant)
            var normalSideZipf = new List<int>();  // positions where Normal has more qualified (zipf-dominant)
            var tumorSideAlfa = new List<int>();    // positions where Tumor has more qualified (alfa-dominant)
            var tumorSideZipf = new List<int>();    // positions where Tumor has more qualified (zipf-dominant)

            for (int center = 0; center < n; center++)
            {
                int start = Math.Max(0, center - SlidingSize / 2);
                int end = Math.Min(n, start + SlidingSize);
                if (end - start < SlidingSize && end == n) start = Math.Max(0, end - SlidingSize);

                int nAlfa = 0, nZipf = 0, tAlfa = 0, tZipf = 0;
                for (int j = start; j < end; j++)
                {
                    if (normalClasses[j] == "alfa") nAlfa++;
                    else if (normalClasses[j] == "zipf") nZipf++;
                    if (tumorClasses[j] == "alfa") tAlfa++;
                    else if (tumorClasses[j] == "zipf") tZipf++;
                }

                int normalQualified = nAlfa + nZipf;
                int tumorQualified = tAlfa + tZipf;
                int delta = normalQualified - tumorQualified;

                if (Math.Abs(delta) >= MinDelta)
                {
                    if (delta > 0)
                    {
                        // Normal has more → tumor decayed
                        if (nAlfa > nZipf)
                            normalSideAlfa.Add(center);
                        else
                            normalSideZipf.Add(center);
                    }
                    else
                    {
                        // Tumor has more → normal decayed
                        if (tAlfa > tZipf)
                            tumorSideAlfa.Add(center);
                        else
                            tumorSideZipf.Add(center);
                    }
                }
            }

            // Also compute combined lists (alfa + zipf on each side)
            var normalSideAll = normalSideAlfa.Concat(normalSideZipf).OrderBy(x => x).ToList();
            var tumorSideAll = tumorSideAlfa.Concat(tumorSideZipf).OrderBy(x => x).ToList();

            // Compute gap series and run DFA for each category
            var results = new List<string>();
            results.Add("Category\tDotCount\tGapCount\tMinGap\tMaxGap\tMeanGap\tDFA_Alpha");
            results.Add("");

            RunCategory("NormalSide_Alfa", normalSideAlfa, results, outDir);
            RunCategory("NormalSide_Zipf", normalSideZipf, results, outDir);
            RunCategory("NormalSide_All", normalSideAll, results, outDir);
            RunCategory("TumorSide_Alfa", tumorSideAlfa, results, outDir);
            RunCategory("TumorSide_Zipf", tumorSideZipf, results, outDir);
            RunCategory("TumorSide_All", tumorSideAll, results, outDir);

            // Write summary
            string summaryPath = Path.Combine(outDir, $"DecayGapDFA_{pairName}_Summary.txt");
            File.WriteAllText(summaryPath, string.Join(GenomeConstants.OutputNewLine, results));
            Console.WriteLine($"  DecayGapDFA done. Summary: {summaryPath}");
        }

        private void RunCategory(string name, List<int> dotPositions,
                                  List<string> results, string outDir)
        {
            if (dotPositions.Count < 2)
            {
                results.Add($"{name}\t{dotPositions.Count}\t0\t-\t-\t-\tN/A (too few dots)");
                Console.WriteLine($"    {name}: {dotPositions.Count} dots – too few for gap analysis");
                return;
            }

            // Compute gaps between consecutive dots
            var gaps = new List<double>();
            for (int i = 1; i < dotPositions.Count; i++)
            {
                gaps.Add(dotPositions[i] - dotPositions[i - 1]);
            }

            // Write gap series to file
            string gapPath = Path.Combine(outDir, $"GapSeries_{name}.tsv");
            var gapLines = new List<string> { "GapIndex\tGapLength" };
            for (int i = 0; i < gaps.Count; i++)
                gapLines.Add($"{i + 1}\t{gaps[i]:F0}");
            File.WriteAllText(gapPath, string.Join(GenomeConstants.OutputNewLine, gapLines));

            // Run DFA
            double alpha = -1000;
            string failureReason = null;
            try
            {
                if (gaps.Count >= DfaThreshold)
                {
                    alpha = WordFluctuations.DfaAlpha(gaps, DfaThreshold);
                }
                else if (gaps.Count >= 15)
                {
                    // Use smaller threshold if we have at least 15 gaps
                    alpha = WordFluctuations.DfaAlpha(gaps, gaps.Count / 2);
                }
            }
            catch (Exception ex)
            {
                // DFA can throw when an internal box-size partition yields a degenerate
                // 1-point regression (e.g. very few or evenly spaced gaps). Treat as N/A
                // so the rest of the categories/run still complete.
                failureReason = ex.Message.Split('\n')[0];
                Console.WriteLine($"    {name}: DFA failed ({failureReason})");
            }

            string alphaStr = failureReason != null
                ? $"N/A (DFA error)"
                : (alpha <= -999 ? "N/A (insufficient data)" : Math.Round(alpha, 4).ToString());
            double minGap = gaps.Min();
            double maxGap = gaps.Max();
            double meanGap = gaps.Average();

            results.Add($"{name}\t{dotPositions.Count}\t{gaps.Count}\t{minGap:F0}\t{maxGap:F0}\t{meanGap:F1}\t{alphaStr}");
            Console.WriteLine($"    {name}: {dotPositions.Count} dots, {gaps.Count} gaps, mean={meanGap:F1}, DFA alpha={alphaStr}");
        }

        private static (List<string> normalClasses, List<string> tumorClasses) ParseGenomeTsv(string path)
        {
            var normalClasses = new List<string>();
            var tumorClasses = new List<string>();

            bool header = true;
            int normalClassCol = -1, tumorClassCol = -1;

            foreach (var line in File.ReadLines(path))
            {
                if (header)
                {
                    header = false;
                    var cols = line.Split('\t');
                    for (int i = 0; i < cols.Length; i++)
                    {
                        if (cols[i] == "Normal_Class") normalClassCol = i;
                        else if (cols[i] == "Tumor_Class") tumorClassCol = i;
                    }
                    if (normalClassCol < 0 || tumorClassCol < 0)
                        throw new Exception($"Missing Normal_Class or Tumor_Class column in {path}");
                    continue;
                }

                var parts = line.Split('\t');
                if (parts.Length > Math.Max(normalClassCol, tumorClassCol))
                {
                    normalClasses.Add(parts[normalClassCol].ToLowerInvariant());
                    tumorClasses.Add(parts[tumorClassCol].ToLowerInvariant());
                }
            }

            return (normalClasses, tumorClasses);
        }
    }
}
