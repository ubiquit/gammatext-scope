using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace GammaTextScope
{
    /// <summary>
    /// Counts stop codons (TAA, TAG, TGA) per 100,000 nt bins for Normal vs Tumor,
    /// outputs per-chromosome and genome-wide summaries, and tests per-bin differences.
    /// Supports feature-type filtering (Intergenic, Gene, Exon, WholeGenome).
    /// </summary>
    public class StopCodonDensity
    {
        private const int BinSize = 100_000;

        public void Run(string speciesNormal, string speciesTumor, string pairName = null,
                        FeatureType featureType = FeatureType.WholeGenome)
        {
            pairName ??= speciesNormal.Replace("\\", "_") + "_VS_" + speciesTumor.Replace("\\", "_");
            string featureLabel = featureType.ToString();

            string outDir = PipelineConfig.Default.TumorOutDir(pairName, "StopCodonDensity", featureLabel) + System.IO.Path.DirectorySeparatorChar;
            EnsureDirectory(outDir);

            // Map: chromosome -> binIndex -> stats
            var perChr = new Dictionary<string, Dictionary<int, BinStats>>(StringComparer.OrdinalIgnoreCase);

            // Gather chromosome file paths
            var normalChrFiles = GetChrFiles(speciesNormal);
            var tumorChrFiles = GetChrFiles(speciesTumor);

            // Load gene annotations if feature-filtering is needed
            Dictionary<string, List<EnsembleGeneEntry>> normalGeneIndex = null;
            Dictionary<string, List<EnsembleGeneEntry>> tumorGeneIndex = null;
            if (featureType != FeatureType.WholeGenome)
            {
                Console.WriteLine($"  Loading gene annotations for {featureLabel}...");
                normalGeneIndex = SequenceHelper.ReadGeneIndexFromStream(PipelineConfig.Default.StreamDir(speciesNormal + "_GeneInfo"), featureType);
                tumorGeneIndex = SequenceHelper.ReadGeneIndexFromStream(PipelineConfig.Default.StreamDir(speciesTumor + "_GeneInfo"), featureType);
            }

            // Chromosome union
            var chromosomes = new HashSet<string>(normalChrFiles.Keys, StringComparer.OrdinalIgnoreCase);
            foreach (var c in tumorChrFiles.Keys) chromosomes.Add(c);

            foreach (var chr in chromosomes.OrderBy(NaturalChrOrderKey))
            {
                if (!PipelineConfig.Default.IncludesChromosome(chr)) continue;
                Console.WriteLine($"  StopCodonDensity ({featureLabel}): {chr}");
                if (!perChr.ContainsKey(chr)) perChr[chr] = new Dictionary<int, BinStats>();

                string normalSeq = LoadSequence(normalChrFiles, chr, normalGeneIndex, featureType);
                string tumorSeq = LoadSequence(tumorChrFiles, chr, tumorGeneIndex, featureType);

                // Normal: bin using Normal coordinates (the reference frame)
                if (normalSeq.Length > 0)
                    CountStopsIntoBins(normalSeq, perChr[chr], isNormal: true);

                // Tumor: align to Normal coordinates via offset map, then bin
                if (tumorSeq.Length > 0 && normalSeq.Length > 0)
                {
                    var offsetMap = AlignmentOffsetMap.Build(normalSeq, tumorSeq);
                    Console.WriteLine($"    alignment: {offsetMap.BreakpointCount} breakpoints, final offset {offsetMap.FinalOffset}");
                    CountStopsIntoBinsAligned(normalSeq, tumorSeq, offsetMap, perChr[chr]);
                }
                else if (tumorSeq.Length > 0)
                {
                    CountStopsIntoBins(tumorSeq, perChr[chr], isNormal: false);
                }
            }

            // ---------- write per-chromosome density ----------
            string perChrPath = Path.Combine(outDir, $"{pairName}_StopCodonDensity_100k_perChromosome.tsv");
            using (var sw = new StreamWriter(perChrPath))
            {
                sw.WriteLine("Chromosome\tBinStart\tBinEnd\tNormal_count\tTumor_count\tNormal_trials\tTumor_trials\tNormal_rate_per100k\tTumor_rate_per100k\tDelta_per100k\tRate_ratio");
                foreach (var chr in perChr.Keys.OrderBy(NaturalChrOrderKey))
                {
                    foreach (var kv in perChr[chr].OrderBy(k => k.Key))
                    {
                        int binIndex = kv.Key;
                        var s = kv.Value;

                        int binStart = binIndex * BinSize;
                        int binEnd = binStart + s.MaxLen();

                        double nRate = s.NormalTrials > 0 ? s.NormalCount * (100000.0 / s.NormalTrials) : 0.0;
                        double tRate = s.TumorTrials > 0 ? s.TumorCount * (100000.0 / s.TumorTrials) : 0.0;
                        double delta = tRate - nRate;
                        double ratio = (s.NormalCount == 0 && s.TumorCount == 0) ? 1.0 :
                                       (s.NormalCount == 0) ? double.PositiveInfinity :
                                       (s.TumorCount == 0) ? 0.0 : (tRate / nRate);

                        sw.WriteLine($"{chr}\t{binStart}\t{binEnd}\t{s.NormalCount}\t{s.TumorCount}\t{s.NormalTrials}\t{s.TumorTrials}\t{nRate:0.####}\t{tRate:0.####}\t{delta:0.####}\t{ratio:0.####}");
                    }
                }
            }

            // ---------- genome-wide bin summary ----------
            string sumPath = Path.Combine(outDir, $"{pairName}_StopCodonDensity_100k_genomeSummary.tsv");
            using (var sw = new StreamWriter(sumPath))
            {
                sw.WriteLine("BinStart\tBinEnd\tNormal_count\tTumor_count\tNormal_trials\tTumor_trials\tNormal_rate_per100k\tTumor_rate_per100k\tDelta_per100k\tRate_ratio");
                var allBinIdx = perChr.Values.SelectMany(d => d.Keys).Distinct().OrderBy(i => i);
                foreach (var i in allBinIdx)
                {
                    long nCount = 0, tCount = 0, nTrials = 0, tTrials = 0, nLen = 0, tLen = 0;
                    foreach (var chr in perChr.Keys)
                    {
                        if (perChr[chr].TryGetValue(i, out var s))
                        {
                            nCount += s.NormalCount;
                            tCount += s.TumorCount;
                            nTrials += s.NormalTrials;
                            tTrials += s.TumorTrials;
                            nLen += s.NormalLen;
                            tLen += s.TumorLen;
                        }
                    }
                    int binStart = i * BinSize;
                    int binEnd = binStart + (int)Math.Max(nLen, tLen);

                    double nRate = nTrials > 0 ? nCount * (100000.0 / nTrials) : 0.0;
                    double tRate = tTrials > 0 ? tCount * (100000.0 / tTrials) : 0.0;
                    double delta = tRate - nRate;
                    double ratio = (nCount == 0 && tCount == 0) ? 1.0 :
                                   (nCount == 0) ? double.PositiveInfinity :
                                   (tCount == 0) ? 0.0 : (tRate / nRate);

                    sw.WriteLine($"{binStart}\t{binEnd}\t{nCount}\t{tCount}\t{nTrials}\t{tTrials}\t{nRate:0.####}\t{tRate:0.####}\t{delta:0.####}\t{ratio:0.####}");
                }
            }

            // ---------- significance testing per bin ----------
            const double AlphaFDR = 0.05;
            const double MinPctDiff = 0.05;

            var tests = new List<BinTestRow>(capacity: 1 << 15);
            foreach (var chr in perChr.Keys)
            {
                foreach (var kv in perChr[chr])
                {
                    int i = kv.Key;
                    var s = kv.Value;

                    var (p, z) = TwoPropZTest(s.NormalCount, s.NormalTrials, s.TumorCount, s.TumorTrials);

                    int binStart = i * BinSize;
                    int binEnd = binStart + s.MaxLen();

                    double nRate = s.NormalTrials > 0 ? s.NormalCount * (100000.0 / s.NormalTrials) : 0.0;
                    double tRate = s.TumorTrials > 0 ? s.TumorCount * (100000.0 / s.TumorTrials) : 0.0;
                    double delta = tRate - nRate;
                    double ratio = (s.NormalCount == 0 && s.TumorCount == 0) ? 1.0 :
                                   (s.NormalCount == 0) ? double.PositiveInfinity :
                                   (s.TumorCount == 0) ? 0.0 : (tRate / nRate);

                    tests.Add(new BinTestRow
                    {
                        Chr = chr, BinIndex = i, BinStart = binStart, BinEnd = binEnd,
                        NormalCount = s.NormalCount, TumorCount = s.TumorCount,
                        NormalTrials = s.NormalTrials, TumorTrials = s.TumorTrials,
                        NormalRate100k = nRate, TumorRate100k = tRate,
                        Delta100k = delta, RateRatio = ratio, Z = z, P = p
                    });
                }
            }

            BenjaminiHochbergInPlace(tests);

            string diffPath = Path.Combine(outDir, $"{pairName}_StopCodonDensity_100k_diff.tsv");
            using (var sw = new StreamWriter(diffPath))
            {
                sw.WriteLine("Chromosome\tBinStart\tBinEnd\tNormal_count\tTumor_count\tNormal_trials\tTumor_trials\tNormal_rate_per100k\tTumor_rate_per100k\tDelta_per100k\tRate_ratio\tZ\tP\tQ");
                foreach (var r in tests.OrderBy(x => x.Q).ThenByDescending(x => Math.Abs(x.Delta100k)))
                    sw.WriteLine($"{r.Chr}\t{r.BinStart}\t{r.BinEnd}\t{r.NormalCount}\t{r.TumorCount}\t{r.NormalTrials}\t{r.TumorTrials}\t{r.NormalRate100k:0.####}\t{r.TumorRate100k:0.####}\t{r.Delta100k:0.####}\t{r.RateRatio:0.####}\t{r.Z:0.###}\t{r.P:0.000E+00}\t{r.Q:0.000E+00}");
            }

            string sigPath = Path.Combine(outDir, $"{pairName}_StopCodonDensity_100k_significant.tsv");
            using (var sw = new StreamWriter(sigPath))
            {
                sw.WriteLine("Chromosome\tBinStart\tBinEnd\tNormal_rate_per100k\tTumor_rate_per100k\tDelta_per100k\tRate_ratio\tNormal_count\tTumor_count\tNormal_trials\tTumor_trials\tZ\tP\tQ");
                foreach (var r in tests.OrderBy(x => x.Chr).ThenBy(x => x.BinStart))
                {
                    if (r.Q > AlphaFDR) continue;
                    double denom = Math.Max(r.NormalRate100k, 1e-9);
                    double relDiff = Math.Abs(r.TumorRate100k - r.NormalRate100k) / denom;
                    if (relDiff < MinPctDiff) continue;
                    sw.WriteLine($"{r.Chr}\t{r.BinStart}\t{r.BinEnd}\t{r.NormalRate100k:0.####}\t{r.TumorRate100k:0.####}\t{r.Delta100k:0.####}\t{r.RateRatio:0.####}\t{r.NormalCount}\t{r.TumorCount}\t{r.NormalTrials}\t{r.TumorTrials}\t{r.Z:0.###}\t{r.P:0.000E+00}\t{r.Q:0.000E+00}");
                }
            }

            Console.WriteLine($"Significant-only bins written:\n  {sigPath}  (FDR≤{AlphaFDR}, min relative diff ≥ {MinPctDiff:P0})");

            // ---------- Charts ----------
            var chartLabels = new List<string>();
            var tumorTotals = new List<double>();
            var normalTotals = new List<double>();
            var diffTotals = new List<double>();

            foreach (var chr in perChr.Keys.OrderBy(NaturalChrOrderKey))
            {
                long nTotal = 0, tTotal = 0;
                foreach (var kv in perChr[chr])
                {
                    nTotal += kv.Value.NormalCount;
                    tTotal += kv.Value.TumorCount;
                }
                chartLabels.Add(chr);
                tumorTotals.Add(tTotal);
                normalTotals.Add(nTotal);
                diffTotals.Add(tTotal - nTotal);
            }

            var seriesTumorNormal = new List<(string Name, List<double> Values, OxyPlot.OxyColor Color)>
            {
                ("Tumor",  tumorTotals,  OxyPlot.OxyColors.IndianRed),
                ("Normal", normalTotals, OxyPlot.OxyColors.CornflowerBlue)
            };

            ChartPrinter.PrintGroupedHorizontalBarChart(
                seriesTumorNormal,
                chartLabels,
                $"Stop codons per chromosome – Tumour vs Normal ({pairName}, {featureLabel})",
                "Number of stop codons",
                "Chromosome",
                Path.Combine(outDir, $"Chart_{pairName}_StopCodon_TumorVsNormal.pdf")
            );

            ChartPrinter.PrintHorizontalBarChart(
                diffTotals,
                chartLabels,
                $"Stop codons per chromosome – Difference (Tumour − Normal) ({pairName}, {featureLabel})",
                "Difference (number of stop codons)",
                "Chromosome",
                Path.Combine(outDir, $"Chart_{pairName}_StopCodon_Difference.pdf")
            );

            Console.WriteLine($"Charts written to {outDir}");
        }

        // ───────── Sequence loading ─────────

        private static string LoadSequence(Dictionary<string, string> chrFiles, string chr,
            Dictionary<string, List<EnsembleGeneEntry>> geneIndex, FeatureType featureType)
        {
            if (!chrFiles.TryGetValue(chr, out var filePath)) return "";

            if (featureType == FeatureType.WholeGenome || geneIndex == null)
            {
                return CleanSequence(File.ReadAllText(filePath));
            }

            // Feature-filtered: use WindowHelper with a cloned gene index
            var clonedIndex = CloneGeneIndex(geneIndex);
            return CleanSequence(WindowHelper.ReadSequenceFromFile(filePath, Path.GetFileName(filePath), ref clonedIndex, featureType));
        }

        private static Dictionary<string, List<EnsembleGeneEntry>> CloneGeneIndex(
            Dictionary<string, List<EnsembleGeneEntry>> src)
        {
            var clone = new Dictionary<string, List<EnsembleGeneEntry>>(src.Count);
            foreach (var kv in src)
                clone[kv.Key] = new List<EnsembleGeneEntry>(kv.Value);
            return clone;
        }

        // ───────── helpers ─────────

        private static Dictionary<string, string> GetChrFiles(string species)
        {
            string dir = PipelineConfig.Default.StreamDir(species);
            var map = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            if (!Directory.Exists(dir)) return map;

            foreach (var path in Directory.GetFiles(dir))
            {
                var name = Path.GetFileName(path);
                var chr = name.Split('_').Last().Split('.').First();
                chr = NormalizeChr(chr);
                map[chr] = path;
            }
            return map;
        }

        private static string NormalizeChr(string chr)
        {
            if (chr.StartsWith("chr", StringComparison.OrdinalIgnoreCase))
                chr = chr.Substring(3);
            return chr.ToUpperInvariant();
        }

        private static void CountStopsIntoBins(string seq, Dictionary<int, BinStats> bins, bool isNormal)
        {
            int n = seq.Length;
            if (n == 0) return;

            int lastBin = Math.Max(0, (n - 1) / BinSize);
            for (int i = 0; i <= lastBin; i++)
                if (!bins.ContainsKey(i)) bins[i] = new BinStats();

            for (int start = 0; start < n; start += BinSize)
            {
                int end = Math.Min(n, start + BinSize);
                int idx = start / BinSize;
                int len = end - start;
                if (isNormal) bins[idx].NormalLen += len;
                else bins[idx].TumorLen += len;
            }

            int maxI = n - 3;
            for (int i = 0; i <= maxI; i++)
            {
                int idx = i / BinSize;
                var b = bins[idx];

                if (isNormal) b.NormalTrials++; else b.TumorTrials++;

                char a = seq[i], d1 = seq[i + 1], d2 = seq[i + 2];
                if (a == 'N' || d1 == 'N' || d2 == 'N') continue;

                bool isStop =
                    (a == 'T' && d1 == 'A' && (d2 == 'A' || d2 == 'G')) ||
                    (a == 'T' && d1 == 'G' && d2 == 'A');
                if (isStop)
                {
                    if (isNormal) b.NormalCount++; else b.TumorCount++;
                }
            }
        }

        private static void CountStopsIntoBinsAligned(string normalSeq, string tumorSeq, AlignmentOffsetMap offsetMap, Dictionary<int, BinStats> bins)
        {
            int normalLen = normalSeq.Length;
            int tumorLen = tumorSeq.Length;
            int lastBin = Math.Max(0, (normalLen - 1) / BinSize);

            for (int binIdx = 0; binIdx <= lastBin; binIdx++)
            {
                if (!bins.ContainsKey(binIdx)) bins[binIdx] = new BinStats();

                int normalStart = binIdx * BinSize;
                int normalEnd = Math.Min(normalLen, normalStart + BinSize);

                var (tStart, tLen) = offsetMap.TranslateRange(normalStart, normalEnd - normalStart);

                if (tStart < 0) { tLen += tStart; tStart = 0; }
                if (tStart + tLen > tumorLen) tLen = tumorLen - tStart;
                if (tLen <= 0) continue;

                bins[binIdx].TumorLen += tLen;

                int maxI = tStart + tLen - 3;
                for (int i = tStart; i <= maxI; i++)
                {
                    bins[binIdx].TumorTrials++;

                    char a = tumorSeq[i], d1 = tumorSeq[i + 1], d2 = tumorSeq[i + 2];
                    if (a == 'N' || d1 == 'N' || d2 == 'N') continue;

                    bool isStop =
                        (a == 'T' && d1 == 'A' && (d2 == 'A' || d2 == 'G')) ||
                        (a == 'T' && d1 == 'G' && d2 == 'A');
                    if (isStop) bins[binIdx].TumorCount++;
                }
            }
        }

        private static string CleanSequence(string s)
        {
            var sb = new System.Text.StringBuilder(s.Length);
            foreach (var ch in s)
            {
                if (char.IsWhiteSpace(ch)) continue;
                char u = char.ToUpperInvariant(ch);
                switch (u)
                {
                    case 'A': case 'C': case 'G': case 'T': case 'N': sb.Append(u); break;
                    default: sb.Append('N'); break;
                }
            }
            return sb.ToString();
        }

        private static int NaturalChrOrderKey(string chr)
        {
            if (int.TryParse(chr, out int num)) return num;
            if (chr == "X") return 23;
            if (chr == "Y") return 24;
            if (chr == "MT" || chr == "M") return 25;
            return int.MaxValue;
        }

        private static (double p, double z) TwoPropZTest(long x1, long n1, long x2, long n2)
        {
            if (n1 == 0 || n2 == 0) return (1.0, 0.0);

            double p1 = (double)x1 / n1;
            double p2 = (double)x2 / n2;
            double pPool = (double)(x1 + x2) / (n1 + n2);

            double denom = Math.Sqrt(pPool * (1.0 - pPool) * (1.0 / n1 + 1.0 / n2));
            if (denom == 0.0) return (1.0, 0.0);

            double z = (p2 - p1) / denom;
            double p = 2.0 * (1.0 - CdfStdNormal(Math.Abs(z)));
            return (p, z);
        }

        private static double CdfStdNormal(double x)
        {
            double t = 1.0 / (1.0 + 0.5 * Math.Abs(x));
            double tau = t * Math.Exp(-x * x - 1.26551223 +
                       1.00002368 * t +
                       0.37409196 * Math.Pow(t, 2) +
                       0.09678418 * Math.Pow(t, 3) -
                       0.18628806 * Math.Pow(t, 4) +
                       0.27886807 * Math.Pow(t, 5) -
                       1.13520398 * Math.Pow(t, 6) +
                       1.48851587 * Math.Pow(t, 7) -
                       0.82215223 * Math.Pow(t, 8) +
                       0.17087277 * Math.Pow(t, 9));
            return (x >= 0) ? 1.0 - tau : tau;
        }

        private static void BenjaminiHochbergInPlace(List<BinTestRow> rows)
        {
            int m = rows.Count;
            var sorted = rows.Select((r, idx) => new { r, idx }).OrderBy(a => a.r.P).ToList();

            double minQ = double.PositiveInfinity;
            for (int i = m; i >= 1; i--)
            {
                double q = sorted[i - 1].r.P * m / i;
                if (q < minQ) minQ = q;
                sorted[i - 1].r.Q = Math.Min(minQ, 1.0);
            }
        }

        private static void EnsureDirectory(string path)
        {
            if (!Directory.Exists(path)) Directory.CreateDirectory(path);
        }

        private class BinStats
        {
            public int NormalCount;
            public int TumorCount;
            public int NormalLen;
            public int TumorLen;
            public long NormalTrials;
            public long TumorTrials;
            public int MaxLen() => Math.Max(NormalLen, TumorLen);
        }

        private class BinTestRow
        {
            public string Chr;
            public int BinIndex;
            public int BinStart;
            public int BinEnd;
            public long NormalCount;
            public long TumorCount;
            public long NormalTrials;
            public long TumorTrials;
            public double NormalRate100k;
            public double TumorRate100k;
            public double Delta100k;
            public double RateRatio;
            public double Z;
            public double P;
            public double Q;
        }
    }
}
