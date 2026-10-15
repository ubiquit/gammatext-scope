using System.Collections.Generic;
using System.IO;
using System.Linq;
using GammaTextScope;

namespace GammaTextScope.Tests
{
    /// <summary>
    /// Re-runs the pipeline and compares its output byte-for-byte with expected/ for
    /// both reference datasets. decay-dfa only needs the tables in expected/ and always
    /// runs; the sequence based steps run when the dataset's samples are present under
    /// the data root (--data-root / GAMMATEXT_DATA).
    /// </summary>
    public class PipelineRegressionTests
    {
        public static IEnumerable<object[]> DfaCases()
        {
            foreach (var v in new[] { "v1", "v2" })
                foreach (var cancer in new[] { 0, 1 })
                    foreach (var mode in new[] { "Intergenic", "Gene", "WholeGenome", "WholeGenomeMerged" })
                        if (!(v == "v1" && mode == "WholeGenome"))   // v1 never had annotation-free whole-genome classifications
                            yield return new object[] { v, cancer, mode };
        }

        [Theory]
        [MemberData(nameof(DfaCases))]
        public void DecayDfaReproducesSummaryAndGapSeries(string version, int cancer, string mode)
        {
            var ds = Dataset.Get(version); string pair = ds.AlphaPair(cancer);
            string outRoot = TestData.StageOutRoot(ds, $"dfa-{pair}-{mode}", new string[0]);
            TestData.UseOutRoot(outRoot);
            string scaRel = $"Tumor/{pair}/StopCodonAlphaComparison_{mode}";
            Directory.CreateDirectory(Path.Combine(outRoot, scaRel));
            File.Copy(ds.Expected($"{scaRel}/StopCodonAlpha_{pair}_AllChromosomes.tsv"), Path.Combine(outRoot, scaRel, $"StopCodonAlpha_{pair}_AllChromosomes.tsv"));

            new DecayGapDfa().Run(pair, mode);

            string dfaRel = $"Tumor/{pair}/DecayGapDFA_{mode}";
            foreach (var expectedFile in Directory.GetFiles(ds.Expected(dfaRel)))
            {
                string produced = Path.Combine(outRoot, dfaRel, Path.GetFileName(expectedFile));
                Assert.True(File.Exists(produced), $"missing {produced}");
                Assert.Equal(File.ReadAllText(expectedFile), File.ReadAllText(produced));
            }
        }

        public static IEnumerable<object[]> DensityCases()
        {
            foreach (var v in new[] { "v1", "v2" })
                foreach (var mode in new[] { "Intergenic", "Gene", "WholeGenome" })
                    yield return new object[] { v, 0, mode };
            yield return new object[] { "v1", 1, "Intergenic" };
            yield return new object[] { "v2", 1, "Intergenic" };
            yield return new object[] { "v2", 1, "WholeGenome" };
        }

        [SkippableTheory]
        [MemberData(nameof(DensityCases))]
        public void DensityReproducesChromosome21(string version, int cancer, string mode)
        {
            var ds = Dataset.Get(version);
            Skip.IfNot(ds.SequencesAvailable, $"{version} sequence data not available under {TestData.DataRoot}");
            var samples = ds.Samples(cancer); string pair = ds.DensityPair(cancer);
            string outRoot = TestData.StageOutRoot(ds, $"density-{pair}-{mode}", new string[0]);
            TestData.UseOutRoot(outRoot, "21");

            new StopCodonDensity().Run(samples[0], samples[1], pair, System.Enum.Parse<FeatureType>(mode));

            string rel = $"Tumor/{pair}/StopCodonDensity_{mode}/{pair}_StopCodonDensity_100k_perChromosome.tsv";
            var produced = TestData.Lines(Path.Combine(outRoot, rel));
            var expected = TestData.Lines(ds.Expected(rel));
            Assert.Equal(expected[0], produced[0]);
            Assert.Equal(TestData.ChromosomeRows(ds.Expected(rel), "21"), produced.Skip(1).ToList());
        }

        public static IEnumerable<object[]> AlphaCases()
        {
            foreach (var v in new[] { "v1", "v2" })
                foreach (var cancer in new[] { 0, 1 })
                    foreach (var mode in new[] { "Intergenic", "Gene", "WholeGenomeMerged", "WholeGenome" })
                        if (!(v == "v1" && mode == "WholeGenome"))
                            yield return new object[] { v, cancer, mode };
        }

        [SkippableTheory]
        [MemberData(nameof(AlphaCases))]
        public void AlphaComparisonReproducesChromosome21(string version, int cancer, string mode)
        {
            var ds = Dataset.Get(version);
            Skip.IfNot(ds.SequencesAvailable, $"{version} sequence data not available under {TestData.DataRoot}");
            var samples = ds.Samples(cancer); string pair = ds.AlphaPair(cancer);
            bool merged = mode.EndsWith("Merged");
            var ft = System.Enum.Parse<FeatureType>(mode.Replace("Merged", ""));
            var staged = version == "v1" ? new[] { FeatureType.Intergenic, FeatureType.Gene } : new[] { FeatureType.Intergenic, FeatureType.Gene, FeatureType.WholeGenome };
            string outRoot = TestData.StageOutRoot(ds, $"alpha-{pair}-{mode}", samples, staged);
            TestData.UseOutRoot(outRoot, "21");

            new StopCodonAlphaComparison().Run(samples[0], samples[1], pair, ft, merged);

            string rel = $"Tumor/{pair}/StopCodonAlphaComparison_{mode}/StopCodonAlpha_{pair}_chr21.tsv";
            Assert.Equal(File.ReadAllText(ds.Expected(rel)), File.ReadAllText(Path.Combine(outRoot, rel)));
            string all = $"Tumor/{pair}/StopCodonAlphaComparison_{mode}/StopCodonAlpha_{pair}_AllChromosomes.tsv";
            Assert.Equal(TestData.ChromosomeRows(ds.Expected(all), "chr21"), TestData.Lines(Path.Combine(outRoot, all)).Skip(1).ToList());
        }

        public static IEnumerable<object[]> ClassifyCases()
        {
            yield return new object[] { "v1", "HomoSapiens_BreastCancer_Normal", "Intergenic" };
            yield return new object[] { "v1", "HomoSapiens_BreastCancer_Tumor", "Gene" };
            yield return new object[] { "v2", "HomoSapiens_BreastCancer_Normal_v2", "Intergenic" };
            yield return new object[] { "v2", "HomoSapiens_BreastCancer_Tumor_v2", "Gene" };
            yield return new object[] { "v2", "HomoSapiens_Melanoma_Normal_v3", "Intergenic" };
            yield return new object[] { "v2", "HomoSapiens_Melanoma_Tumor_v3", "WholeGenome" };
        }

        [SkippableTheory]
        [MemberData(nameof(ClassifyCases))]
        public void ClassifyReproducesChromosome21(string version, string sample, string mode)
        {
            var ds = Dataset.Get(version);
            Skip.IfNot(ds.SequencesAvailable, $"{version} sequence data not available under {TestData.DataRoot}");
            var ft = System.Enum.Parse<FeatureType>(mode);
            string outRoot = TestData.StageOutRoot(ds, $"classify-{sample}-{mode}", new string[0]);
            var cfg = TestData.UseOutRoot(outRoot, "21");

            new WindowClassifier().RunSimplified(sample, cfg.WindowWidth, cfg.WordsPerWindow, ft, cfg.ClassificationDir(sample, ft));

            string rel = $"{sample}/Experiment_{sample}_{mode}_7000/WindowWidth1/chr21WW1_window_positions_indexLength.txt";
            // columns 1–5 (window, start, length, status, isAlpha) are what the experiments consume
            List<string> Key(string path) => TestData.Lines(path).Select(l => string.Join('\t', l.Split('\t').Take(5))).ToList();
            Assert.Equal(Key(ds.Expected(rel)), Key(Path.Combine(outRoot, rel)));
        }
    }
}
