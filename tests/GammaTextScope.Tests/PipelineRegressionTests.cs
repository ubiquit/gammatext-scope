using System.Collections.Generic;
using System.IO;
using System.Linq;
using GammaTextScope;

namespace GammaTextScope.Tests
{
    /// <summary>
    /// Re-runs the pipeline and compares its output byte-for-byte with expected/.
    /// decay-dfa only needs the tables in expected/ and always runs; the sequence
    /// based steps run when the four samples are present under the data root.
    /// </summary>
    public class PipelineRegressionTests
    {
        private static readonly string[] BreastPair = { "HomoSapiens_BreastCancer_Normal", "HomoSapiens_BreastCancer_Tumor" };
        private static readonly string[] MelanomaPair = { "HomoSapiens_Melanoma_Normal", "HomoSapiens_Melanoma_Tumor" };

        [Theory]
        [InlineData("BreastCancer", "Intergenic")]
        [InlineData("BreastCancer", "Gene")]
        [InlineData("BreastCancer", "WholeGenome")]
        [InlineData("Melanoma", "Intergenic")]
        [InlineData("Melanoma", "Gene")]
        [InlineData("Melanoma", "WholeGenome")]
        public void DecayDfaReproducesSummaryAndGapSeries(string pair, string mode)
        {
            string outRoot = TestData.StageOutRoot($"dfa-{pair}-{mode}", new string[0]);
            TestData.UseOutRoot(outRoot);
            string scaRel = $"Tumor/{pair}/StopCodonAlphaComparison_{mode}";
            Directory.CreateDirectory(Path.Combine(outRoot, scaRel));
            File.Copy(TestData.Expected($"{scaRel}/StopCodonAlpha_{pair}_AllChromosomes.tsv"),
                      Path.Combine(outRoot, scaRel, $"StopCodonAlpha_{pair}_AllChromosomes.tsv"));

            new DecayGapDfa().Run(pair, mode);

            string dfaRel = $"Tumor/{pair}/DecayGapDFA_{mode}";
            foreach (var expectedFile in Directory.GetFiles(TestData.Expected(dfaRel)))
            {
                string produced = Path.Combine(outRoot, dfaRel, Path.GetFileName(expectedFile));
                Assert.True(File.Exists(produced), $"missing {produced}");
                Assert.Equal(File.ReadAllText(expectedFile), File.ReadAllText(produced));
            }
        }

        [SkippableTheory]
        [InlineData("BreastCancer_Normal_vs_Tumor", 0, "Intergenic")]
        [InlineData("BreastCancer_Normal_vs_Tumor", 0, "Gene")]
        [InlineData("BreastCancer_Normal_vs_Tumor", 0, "WholeGenome")]
        [InlineData("Melanoma_Normal_vs_Tumor", 1, "Intergenic")]
        public void DensityReproducesChromosome21(string pair, int which, string mode)
        {
            Skip.IfNot(TestData.SequencesAvailable, "sequence data not available under " + TestData.DataRoot);
            var samples = which == 0 ? BreastPair : MelanomaPair;
            var ft = System.Enum.Parse<FeatureType>(mode);
            string outRoot = TestData.StageOutRoot($"density-{pair}-{mode}", new string[0]);
            TestData.UseOutRoot(outRoot, "21");

            new StopCodonDensity().Run(samples[0], samples[1], pair, ft);

            string rel = $"Tumor/{pair}/StopCodonDensity_{mode}/{pair}_StopCodonDensity_100k_perChromosome.tsv";
            var produced = TestData.Lines(Path.Combine(outRoot, rel));
            var expected = TestData.Lines(TestData.Expected(rel));
            Assert.Equal(expected[0], produced[0]);
            Assert.Equal(TestData.ChromosomeRows(TestData.Expected(rel), "21"), produced.Skip(1).ToList());
        }

        [SkippableTheory]
        [InlineData("BreastCancer", 0, "Intergenic")]
        [InlineData("BreastCancer", 0, "Gene")]
        [InlineData("BreastCancer", 0, "WholeGenome")]
        [InlineData("Melanoma", 1, "Intergenic")]
        [InlineData("Melanoma", 1, "Gene")]
        [InlineData("Melanoma", 1, "WholeGenome")]
        public void AlphaComparisonReproducesChromosome21(string pair, int which, string mode)
        {
            Skip.IfNot(TestData.SequencesAvailable, "sequence data not available under " + TestData.DataRoot);
            var samples = which == 0 ? BreastPair : MelanomaPair;
            var ft = System.Enum.Parse<FeatureType>(mode);
            string outRoot = TestData.StageOutRoot($"alpha-{pair}-{mode}", samples, FeatureType.Intergenic, FeatureType.Gene);
            TestData.UseOutRoot(outRoot, "21");

            new StopCodonAlphaComparison().Run(samples[0], samples[1], pair, ft);

            string rel = $"Tumor/{pair}/StopCodonAlphaComparison_{mode}/StopCodonAlpha_{pair}_chr21.tsv";
            Assert.Equal(File.ReadAllText(TestData.Expected(rel)), File.ReadAllText(Path.Combine(outRoot, rel)));
            string all = $"Tumor/{pair}/StopCodonAlphaComparison_{mode}/StopCodonAlpha_{pair}_AllChromosomes.tsv";
            Assert.Equal(TestData.ChromosomeRows(TestData.Expected(all), "chr21"), TestData.Lines(Path.Combine(outRoot, all)).Skip(1).ToList());
        }

        [SkippableTheory]
        [InlineData("HomoSapiens_BreastCancer_Normal")]
        [InlineData("HomoSapiens_BreastCancer_Tumor")]
        [Trait("Category", "Slow")]
        public void ClassifyReproducesChromosome21Intergenic(string sample)
        {
            Skip.IfNot(TestData.SequencesAvailable, "sequence data not available under " + TestData.DataRoot);
            string outRoot = TestData.StageOutRoot($"classify-{sample}", new string[0]);
            var cfg = TestData.UseOutRoot(outRoot, "21");

            new WindowClassifier().RunSimplified(sample, cfg.WindowWidth, cfg.WordsPerWindow, FeatureType.Intergenic,
                cfg.ClassificationDir(sample, FeatureType.Intergenic));

            string rel = $"{sample}/Experiment_{sample}_Intergenic_7000/WindowWidth1/chr21WW1_window_positions_indexLength.txt";
            // columns 1–5 (window, start, length, status, isAlpha) are what the experiments consume;
            // the legacy writer concatenated the two trailing diagnostic columns
            List<string> Key(string path) => TestData.Lines(path).Select(l => string.Join('\t', l.Split('\t').Take(5))).ToList();
            Assert.Equal(Key(TestData.Expected(rel)), Key(Path.Combine(outRoot, rel)));
        }
    }
}
