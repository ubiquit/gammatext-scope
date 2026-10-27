using System.Collections.Generic;
using System.IO;
using System.Linq;
using GammaTextScope;

namespace GammaTextScope.Tests
{
    public class SequenceMaskingTests
    {
        private static Dictionary<string, List<EnsembleGeneEntry>> Genes(params (int Start, int End)[] genes)
        {
            var index = new Dictionary<string, List<EnsembleGeneEntry>>();
            for (int i = 0; i < genes.Length; i++)
                index["g" + i] = new List<EnsembleGeneEntry> { new EnsembleGeneEntry { Chromosome = "chr1", StartPosition = genes[i].Start, EndPosition = genes[i].End } };
            return index;
        }

        private static string Masked(FeatureType mode, string line1, string line2, params (int Start, int End)[] genes)
        {
            string dir = Path.Combine(Path.GetTempPath(), "gts-tests", "mask");
            Directory.CreateDirectory(dir);
            string path = Path.Combine(dir, "sample_chr1.txt");
            File.WriteAllText(path, line1 + "\n" + line2 + "\n");
            return WindowHelper.ReadSequenceFromFile(path, "sample_chr1.txt", Genes(genes), mode);
        }

        // Legacy gene-mask conventions pinned by expected/ (see WindowHelper): a gene (S, E) covers the
        // 0-based positions S-2..E-2, and the first gene of a chromosome is skipped.
        [Fact]
        public void GeneMaskKeepsTheLegacyOffsetAndSkipsTheFirstGene()
        {
            // positions: 0123456789 0123456789
            //            ACGTACGTAC GTACGTACGT
            var genes = new[] { (1, 3), (4, 6), (10, 12) };   // (1,3) skipped; (4,6) -> 2..4 "GTA"; (10,12) -> 8..10 "ACG"
            Assert.Equal("GTAACG", Masked(FeatureType.Gene, "ACGTACGTAC", "GTACGTACGT", genes));
            Assert.Equal("ACCGTTACGTACGT", Masked(FeatureType.Intergenic, "ACGTACGTAC", "GTACGTACGT", genes));
            Assert.Equal("ACGTACGTACGTACGTACGT", Masked(FeatureType.WholeGenome, "ACGTACGTAC", "GTACGTACGT", genes));
        }

        [Fact]
        public void TokenizerKeepsEmptyTokensBetweenAdjacentStopCodons()
        {
            var tokens = Tokenizer.TokenizeOnStopCodons("ACGTAATAGCCC", GenomeConstants.stopWords, removeEmptyEntries: false);
            Assert.Equal(new[] { "ACG", "", "CCC" }, tokens);
            Assert.Equal(12, SequenceHelper.GetTokenizedSequenceLength(tokens, isEnd: true));
        }
    }

    public class ClassificationTests
    {
        [Fact]
        public void AlphaHistogramHasExactBinKeysUsedByTheClassifier()
        {
            // ClassifyDistribution indexes the histogram with the literal keys 1.0, 0.6 and 0.7;
            // the bins are built by repeated += 0.1, so this guards the floating-point layout.
            var bins = DistributionHelper.CountNumberGroupings(new List<double> { 0.62, 0.68, 1.02 }, 0.1);
            Assert.Equal(35, bins.Count);
            Assert.True(bins.ContainsKey(1.0));
            Assert.Equal(1, bins[1.0]);
            Assert.Equal(2, bins.Single(kv => System.Math.Round(kv.Key, 2) == 0.6).Value + bins.Single(kv => System.Math.Round(kv.Key, 2) == 0.7).Value);
        }

        [Fact]
        public void ZipfMidpointTestAcceptsAZipfianRankFrequencyList()
        {
            var counts = new Dictionary<string, int>();
            for (int rank = 1; rank <= 400; rank++) counts["w" + rank] = (int)System.Math.Round(4000.0 / rank);
            Assert.True(ZipfHelper.ClassifyDistributionDeltaCenter(counts).Item1);
        }
    }

    public class DfaRegressionTests
    {
        private static string ExpectedDir =>
            Path.GetFullPath(Path.Combine(System.AppContext.BaseDirectory, "..", "..", "..", "..", "..", "expected", "v1"));

        [Theory]
        [InlineData("Tumor/BreastCancer/DecayGapDFA_Intergenic/GapSeries_NormalSide_Zipf.tsv", 0.4487)]
        [InlineData("Tumor/BreastCancer/DecayGapDFA_Intergenic/GapSeries_TumorSide_Zipf.tsv", 0.5131)]
        [InlineData("Tumor/Melanoma/DecayGapDFA_Intergenic/GapSeries_NormalSide_Zipf.tsv", 0.4787)]
        [InlineData("Tumor/Melanoma/DecayGapDFA_Intergenic/GapSeries_TumorSide_Zipf.tsv", 0.4982)]
        public void DfaAlphaReproducesPublishedGapSeriesExponents(string relativePath, double expectedAlpha)
        {
            var gaps = File.ReadLines(Path.Combine(ExpectedDir, relativePath)).Skip(1)
                .Select(l => double.Parse(l.Split('\t')[1], System.Globalization.CultureInfo.InvariantCulture)).ToList();
            double alpha = WordFluctuations.DfaAlpha(gaps, 30);
            Assert.Equal(expectedAlpha, System.Math.Round(alpha, 4));
        }

        [Fact]
        public void DfaReturnsSentinelBelowThreshold()
        {
            Assert.Equal(-1000, WordFluctuations.DfaAlpha(new List<double> { 1, 2, 3 }, 30));
        }
    }
}
