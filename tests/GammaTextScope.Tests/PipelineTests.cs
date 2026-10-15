using System.Collections.Generic;
using System.IO;
using System.Linq;
using GammaTextScope;

namespace GammaTextScope.Tests
{
    public class SequenceMaskingTests
    {
        [Fact]
        public void IntergenicDistanceMergingMatchesLegacySelfTest() => Assert.True(SequenceHelper.TestCreateGeneticDistances());

        [Fact]
        public void SequenceInverterMatchesLegacySelfTest() => Assert.True(SequenceHelper.TestCreateGeneticTextNew());

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
