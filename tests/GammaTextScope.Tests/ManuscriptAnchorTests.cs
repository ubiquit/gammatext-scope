using System.Linq;
using GammaTextScope;

namespace GammaTextScope.Tests
{
    /// <summary>
    /// The numbers printed in Glattre, Glattre &amp; Moe must be exactly the ones in
    /// expected/. These tests pin the reference outputs to the manuscript; the
    /// pipeline tests pin the code to the reference outputs.
    /// </summary>
    public class ManuscriptAnchorTests
    {
        [Theory]
        // Figure 1 / Results: per-chromosome stop-codon difference (Tumour − Normal), intergenic mode
        [InlineData("BreastCancer_Normal_vs_Tumor", "5", -425)]
        [InlineData("BreastCancer_Normal_vs_Tumor", "10", 528)]
        [InlineData("BreastCancer_Normal_vs_Tumor", "12", 471)]
        public void Figure1PerChromosomeDifferences(string pair, string chr, long expectedDelta)
        {
            var rows = TestData.ChromosomeRows(TestData.Expected($"Tumor/{pair}/StopCodonDensity_Intergenic/{pair}_StopCodonDensity_100k_perChromosome.tsv"), chr);
            long delta = rows.Sum(r => long.Parse(r.Split('\t')[4]) - long.Parse(r.Split('\t')[3]));
            Assert.Equal(expectedDelta, delta);
        }

        [Theory]
        [InlineData("BreastCancer_Normal_vs_Tumor", -836)]
        [InlineData("Melanoma_Normal_vs_Tumor", 1099)]
        public void GenomeWideStopCodonDifferenceIntergenic(string pair, long expectedTotal)
        {
            var rows = TestData.Lines(TestData.Expected($"Tumor/{pair}/StopCodonDensity_Intergenic/{pair}_StopCodonDensity_100k_perChromosome.tsv")).Skip(1);
            long total = rows.Sum(r => long.Parse(r.Split('\t')[4]) - long.Parse(r.Split('\t')[3]));
            Assert.Equal(expectedTotal, total);
        }

        [Theory]
        // Table 2: npc-DNA (intergenic), Zipf-only decay-marker gap series
        [InlineData("BreastCancer", "Intergenic", "NormalSide_Zipf", "0.4487")]
        [InlineData("BreastCancer", "Intergenic", "TumorSide_Zipf", "0.5131")]
        [InlineData("Melanoma", "Intergenic", "NormalSide_Zipf", "0.4787")]
        [InlineData("Melanoma", "Intergenic", "TumorSide_Zipf", "0.4982")]
        // Table 3: pc-DNA (gene) — "0.4892 → 0.5690 in breast cancer … 0.5863 → 0.5250 in melanoma"
        [InlineData("BreastCancer", "Gene", "NormalSide_Zipf", "0.4892")]
        [InlineData("BreastCancer", "Gene", "TumorSide_Zipf", "0.569")]
        [InlineData("Melanoma", "Gene", "NormalSide_Zipf", "0.5863")]
        [InlineData("Melanoma", "Gene", "TumorSide_Zipf", "0.525")]
        public void DfaAlphaValuesInTables2And3(string pair, string mode, string category, string expectedAlpha)
        {
            var summary = TestData.DfaSummary(TestData.Expected($"Tumor/{pair}/DecayGapDFA_{mode}/DecayGapDFA_{pair}_Summary.txt"));
            Assert.Equal(expectedAlpha, summary[category].Alpha);
        }

        [Theory]
        // Table 3: "alfa-qualified windows in pc-DNA are scarce (Normal_Alfa: 1 in BC, 5 in Mel; Tumor_Alfa: 0 in BC, 12 in Mel)"
        [InlineData("BreastCancer", "NormalSide_Alfa", 1)]
        [InlineData("BreastCancer", "TumorSide_Alfa", 0)]
        [InlineData("Melanoma", "NormalSide_Alfa", 5)]
        [InlineData("Melanoma", "TumorSide_Alfa", 12)]
        public void Table3AlfaDecayDotCountsInGeneMode(string pair, string category, int expectedDots)
        {
            var summary = TestData.DfaSummary(TestData.Expected($"Tumor/{pair}/DecayGapDFA_Gene/DecayGapDFA_{pair}_Summary.txt"));
            Assert.Equal(expectedDots, summary[category].Dots);
        }
    }
}
