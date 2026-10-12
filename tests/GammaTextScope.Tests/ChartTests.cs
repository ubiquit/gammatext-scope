using System.Collections.Generic;
using System.IO;
using System.Linq;
using GammaTextScope;

namespace GammaTextScope.Tests
{
    /// <summary>
    /// The charts are rendered from the article's data in expected/. OxyPlot's PDF
    /// output is not byte-stable across versions, so these tests check that the
    /// renderers accept the real data and emit well-formed, non-trivial PDFs; the
    /// numbers behind the charts are covered by the pipeline and anchor tests.
    /// </summary>
    public class ChartTests
    {
        private static string Out(string name)
        {
            string dir = Path.Combine(Path.GetTempPath(), "gts-tests", "charts");
            Directory.CreateDirectory(dir);
            return Path.Combine(dir, name);
        }

        private static void AssertPdf(string path, int minBytes)
        {
            var bytes = File.ReadAllBytes(path);
            Assert.True(bytes.Length >= minBytes, $"{path}: {bytes.Length} bytes");
            Assert.Equal("%PDF-", System.Text.Encoding.ASCII.GetString(bytes, 0, 5));
            Assert.Contains("%%EOF", System.Text.Encoding.ASCII.GetString(bytes, bytes.Length - 64, 64));
        }

        [Fact]
        public void Figure1DifferenceBarChartRendersFromExpectedDensityTable()
        {
            string table = TestData.Expected("Tumor/BreastCancer_Normal_vs_Tumor/StopCodonDensity_Intergenic/BreastCancer_Normal_vs_Tumor_StopCodonDensity_100k_perChromosome.tsv");
            var perChr = TestData.Lines(table).Skip(1).Select(l => l.Split('\t'))
                .GroupBy(p => p[0]).ToDictionary(g => g.Key, g => g.Sum(p => long.Parse(p[4]) - long.Parse(p[3])));
            var labels = perChr.Keys.ToList();
            var values = labels.Select(k => (double)perChr[k]).ToList();
            Assert.Equal(-425, perChr["5"]);

            string pdf = Out("fig1_difference.pdf");
            ChartPrinter.PrintHorizontalBarChart(values, labels, "Stop codons per chromosome – Difference (Tumour − Normal)", "Difference (number of stop codons)", "Chromosome", pdf);
            AssertPdf(pdf, 4_000);
        }

        [Theory]
        [InlineData("BreastCancer")]
        [InlineData("Melanoma")]
        public void DnaScopeChartRendersFromExpectedChromosome1Table(string pair)
        {
            string table = TestData.Expected($"Tumor/{pair}/StopCodonAlphaComparison_Intergenic/StopCodonAlpha_{pair}_chr1.tsv");
            var rows = TestData.Lines(table).Skip(1).Select(l => l.Split('\t')).ToList();
            Assert.True(rows.Count > 500, "chr1 should have hundreds of windows");

            string pdf = Out($"dnascope_{pair}_chr1.pdf");
            ChartPrinter.PrintStopCodonAlphaLineChart(
                rows.Select(r => long.Parse(r[7])).ToList(),
                rows.Select(r => r[3]).ToList(),
                rows.Select(r => r[4]).ToList(),
                $"Stop-codon delta and qualification – chr1 ({pair})", "Window (sequential)", "Stop-codon delta (Tumour − Normal)",
                pdf, chromosomeBoundaries: null, subtitle: "Intergenic", centromereRanges: new List<(int, int)> { (640, 646) });
            AssertPdf(pdf, 20_000);
        }
    }
}
