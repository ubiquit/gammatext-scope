using System;
using System.Collections.Generic;
using System.IO;

namespace GammaTextScope
{
    /// <summary>
    /// Locations of inputs and outputs. The legacy pipeline hard-coded
    /// D:\Files\Stream (inputs) and D:\Files\Harvests (outputs); those remain the
    /// defaults so existing data layouts keep working, but both can be overridden
    /// with --data-root / --out-root or the GAMMATEXT_DATA / GAMMATEXT_OUT variables.
    /// </summary>
    public sealed class PipelineConfig
    {
        public string DataRoot { get; init; } = Environment.GetEnvironmentVariable("GAMMATEXT_DATA") ?? @"D:\Files\Stream";
        public string OutRoot  { get; init; } = Environment.GetEnvironmentVariable("GAMMATEXT_OUT")  ?? @"D:\Files\Harvests";

        /// <summary>Optional chromosome filter ("21", "X", ...; without the chr prefix). Null = all.</summary>
        public HashSet<string> Chromosomes { get; init; }

        public int WordsPerWindow { get; init; } = 7000;
        public int WindowWidth { get; init; } = 1;

        public static PipelineConfig Default { get; set; } = new PipelineConfig();

        /// <summary>Folder holding one text file per chromosome for a sample, e.g. &lt;DataRoot&gt;/HomoSapiens_BreastCancer_Normal.</summary>
        public string StreamDir(string sample) => Path.Combine(DataRoot, sample);

        /// <summary>Folder with the per-window classification output of <c>classify</c> for a sample and feature type.</summary>
        public string ClassificationDir(string sample, FeatureType featureType) =>
            Path.Combine(OutRoot, sample, $"Experiment_{sample}_{featureType}_{WordsPerWindow}");

        public string TumorOutDir(string pairName, string experiment, string featureLabel) =>
            Path.Combine(OutRoot, "Tumor", pairName, $"{experiment}_{featureLabel}");

        public bool IncludesChromosome(string chr)
        {
            if (Chromosomes == null || Chromosomes.Count == 0) return true;
            if (chr.StartsWith("chr", StringComparison.OrdinalIgnoreCase)) chr = chr.Substring(3);
            return Chromosomes.Contains(chr.ToUpperInvariant());
        }

        public static HashSet<string> ParseChromosomes(string list)
        {
            if (string.IsNullOrWhiteSpace(list)) return null;
            var set = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
            foreach (var raw in list.Split(',', StringSplitOptions.RemoveEmptyEntries))
            {
                var c = raw.Trim();
                if (c.StartsWith("chr", StringComparison.OrdinalIgnoreCase)) c = c.Substring(3);
                set.Add(c.ToUpperInvariant());
            }
            return set;
        }
    }
}
