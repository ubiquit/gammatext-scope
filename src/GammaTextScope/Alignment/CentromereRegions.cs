using System;
using System.Collections.Generic;

namespace GammaTextScope
{
    /// <summary>
    /// GRCh38 (hg38) centromere coordinates from UCSC cytoband data.
    /// These are the approximate start and end positions of the centromeric
    /// region on each chromosome in the reference assembly.
    ///
    /// Used to annotate charts with centromere locations so that assembly
    /// artifacts in these regions can be visually identified.
    /// </summary>
    public static class CentromereRegions
    {
        // (start, end) in raw chromosome coordinates (0-based, bp)
        // Source: UCSC Table Browser → cytoBand → GRCh38, filtered for acen bands
        private static readonly Dictionary<string, (long Start, long End)> _regions =
            new Dictionary<string, (long, long)>(StringComparer.OrdinalIgnoreCase)
            {
                { "chr1",  (121_700_000, 125_100_000) },
                { "chr2",  ( 91_800_000,  96_000_000) },
                { "chr3",  ( 87_800_000,  94_000_000) },
                { "chr4",  ( 48_200_000,  51_800_000) },
                { "chr5",  ( 46_100_000,  51_400_000) },
                { "chr6",  ( 58_500_000,  62_600_000) },
                { "chr7",  ( 58_100_000,  62_100_000) },
                { "chr8",  ( 43_200_000,  47_200_000) },
                { "chr9",  ( 42_200_000,  45_500_000) },
                { "chr10", ( 38_000_000,  41_600_000) },
                { "chr11", ( 51_000_000,  55_800_000) },
                { "chr12", ( 33_200_000,  37_800_000) },
                { "chr13", (  0,          18_500_000) }, // acrocentric — short arm is satellite
                { "chr14", (  0,          18_200_000) }, // acrocentric
                { "chr15", (  0,          19_000_000) }, // acrocentric
                { "chr16", ( 35_300_000,  38_400_000) },
                { "chr17", ( 22_700_000,  27_400_000) },
                { "chr18", ( 15_400_000,  21_500_000) },
                { "chr19", ( 24_200_000,  28_100_000) },
                { "chr20", ( 26_200_000,  30_400_000) },
                { "chr21", (  0,          13_000_000) }, // acrocentric
                { "chr22", (  0,          17_400_000) }, // acrocentric
                { "chrX",  ( 58_100_000,  63_800_000) },
                { "chrY",  ( 10_300_000,  10_600_000) },
            };

        /// <summary>
        /// Get centromere region for a chromosome. Returns null if not found.
        /// Accepts formats: "chr1", "1", "CHR1", etc.
        /// </summary>
        public static (long Start, long End)? GetRegion(string chr)
        {
            string key = chr;
            if (!key.StartsWith("chr", StringComparison.OrdinalIgnoreCase))
                key = "chr" + key;

            if (_regions.TryGetValue(key, out var region))
                return region;
            return null;
        }

        /// <summary>
        /// Check if a genomic position falls within the centromere.
        /// </summary>
        public static bool IsInCentromere(string chr, long position)
        {
            var region = GetRegion(chr);
            if (region == null) return false;
            return position >= region.Value.Start && position <= region.Value.End;
        }

        /// <summary>
        /// Check if a range overlaps the centromere.
        /// </summary>
        public static bool OverlapsCentromere(string chr, long start, long length)
        {
            var region = GetRegion(chr);
            if (region == null) return false;
            long end = start + length;
            return start < region.Value.End && end > region.Value.Start;
        }

        /// <summary>
        /// Get all centromere regions.
        /// </summary>
        public static IReadOnlyDictionary<string, (long Start, long End)> All => _regions;
    }
}
