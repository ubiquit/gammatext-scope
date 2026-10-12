using System.Collections.Generic;

namespace GammaTextScope
{
    public static class GenomeConstants
    {
        /// <summary>The three canonical stop codons act as γword delimiters.</summary>
        public static readonly List<string> stopWords = new List<string> { "TAG", "TAA", "TGA" };

        /// <summary>Minimum number of occurrences a word needs before its occurrence-gap series is submitted to DFA.</summary>
        public const int AlphaWordThreshold = 20;

        public const int SlidingWindowKWordThreshold = 50;
        public const int WindowSize = 5000000;
    }

    public enum FeatureType
    {
        Intergenic,
        Gene,
        Exon,
        WholeGenome
    }
}
