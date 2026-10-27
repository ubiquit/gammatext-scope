using System.Collections.Generic;

namespace GammaTextScope
{
    public static class GenomeConstants
    {
        /// <summary>The three canonical stop codons act as γword delimiters.</summary>
        public static readonly List<string> stopWords = new List<string> { "TAG", "TAA", "TGA" };

        /// <summary>Minimum number of occurrences a word needs before its occurrence-gap series is submitted to DFA.</summary>
        public const int AlphaWordThreshold = 20;

        /// <summary>
        /// Line terminator of every file the pipeline writes. The reference outputs in
        /// expected/ were produced on Windows and are compared byte-for-byte, so the
        /// terminator is fixed to CRLF instead of following Environment.NewLine.
        /// </summary>
        public const string OutputNewLine = "\r\n";
    }

    public enum FeatureType
    {
        Intergenic,
        Gene,
        Exon,
        WholeGenome
    }
}
