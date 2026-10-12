using System;
using System.Collections.Generic;
using System.Linq;

namespace GammaTextScope
{
    public static class Tokenizer
    {
        /// <summary>
        /// Splits a nucleotide string on the stop codons. With removeEmptyEntries = false,
        /// adjacent stop codons yield empty tokens; the window-length bookkeeping in
        /// SequenceHelper.GetTokenizedSequenceLength relies on those empties.
        /// </summary>
        public static List<string> TokenizeOnStopCodons(string text, List<string> stopWords, bool removeEmptyEntries = true)
        {
            return removeEmptyEntries
                ? text.Split(stopWords.ToArray(), StringSplitOptions.RemoveEmptyEntries).ToList()
                : text.Split(stopWords.ToArray(), StringSplitOptions.None).ToList();
        }
    }
}
