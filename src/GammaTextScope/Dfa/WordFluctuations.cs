using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace GammaTextScope
{
    /// <summary>
    /// Per-word DFA exponents for one window of γtext: the occurrence-gap series of
    /// every word and the DFA α of each series that has at least
    /// <see cref="GenomeConstants.AlphaWordThreshold"/> occurrences. The gap bookkeeping
    /// (the first gap is the position of the first occurrence, counted from 1) follows
    /// the word index and is pinned by expected/.
    /// </summary>
    public static class WordFluctuations
    {
        public static Dictionary<string, double> Compute(string[] tokenizedText, int threshold = GenomeConstants.AlphaWordThreshold)
        {
            var lexicon = new HashSet<string>(tokenizedText);
            var flucts = new ConcurrentDictionary<string, double>();

            Parallel.ForEach(lexicon, currentWord =>
            {
                int distCounter = 0;
                var distList = new List<int>();
                foreach (string word in tokenizedText)
                {
                    if (currentWord.Equals(word, StringComparison.OrdinalIgnoreCase))
                    {
                        distList.Add(distCounter);
                        distCounter = 0;
                    }
                    distCounter++;
                }

                if (distList.Count >= threshold)
                {
                    double alpha = new Dfa().CalculateDFABoxController(distList.Select(d => (double)d).ToList());
                    if (!double.IsNaN(alpha))
                        flucts.TryAdd(currentWord, alpha);
                }
            });

            return flucts.ToDictionary(kvp => kvp.Key, kvp => kvp.Value);
        }

        /// <summary>DFA α of an arbitrary series; -1000 when the series is shorter than the threshold (sentinel).</summary>
        public static double DfaAlpha(List<double> series, int threshold)
        {
            if (series.Count < threshold) return -1000;
            return new Dfa().CalculateDFABoxController(series);
        }
    }
}
