using System.Collections.Generic;

namespace GammaTextScope
{
    public static class DistributionHelper
    {
        /// <summary>
        /// Histogram of the per-word DFA exponents: bins of width <paramref name="intervals"/>
        /// from −1.25 upwards to 2.25, keyed by the bin centre (step + intervals/2). The
        /// classifier reads the result with the literal key 1.0 and compares the rounded peak
        /// key with 0.6 and 0.7; ClassificationTests guards that the repeated += 0.1 keeps
        /// producing those keys.
        /// </summary>
        public static Dictionary<double, double> CountNumberGroupings(List<double> sortedInputList, double intervals)
        {
            List<double> steps = new List<double>();
            Dictionary<double, double> stepCounts = new Dictionary<double, double>();
            for (double i = -1.25; i <= 2.25; i += intervals)
            {
                steps.Add(i);
            }
            foreach (double step in steps)
            {
                int counter = 0;
                foreach (double d in sortedInputList)
                {
                    if (d >= step && d < step + intervals)
                    {
                        counter++;
                    }
                }
                stepCounts.Add(step + (intervals / 2), counter);
            }
            return stepCounts;
        }
    }
}
