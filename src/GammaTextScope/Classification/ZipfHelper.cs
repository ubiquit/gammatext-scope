using System;
using System.Collections.Generic;
using System.Linq;

namespace GammaTextScope
{
    /// <summary>
    /// Zipf test of one window. The rank–frequency list of the γwords seen more than three
    /// times must decline monotonically on log–log axes, and the log-frequency of the point
    /// halfway along the list must lie within 5 % of the mean of the points at 30 % and 70 %.
    /// </summary>
    public static class ZipfHelper
    {
        /// <returns>(isZipfian, log–log rank–frequency list, |midpoint − interpolated| , slope between the 30 % and 70 % points)</returns>
        public static Tuple<bool, List<KeyValuePair<double, double>>, double, double> ClassifyDistributionDeltaCenter(Dictionary<string, int> distribution)
        {
            var rankFreqList = CreateZipfDistribution(distribution);
            var logRankFreqList = LogTransformDistribution(rankFreqList);

            if (!IsMonotonicDecline(logRankFreqList))
                return new Tuple<bool, List<KeyValuePair<double, double>>, double, double>(false, null, 0.0, 0.0);

            var p1Real = FindPointInDistribution(logRankFreqList, 0.5);
            var p2Inter = FindPointInDistribution(logRankFreqList, 0.3);
            var p3Inter = FindPointInDistribution(logRankFreqList, 0.7);

            if (p2Inter == null)
                return new Tuple<bool, List<KeyValuePair<double, double>>, double, double>(false, null, 0.0, 0.0);

            double interpolatedZ = Math.Abs(p2Inter.Item2 + p3Inter.Item2) * 0.5;

            double threshold = 0.05 * p1Real.Item2;
            double interpolatedDiff = Math.Abs(interpolatedZ - p1Real.Item2);
            bool isZipfian = false;

            if (interpolatedDiff <= threshold)
                isZipfian = true;

            //for metadata print
            double coeffientSlope = FindSlope(p2Inter.Item1, p2Inter.Item2, p3Inter.Item1, p3Inter.Item2);

            return new Tuple<bool, List<KeyValuePair<double, double>>, double, double>(isZipfian, logRankFreqList, interpolatedDiff, coeffientSlope);
        }

        public static List<KeyValuePair<int, int>> CreateZipfDistribution(Dictionary<string, int> distribution)
        {
            var sortedList = distribution.ToList().OrderByDescending(x => x.Value).ToList();
            List<KeyValuePair<int, int>> rankFreqList = new List<KeyValuePair<int, int>>();
            int rank = 1;

            foreach (var kvp in sortedList)
            {
                rankFreqList.Add(new KeyValuePair<int, int>(rank, kvp.Value));
                rank++;
            }
            rankFreqList = rankFreqList.OrderByDescending(x => x.Value).ThenBy(x => x.Key).ToList();

            return rankFreqList;
        }

        private static List<KeyValuePair<double, double>> LogTransformDistribution(List<KeyValuePair<int, int>> distribution)
        {
            var loggedDistribution = new List<KeyValuePair<double, double>>();
            foreach (var kvp in distribution)
            {
                // words seen three times or fewer are left out of the fit
                if (kvp.Value > 3)
                    loggedDistribution.Add(new KeyValuePair<double, double>(Math.Log(kvp.Key), Math.Log(kvp.Value)));
            }

            return loggedDistribution;
        }

        private static bool IsMonotonicDecline(List<KeyValuePair<double, double>> distribution)
        {
            double y = -1;
            foreach (var kvp in distribution)
            {
                if (y == -1)
                    y = kvp.Value;

                if (y < kvp.Value)
                    return false;

                y = kvp.Value;
            }

            return true;
        }

        private static Tuple<double, double> FindPointInDistribution(List<KeyValuePair<double, double>> distribution, double partOfDistribution)
        {
            int indexOfPoint = (int)(distribution.Count * partOfDistribution) - 1;

            if (indexOfPoint + 1 > distribution.Count)
                return null;

            if (indexOfPoint < 0)
                return null;

            return new Tuple<double, double>(distribution[indexOfPoint].Key, distribution[indexOfPoint].Value);
        }

        private static double FindSlope(double x1, double y1, double x2, double y2)
        {
            return (y2 - y1) / (x2 - x1);
        }
    }
}
