using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GammaTextScope
{
    public static class ZipfHelper
    {
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

        public static Tuple<bool, List<KeyValuePair<double, double>>> ClassifyDistributionSlope(Dictionary<string, int> distribution)
        {
            var rankFreqList = CreateZipfDistribution(distribution);
            var logRankFreqList = LogTransformDistribution(rankFreqList);

            if (!IsMonotonicDecline(logRankFreqList))
                return new Tuple<bool, List<KeyValuePair<double, double>>>(false, null);

            var p1top = FindPointInDistribution(logRankFreqList, 0.1);
            var p2top = FindPointInDistribution(logRankFreqList, 0.4);

            if (p1top == null || p2top == null)
                return new Tuple<bool, List<KeyValuePair<double, double>>>(false, null);

            double slopeTop = FindSlope(p1top.Item1, p1top.Item2, p2top.Item1, p2top.Item2);

            var p1bot = FindPointInDistribution(logRankFreqList, 0.4);
            var p2bot = FindPointInDistribution(logRankFreqList, 0.7);

            if (p1bot == null || p2bot == null)
                return new Tuple<bool, List<KeyValuePair<double, double>>>(false, null);

            double slopeBot = FindSlope(p1bot.Item1, p1bot.Item2, p2bot.Item1, p2bot.Item2);

            return new Tuple<bool, List<KeyValuePair<double, double>>>(IsSimilarSlope(slopeTop, slopeBot), logRankFreqList);
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
            List<KeyValuePair<double, double>> loggedDistribution = new List<KeyValuePair<double, double>>();
            int counter = 0;
            double lastVal = 0;
            foreach (var kvp in distribution)
            {
                //commented out to make full running
                if (kvp.Value > 3)
                {
                    double x = Math.Log(kvp.Key);
                    if (x < lastVal)
                    {
                        int i = 0;
                    }
                    counter++;
                    lastVal = x;
                    loggedDistribution.Add(new KeyValuePair<double, double>(x, Math.Log(kvp.Value)));
                }
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

        private static bool IsSimilarSlope(double slope1, double slope2)
        {
            double delta = Math.Abs(slope1 - slope2);

            double numerator = slope1 < slope2 ? slope2 : slope1;
            double denominator = slope1 < slope2 ? slope1 : slope2;

            double pctOfSlope1 = 1 - (numerator / denominator);

            if (pctOfSlope1 <= 0.1)
                return true;

            return false;
        }
    }
}
