using System;
using System.Collections.Generic;
using System.Text;

namespace GammaTextScope
{
    public static class DistributionHelper
    {
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

        public static bool DeterminePeakFromAplha(double alpha, Color color)
        {
            if (IdentifyAlpha(alpha) == color)
                return true;

            return false;
        }

        private static Color IdentifyAlpha(double alpha)
        {
            if (alpha > 0 && alpha < 0.4)
            {
                return Color.tail;
            }
            else if (alpha >= 0.4 && alpha < 0.7)
            {
                return Color.white;

            }
            else if (alpha >= 0.7 && alpha < 0.9)
            {
                return Color.bulge;

            }
            else if (alpha >= 0.9 && alpha < 1.1)
            {
                return Color.pink;

            }
            else if (alpha >= 1.1 && alpha < 2.0)
            {
                return Color.brown;

            }

            return Color.na;
        }
    }

    public enum Color
    {
        tail,
        white,
        bulge,
        pink,
        brown,
        na
    }
}
