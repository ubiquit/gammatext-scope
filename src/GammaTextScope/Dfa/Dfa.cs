using System;
using System.Collections.Generic;
using System.Linq;

namespace GammaTextScope
{
    /// <summary>
    /// Detrended fluctuation analysis exponent α: the PhysioNet dfa.c port produces
    /// log10(box size) and log10(F(n)) for a geometric series of box sizes (ratio
    /// 2^(1/8), linear detrending, non-overlapping boxes, min box 4, max box N/4);
    /// α is the least-squares slope of that log–log relation.
    /// </summary>
    public class Dfa
    {
        public double CalculateDFABoxController(List<double> values)
        {
            Tuple<double[], double[]> logFns = new PhysioNetDfa().main(0, new string[] { "TheCalculator" }, values.ToArray());
            return FitLocalTrend(logFns.Item1, logFns.Item2);
        }

        private static double FitLocalTrend(double[] xvals, double[] yvals)
        {
            Tuple<double, double> p = MathNet.Numerics.LinearRegression.SimpleRegression.Fit(xvals, yvals);
            return p.Item2; // slope
        }
    }
}
