using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GammaTextScope
{
    internal static class DfaArrayExtensions
    {
        public static T[] SubArray<T>(this T[] data, long index, long length)
        {
            T[] result = new T[length];
            if (index + length > data.Length)
            { length--; }
            Array.Copy(data, index, result, 0, length);
            return result;
        }
    }
    public class PhysioNetDfa
    {
        /* file: dfa.c	  J. Mietus, C-K Peng, and G. Moody	8 February 2001
Last revised:			        25 January 2005  v4.9

-------------------------------------------------------------------------------
dfa: Detrended Fluctuation Analysis (translated from C-K Peng's Fortran code)
Copyright (C) 2001-2005 Joe Mietus, C-K Peng, and George B. Moody

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place - Suite 330, Boston, MA 02111-1307, USA.

You may contact the authors by e-mail (peng@physionet.org) or postal mail
(Beth Israel Deaconess Medical Center, Room KS-B26, 330 Brookline Ave., Boston,
MA 02215 USA).  For updates to this software, please visit PhysioNet
(http://www.physionet.org/).
_______________________________________________________________________________

This method was first proposed in:
Peng C-K, Buldyrev SV, Havlin S, Simons M, Stanley HE, Goldberger AL. Mosaic
organization of DNA nucleotides. Phys Rev E 1994;49:1685-1689.  [Available
on-line at http://prola.aps.org/abstract/PRE/v49/i2/p1685_1]

A detailed description of the algorithm and its application to physiologic
signals can be found in:
Peng C-K, Havlin S, Stanley HE, Goldberger AL. Quantification of scaling
exponents and crossover phenomena in nonstationary heartbeat time series.
Chaos 1995;5:82-87. [Abstract online at http://www.ncbi.nlm.nih.gov/entrez/-
query.fcgi?cmd=Retrieve&db=PubMed&list_uids=11538314&dopt=Abstract]

If you use this program in support of published research, please include a
citation of at least one of the two references above, as well as the standard
citation for PhysioNet:
Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG,
Mietus JE, Moody GB, Peng CK, Stanley HE.  PhysioBank, PhysioToolkit, and
Physionet: Components of a New Research Resource for Complex Physiologic
Signals. Circulation 101(23):e215-e220 [Circulation Electronic Pages;
http://circ.ahajournals.org/cgi/content/full/101/23/e215]; 2000 (June 13).
*/

        static void SWAP<T>(ref T a, ref T b) { T temp = a; a = b; b = temp; }
        /* Global variables. */
        double[] seq;    /* input data buffer; allocated and filled by input() */
        long[] rs;   /* box size array; allocated and filled by rscale() */
        double[] mse;    /* fluctuation array; allocated by setup(), filled by dfa() */
        bool iflag = true;  /* integrate the input data if non-zero */
        int nfit = 2;   /* order of the regression fit, plus 1 */
        int nr;		/* number of box sizes */

        public Tuple<double[], double[]> main(int argc, string[] argv, double[] pointsIn)
        {
            int i, sw = 0;
            long minbox = 0L, maxbox = 0L, npts;

            /* dfa.c reads its options here; the port always uses the defaults: linear
               detrending (nfit = 2), integrate the input, minbox = 2*nfit, maxbox = npts/4,
               no sliding window. */

            /* Allocate and fill the input data array seq[]. */
            npts = input(pointsIn);

            /* Set minimum and maximum box sizes. */
            if (minbox < 2 * nfit) minbox = 2 * nfit;
            if (maxbox == 0 || maxbox > npts / 4) maxbox = npts / 4;
            if (minbox > maxbox)
            {
                // dfa.c swaps minbox and maxbox here; the port does not, which only matters for series shorter than 16 points
                if (minbox < 2 * nfit) minbox = 2 * nfit;
            }

            /* Allocate and fill the box size array rs[].  rscale's third argument
            specifies that the ratio between successive box sizes is 2^(1/8). */
            nr = rscale(minbox, maxbox, Math.Pow(2.0, 1.0 / 8.0));

            /* Allocate memory for dfa() and the functions it calls. */
            setup();

            /* Measure the fluctuations of the detrended input data at each box size
            using the DFA algorithm; fill mse[] with these results. */
            dfa(seq, npts, nfit, rs, nr, sw);

            string[] outpoints = new string[nr + 1];
            List<double> rsOut = new List<double>();
            List<double> mseOut = new List<double>();
            /* Output the results. */
            for (i = 1; i <= nr; i++)
            {
                outpoints[i] = string.Format("{0} {1}", Math.Log10((double)rs[i]), Math.Log10(mse[i]) / 2.0);
                rsOut.Add(Math.Log10((double)rs[i]));
                mseOut.Add(Math.Log10(mse[i]) / 2.0);
            }

            Tuple<double[], double[]> outTuple = new Tuple<double[], double[]>(rsOut.ToArray(), mseOut.ToArray());
            return outTuple;//exit(0);
        }

        /* Read input data, allocating and filling seq[], integrating if iflag != 0.
        Following the convention used for other arrays in this program, seq[0] is
        unused, and the first point is stored in seq[1].  The return value is the
        number of points read.

        This function allows the input buffer to grow as large as necessary, up to
        the available memory (assuming that a long int is large enough to address
        any memory location).  Note that the integration is done using double
        precision arithmetic to avoid complete loss of precision when the integrated
        data reach large amplitudes.  */

        long input(double[] timeSeries)
        {
            long lnpts;
            double y, yp = 0.0;
            seq = new double[timeSeries.Count() + 1];
            for (lnpts = 1L; lnpts <= timeSeries.Count(); lnpts++)//  while (scanf("%lf", &y) == 1)
            {
                y = timeSeries[lnpts - 1];
                seq[lnpts] = iflag ? (yp += y) : y;
            }

            return (timeSeries.Count());
        }
        int rslen;  /* length of rs[] */

        /* rscale() allocates and fills rs[], the array of box sizes used by dfa()
        below.  The box sizes range from (exactly) minbox to (approximately) maxbox,
        and are arranged in a geometric series such that the ratio between
        consecutive box sizes is (approximately) boxratio.  The return value is
        the number of box sizes in rs[].
        */
        int rscale(long minbox, long maxbox, double boxratio)
        {
            int ir, n;
            long rw;

            /* Determine how many scales are needed. */
            rslen = (int)(Math.Log10(maxbox / (double)minbox) / Math.Log10(boxratio) + 1.5);
            /* Thanks to Peter Domitrovich for pointing out that a previous version
            of the above calculation undercounted the number of scales in some
            situations. */
            rs = lvector(1, rslen);
            for (ir = 1, n = 2, rs[1] = minbox; n <= rslen && rs[n - 1] < maxbox; ir++)
                if ((rw = (long)(minbox * Math.Pow(boxratio, ir) + 0.5)) > rs[n - 1])
                    rs[n++] = rw;
            if (rs[--n] > maxbox) --n;
            return (n);
        }

        double[][] x; /* matrix of abscissas and their powers, for polyfit(). */

        /* Detrended fluctuation analysis
        seq:	input data array
        npts:	number of input points
        nfit:	order of detrending (2: linear, 3: quadratic, etc.)
        rs:		array of box sizes (uniformly distributed on log scale)
        nr:		number of entries in rs[] and mse[]
        sw:		mode (0: non-overlapping windows, 1: sliding window)
        This function returns the mean squared fluctuations in mse[].
        */
        private void dfa(double[] seq, long npts, int nfit, long[] rs, int nr, int sw)
        {
            long i, boxsize, inc, j;
            double stat;

            for (i = 1; i <= nr; i++) //changed
            {
                boxsize = rs[i];
                int c = 0;
                if (sw > 0) { inc = 1; stat = (int)(npts - boxsize + 1) * boxsize; }
                else { inc = boxsize; stat = (int)(npts / boxsize) * boxsize; }
                for (mse[i] = 0.0, j = 0; j <= npts - boxsize; j += inc)
                {
                    mse[i] += polyfit(x, seq.SubArray(j, boxsize + 1)/*seq + j*/, boxsize, nfit);
                    c++;
                    if (c == 9)
                    {

                    }
                }
                c = 0;
                mse[i] /= stat;
            }
        }

        /* workspace for polyfit() */
        double[] beta;
        double[][] covar, covar0;
        int[] indxc, indxr, ipiv;

        /* This function allocates workspace for dfa() and polyfit(), and sets
        x[i][j] = i**(j-1), in preparation for polyfit(). */
        private void setup()
        {
            long i;
            int j;

            beta = vector(1, nfit);
            covar = matrix(1, nfit, 1, nfit);
            covar0 = matrix(1, nfit, 1, nfit);
            indxc = ivector(1, nfit);
            indxr = ivector(1, nfit);
            ipiv = ivector(1, nfit);
            mse = vector(1, nr);
            x = matrix(1, rs[nr], 1, nfit);
            for (i = 1; i <= rs[nr]; i++)
            {
                x[i][1] = 1.0;
                x[i][2] = i;
                for (j = 3; j <= nfit; j++)
                    x[i][j] = x[i][j - 1] * i;
            }
        }

        /* polyfit() is based on lfit() and gaussj() from Numerical Recipes in C
        (Press, Teukolsky, Vetterling, and Flannery; Cambridge U. Press, 1992).  It
        fits a polynomial of degree (nfit-1) to a set of boxsize points given by
        x[1...boxsize][2] and y[1...boxsize].  The return value is the sum of the
        squared errors (chisq) between the (x,y) pairs and the fitted polynomial.
        */
        double polyfit(double[][] x, double[] y, long boxsize, int nfit)
        {
            int icol, irow, j, k;
            double big, chisq, pivinv, temp;
            long i;
            long pboxsize = 0L; //static
            icol = -1;
            irow = -1;
            /* This block sets up the covariance matrix.  Provided that boxsize
            never decreases (which is true in this case), covar0 can be calculated
            incrementally from the previous value. */
            if (pboxsize != boxsize)
            {   /* this will be false most of the time */
                if (pboxsize > boxsize) /* this should never happen */
                    pboxsize = 0L;
                if (pboxsize == 0L) /* this should be true the first time only */
                    for (j = 1; j <= nfit; j++)
                        for (k = 1; k <= nfit; k++)
                        {
                            covar0[j][k] = 0.0;
                        }
                for (i = pboxsize + 1; i <= boxsize; i++)
                    for (j = 1; j <= nfit; j++)
                        for (k = 1, temp = x[i][j]; k <= j; k++)
                            covar0[j][k] += temp * x[i][k];
                for (j = 2; j <= nfit; j++)
                    for (k = 1; k < j; k++)
                        covar0[k][j] = covar0[j][k];
                pboxsize = boxsize;
            }
            for (j = 1; j <= nfit; j++)
            {
                beta[j] = ipiv[j] = 0;
                for (k = 1; k <= nfit; k++)
                    covar[j][k] = covar0[j][k];
            }
            for (i = 1; i <= boxsize; i++)//changed
            {
                beta[1] += (temp = y[i]);
                beta[2] += temp * i;
            }
            if (nfit > 2)
                for (i = 1; i <= boxsize; i++)
                    for (j = 3, temp = y[i]; j <= nfit; j++)
                        beta[j] += temp * x[i][j];
            for (i = 1; i <= nfit; i++)
            {
                big = 0.0;
                for (j = 1; j <= nfit; j++)
                {
                    if (ipiv[j] != 1)
                    {
                        for (k = 1; k <= nfit; k++)
                        {
                            if (ipiv[k] == 0)
                            {
                                if ((temp = covar[j][k]) >= big ||
                                    (temp = -temp) >= big)
                                {
                                    big = temp;
                                    irow = j;
                                    icol = k;
                                }
                            }
                            else if (ipiv[k] > 1)
                            {
                                // dfa.c calls error("singular matrix") and exits here; the port continues
                            }
                        }
                    }
                }
                ++(ipiv[icol]);
                if (irow != icol)
                {
                    for (j = 1; j <= nfit; j++) SWAP(ref covar[irow][j], ref covar[icol][j]);
                    SWAP(ref beta[irow], ref beta[icol]);
                }
                indxr[i] = irow;
                indxc[i] = icol;
                if (covar[icol][icol] == 0.0)
                {
                    // dfa.c calls error("singular matrix") and exits here; the port continues
                }
                pivinv = 1.0 / covar[icol][icol];
                covar[icol][icol] = 1.0;
                for (j = 1; j <= nfit; j++) covar[icol][j] *= pivinv;
                beta[icol] *= pivinv;
                for (j = 1; j <= nfit; j++)
                    if (j != icol)
                    {
                        temp = covar[j][icol];
                        covar[j][icol] = 0.0;
                        for (k = 1; k <= nfit; k++) covar[j][k] -= covar[icol][k] * temp;
                        beta[j] -= beta[icol] * temp;
                    }
            }
            chisq = 0.0;
            if (nfit <= 2)
            {
                for (i = 1; i <= boxsize; i++) //changed index
                {
                    temp = beta[1] + beta[2] * i - y[i];
                    chisq += temp * temp;
                }
            }
            else
            {
                for (i = 1; i <= boxsize; i++)
                {
                    temp = beta[1] + beta[2] * i - y[i];
                    for (j = 3; j <= nfit; j++) temp += beta[j] * x[i][j];
                    chisq += temp * temp;
                }
            }
            return (chisq);
        }

        /* The allocators below correspond to vector(), ivector(), lvector() and matrix() of
        Numerical Recipes (see above). The C versions malloc and offset the pointer so that
        v[nl..nh] is addressable; here the arrays are simply one element longer than needed
        and indexed directly, and there is nothing to free. */
        double[] vector(long nl, long nh) => new double[nh - nl + 2];

        int[] ivector(long nl, long nh) => new int[nh - nl + 2];

        long[] lvector(long nl, long nh) => new long[nh - nl + 2];

        /* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
        double[][] matrix(long nrl, long nrh, long ncl, long nch)
        {
            long nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
            var m = new double[nrow + 1][];
            for (int jj = 0; jj < m.Length; jj++)
                m[jj] = new double[nrow * ncol + 1];
            return m;
        }
    }
}
