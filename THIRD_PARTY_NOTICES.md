# Third-party notices

## PhysioNet `dfa.c` (GPL-2.0-or-later)

`src/GammaTextScope/Dfa/PhysioNetDfa.cs` is a C# translation of `dfa.c` v4.9
(J. Mietus, C-K Peng and G. B. Moody, 2001–2005, PhysioNet), distributed under
the GNU General Public License version 2 or later. The original copyright and
licence header is retained in the file. Because of this component the whole
repository is licensed GPL-3.0-or-later.

Known deviations of the port from `dfa.c`:

* `dfa.c` swaps `minbox` and `maxbox` when `minbox > maxbox`; the port resets
  `minbox` to `2 * nfit` instead. This only matters for series shorter than 16
  points, which the pipeline never submits (DFA is reported for 30 or more
  decay-marker gaps).
* `dfa.c` exits on a singular matrix in the Gauss–Jordan step; the port
  continues. With linear detrending (`nfit = 2`) the pivot is always on the
  diagonal, so the case does not arise.
* The port allocates its work matrix for the largest box size up front
  (O(maxbox²) doubles per call), where `dfa.c` allocates per box size.

If you use this software in published research please cite, as requested by the
original authors:

- Peng C-K, Buldyrev SV, Havlin S, Simons M, Stanley HE, Goldberger AL.
  Mosaic organization of DNA nucleotides. Phys Rev E 1994;49:1685–1689.
- Peng C-K, Havlin S, Stanley HE, Goldberger AL. Quantification of scaling
  exponents and crossover phenomena in nonstationary heartbeat time series.
  Chaos 1995;5:82–87.
- Goldberger AL, Amaral LAN, Glass L, et al. PhysioBank, PhysioToolkit, and
  PhysioNet. Circulation 2000;101(23):e215–e220.

## Math.NET Numerics (MIT)

Used for the linear regression of the DFA log–log fit
(`MathNet.Numerics` NuGet package).

    Copyright (c) 2002-2022 Math.NET

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.

## OxyPlot (MIT)

Used for the PDF chart export of the experiments (`OxyPlot.Core` NuGet package).

    Copyright (c) 2014 OxyPlot contributors

    Licensed under the MIT licence; the licence text is the one reproduced
    above for Math.NET Numerics.

## Data and reference tables

* **GRCh38 centromere coordinates** (`CentromereRegions.cs`)
  are the `acen` bands of the UCSC Genome Browser `cytoBand` table for hg38
  (https://genome.ucsc.edu), rounded to 100 kb; for the acrocentric
  chromosomes 13, 14, 15, 21 and 22 the whole short arm up to the centromere is
  marked.
* **Gene annotation**: GENCODE release 44 (Ensembl 110), https://www.gencodegenes.org,
  lifted onto each consensus genome as described in the README (*Input provenance*).
* **HCC1395 / HCC1395BL**: SEQC2 somatic reference sample, high-confidence
  sSNV/sINDEL v1.2.1 — Fang L.T. et al., *Establishing community reference
  samples, data and call sets for benchmarking cancer mutation detection using
  whole-genome sequencing*, Nature Biotechnology 39, 1151–1160 (2021); normal
  assembly GCA_021234545.1 (NCBI).
* **COLO829 / COLO829BL**: Oxford Nanopore Technologies open data (Clair3 and
  ClairS calls); cell line characterised in Pleasance E.D. et al., *A
  comprehensive catalogue of somatic mutations from a human cancer genome*,
  Nature 463, 191–196 (2010).
* The sequences themselves are not redistributed in this repository; the
  README (*Input provenance*) describes how to rebuild them from the public
  sources.
