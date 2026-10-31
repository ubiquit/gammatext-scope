using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace GammaTextScope
{
    /// <summary>
    /// Combines alpha/zipf window classification with stop codon counts per window.
    /// For each WW1 window in each chromosome, outputs:
    ///   - Classification (alfa / zipf / ikke) for Normal and Tumor
    ///   - Stop codon counts and delta for that genomic region
    /// Produces per-chromosome TSV + chart, and a genome-wide TSV + chart.
    ///
    /// When FeatureType.WholeGenome is requested, there are no WholeGenome alpha index
    /// files (too expensive to compute). Instead, Intergenic and Gene alpha results are
    /// combined: windows from both are mapped to raw chromosome positions, merged, and
    /// sorted. Stop codons are counted on the respective filtered sequences (same DNA
    /// that was classified).
    /// </summary>
    public class StopCodonAlphaComparison
    {
        private const int WindowWidth = 1; // WW1 ≈ 133k bp, closest to 100k

        /// <param name="merged">WholeGenome only: use the union of the Intergenic and Gene window
        /// classifications (legacy behaviour, output suffix WholeGenomeMerged) instead of the
        /// annotation-free whole-genome classification.</param>
        public void Run(string speciesNormal, string speciesTumor, string pairName = null,
                        FeatureType featureType = FeatureType.Intergenic, bool merged = false)
        {
            pairName ??= speciesNormal + "_VS_" + speciesTumor;

            if (featureType == FeatureType.WholeGenome && merged)
            {
                RunWholeGenome(speciesNormal, speciesTumor, pairName);
                return;
            }
            if (featureType == FeatureType.WholeGenome)
            {
                foreach (var s in new[] { speciesNormal, speciesTumor })
                    if (!Directory.Exists(GetIndexDir(s, FeatureType.WholeGenome)))
                        throw new ArgumentException($"no WholeGenome classification for {s} ({GetIndexDir(s, FeatureType.WholeGenome)}): run 'classify --mode WholeGenome' first, or pass --merged for the legacy union of Intergenic and Gene");
            }

            RunSingleFeature(speciesNormal, speciesTumor, pairName, featureType);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Single feature type path (Intergenic, Gene, Exon)
        // ═══════════════════════════════════════════════════════════════

        private void RunSingleFeature(string speciesNormal, string speciesTumor,
                                      string pairName, FeatureType featureType)
        {
            string featureLabel = featureType.ToString();
            string outDir = PipelineConfig.Default.TumorOutDir(pairName, "StopCodonAlphaComparison", featureLabel) + System.IO.Path.DirectorySeparatorChar;
            EnsureDirectory(outDir);

            string normalIndexDir = GetIndexDir(speciesNormal, featureType);
            string tumorIndexDir = GetIndexDir(speciesTumor, featureType);

            Console.WriteLine($"  Loading gene annotations for {featureLabel}...");
            var normalGeneIndex = SequenceHelper.ReadGeneIndexFromStream(PipelineConfig.Default.StreamDir(speciesNormal + "_GeneInfo"), featureType);
            var tumorGeneIndex = SequenceHelper.ReadGeneIndexFromStream(PipelineConfig.Default.StreamDir(speciesTumor + "_GeneInfo"), featureType);

            var normalChrFiles = GetChrSequenceFiles(speciesNormal);
            var tumorChrFiles = GetChrSequenceFiles(speciesTumor);

            var chrNames = DiscoverChromosomes(normalIndexDir);

            var allRows = new List<WindowRow>();
            var chrOffsetInfo = new Dictionary<string, (int Breakpoints, int FinalOffset)>();
            int totalBreakpoints = 0;
            int totalFinalOffset = 0;

            foreach (var chr in chrNames.OrderBy(NaturalChrOrderKey))
            {
                Console.WriteLine($"  StopCodonAlpha ({featureLabel}): {chr}");

                var result = ProcessChrFeature(chr, normalIndexDir, tumorIndexDir,
                    normalChrFiles, tumorChrFiles, normalGeneIndex, tumorGeneIndex, featureType);

                if (result.Breakpoints >= 0)
                {
                    chrOffsetInfo[chr] = (result.Breakpoints, result.FinalOffset);
                    totalBreakpoints += result.Breakpoints;
                    totalFinalOffset += result.FinalOffset;
                }

                // Per-chromosome TSV + chart
                string chrTsvPath = Path.Combine(outDir, $"StopCodonAlpha_{pairName}_{chr}.tsv");
                WriteTsv(result.Rows, chrTsvPath, includeChromosomeColumn: false);

                // Centromere range for this chromosome
                var chrCentro = FindCentromereWindowRange(result.Rows);
                var chrCentroList = chrCentro.HasValue
                    ? new List<(int, int)> { chrCentro.Value }
                    : null;

                string chrChartPath = Path.Combine(outDir, $"Chart_StopCodonAlpha_{pairName}_{chr}.pdf");
                string chrSubtitle = chrOffsetInfo.TryGetValue(chr, out var oi)
                    ? $"Alignment: {oi.Breakpoints} breakpoints, final offset {oi.FinalOffset:+#;-#;0} bp | {featureLabel}"
                    : featureLabel;
                ChartPrinter.PrintStopCodonAlphaLineChart(
                    result.Rows.Select(r => r.Delta).ToList(),
                    result.Rows.Select(r => r.NormalClass).ToList(),
                    result.Rows.Select(r => r.TumorClass).ToList(),
                    $"Stop-codon delta and qualification – {chr} ({pairName})",
                    "Window (sequential)",
                    "Stop-codon delta (Tumour − Normal)",
                    chrChartPath,
                    chromosomeBoundaries: null,
                    subtitle: chrSubtitle,
                    centromereRanges: chrCentroList);

                allRows.AddRange(result.Rows);
            }

            WriteGenomeWideOutput(allRows, chrNames, chrOffsetInfo,
                totalBreakpoints, totalFinalOffset, featureLabel, pairName, outDir);

            Console.WriteLine($"StopCodonAlphaComparison ({featureLabel}) done. Output in:\n  {outDir}");
        }

        // ═══════════════════════════════════════════════════════════════
        //  WholeGenomeMerged path (legacy): merge Intergenic + Gene alpha results
        // ═══════════════════════════════════════════════════════════════

        private void RunWholeGenome(string speciesNormal, string speciesTumor, string pairName)
        {
            const string featureLabel = "WholeGenomeMerged";
            string outDir = PipelineConfig.Default.TumorOutDir(pairName, "StopCodonAlphaComparison", featureLabel) + System.IO.Path.DirectorySeparatorChar;
            EnsureDirectory(outDir);

            // Gene annotations are the same for Gene and Intergenic
            Console.WriteLine("  Loading gene annotations...");
            var normalGeneIndex = SequenceHelper.ReadGeneIndexFromStream(PipelineConfig.Default.StreamDir(speciesNormal + "_GeneInfo"), FeatureType.Gene);
            var tumorGeneIndex = SequenceHelper.ReadGeneIndexFromStream(PipelineConfig.Default.StreamDir(speciesTumor + "_GeneInfo"), FeatureType.Gene);

            var normalChrFiles = GetChrSequenceFiles(speciesNormal);
            var tumorChrFiles = GetChrSequenceFiles(speciesTumor);

            // Index dirs for both feature types
            string normalInterDir = GetIndexDir(speciesNormal, FeatureType.Intergenic);
            string tumorInterDir = GetIndexDir(speciesTumor, FeatureType.Intergenic);
            string normalGeneDir = GetIndexDir(speciesNormal, FeatureType.Gene);
            string tumorGeneDir = GetIndexDir(speciesTumor, FeatureType.Gene);

            // Discover chromosomes (intergenic has the largest set)
            var chrNames = DiscoverChromosomes(normalInterDir);
            var geneChrNames = DiscoverChromosomes(normalGeneDir);
            foreach (var c in geneChrNames)
                if (!chrNames.Contains(c)) chrNames.Add(c);

            var allRows = new List<WindowRow>();
            var chrOffsetInfo = new Dictionary<string, (int Breakpoints, int FinalOffset)>();
            int totalBreakpoints = 0;
            int totalFinalOffset = 0;

            foreach (var chr in chrNames.OrderBy(NaturalChrOrderKey))
            {
                Console.WriteLine($"  StopCodonAlpha (WholeGenome merge): {chr}");

                // Process both feature types
                var interResult = ProcessChrFeature(chr, normalInterDir, tumorInterDir,
                    normalChrFiles, tumorChrFiles,
                    CloneGeneIndex(normalGeneIndex), CloneGeneIndex(tumorGeneIndex),
                    FeatureType.Intergenic);

                var geneResult = ProcessChrFeature(chr, normalGeneDir, tumorGeneDir,
                    normalChrFiles, tumorChrFiles,
                    CloneGeneIndex(normalGeneIndex), CloneGeneIndex(tumorGeneIndex),
                    FeatureType.Gene);

                // Build filtered-to-raw mappers to sort windows by genomic position
                var chrGeneList = GetChrGeneList(chr, normalGeneIndex);
                long rawChrLen = CountRawChars(normalChrFiles, chr);

                if (rawChrLen > 0 && chrGeneList.Count > 0)
                {
                    var interMapper = FilteredToRawMapper.Build(chrGeneList, rawChrLen, isGene: false);
                    var geneMapper = FilteredToRawMapper.Build(chrGeneList, rawChrLen, isGene: true);

                    Console.WriteLine($"    mapper filtered lengths: inter={interMapper.TotalFilteredLength:N0}, gene={geneMapper.TotalFilteredLength:N0}, raw={rawChrLen:N0}");

                    foreach (var r in interResult.Rows)
                        r.RawMidpoint = interMapper.GetRawMidpoint(r.StartPos, r.Length);
                    foreach (var r in geneResult.Rows)
                        r.RawMidpoint = geneMapper.GetRawMidpoint(r.StartPos, r.Length);
                }
                else
                {
                    // Fallback: no gene data, just concatenate
                    int offset = interResult.Rows.Count > 0
                        ? interResult.Rows.Last().RawMidpoint + 1
                        : 0;
                    foreach (var r in geneResult.Rows)
                        r.RawMidpoint = offset++;
                }

                // Merge and sort by raw genomic position
                var merged = interResult.Rows.Concat(geneResult.Rows)
                    .OrderBy(r => r.RawMidpoint)
                    .ToList();

                // Re-number windows sequentially
                for (int i = 0; i < merged.Count; i++)
                    merged[i].WindowIndex = i + 1;

                // Track combined alignment stats
                int chrBp = 0, chrOffset = 0;
                if (interResult.Breakpoints >= 0) { chrBp += interResult.Breakpoints; chrOffset += interResult.FinalOffset; }
                if (geneResult.Breakpoints >= 0) { chrBp += geneResult.Breakpoints; chrOffset += geneResult.FinalOffset; }
                if (interResult.Breakpoints >= 0 || geneResult.Breakpoints >= 0)
                {
                    chrOffsetInfo[chr] = (chrBp, chrOffset);
                    totalBreakpoints += chrBp;
                    totalFinalOffset += chrOffset;
                }

                Console.WriteLine($"    merged: {merged.Count} windows ({interResult.Rows.Count} intergenic + {geneResult.Rows.Count} gene)");

                // Per-chromosome TSV + chart
                string chrTsvPath = Path.Combine(outDir, $"StopCodonAlpha_{pairName}_{chr}.tsv");
                WriteTsv(merged, chrTsvPath, includeChromosomeColumn: false);

                // Centromere range for this chromosome
                var chrCentroWG = FindCentromereWindowRange(merged);
                var chrCentroListWG = chrCentroWG.HasValue
                    ? new List<(int, int)> { chrCentroWG.Value }
                    : null;

                string chrChartPath = Path.Combine(outDir, $"Chart_StopCodonAlpha_{pairName}_{chr}.pdf");
                string chrSubtitle = chrOffsetInfo.TryGetValue(chr, out var oi)
                    ? $"Alignment: {oi.Breakpoints} breakpoints, final offset {oi.FinalOffset:+#;-#;0} bp | {featureLabel} (Intergenic+Gene)"
                    : $"{featureLabel} (Intergenic+Gene)";
                ChartPrinter.PrintStopCodonAlphaLineChart(
                    merged.Select(r => r.Delta).ToList(),
                    merged.Select(r => r.NormalClass).ToList(),
                    merged.Select(r => r.TumorClass).ToList(),
                    $"Stop-codon delta and qualification – {chr} ({pairName})",
                    "Window (sequential)",
                    "Stop-codon delta (Tumour − Normal)",
                    chrChartPath,
                    chromosomeBoundaries: null,
                    subtitle: chrSubtitle,
                    centromereRanges: chrCentroListWG);

                allRows.AddRange(merged);
            }

            WriteGenomeWideOutput(allRows, chrNames, chrOffsetInfo,
                totalBreakpoints, totalFinalOffset, featureLabel, pairName, outDir);

            Console.WriteLine($"StopCodonAlphaComparison ({featureLabel}) done. Output in:\n  {outDir}");
        }

        // ═══════════════════════════════════════════════════════════════
        //  Shared: process one chromosome for one feature type
        // ═══════════════════════════════════════════════════════════════

        private static ChrFeatureResult ProcessChrFeature(
            string chr,
            string normalIndexDir, string tumorIndexDir,
            Dictionary<string, string> normalChrFiles, Dictionary<string, string> tumorChrFiles,
            Dictionary<string, List<EnsembleGeneEntry>> normalGeneIndex,
            Dictionary<string, List<EnsembleGeneEntry>> tumorGeneIndex,
            FeatureType featureType)
        {
            var normalWindows = ParseIndexFile(normalIndexDir, chr);
            var tumorWindows = ParseIndexFile(tumorIndexDir, chr);

            string normalSeq = "";
            string tumorSeq = "";

            if (normalChrFiles.TryGetValue(chr, out var nf))
            {
                var nGeneIdx = CloneGeneIndex(normalGeneIndex);
                normalSeq = CleanSequence(WindowHelper.ReadSequenceFromFile(nf, Path.GetFileName(nf), nGeneIdx, featureType));
            }
            if (tumorChrFiles.TryGetValue(chr, out var tf))
            {
                var tGeneIdx = CloneGeneIndex(tumorGeneIndex);
                tumorSeq = CleanSequence(WindowHelper.ReadSequenceFromFile(tf, Path.GetFileName(tf), tGeneIdx, featureType));
            }

            // Build alignment offset map
            AlignmentOffsetMap offsetMap = null;
            int breakpoints = -1;
            int finalOffset = 0;

            if (normalSeq.Length > 0 && tumorSeq.Length > 0)
            {
                offsetMap = AlignmentOffsetMap.Build(normalSeq, tumorSeq);
                breakpoints = offsetMap.BreakpointCount;
                finalOffset = offsetMap.FinalOffset;
            }

            // Build window rows
            int count = Math.Min(normalWindows.Count, tumorWindows.Count);
            var rows = new List<WindowRow>(count);

            for (int i = 0; i < count; i++)
            {
                var nw = normalWindows[i];
                var tw = tumorWindows[i];

                int startPos = nw.StartPos;
                int length = nw.Length;

                long normalStops = CountStopsInRange(normalSeq, startPos, length);

                long tumorStops;
                if (offsetMap != null)
                {
                    tumorStops = CountStopsAlignedChunked(tumorSeq, offsetMap, startPos, length);
                }
                else
                {
                    tumorStops = CountStopsInRange(tumorSeq, startPos, length);
                }

                long delta = tumorStops - normalStops;

                rows.Add(new WindowRow
                {
                    Chr = chr,
                    WindowIndex = i + 1,
                    StartPos = startPos,
                    Length = length,
                    NormalClass = nw.Classification,
                    TumorClass = tw.Classification,
                    NormalStopCodons = normalStops,
                    TumorStopCodons = tumorStops,
                    Delta = delta
                });
            }

            return new ChrFeatureResult { Rows = rows, Breakpoints = breakpoints, FinalOffset = finalOffset };
        }

        // ═══════════════════════════════════════════════════════════════
        //  Shared: genome-wide output
        // ═══════════════════════════════════════════════════════════════

        private static void WriteGenomeWideOutput(
            List<WindowRow> allRows, List<string> chrNames,
            Dictionary<string, (int Breakpoints, int FinalOffset)> chrOffsetInfo,
            int totalBreakpoints, int totalFinalOffset,
            string featureLabel, string pairName, string outDir)
        {
            // Genome-wide TSV
            string genomeTsvPath = Path.Combine(outDir, $"StopCodonAlpha_{pairName}_AllChromosomes.tsv");
            WriteTsv(allRows, genomeTsvPath, includeChromosomeColumn: true);

            // Chromosome boundary indices for the genome-wide chart
            var boundaries = new List<(int Index, string Label)>();
            int runningIdx = 0;
            foreach (var chr in chrNames.OrderBy(NaturalChrOrderKey))
            {
                int chrCount = allRows.Count(r => r.Chr == chr);
                if (chrCount > 0)
                {
                    boundaries.Add((runningIdx, chr));
                    runningIdx += chrCount;
                }
            }

            // Compute centromere ranges for the genome-wide chart
            var centromereRanges = ComputeGenomeWideCentromereRanges(allRows, chrNames);

            string genomeChartPath = Path.Combine(outDir, $"Chart_StopCodonAlpha_{pairName}_AllChromosomes.pdf");
            string genomeSubtitle = $"Alignment totals: {totalBreakpoints} breakpoints across {chrOffsetInfo.Count} chromosomes, net offset {totalFinalOffset:+#;-#;0} bp | {featureLabel}";
            ChartPrinter.PrintStopCodonAlphaLineChart(
                allRows.Select(r => r.Delta).ToList(),
                allRows.Select(r => r.NormalClass).ToList(),
                allRows.Select(r => r.TumorClass).ToList(),
                $"Stop-codon delta and qualification – whole genome ({pairName})",
                "Window (sequential)",
                "Stop-codon delta (Tumour − Normal)",
                genomeChartPath,
                chromosomeBoundaries: boundaries,
                subtitle: genomeSubtitle,
                centromereRanges: centromereRanges);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Gene list helpers (for FilteredToRawMapper)
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Get sorted gene list for a chromosome, replicating WindowHelper.GetChrGenes logic.
        /// </summary>
        private static List<EnsembleGeneEntry> GetChrGeneList(
            string chr, Dictionary<string, List<EnsembleGeneEntry>> geneIndex)
        {
            // chr comes in as "chr1", "chrX", etc.
            var result = new List<EnsembleGeneEntry>();
            foreach (var kv in geneIndex)
                foreach (var gene in kv.Value)
                    if (gene.Chromosome != null && gene.Chromosome.Equals(chr, StringComparison.OrdinalIgnoreCase))
                        result.Add(gene);

            return result.OrderBy(g => g.StartPosition).ThenBy(g => g.EndPosition).ToList();
        }

        /// <summary>
        /// Count total characters in the raw chromosome file (excluding newlines).
        /// </summary>
        private static long CountRawChars(Dictionary<string, string> chrFiles, string chr)
        {
            if (!chrFiles.TryGetValue(chr, out var path)) return 0;

            long count = 0;
            foreach (var line in File.ReadLines(path))
                count += line.Length;
            return count;
        }

        // ───────── Index file parsing ─────────

        private static string GetIndexDir(string species, FeatureType featureType)
        {
            return Path.Combine(PipelineConfig.Default.ClassificationDir(species, featureType),
                $"WindowWidth{WindowWidth}");
        }

        private static List<string> DiscoverChromosomes(string indexDir)
        {
            var result = new List<string>();
            if (!Directory.Exists(indexDir)) return result;

            foreach (var f in Directory.GetFiles(indexDir, "*_window_positions_indexLength.txt"))
            {
                string name = Path.GetFileName(f);
                int wwPos = name.IndexOf("WW", StringComparison.OrdinalIgnoreCase);
                if (wwPos > 0 && PipelineConfig.Default.IncludesChromosome(name.Substring(0, wwPos)))
                    result.Add(name.Substring(0, wwPos));
            }
            return result;
        }

        private static List<WindowInfo> ParseIndexFile(string indexDir, string chr)
        {
            string path = Path.Combine(indexDir, $"{chr}WW{WindowWidth}_window_positions_indexLength.txt");
            var windows = new List<WindowInfo>();
            if (!File.Exists(path)) return windows;

            foreach (var line in File.ReadAllLines(path))
            {
                if (string.IsNullOrWhiteSpace(line)) continue;
                var parts = line.Split('\t');
                if (parts.Length < 4) continue;

                string status = parts[3].Trim();
                string classification;

                if (status == "QUALIFIED")
                {
                    bool isAlpha = parts.Length > 4 && parts[4].Trim().Equals("True", StringComparison.OrdinalIgnoreCase);
                    classification = isAlpha ? "alfa" : "zipf";
                }
                else
                {
                    classification = "ikke";
                }

                int startPos = int.TryParse(parts[1].Trim(), out var sp) ? sp : 0;
                int length = int.TryParse(parts[2].Trim(), out var len) ? len : 0;

                windows.Add(new WindowInfo
                {
                    StartPos = startPos,
                    Length = length,
                    Classification = classification
                });
            }
            return windows;
        }

        // ───────── Stop codon counting ─────────

        /// <summary>
        /// Count stop codons in the tumour range that corresponds to a normal window. A window
        /// of up to 200 kb is translated through the offset map in one go. A longer window
        /// (one that spans an assembly gap, a centromere or an acrocentric arm) is translated
        /// in 100-kb chunks, so every chunk is re-anchored at its own translated start and an
        /// offset that drifts inside the window cannot carry the count into the wrong region.
        /// Every chunk is counted; a stop codon that straddles a chunk boundary is counted by
        /// the chunk its first base falls in, so the sum equals a count of the whole range.
        ///
        /// Until 29 October 2026 a chunk whose stop-codon density differed by more than 2x
        /// from the window mean was dropped and the sum of the remaining chunks was scaled up
        /// by chunks/accepted. An N-run has density 0, so its chunks were always dropped and
        /// the real chunks of the window were multiplied by 1.5, 2, 4 ... — every |Delta| above
        /// about 100 in the dataset v1 and v2 tables, and every "Tumor_StopCodons = 0" row,
        /// came from that rule and not from the sequences (see README, Reproducing the article).
        /// </summary>
        private static long CountStopsAlignedChunked(string tumorSeq, AlignmentOffsetMap offsetMap,
            int normalStart, int normalLength)
        {
            const int ChunkSize = 100_000;

            if (normalLength <= ChunkSize * 2)
            {
                var (tStart, tLen) = offsetMap.TranslateRange(normalStart, normalLength);
                return CountStopsInRange(tumorSeq, tStart, tLen);
            }

            long total = 0;
            for (int offset = 0; offset < normalLength; offset += ChunkSize)
            {
                int chunkLen = Math.Min(ChunkSize, normalLength - offset);
                bool lastChunk = offset + chunkLen >= normalLength;
                var (tStart, tLen) = offsetMap.TranslateRange(normalStart + offset, chunkLen);
                // CountStopsInRange counts codons that start up to two bases before the end of
                // the range, so a chunk is given two extra bases to catch the codons that start
                // in its last two bases; the last chunk ends where the window ends.
                total += CountStopsInRange(tumorSeq, tStart, lastChunk ? tLen : tLen + 2);
            }
            return total;
        }

        private static long CountStopsInRange(string seq, int startPos, int length)
        {
            if (string.IsNullOrEmpty(seq) || length <= 0) return 0;
            if (startPos < 0) { length += startPos; startPos = 0; }
            if (length <= 0) return 0;

            int end = Math.Min(seq.Length, startPos + length);
            int maxI = end - 3;
            long count = 0;

            for (int i = startPos; i <= maxI; i++)
            {
                char a = seq[i], b = seq[i + 1], c = seq[i + 2];
                if (a == 'N' || b == 'N' || c == 'N') continue;

                bool isStop =
                    (a == 'T' && b == 'A' && (c == 'A' || c == 'G')) ||
                    (a == 'T' && b == 'G' && c == 'A');
                if (isStop) count++;
            }
            return count;
        }

        // ───────── Sequence files ─────────

        private static Dictionary<string, string> GetChrSequenceFiles(string species)
        {
            string dir = PipelineConfig.Default.StreamDir(species);
            var map = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            if (!Directory.Exists(dir)) return map;

            foreach (var path in Directory.GetFiles(dir))
            {
                var name = Path.GetFileName(path);
                var chr = name.Split('_').Last().Split('.').First();
                if (!chr.StartsWith("chr", StringComparison.Ordinal))
                {
                    Console.Error.WriteLine($"  skipping {name}: not a chromosome file (expected chrN.<anything>)");
                    continue;
                }
                chr = NormalizeChr(chr);
                map["chr" + chr] = path;
            }
            return map;
        }

        private static Dictionary<string, List<EnsembleGeneEntry>> CloneGeneIndex(
            Dictionary<string, List<EnsembleGeneEntry>> src)
        {
            var clone = new Dictionary<string, List<EnsembleGeneEntry>>(src.Count);
            foreach (var kv in src)
                clone[kv.Key] = new List<EnsembleGeneEntry>(kv.Value);
            return clone;
        }

        private static string NormalizeChr(string chr)
        {
            if (chr.StartsWith("chr", StringComparison.OrdinalIgnoreCase))
                chr = chr.Substring(3);
            return chr.ToUpperInvariant();
        }

        // ───────── TSV output ─────────

        private static void WriteTsv(List<WindowRow> rows, string path, bool includeChromosomeColumn)
        {
            using var sw = new StreamWriter(path) { NewLine = GenomeConstants.OutputNewLine };

            if (includeChromosomeColumn)
                sw.WriteLine("Chr\tWindow\tStartPos\tLength\tNormal_Class\tTumor_Class\tNormal_StopCodons\tTumor_StopCodons\tDelta");
            else
                sw.WriteLine("Window\tStartPos\tLength\tNormal_Class\tTumor_Class\tNormal_StopCodons\tTumor_StopCodons\tDelta");

            foreach (var r in rows)
            {
                if (includeChromosomeColumn)
                    sw.WriteLine($"{r.Chr}\t{r.WindowIndex}\t{r.StartPos}\t{r.Length}\t{r.NormalClass}\t{r.TumorClass}\t{r.NormalStopCodons}\t{r.TumorStopCodons}\t{r.Delta}");
                else
                    sw.WriteLine($"{r.WindowIndex}\t{r.StartPos}\t{r.Length}\t{r.NormalClass}\t{r.TumorClass}\t{r.NormalStopCodons}\t{r.TumorStopCodons}\t{r.Delta}");
            }
        }

        // ───────── Helpers ─────────

        private static string CleanSequence(string s)
        {
            var sb = new System.Text.StringBuilder(s.Length);
            foreach (var ch in s)
            {
                if (char.IsWhiteSpace(ch)) continue;
                char u = char.ToUpperInvariant(ch);
                switch (u)
                {
                    case 'A': case 'C': case 'G': case 'T': case 'N': sb.Append(u); break;
                    default: sb.Append('N'); break;
                }
            }
            return sb.ToString();
        }

        private static int NaturalChrOrderKey(string chr)
        {
            string num = chr.StartsWith("chr", StringComparison.OrdinalIgnoreCase) ? chr.Substring(3) : chr;
            if (int.TryParse(num, out int n)) return n;
            if (num.Equals("X", StringComparison.OrdinalIgnoreCase)) return 23;
            if (num.Equals("Y", StringComparison.OrdinalIgnoreCase)) return 24;
            if (num.Equals("MT", StringComparison.OrdinalIgnoreCase) || num.Equals("M", StringComparison.OrdinalIgnoreCase)) return 25;
            return int.MaxValue;
        }

        private static void EnsureDirectory(string path)
        {
            if (!Directory.Exists(path)) Directory.CreateDirectory(path);
        }

        // ───────── Centromere helpers ─────────

        /// <summary>
        /// For a list of window rows from a single chromosome, find the window index range
        /// that overlaps the centromere. Uses a two-pass approach:
        ///
        ///  1. Primary signal: windows with Tumor_StopCodons = 0. The counting rule in force
        ///     until 29 October 2026 wrote such rows where no chunk of a long window passed its
        ///     density check (windows spanning an assembly gap); the current rule never does, so
        ///     this pass finds nothing and the fallback below applies. The run from the first
        ///     to the last such window is marked, widened by adjacent windows whose
        ///     tumour/normal ratio is outside 0.5–2.
        ///
        ///  2. If there is no such window, fall back to the UCSC centromere coordinates mapped
        ///     proportionally onto the filtered coordinate (approximate; the figure script in
        ///     figures/ maps them exactly through the gene mask instead).
        ///
        /// Returns null if no centromere region can be identified.
        /// </summary>
        private static (int StartIndex, int EndIndex)? FindCentromereWindowRange(List<WindowRow> rows)
        {
            if (rows.Count == 0) return null;
            string chr = rows[0].Chr;

            // Pass 1: find all tStops=0 windows and their neighbourhood
            int firstZero = -1, lastZero = -1;
            for (int i = 0; i < rows.Count; i++)
            {
                if (rows[i].TumorStopCodons == 0)
                {
                    if (firstZero < 0) firstZero = i;
                    lastZero = i;
                }
            }

            if (firstZero >= 0)
            {
                // Expand range to include adjacent windows with very low or very high tStops
                // (near-centromeric alignment artifacts)
                int startIdx = firstZero;
                int endIdx = lastZero;

                // Expand left
                while (startIdx > 0)
                {
                    var r = rows[startIdx - 1];
                    if (r.NormalStopCodons > 0)
                    {
                        double ratio = (double)r.TumorStopCodons / r.NormalStopCodons;
                        if (ratio < 0.5 || ratio > 2.0)
                            startIdx--;
                        else break;
                    }
                    else break;
                }

                // Expand right
                while (endIdx < rows.Count - 1)
                {
                    var r = rows[endIdx + 1];
                    if (r.NormalStopCodons > 0)
                    {
                        double ratio = (double)r.TumorStopCodons / r.NormalStopCodons;
                        if (ratio < 0.5 || ratio > 2.0)
                            endIdx++;
                        else break;
                    }
                    else break;
                }

                return (startIdx, endIdx);
            }

            // Pass 2: fallback to UCSC centromere coordinates with proportional mapping
            var centroRegion = CentromereRegions.GetRegion(chr);
            if (centroRegion == null) return null;

            var lastRow = rows[rows.Count - 1];
            long totalFiltered = lastRow.StartPos + lastRow.Length;
            if (totalFiltered <= 0) return null;

            long centroStart = centroRegion.Value.Start;
            long centroEnd = centroRegion.Value.End;

            // Approximate: rawLen ≈ totalFiltered * 1.6 for intergenic
            double rawLenEstimate = totalFiltered * 1.6;
            double mapRatio = totalFiltered / rawLenEstimate;

            long filteredCentroStart = (long)(centroStart * mapRatio);
            long filteredCentroEnd = (long)(centroEnd * mapRatio);

            // Generous margin (30% of centromere size)
            long margin = (long)((filteredCentroEnd - filteredCentroStart) * 0.3);
            filteredCentroStart = Math.Max(0, filteredCentroStart - margin);
            filteredCentroEnd = Math.Min(totalFiltered, filteredCentroEnd + margin);

            int sIdx = -1, eIdx = -1;
            for (int i = 0; i < rows.Count; i++)
            {
                long wStart = rows[i].StartPos;
                long wEnd = wStart + rows[i].Length;

                if (wStart < filteredCentroEnd && wEnd > filteredCentroStart)
                {
                    if (sIdx < 0) sIdx = i;
                    eIdx = i;
                }
            }

            if (sIdx >= 0 && eIdx >= sIdx)
                return (sIdx, eIdx);

            return null;
        }

        /// <summary>
        /// Compute centromere window ranges for the genome-wide chart,
        /// offset by each chromosome's starting window index.
        /// </summary>
        private static List<(int StartIndex, int EndIndex)> ComputeGenomeWideCentromereRanges(
            List<WindowRow> allRows, List<string> chrNames)
        {
            var ranges = new List<(int StartIndex, int EndIndex)>();
            int offset = 0;

            foreach (var chr in chrNames.OrderBy(NaturalChrOrderKey))
            {
                var chrRows = allRows.Where(r => r.Chr == chr).ToList();
                if (chrRows.Count == 0) continue;

                var centroRange = FindCentromereWindowRange(chrRows);
                if (centroRange.HasValue)
                {
                    ranges.Add((offset + centroRange.Value.StartIndex,
                                offset + centroRange.Value.EndIndex));
                }

                offset += chrRows.Count;
            }

            return ranges;
        }

        // ───────── Data classes ─────────

        private class WindowInfo
        {
            public int StartPos;
            public int Length;
            public string Classification;
        }

        private class WindowRow
        {
            public string Chr;
            public int WindowIndex;
            public int StartPos;
            public int Length;
            public string NormalClass;
            public string TumorClass;
            public long NormalStopCodons;
            public long TumorStopCodons;
            public long Delta;
            public int RawMidpoint; // only used for WholeGenome sorting
        }

        private class ChrFeatureResult
        {
            public List<WindowRow> Rows;
            public int Breakpoints; // -1 = no alignment built
            public int FinalOffset;
        }
    }
}
