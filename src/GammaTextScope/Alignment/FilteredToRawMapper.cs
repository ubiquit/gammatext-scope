using System;
using System.Collections.Generic;
using System.Linq;

namespace GammaTextScope
{
    /// <summary>
    /// Maps positions in a feature-filtered sequence (Intergenic or Gene) back to
    /// positions in the whole-genome (unfiltered) sequence for the same chromosome.
    ///
    /// Uses sparse intervals computed from gene boundaries. Each interval represents
    /// a contiguous run of filtered characters that maps to a contiguous run of raw
    /// characters, enabling O(log N) position lookups via binary search.
    ///
    /// Replicates the exact gene-consumption behaviour of WindowHelper:
    ///   - index starts at 1, incremented before use (first char = index 2)
    ///   - Gene included when: index >= gene.StartPosition && index <= gene.EndPosition
    ///   - FindNextGene always skips (removes) the first entry in the sorted list
    ///   - At init, FindNextGene is called once, which discards the first gene
    /// </summary>
    public class FilteredToRawMapper
    {
        // Each interval: filtered positions [FilteredStart, FilteredStart+Length)
        // correspond to raw whole-genome positions [RawStart, RawStart+Length).
        private readonly List<Interval> _intervals;
        private readonly int _totalFilteredLength;

        private FilteredToRawMapper(List<Interval> intervals, int totalLen)
        {
            _intervals = intervals;
            _totalFilteredLength = totalLen;
        }

        public int TotalFilteredLength => _totalFilteredLength;

        /// <summary>
        /// Build mapper from gene boundaries.
        /// </summary>
        /// <param name="sortedGenes">Genes for this chromosome, sorted by StartPosition (will be cloned).</param>
        /// <param name="chrRawCharCount">Total character count in the raw chromosome file (excluding newlines).</param>
        /// <param name="isGene">true = Gene feature type, false = Intergenic feature type.</param>
        public static FilteredToRawMapper Build(List<GeneEntry> sortedGenes, long chrRawCharCount, bool isGene)
        {
            // Replicate initial FindNextGene: removes first gene, returns second
            var genes = new List<GeneEntry>(sortedGenes);
            if (genes.Count > 0)
                genes.RemoveAt(0);

            // Merge overlapping gene regions into non-overlapping sorted intervals
            // Gene coordinates are in "index space" (index = raw_0based_pos + 2)
            var merged = MergeGeneRegions(genes);

            // index range for the whole file: [2, chrRawCharCount + 1]
            // rawPos (0-based in whole-genome string) = index - 2
            var intervals = new List<Interval>();
            int filteredPos = 0;

            if (isGene)
            {
                // Gene: include characters within gene regions
                foreach (var (gStart, gEnd) in merged)
                {
                    int rawStart = (int)(gStart - 2);
                    int length = (int)(gEnd - gStart + 1);

                    // Clamp to valid range
                    if (rawStart < 0) { length += rawStart; rawStart = 0; }
                    if (rawStart + length > chrRawCharCount) length = (int)(chrRawCharCount - rawStart);
                    if (length <= 0) continue;

                    intervals.Add(new Interval { FilteredStart = filteredPos, RawStart = rawStart, Length = length });
                    filteredPos += length;
                }
            }
            else
            {
                // Intergenic: include characters between/outside gene regions
                long prevEndIndex = 1; // index before first character

                foreach (var (gStart, gEnd) in merged)
                {
                    // Gap before this gene: index [prevEndIndex+1, gStart-1]
                    long gapStartIdx = prevEndIndex + 1;
                    long gapEndIdx = gStart - 1;

                    if (gapEndIdx >= gapStartIdx)
                    {
                        int rawStart = (int)(gapStartIdx - 2);
                        int length = (int)(gapEndIdx - gapStartIdx + 1);

                        if (rawStart < 0) { length += rawStart; rawStart = 0; }
                        if (length > 0)
                        {
                            intervals.Add(new Interval { FilteredStart = filteredPos, RawStart = rawStart, Length = length });
                            filteredPos += length;
                        }
                    }

                    prevEndIndex = gEnd;
                }

                // Tail after last gene: index [prevEndIndex+1, chrRawCharCount+1]
                long tailStartIdx = prevEndIndex + 1;
                long tailEndIdx = chrRawCharCount + 1;

                if (tailEndIdx >= tailStartIdx)
                {
                    int rawStart = (int)(tailStartIdx - 2);
                    int length = (int)(tailEndIdx - tailStartIdx + 1);

                    if (rawStart < 0) { length += rawStart; rawStart = 0; }
                    if (rawStart + length > chrRawCharCount) length = (int)(chrRawCharCount - rawStart);
                    if (length > 0)
                    {
                        intervals.Add(new Interval { FilteredStart = filteredPos, RawStart = rawStart, Length = length });
                        filteredPos += length;
                    }
                }
            }

            return new FilteredToRawMapper(intervals, filteredPos);
        }

        /// <summary>
        /// Convert a filtered-sequence position to a raw whole-genome position.
        /// </summary>
        public int ToRaw(int filteredPos)
        {
            if (_intervals.Count == 0) return filteredPos; // fallback

            // Binary search for the interval containing filteredPos
            int lo = 0, hi = _intervals.Count - 1;
            while (lo < hi)
            {
                int mid = (lo + hi + 1) / 2;
                if (_intervals[mid].FilteredStart <= filteredPos)
                    lo = mid;
                else
                    hi = mid - 1;
            }

            var iv = _intervals[lo];
            return iv.RawStart + (filteredPos - iv.FilteredStart);
        }

        /// <summary>
        /// Get the raw midpoint for a filtered window, useful for sorting windows
        /// by genomic position.
        /// </summary>
        public int GetRawMidpoint(int filteredStart, int filteredLength)
        {
            return ToRaw(filteredStart + filteredLength / 2);
        }

        // ───────── internals ─────────

        private static List<(long Start, long End)> MergeGeneRegions(List<GeneEntry> genes)
        {
            if (genes.Count == 0) return new List<(long, long)>();

            var sorted = genes.OrderBy(g => g.StartPosition).ThenBy(g => g.EndPosition).ToList();
            var merged = new List<(long Start, long End)>();

            long curStart = sorted[0].StartPosition;
            long curEnd = sorted[0].EndPosition;

            for (int i = 1; i < sorted.Count; i++)
            {
                if (sorted[i].StartPosition <= curEnd + 1)
                {
                    curEnd = Math.Max(curEnd, sorted[i].EndPosition);
                }
                else
                {
                    merged.Add((curStart, curEnd));
                    curStart = sorted[i].StartPosition;
                    curEnd = sorted[i].EndPosition;
                }
            }
            merged.Add((curStart, curEnd));

            return merged;
        }

        private struct Interval
        {
            public int FilteredStart;
            public int RawStart;
            public int Length;
        }
    }
}
