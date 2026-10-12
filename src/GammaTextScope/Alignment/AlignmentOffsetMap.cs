using System;
using System.Collections.Generic;

namespace GammaTextScope
{
    /// <summary>
    /// Scans two nearly-identical sequences (Normal vs Tumor) and builds a
    /// cumulative-offset map that tracks indels.  Given a position in the
    /// Normal sequence the map returns the corresponding position in the
    /// Tumor sequence, so that the same biological region is compared.
    ///
    /// Algorithm:
    ///   Walk both sequences in lock-step using an anchor window (e.g. 32 bp).
    ///   When the anchor in Normal no longer matches at the expected Tumor
    ///   position, do a small local search (±MaxSearchRadius) to find where
    ///   it went.  The shift that re-aligns them is recorded as a new
    ///   offset breakpoint.  Between breakpoints the offset is constant.
    /// </summary>
    public class AlignmentOffsetMap
    {
        private const int AnchorLen = 32;       // length of the anchor k-mer
        private const int StepSize = 500;       // check every N bp (speed vs precision)
        private const int MaxSearchRadius = 200; // max indel size we can detect at once

        // Re-anchoring parameters: at each segment boundary, do a fresh wide search
        // to prevent accumulated drift across the chromosome.
        private const int SegmentSize = 10_000_000;       // independent alignment every 10M bp
        private const int ReanchorLen = 64;                // longer anchor for re-anchoring (more specific)
        private const int ReanchorSearchRadius = 5000;     // wider search at segment boundary
        private const int ReanchorProbes = 5;              // number of probe positions to vote on offset

        // Sorted list of (normalPos, cumulativeOffset) breakpoints.
        // Between breakpoints the offset is constant.
        // tumorPos = normalPos + offset
        private readonly List<(int NormalPos, int Offset)> _breakpoints = new List<(int, int)>();

        /// <summary>
        /// Build the offset map by scanning both cleaned sequences.
        /// The chromosome is divided into independent segments (default 10M bp).
        /// At each segment boundary, a fresh wide-radius re-anchor search establishes
        /// the correct offset from scratch, preventing accumulated drift.
        /// Within each segment, fine-grained scanning (32bp anchors every 500bp)
        /// tracks small local indels as before.
        /// </summary>
        public static AlignmentOffsetMap Build(string normalSeq, string tumorSeq)
        {
            var map = new AlignmentOffsetMap();
            int maxNormal = normalSeq.Length - AnchorLen;
            int tumorLen = tumorSeq.Length;

            // Process chromosome in independent segments
            for (int segStart = 0; segStart < maxNormal; segStart += SegmentSize)
            {
                int segEnd = Math.Min(segStart + SegmentSize, maxNormal);

                // --- Re-anchor at segment start ---
                // Use a fresh wide search to find the true offset, independent of
                // any accumulated drift from previous segments.
                int segOffset = FindOffsetAtPosition(normalSeq, tumorSeq, segStart, tumorLen);

                // Record the segment start offset
                if (map._breakpoints.Count == 0 || map._breakpoints[map._breakpoints.Count - 1].Offset != segOffset)
                    map._breakpoints.Add((segStart, segOffset));

                // --- Fine-grained scanning within the segment ---
                int currentOffset = segOffset;

                for (int nPos = segStart + StepSize; nPos < segEnd; nPos += StepSize)
                {
                    int expectedTPos = nPos + currentOffset;

                    // Check if anchor still matches at the expected position
                    if (expectedTPos >= 0 && expectedTPos + AnchorLen <= tumorLen &&
                        SequenceMatch(normalSeq, nPos, tumorSeq, expectedTPos, AnchorLen))
                    {
                        continue; // still aligned
                    }

                    // Lost alignment — search nearby for the anchor
                    bool found = false;
                    for (int delta = 1; delta <= MaxSearchRadius; delta++)
                    {
                        // Try shift right (tumor has insertion)
                        int tryPos = expectedTPos + delta;
                        if (tryPos >= 0 && tryPos + AnchorLen <= tumorLen &&
                            SequenceMatch(normalSeq, nPos, tumorSeq, tryPos, AnchorLen))
                        {
                            currentOffset += delta;
                            map._breakpoints.Add((nPos, currentOffset));
                            found = true;
                            break;
                        }

                        // Try shift left (tumor has deletion)
                        tryPos = expectedTPos - delta;
                        if (tryPos >= 0 && tryPos + AnchorLen <= tumorLen &&
                            SequenceMatch(normalSeq, nPos, tumorSeq, tryPos, AnchorLen))
                        {
                            currentOffset -= delta;
                            map._breakpoints.Add((nPos, currentOffset));
                            found = true;
                            break;
                        }
                    }

                    // If not found, keep current offset and try again at next step.
                }
            }

            // Ensure at least one breakpoint exists
            if (map._breakpoints.Count == 0)
                map._breakpoints.Add((0, 0));

            return map;
        }

        /// <summary>
        /// Find the correct offset at a given Normal position using a fresh wide-radius search.
        /// Uses multiple probe positions around the target and longer anchors for reliability.
        /// Returns the consensus offset via majority vote among successful probes.
        /// </summary>
        private static int FindOffsetAtPosition(string normalSeq, string tumorSeq, int normalPos, int tumorLen)
        {
            // Try multiple probe positions near normalPos to get a robust offset estimate.
            // Each probe uses a longer anchor (64bp) and wider search radius (5000bp).
            var offsets = new List<int>();
            int probeSpacing = StepSize * 4; // 2000bp between probes

            for (int p = 0; p < ReanchorProbes; p++)
            {
                int probePos = normalPos + p * probeSpacing;
                if (probePos + ReanchorLen > normalSeq.Length) break;

                // First check offset = 0 (identity)
                if (probePos + ReanchorLen <= tumorLen &&
                    SequenceMatch(normalSeq, probePos, tumorSeq, probePos, ReanchorLen))
                {
                    offsets.Add(0);
                    continue;
                }

                // Search outward from identity position
                bool found = false;
                for (int delta = 1; delta <= ReanchorSearchRadius; delta++)
                {
                    int tryPos = probePos + delta;
                    if (tryPos >= 0 && tryPos + ReanchorLen <= tumorLen &&
                        SequenceMatch(normalSeq, probePos, tumorSeq, tryPos, ReanchorLen))
                    {
                        offsets.Add(delta);
                        found = true;
                        break;
                    }

                    tryPos = probePos - delta;
                    if (tryPos >= 0 && tryPos + ReanchorLen <= tumorLen &&
                        SequenceMatch(normalSeq, probePos, tumorSeq, tryPos, ReanchorLen))
                    {
                        offsets.Add(-delta);
                        found = true;
                        break;
                    }
                }
            }

            if (offsets.Count == 0)
                return 0; // can't determine — assume identity

            // Return median offset (robust against outliers)
            offsets.Sort();
            return offsets[offsets.Count / 2];
        }

        /// <summary>
        /// Given a position in the Normal sequence, return the corresponding
        /// position in the Tumor sequence.
        /// </summary>
        public int NormalToTumor(int normalPos)
        {
            return normalPos + OffsetAt(normalPos);
        }

        /// <summary>
        /// Translate a Normal range [start, start+length) to the corresponding
        /// Tumor range. Returns (tumorStart, tumorLength).
        ///
        /// When the offset drift within the range is small (a few hundred bp from
        /// real indels), both endpoints are translated independently so the tumor
        /// length accurately reflects insertions/deletions.
        ///
        /// When the drift is large (e.g. centromeric deserts where alignment is
        /// lost), the start position is translated but the original length is
        /// preserved. This prevents pathological cases where a misaligned offset
        /// produces a tumor range many times larger than the normal range.
        /// </summary>
        public (int TumorStart, int TumorLength) TranslateRange(int normalStart, int normalLength)
        {
            int startOffset = OffsetAt(normalStart);
            int endOffset = OffsetAt(normalStart + normalLength);
            int drift = Math.Abs(endOffset - startOffset);

            int tumorStart = normalStart + startOffset;

            // Fast path: offset is stable across the range — use accurate end translation
            const int MaxDrift = 5000;
            if (drift <= MaxDrift)
            {
                int tumorEnd = normalStart + normalLength + endOffset;
                return (tumorStart, tumorEnd - tumorStart);
            }

            // Large drift: alignment is unreliable over this range (e.g. centromere,
            // large structural variant). Use the start offset with the original length
            // so we stay in the correct neighbourhood and count a comparable amount of DNA.
            return (tumorStart, normalLength);
        }

        /// <summary>
        /// Get the cumulative offset at a given Normal position (without adding
        /// the position itself, unlike NormalToTumor which returns pos + offset).
        /// </summary>
        private int OffsetAt(int normalPos)
        {
            int lo = 0, hi = _breakpoints.Count - 1;
            while (lo < hi)
            {
                int mid = (lo + hi + 1) / 2;
                if (_breakpoints[mid].NormalPos <= normalPos)
                    lo = mid;
                else
                    hi = mid - 1;
            }
            return _breakpoints[lo].Offset;
        }

        /// <summary>
        /// Number of offset breakpoints detected (= number of indel events).
        /// </summary>
        public int BreakpointCount => _breakpoints.Count;

        /// <summary>
        /// Final cumulative offset at the end of the sequences.
        /// </summary>
        public int FinalOffset => _breakpoints.Count > 0 ? _breakpoints[_breakpoints.Count - 1].Offset : 0;

        private static bool SequenceMatch(string a, int aStart, string b, int bStart, int length)
        {
            for (int i = 0; i < length; i++)
            {
                char ca = a[aStart + i];
                char cb = b[bStart + i];
                // Skip N positions — they shouldn't cause a false mismatch
                if (ca == 'N' || cb == 'N') continue;
                if (ca != cb) return false;
            }
            return true;
        }
    }
}
