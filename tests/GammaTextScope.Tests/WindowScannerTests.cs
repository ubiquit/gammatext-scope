using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using GammaTextScope;

namespace GammaTextScope.Tests
{
    /// <summary>
    /// The incremental window scanner must return exactly what the reference
    /// split-the-whole-remainder scanner returned, window after window.
    /// </summary>
    public class WindowScannerTests
    {
        private static string RandomSequence(int seed, int length, double stopDensity)
        {
            var rng = new Random(seed);
            var sb = new StringBuilder(length);
            string[] stops = { "TAG", "TAA", "TGA" };
            while (sb.Length < length)
            {
                if (rng.NextDouble() < stopDensity) sb.Append(stops[rng.Next(3)]);
                else sb.Append("ACGTN"[rng.Next(5)]);
            }
            return sb.ToString(0, length);
        }

        private static void AssertWindowsIdentical(string sequence, int wordsPerWindow)
        {
            int startA = 0, startB = 0, finishA = 0, finishB = 0, retryA = 1, retryB = 1, windows = 0;
            while (startA < sequence.Length)
            {
                var a = SequenceHelper.FindNextLengthWindow(sequence, startA, ref finishA, wordsPerWindow, out bool endA, ref retryA, 1);
                var b = SequenceHelper.FindNextLengthWindowReference(sequence, startB, ref finishB, wordsPerWindow, out bool endB, ref retryB, 1);
                Assert.Equal(b, a);
                Assert.Equal(endB, endA);
                Assert.Equal(finishB, finishA);
                int len = SequenceHelper.GetTokenizedSequenceLength(a, endA);
                startA += len; startB += len;
                finishA = 0; finishB = 0;
                windows++;
                if (endA) break;
            }
            Assert.True(windows >= 2, "test sequence too short to form several windows");
        }

        [Theory]
        [InlineData(1, 200_000, 0.05, 300)]
        [InlineData(2, 200_000, 0.30, 300)]   // many adjacent stop codons -> many empty tokens
        [InlineData(3, 50_000, 0.02, 300)]    // few tokens: exercises the end-of-sequence path
        [InlineData(4, 120_000, 0.08, 7000)]  // fewer tokens than one window
        public void IncrementalScannerMatchesReferenceOnSyntheticSequences(int seed, int length, double density, int words)
        {
            AssertWindowsIdentical(RandomSequence(seed, length, density), words);
        }

        [Fact]
        public void IncrementalScannerHandlesStopCodonAtVeryEnd()
        {
            AssertWindowsIdentical(RandomSequence(5, 20_000, 0.1) + "TAA", 50);
            AssertWindowsIdentical("TAATAGTGA" + RandomSequence(6, 20_000, 0.1), 50);
        }

        [SkippableFact]
        public void IncrementalScannerMatchesReferenceOnRealChromosomePrefix()
        {
            Skip.IfNot(Dataset.V1.SequencesAvailable, "v1 sequence data not available under " + TestData.DataRoot);
            // 4 Mb from the middle of chr21 — the start of the chromosome is an unassembled N run without stop codons
            string file = Directory.GetFiles(Path.Combine(TestData.DataRoot, "HomoSapiens_BreastCancer_Normal"), "chr21.*").Single();
            var sb = new StringBuilder();
            foreach (var line in File.ReadLines(file)) { sb.Append(line); if (sb.Length > 24_000_000) break; }
            AssertWindowsIdentical(sb.ToString(20_000_000, 4_000_000), 7000);
        }
    }
}
