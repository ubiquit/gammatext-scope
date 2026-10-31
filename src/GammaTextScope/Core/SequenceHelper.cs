using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace GammaTextScope
{
    /// <summary>
    /// Gene-annotation loading and the γword window scanner shared by classify and the experiments.
    /// </summary>
    public static class SequenceHelper
    {
        /// <summary>
        /// Reads the protein-coding genes (or, for <see cref="FeatureType.Exon"/>, their exons) from
        /// every GFF3 file in <paramref name="geneInfoDir"/>, the <c>&lt;sample&gt;_GeneInfo</c> folder
        /// written by <c>prepare-annotation</c>. GFF3 coordinates are 1-based and inclusive; both start
        /// and end are decremented by one. The index is keyed by the line's first attribute (ID).
        /// </summary>
        public static Dictionary<string, List<GeneEntry>> ReadGeneIndexFromStream(string geneInfoDir, FeatureType featureType)
        {
            var geneIndex = new Dictionary<string, List<GeneEntry>>();

            // the unmasked whole genome needs no annotation; do not require the folder
            if (featureType == FeatureType.WholeGenome)
                return geneIndex;

            if (!Directory.Exists(geneInfoDir))
                throw new DirectoryNotFoundException($"{geneInfoDir}: gene annotation folder not found (run prepare-annotation first, or check --data-root)");
            string[] paths = Directory.GetFiles(geneInfoDir);

            // exons are attributed to their parent gene, so the gene index is needed first
            var geneIndexForExons = featureType == FeatureType.Exon ? ReadGeneIndexFromStream(geneInfoDir, FeatureType.Gene) : null;

            foreach (string path in paths)
            {
                using var reader = new StreamReader(path);
                while (!reader.EndOfStream)
                {
                    var gene = ReadGeneFeatureLine(reader.ReadLine(), featureType, geneIndexForExons);
                    if (gene != null)
                        UpdateGeneIndex(gene, geneIndex);
                }
            }

            if (geneIndex.Count == 0)
                Console.Error.WriteLine($"  WARNING: no protein-coding {featureType} features read from {geneInfoDir}; the {featureType} mask will be empty");

            return geneIndex;
        }

        // ───────── Window scanner ─────────

        /// <summary>
        /// Returns the tokens of the next window: the shortest prefix of the stop-codon
        /// tokenisation of sequence[lastStartPosition..] that has more than length*retry
        /// tokens, ends in a non-empty token and extends past lastFinishIndex. If no such
        /// prefix exists the whole remainder is returned with isEnd = true.
        ///
        /// Incremental form of <see cref="FindNextLengthWindowReference"/>: it tokenises only as
        /// far as the window needs instead of splitting the entire remainder of the chromosome
        /// for every window (which makes the reference scanner quadratic in chromosome length). The tokenisation rule is that of
        /// String.Split(new[]{"TAG","TAA","TGA"}, StringSplitOptions.None): scan left to
        /// right, split at the first position where a stop codon starts, resume after it.
        /// </summary>
        public static List<string> FindNextLengthWindow(string sequence, int lastStartPosition, ref int lastFinishIndex, int length, out bool isEnd, ref int retry, int windowWidth)
        {
            var tokens = new List<string>();
            int minIndex = length * retry;
            int tokenStart = lastStartPosition;
            int n = sequence.Length;
            int i = lastStartPosition;

            while (i <= n - 3)
            {
                char a = sequence[i], b = sequence[i + 1], c = sequence[i + 2];
                bool isStop = a == 'T' && ((b == 'A' && (c == 'A' || c == 'G')) || (b == 'G' && c == 'A'));
                if (!isStop) { i++; continue; }

                string token = sequence.Substring(tokenStart, i - tokenStart);
                tokens.Add(token);
                int index = tokens.Count - 1;
                if (index >= minIndex && index > lastFinishIndex && token.Length > 0)
                {
                    lastFinishIndex += length; // callers reset it after every window
                    isEnd = false;
                    return tokens;
                }

                i += 3;
                tokenStart = i;
            }

            // end of sequence: the trailing token is everything after the last stop codon
            tokens.Add(sequence.Substring(tokenStart));
            int last = tokens.Count - 1;
            if (last >= minIndex && last > lastFinishIndex && tokens[last].Length > 0)
            {
                lastFinishIndex += length;
                isEnd = false;
                return tokens;
            }

            isEnd = true;
            return tokens;
        }

        /// <summary>The original scanner, kept only for the equivalence test in WindowScannerTests.</summary>
        public static List<string> FindNextLengthWindowReference(string sequence, int lastStartPosition, ref int lastFinishIndex, int length, out bool isEnd, ref int retry, int windowWidth)
        {
            string endString = sequence.Substring(lastStartPosition);
            var tokenizedProcessingString = Tokenizer.TokenizeOnStopCodons(endString, GenomeConstants.stopWords, false);
            isEnd = false;

            for (int i = 0; i < tokenizedProcessingString.Count; i++)
            {
                if (i < length * retry)
                    continue;

                if (length * retry > tokenizedProcessingString.Count)
                    break;

                if (i <= lastFinishIndex)
                    continue;

                string word = tokenizedProcessingString[i];
                if (word.Equals(string.Empty))
                    continue;

                var tokenizedList = tokenizedProcessingString.GetRange(0, i + 1);
                lastFinishIndex += length;

                return tokenizedList;
            }

            //for unqualified end string
            isEnd = true;
            return tokenizedProcessingString;
        }

        /// <summary>
        /// Number of nucleotides a token list covers: every token plus 3 for the stop codon that
        /// ended it (empty tokens stand for adjacent stop codons). The window that reaches the end
        /// of the chromosome (<paramref name="isEnd"/>) has no closing stop codon, so 3 is taken off
        /// again. Convention pinned by expected/.
        /// </summary>
        public static int GetTokenizedSequenceLength(List<string> tokenized, bool isEnd)
        {
            int length = 0;
            int emptyCounter = 1; // emptyCounter is to compensate for areas where there are several stop codons in a row
            for (int i = 0; i < tokenized.Count; i++)
            {
                string w = tokenized[i];

                if (w == string.Empty)
                {
                    emptyCounter++;

                    if (i == tokenized.Count - 1)
                        length += (emptyCounter * 3); // addition here to add possible empty window lengths
                }
                else
                {
                    length += (w.Length + (emptyCounter * 3)); // + 3 to compensate for lost stop codon length
                    emptyCounter = 1;
                }
            }

            if (isEnd)
                length -= 3; //subtract here if this is the last window. there is no end element

            return length;
        }

        public static List<string> ConvertTokenizedWithEptyPositions(List<string> tokenizedwithEmptyPositions)
        {
            List<string> tokenized = new List<string>();

            foreach (var w in tokenizedwithEmptyPositions)
                if (w != string.Empty)
                    tokenized.Add(w);

            return tokenized;
        }

        public static Dictionary<string, int> CreateWordCounts(List<string> tokenizedSequence)
        {
            Dictionary<string, int> countDict = new Dictionary<string, int>();

            foreach (string word in tokenizedSequence)
                if (countDict.ContainsKey(word))
                    countDict[word]++;
                else
                    countDict.Add(word, 1);

            return countDict;
        }

        // ───────── GFF3 parsing ─────────

        private static GeneEntry ReadGeneFeatureLine(string lineIn, FeatureType featureType, Dictionary<string, List<GeneEntry>> geneIndex)
        {
            string[] line = lineIn.Split('\t');

            if (line.Length <= 1)
                return null;

            return featureType == FeatureType.Exon ? CreateExonFromLine(line, geneIndex) : CreateGeneFromLine(line);
        }

        private static GeneEntry CreateGeneFromLine(string[] line)
        {
            int startPos = -1;
            int endPos = -1;

            if (line.Length > 8 && line[2].Equals("gene") && int.TryParse(line[3], out startPos) && int.TryParse(line[4], out endPos))
            {
                var geneData = line[8].Split(';');
                var geneType = AttributeValue(geneData, "gene_biotype") ?? AttributeValue(geneData, "gene_type") ?? AttributeValue(geneData, "biotype");
                var transcriptType = line[1];
                string id = geneData[0].Contains('=') ? geneData[0].Split('=')[1] : null;

                if (id != null && geneType != null && geneType.Equals("protein_coding"))
                {
                    return new GeneEntry
                    {
                        ID = id,
                        GeneName = id,
                        GeneType = geneType,
                        TranscriptType = transcriptType,
                        StartPosition = startPos - 1,
                        EndPosition = endPos - 1,
                        Chromosome = line[0]
                    };
                }
            }

            return null;
        }

        /// <summary>Value of the first GFF attribute starting with <paramref name="key"/>, or null if absent or malformed.</summary>
        private static string AttributeValue(string[] attributes, string key)
        {
            var a = attributes.FirstOrDefault(x => x.StartsWith(key));
            if (a == null) return null;
            int eq = a.IndexOf('=');
            return eq < 0 ? null : a.Substring(eq + 1);
        }

        private static GeneEntry CreateExonFromLine(string[] line, Dictionary<string, List<GeneEntry>> workingGeneIndex)
        {
            int startPos = -1;
            int endPos = -1;

            if (line.Length > 2 && line[2].Equals("exon") && int.TryParse(line[3], out startPos) && int.TryParse(line[4], out endPos))
            {
                var geneData = line[8].Split(';');
                var geneDataField = geneData.FirstOrDefault(x => x.StartsWith("gene="));
                string parentName = null;
                if (geneDataField != null)
                {
                    parentName = geneDataField.Split("=")[1];
                }

                // exons of some annotations name their gene through a locus tag
                string locus_tag = null;
                if (geneData.Any(x => x.StartsWith("locus_tag")))
                    locus_tag = geneData.FirstOrDefault(x => x.StartsWith("locus_tag")).Split("=")[1];

                var transcriptType = line[1];

                bool hasLocusTagGene = locus_tag == null ? false : workingGeneIndex.Any(x => x.Value.Any(y => y.GeneName.Equals("gene-" + locus_tag)));

                if (parentName != null || (workingGeneIndex.ContainsKey("gene-" + parentName) || hasLocusTagGene))
                {
                    try
                    {
                        var gene = hasLocusTagGene ? workingGeneIndex["gene-" + locus_tag].First() : workingGeneIndex["gene-" + parentName].First();
                        return new GeneEntry
                        {
                            ID = geneData[0].Split('=')[1],
                            GeneName = geneData[0].Split('=')[1],
                            GeneType = gene.GeneType,
                            TranscriptType = transcriptType,
                            StartPosition = startPos - 1,
                            EndPosition = endPos - 1,
                            Chromosome = line[0]
                        };
                    }
                    catch (Exception)
                    {
                        return null;
                    }
                }

                // Ensembl format: exons have Parent=transcript:... instead of gene=...
                // Match exon to a gene by checking if it falls within a gene's position range on the same chromosome
                if (parentName == null && !hasLocusTagGene)
                {
                    string chromosome = line[0];
                    int exonStart = startPos - 1;
                    int exonEnd = endPos - 1;

                    foreach (var kvp in workingGeneIndex)
                    {
                        var gene = kvp.Value.FirstOrDefault(g => g.Chromosome == chromosome && exonStart >= g.StartPosition && exonEnd <= g.EndPosition);
                        if (gene != null)
                        {
                            return new GeneEntry
                            {
                                ID = geneData[0].Split('=')[1],
                                GeneName = geneData[0].Split('=')[1],
                                GeneType = gene.GeneType,
                                TranscriptType = transcriptType,
                                StartPosition = exonStart,
                                EndPosition = exonEnd,
                                Chromosome = chromosome
                            };
                        }
                    }
                }
            }

            return null;
        }

        private static void UpdateGeneIndex(GeneEntry entry, Dictionary<string, List<GeneEntry>> index)
        {
            if (index.TryGetValue(entry.GeneName, out var entries))
                entries.Add(entry);
            else
                index.Add(entry.GeneName, new List<GeneEntry> { entry });
        }
    }
}
