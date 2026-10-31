using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace GammaTextScope
{
    /// <summary>
    /// Reads one chromosome file and returns the text a feature type keeps: the protein-coding
    /// gene stretches (Gene, Exon), everything outside them (Intergenic) or the whole sequence.
    ///
    /// Conventions pinned by the reference outputs in expected/ and replicated by
    /// FilteredToRawMapper and figures/make_figures.py — do not "fix":
    ///  * the running index starts at 1 and is incremented before each character is tested,
    ///    so a gene (Start, End) covers the 0-based file positions Start-2 .. End-2;
    ///  * the first gene of a chromosome (lowest start) is skipped, because the gene list is
    ///    advanced once before the first character is read;
    ///  * in Intergenic mode the gene list only advances one character after a gene has ended,
    ///    so a gene that starts exactly at the previous gene's End+1 keeps that one character.
    /// </summary>
    public static class WindowHelper
    {
        public static string ReadSequenceFromFile(string path, string fileName, Dictionary<string, List<GeneEntry>> geneIndex, FeatureType featureType)
        {
            switch (featureType)
            {
                case FeatureType.Intergenic:
                    return ReadIntergenicSequenceFromFile(path, fileName, geneIndex);
                case FeatureType.WholeGenome:
                    return ReadWholeSequenceFromFile(path);
                default: // Gene, Exon
                    return ReadGeneticSequenceFromFile(path, fileName, geneIndex);
            }
        }

        private static string ReadGeneticSequenceFromFile(string path, string fileName, Dictionary<string, List<GeneEntry>> geneIndex)
        {
            var geneList = GetChrGenes(ChromosomeOf(fileName), geneIndex);
            long index = 1;
            var sb = new StringBuilder();

            using var reader = new StreamReader(path);
            GeneEntry currentGene = FindNextGene(geneList); // skips the first gene (see the class summary)

            while (!reader.EndOfStream)
            {
                foreach (char c in reader.ReadLine())
                {
                    index++;

                    if (currentGene == null)
                        continue;

                    if (index >= currentGene.StartPosition && index <= currentGene.EndPosition)
                        sb.Append(c);

                    if (index >= currentGene.EndPosition)
                        currentGene = FindNextGene(geneList);
                }
            }

            return sb.ToString();
        }

        private static string ReadIntergenicSequenceFromFile(string path, string fileName, Dictionary<string, List<GeneEntry>> geneIndex)
        {
            var geneList = GetChrGenes(ChromosomeOf(fileName), geneIndex);
            long index = 1;
            var sb = new StringBuilder();

            using var reader = new StreamReader(path);
            GeneEntry currentGene = FindNextGene(geneList); // skips the first gene (see the class summary)

            while (!reader.EndOfStream)
            {
                foreach (char c in reader.ReadLine())
                {
                    index++;

                    if (currentGene == null)
                    {
                        sb.Append(c);
                        continue;
                    }

                    if (index >= currentGene.StartPosition && index <= currentGene.EndPosition)
                        continue; // inside a gene: dropped; the gene list is not advanced here (see the class summary)

                    sb.Append(c);

                    if (index >= currentGene.EndPosition)
                        currentGene = FindNextGene(geneList);
                }
            }

            return sb.ToString();
        }

        private static string ReadWholeSequenceFromFile(string path)
        {
            var sb = new StringBuilder();
            foreach (var line in File.ReadLines(path))
                sb.Append(line);
            return sb.ToString();
        }

        /// <summary>Chromosome name from a sequence file name such as HomoSapiens_BreastCancer_Normal_chr21.txt.</summary>
        private static string ChromosomeOf(string fileName) => fileName.Split('_').Last().Split('.').First();

        /// <summary>The genes of one chromosome, ordered by start and then end position (a fresh list; the index is not modified).</summary>
        private static List<GeneEntry> GetChrGenes(string chromosome, Dictionary<string, List<GeneEntry>> geneIndex)
        {
            var retlist = new List<GeneEntry>();
            foreach (var gene in geneIndex)
                foreach (var occ in gene.Value)
                    if (occ.Chromosome.Equals(chromosome))
                        retlist.Add(occ);

            return retlist.OrderBy(x => x.StartPosition).ThenBy(x => x.EndPosition).ToList();
        }

        /// <summary>Removes the head of the list and returns the new head (null when the list is exhausted).</summary>
        private static GeneEntry FindNextGene(List<GeneEntry> chrGenes)
        {
            if (chrGenes.Count == 0)
                return null;

            chrGenes.RemoveAt(0);

            return chrGenes.Count == 0 ? null : chrGenes[0];
        }
    }
}
