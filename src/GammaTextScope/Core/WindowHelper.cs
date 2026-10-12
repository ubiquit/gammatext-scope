using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GammaTextScope
{
    public static class WindowHelper
    {
        public static List<Tuple<string, int, string>> ReadStreamFromFileIntragenic(string path, string fileName, ref Dictionary<string, List<EnsembleGeneEntry>> geneIndex, int windowSize)
        {
            double numberOfWindows = 0;
            double numberOfLines = 0;
            int numberOfSymbols = 0;
            int prevSymbolPosition = 0;
            int symbolPosition = 0;

            //Return list
            List<Tuple<string, int, string>> windowList = new List<Tuple<string, int, string>>();

            //GFF gene intergenic experiment
            string chromosome = fileName.Split('_').Last().Split('.').First();

            //calculating totals
            using (var reader = new StreamReader(path))
            {
                while (!reader.EndOfStream)
                {
                    string line = reader.ReadLine();
                    numberOfSymbols += line.Length;
                    numberOfLines++;
                }
            }

            numberOfWindows = numberOfLines / windowSize;

            using (var reader = new StreamReader(path))
            {
                StringBuilder sb = new StringBuilder();
                int lineNum = 1;
                int windownum = 1;

                while (!reader.EndOfStream)
                {
                    string line = reader.ReadLine();

                    //first line throwaway in GFF
                    if (line.StartsWith(">"))
                        line = reader.ReadLine();

                    sb.Append(line);
                    List<string> lineList = line.Select(c => c.ToString()).ToList();
                    lineNum++;
                    symbolPosition += line.Length;

                    if (lineNum >= windowSize)
                    {
                        //GFF gene intergenic experiment
                        var geneDistances = SequenceHelper.SequenceInverter(geneIndex, sb.ToString(), prevSymbolPosition, symbolPosition, chromosome);
                        string intergenicText = string.Join(string.Empty, geneDistances);
                        windowList.Add(new Tuple<string, int, string>(chromosome, windownum, intergenicText));

                        //reset objects and counters
                        sb = new StringBuilder();
                        windownum++;
                        lineNum = 1;
                        prevSymbolPosition = symbolPosition;
                    }
                }

                //last window (less than window width)
                if (sb.Length >= 30)
                {
                    var geneDistances = SequenceHelper.SequenceInverter(geneIndex, sb.ToString(), prevSymbolPosition, symbolPosition, chromosome);
                    string intergenicText = string.Join(string.Empty, geneDistances);

                    windowList.Add(new Tuple<string, int, string>(chromosome, windownum, intergenicText));
                }
            }

            int intergenicLength = windowList.Sum(x => x.Item3.Length);

            return windowList;
        }

        public static string ReadSequenceFromFile(string path, string fileName, ref Dictionary<string, List<EnsembleGeneEntry>> geneIndex, FeatureType featureType)
        {
            string sequence;
            switch (featureType)
            {
                case FeatureType.Gene:
                    sequence = ReadGeneticSequenceFromFile(path, fileName, ref geneIndex);
                    break;
                case FeatureType.Exon:
                    sequence = ReadGeneticSequenceFromFile(path, fileName, ref geneIndex);
                    break;
                case FeatureType.Intergenic:
                    sequence = ReadIntergenicSequenceFromFile(path, fileName, ref geneIndex);
                    break;
                case FeatureType.WholeGenome:
                    sequence = ReadWholeSequenceFromFile(path, fileName);
                    break;
                default:
                    sequence = ReadGeneticSequenceFromFile(path, fileName, ref geneIndex);
                    break;
            }

            return sequence;
        }

        private static string ReadGeneticSequenceFromFile(string path, string fileName, ref Dictionary<string, List<EnsembleGeneEntry>> geneIndex)
        {
            string chromosome = fileName.Split('_').Last().Split('.').First();
            var geneList = GetChrGenes(chromosome, ref geneIndex);
            long index = 1;
            long lengthOfGenes = 0;
            long lengthOfGeneticMaterial = 0;
            long numberOfProcessedGenes = 0;

            foreach (var g in geneList) lengthOfGenes += (g.EndPosition - g.StartPosition);

            StringBuilder sb = new StringBuilder();
            using (var reader = new StreamReader(path))
            {
                EnsembleGeneEntry currentGene = FindNextGene(index, geneList);
                numberOfProcessedGenes++;

                while (!reader.EndOfStream)
                {
                    string line = reader.ReadLine();
                    char[] characters = line.ToCharArray();
                    for (int i = 0; i < characters.Length; i++)
                    {
                        index++;

                        if (currentGene == null)
                            continue;

                        if (index >= currentGene.StartPosition && index <= currentGene.EndPosition)
                        {
                            lengthOfGeneticMaterial++;
                            sb.Append(characters[i]);
                        }

                        if (index >= currentGene.EndPosition)
                        {
                            currentGene = FindNextGene(index, geneList);
                            if (currentGene != null)
                                numberOfProcessedGenes++;
                        }
                    }
                }
            }

            return sb.ToString();
        }

        private static string ReadIntergenicSequenceFromFile(string path, string fileName, ref Dictionary<string, List<EnsembleGeneEntry>> geneIndex)
        {
            string chromosome = fileName.Split('_').Last().Split('.').First();
            var geneList = GetChrGenes(chromosome, ref geneIndex);
            long index = 1;
            long lengthOfGenes = 0;
            long lengthOfGeneticMaterial = 0;
            long numberOfProcessedGenes = 0;

            foreach (var g in geneList) lengthOfGenes += (g.EndPosition - g.StartPosition);

            StringBuilder sb = new StringBuilder();
            using (var reader = new StreamReader(path))
            {
                EnsembleGeneEntry currentGene = FindNextGene(index, geneList);
                numberOfProcessedGenes++;

                while (!reader.EndOfStream)
                {
                    string line = reader.ReadLine();
                    char[] characters = line.ToCharArray();
                    for (int i = 0; i < characters.Length; i++)
                    {
                        index++;

                        if (currentGene == null)
                        {
                            sb.Append(characters[i]);
                            continue;
                        }

                        if (index >= currentGene.StartPosition && index <= currentGene.EndPosition)
                        {
                            lengthOfGeneticMaterial++;
                            continue;
                        }
                        else
                            sb.Append(characters[i]);

                        if (index >= currentGene.EndPosition)
                        {
                            currentGene = FindNextGene(index, geneList);
                            if (currentGene != null)
                                numberOfProcessedGenes++;
                        }
                    }
                }
            }

            return sb.ToString();
        }

        private static string ReadWholeSequenceFromFile(string path, string fileName)
        {
            string chromosome = fileName.Split('_').Last().Split('.').First();
            long index = 1;

            StringBuilder sb = new StringBuilder();
            using (var reader = new StreamReader(path))
            {
                while (!reader.EndOfStream)
                {
                    string line = reader.ReadLine();
                    char[] characters = line.ToCharArray();
                    for (int i = 0; i < characters.Length; i++)
                    {
                        index++;

                        sb.Append(characters[i]);
                    }
                }
            }

            return sb.ToString();
        }

        public static string ReadStreamFromFileIntragenic(string path, string fileName, ref Dictionary<string, List<EnsembleGeneEntry>> geneIndex)
        {
            string chromosome = fileName.Split('_').Last().Split('.').First();
            var geneList = GetChrGenes(chromosome, ref geneIndex);
            long index = 1;
            long lengthOfGenes = 0;
            long lengthOfGeneticMaterial = 0;
            long numberOfProcessedGenes = 0;

            foreach (var g in geneList) lengthOfGenes += (g.EndPosition - g.StartPosition);

            StringBuilder sb = new StringBuilder();
            using (var reader = new StreamReader(path))
            {
                EnsembleGeneEntry currentGene = FindNextGene(index, geneList);
                numberOfProcessedGenes++;

                while (!reader.EndOfStream)
                {
                    string line = reader.ReadLine();
                    char[] characters = line.ToCharArray();
                    for (int i = 0; i < characters.Length; i++)
                    {
                        index++;

                        if (currentGene == null)
                        {
                            sb.Append(characters[i]);
                            continue;
                        }

                        if (index >= currentGene.StartPosition && index <= currentGene.EndPosition)
                        {
                            lengthOfGeneticMaterial++;
                            continue;
                        }
                        else
                            sb.Append(characters[i]);

                        if (index >= currentGene.EndPosition)
                        {
                            currentGene = FindNextGene(index, geneList);
                            if (currentGene != null)
                                numberOfProcessedGenes++;
                        }
                    }
                }
            }

            return sb.ToString();
        }

        private static List<EnsembleGeneEntry> GetChrGenes(string chromosome, ref Dictionary<string, List<EnsembleGeneEntry>> geneIndex)
        {
            List<EnsembleGeneEntry> retlist = new List<EnsembleGeneEntry>();
            foreach (var gene in geneIndex)
                foreach (var occ in gene.Value)
                    if (occ.Chromosome.Equals(chromosome))
                        retlist.Add(occ);

            retlist = retlist.OrderBy(x => x.StartPosition).ThenBy(x => x.EndPosition).ToList();

            return retlist;
        }

        private static EnsembleGeneEntry FindNextGene(long index, List<EnsembleGeneEntry> chrGenes)
        {
            if (chrGenes.Count == 0)
                return null;

            chrGenes.RemoveAt(0);

            if (chrGenes.Count == 0)
                return null;

            var gene = chrGenes.First();

            return gene;
        }

        public static List<Tuple<string, int, int>> GetGeneRelativePositions(string pathStream, string chromosomeName, ref Dictionary<string, List<EnsembleGeneEntry>> geneIndex, bool isIntergenic)
        {
            int index = 1;
            List<Tuple<string, int, int>> retGenePosList = new List<Tuple<string, int, int>>();

            foreach (var g in geneIndex)
            {
                var geneEntry = g.Value.First();
                int genLen = geneEntry.EndPosition - geneEntry.StartPosition;
                retGenePosList.Add(new Tuple<string, int, int>(geneEntry.ID, index, index + genLen));
                index += genLen;
            }

            return retGenePosList;
        }
    }
}
