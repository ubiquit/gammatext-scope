using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;

namespace GammaTextScope
{
    public class AnnotationCleaner
    {

        public void PrepareAnnotationFile(string species)
        {
            string filenameSpecies = species.Split("\\").Last();
            string AnnotationChromosomesPath = Path.Combine(PipelineConfig.Default.StreamDir(species + "_AccessionProcessing"), filenameSpecies + "_Chromosomes.json");
            string AnnotationFilesPath = Path.Combine(PipelineConfig.Default.StreamDir(species + "_AccessionProcessing"), "genomic.gff");
            string PrintPath = Path.Combine(PipelineConfig.Default.StreamDir(species + "_GeneInfo"), filenameSpecies + "_Annotation_Converted.txt");
            Directory.CreateDirectory(Path.GetDirectoryName(PrintPath));

            string jsonString = File.ReadAllText(AnnotationChromosomesPath);
            var dsList = JsonSerializer.Deserialize<AnnotationItem>(jsonString);
            StringBuilder sb = new StringBuilder();
            int originalCount = 0;
            int newCount = 0;
            HashSet<string> chromosomes = new HashSet<string>();
            File.WriteAllText(PrintPath, "");

            HashSet<string> annotations = new HashSet<string>();
            using (var reader = new StreamReader(AnnotationFilesPath))
            {
                while (!reader.EndOfStream)
                {
                    string line = reader.ReadLine();
                    var values = line.Split("\t");
                    originalCount++;

                    //for beggining of file
                    if (values[0].StartsWith("#"))
                    {
                        sb.Append(line + GenomeConstants.OutputNewLine);
                        continue;
                    }
                    var acc = dsList.seq.FirstOrDefault(x => x.acc.Equals(values[0]));
                    annotations.Add(values[0]);

                    //some accessions are missing
                    if (acc == null)
                        continue;

                    string chromosome = "chr" + acc.name;

                    if (chromosome != null)
                    {
                        string[] newValues = new string[values.Length];
                        newValues[0] = chromosome;
                        for (int i = 1; i < values.Length; i++)
                            newValues[i] = values[i];

                        newCount++;
                        sb.Append(CreatePrintLine(newValues));
                        chromosomes.Add(chromosome);
                    }

                    //reset this since big files causes out of memory
                    if (originalCount % 50000 == 0)
                    {
                        File.AppendAllText(PrintPath, sb.ToString());
                        sb = new StringBuilder();
                    }
                }
            }

            File.AppendAllText(PrintPath, sb.ToString());
        }

        private string CreatePrintLine(string[] values)
        {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < values.Length; i++)
            {
                sb.Append(values[i]);
                if (i != values.Length - 1)
                    sb.Append("\t");
            }

            sb.Append(GenomeConstants.OutputNewLine);

            return sb.ToString();
        }
    }
}

