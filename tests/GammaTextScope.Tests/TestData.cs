using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using GammaTextScope;

[assembly: CollectionBehavior(DisableTestParallelization = true)]

namespace GammaTextScope.Tests
{
    /// <summary>Locations of the reference outputs (expected/) and, when present, the sequence data.</summary>
    public static class TestData
    {
        public static string ExpectedDir { get; } =
            Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", "..", "expected"));

        public static string DataRoot { get; } = new PipelineConfig().DataRoot;

        /// <summary>True when the four human samples' chromosome files are available locally.</summary>
        public static bool SequencesAvailable =>
            Directory.Exists(Path.Combine(DataRoot, "HomoSapiens_BreastCancer_Normal")) &&
            Directory.Exists(Path.Combine(DataRoot, "HomoSapiens_Melanoma_Tumor_GeneInfo"));

        public static string Expected(string relative) => Path.Combine(ExpectedDir, relative.Replace('/', Path.DirectorySeparatorChar));

        /// <summary>Fresh output root for one test, with the expected classifications staged for the given samples.</summary>
        public static string StageOutRoot(string name, IEnumerable<string> samples, params FeatureType[] modes)
        {
            string root = Path.Combine(Path.GetTempPath(), "gts-tests", name);
            if (Directory.Exists(root)) Directory.Delete(root, recursive: true);
            foreach (var s in samples)
                foreach (var m in modes)
                {
                    string rel = $"{s}/Experiment_{s}_{m}_7000/WindowWidth1";
                    Directory.CreateDirectory(Path.Combine(root, rel));
                    foreach (var f in Directory.GetFiles(Expected(rel), "*_window_positions_indexLength.txt"))
                        File.Copy(f, Path.Combine(root, rel, Path.GetFileName(f)));
                }
            return root;
        }

        public static PipelineConfig UseOutRoot(string outRoot, string chromosomes = null)
        {
            var cfg = new PipelineConfig { OutRoot = outRoot, Chromosomes = PipelineConfig.ParseChromosomes(chromosomes) };
            PipelineConfig.Default = cfg;
            return cfg;
        }

        public static List<string> Lines(string path) => File.ReadAllLines(path).Where(l => l.Length > 0).ToList();

        public static List<string> ChromosomeRows(string path, string chrValue, int chrColumn = 0) =>
            Lines(path).Skip(1).Where(l => l.Split('\t')[chrColumn] == chrValue).ToList();

        /// <summary>Summary.txt rows keyed by category: (DotCount, DFA_Alpha).</summary>
        public static Dictionary<string, (int Dots, string Alpha)> DfaSummary(string path)
        {
            var result = new Dictionary<string, (int, string)>();
            foreach (var line in Lines(path).Skip(1))
            {
                var p = line.Split('\t');
                if (p.Length >= 7) result[p[0]] = (int.Parse(p[1]), p[6]);
            }
            return result;
        }
    }
}
