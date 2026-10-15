using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using GammaTextScope;

[assembly: CollectionBehavior(DisableTestParallelization = true)]

namespace GammaTextScope.Tests
{
    /// <summary>
    /// A reference dataset in expected/: v1 is what the submitted manuscript was built
    /// from, v2 is the regenerated dataset (breast cancer v2, melanoma v3 inputs).
    /// </summary>
    public sealed record Dataset(string Version, string BreastNormal, string BreastTumor, string MelanomaNormal, string MelanomaTumor,
                                 string BreastDensityPair, string BreastAlphaPair, string MelanomaDensityPair, string MelanomaAlphaPair)
    {
        public static readonly Dataset V1 = new("v1",
            "HomoSapiens_BreastCancer_Normal", "HomoSapiens_BreastCancer_Tumor", "HomoSapiens_Melanoma_Normal", "HomoSapiens_Melanoma_Tumor",
            "BreastCancer_Normal_vs_Tumor", "BreastCancer", "Melanoma_Normal_vs_Tumor", "Melanoma");
        public static readonly Dataset V2 = new("v2",
            "HomoSapiens_BreastCancer_Normal_v2", "HomoSapiens_BreastCancer_Tumor_v2", "HomoSapiens_Melanoma_Normal_v3", "HomoSapiens_Melanoma_Tumor_v3",
            "BreastCancer_v2_Normal_vs_Tumor", "BreastCancer_v2", "Melanoma_v3_Normal_vs_Tumor", "Melanoma_v3");
        public static Dataset Get(string version) => version == "v1" ? V1 : V2;

        public string[] Samples(int cancer) => cancer == 0 ? new[] { BreastNormal, BreastTumor } : new[] { MelanomaNormal, MelanomaTumor };
        public string DensityPair(int cancer) => cancer == 0 ? BreastDensityPair : MelanomaDensityPair;
        public string AlphaPair(int cancer) => cancer == 0 ? BreastAlphaPair : MelanomaAlphaPair;
        public string Expected(string relative) => TestData.Expected(Version, relative);
        public bool SequencesAvailable => Directory.Exists(Path.Combine(TestData.DataRoot, BreastNormal)) && Directory.Exists(Path.Combine(TestData.DataRoot, MelanomaTumor + "_GeneInfo"));
    }

    public static class TestData
    {
        public static string ExpectedRoot { get; } =
            Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", "..", "expected"));

        public static string DataRoot { get; } = new PipelineConfig().DataRoot;

        public static string Expected(string version, string relative) =>
            Path.Combine(ExpectedRoot, version, relative.Replace('/', Path.DirectorySeparatorChar));

        /// <summary>Reference outputs of the submitted manuscript (dataset v1).</summary>
        public static string ExpectedV1(string relative) => Expected("v1", relative);

        /// <summary>Fresh output root for one test, with the expected classifications staged for the given samples.</summary>
        public static string StageOutRoot(Dataset ds, string name, IEnumerable<string> samples, params FeatureType[] modes)
        {
            string root = Path.Combine(Path.GetTempPath(), "gts-tests", ds.Version, name);
            if (Directory.Exists(root)) Directory.Delete(root, recursive: true);
            foreach (var s in samples)
                foreach (var m in modes)
                {
                    string rel = $"{s}/Experiment_{s}_{m}_7000/WindowWidth1";
                    Directory.CreateDirectory(Path.Combine(root, rel));
                    foreach (var f in Directory.GetFiles(ds.Expected(rel), "*_window_positions_indexLength.txt"))
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
