using System;
using System.Collections.Generic;
using System.Linq;

namespace GammaTextScope
{
    internal static class Program
    {
        private const string Usage = @"gammatext-scope — stop-codon γtext comparison of tumour and matched-normal genomes

usage: gammatext-scope <command> [options]

commands
  prepare-annotation  --sample S                         relabel <S>_AccessionProcessing/genomic.gff into <S>_GeneInfo
  classify            --sample S --mode M                per-window Zipf/alfa classification (writes Experiment_<S>_<M>_7000/WindowWidth1)
  density             --normal N --tumor T --pair P --mode M   stop codons per 100 kb bin and per chromosome
  alpha-compare       --normal N --tumor T --pair P --mode M [--merged]   per-window stop-codon delta + qualification (DNA-scope data)
  decay-dfa           --pair P --mode M [--merged]       DFA of the qualification-decay marker gaps (needs alpha-compare)

options
  --mode M            Intergenic | Gene | Exon | WholeGenome
  --merged            with --mode WholeGenome: use the union of the Intergenic and Gene classifications
                      (legacy behaviour, outputs *_WholeGenomeMerged) instead of whole-genome classifications
  --data-root DIR     input root (default D:\Files\Stream or $GAMMATEXT_DATA)
  --out-root DIR      output root (default D:\Files\Harvests or $GAMMATEXT_OUT)
  --chromosomes LIST  restrict to chromosomes, e.g. 21,22,X
  --words N           γwords per window for classify (default 7000)
";

        private static int Main(string[] args)
        {
            if (args.Length == 0 || args[0] == "-h" || args[0] == "--help")
            {
                Console.Write(Usage);
                return args.Length == 0 ? 1 : 0;
            }

            var opts = ParseOptions(args.Skip(1));
            var cfg = new PipelineConfig
            {
                DataRoot = opts.GetValueOrDefault("data-root") ?? new PipelineConfig().DataRoot,
                OutRoot = opts.GetValueOrDefault("out-root") ?? new PipelineConfig().OutRoot,
                Chromosomes = PipelineConfig.ParseChromosomes(opts.GetValueOrDefault("chromosomes")),
                WordsPerWindow = int.Parse(opts.GetValueOrDefault("words") ?? "7000"),
            };
            PipelineConfig.Default = cfg;

            try
            {
                switch (args[0])
                {
                    case "prepare-annotation":
                        new AnnotationCleaner().PrepareAnnotationFile(Require(opts, "sample"));
                        return 0;

                    case "classify":
                        new WindowClassifier().RunSimplified(Require(opts, "sample"), cfg.WindowWidth, cfg.WordsPerWindow,
                            Mode(opts), cfg.ClassificationDir(Require(opts, "sample"), Mode(opts)));
                        return 0;

                    case "density":
                        new StopCodonDensity().Run(Require(opts, "normal"), Require(opts, "tumor"), Require(opts, "pair"), Mode(opts));
                        return 0;

                    case "alpha-compare":
                        new StopCodonAlphaComparison().Run(Require(opts, "normal"), Require(opts, "tumor"), Require(opts, "pair"), Mode(opts), opts.ContainsKey("merged"));
                        return 0;

                    case "decay-dfa":
                        new DecayGapDfa().Run(Require(opts, "pair"), Mode(opts).ToString() + (opts.ContainsKey("merged") ? "Merged" : ""));
                        return 0;

                    default:
                        Console.Error.WriteLine($"unknown command '{args[0]}'");
                        Console.Write(Usage);
                        return 1;
                }
            }
            catch (ArgumentException ex)
            {
                Console.Error.WriteLine(ex.Message);
                return 2;
            }
        }

        private static FeatureType Mode(Dictionary<string, string> opts) =>
            Enum.Parse<FeatureType>(Require(opts, "mode"), ignoreCase: true);

        private static string Require(Dictionary<string, string> opts, string key) =>
            opts.TryGetValue(key, out var v) ? v : throw new ArgumentException($"missing --{key}");

        private static Dictionary<string, string> ParseOptions(IEnumerable<string> args)
        {
            var result = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            string key = null;
            foreach (var a in args)
            {
                if (a.StartsWith("--"))
                {
                    key = a.Substring(2);
                    result[key] = "true";
                }
                else if (key != null)
                {
                    result[key] = a;
                    key = null;
                }
            }
            return result;
        }
    }
}
