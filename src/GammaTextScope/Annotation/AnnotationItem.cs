using System.Collections.Generic;

namespace GammaTextScope
{
    /// <summary>Shape of &lt;sample&gt;_Chromosomes.json: maps sequence accessions to chromosome names.</summary>
    public class AnnotationItem
    {
        public string name { get; set; }
        public string acc { get; set; }
        public int taxid { get; set; }
        public string type { get; set; }
        public List<AccessionItem> seq { get; set; }
    }

    public class AccessionItem
    {
        public string name { get; set; }
        public int chrom { get; set; }
        public int topl { get; set; }
        public long length { get; set; }
        public string acc { get; set; }
        public int? gi { get; set; }
    }
}
