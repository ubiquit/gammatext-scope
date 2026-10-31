namespace GammaTextScope
{
    /// <summary>
    /// One gene (or exon) feature from the GFF3 annotation. Positions are the GFF3
    /// coordinates minus one (see SequenceHelper.ReadGeneIndexFromStream); how they are
    /// applied to the sequence is described in WindowHelper.
    /// </summary>
    public class GeneEntry
    {
        public string ID { get; set; }
        public string GeneName { get; set; }
        public string Chromosome { get; set; }
        public string GeneType { get; set; }
        public string TranscriptType { get; set; }
        public int StartPosition { get; set; }
        public int EndPosition { get; set; }
    }
}
