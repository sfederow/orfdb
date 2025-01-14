import Bio.Seq

from orf_finding import sequence_tools, orf_classes


def test_split_sequence_into_codons():
    # Test normal case
    result = sequence_tools.split_sequence_into_codons("ATGCGT", 0)
    assert result == ["ATG", "CGT"]

    # Test edge cases
    result = sequence_tools.split_sequence_into_codons("ATGCGT", 1)
    assert result == ["TGC"]

    result = sequence_tools.split_sequence_into_codons("ATGCGT", 2)
    assert result == ["GCG"]


def test_extract_orf_start_context_positive_strand():
    """
    Test the extraction of ORF start context on the positive strand.
    """
    orf_object = orf_classes.OrfBase(
        parent_sequence_id='chrom1', start_pos=6, strand='+', codons=['ATG', 'CAT', 'CAT'])

    parent_transcript_seq = 'GAGAGAATGCATCATCATATATATAT'

    upstream_context, downstream_context = sequence_tools.extract_orf_start_context(orf_object=orf_object,
                                                                                    parent_transcript_seq=parent_transcript_seq)

    assert str(
        upstream_context) == 'GAGA', "The extracted ORF upstream start context on the positive strand does not match expected."
    assert str(
        downstream_context) == 'C', "The extracted ORF downstream start context on the positive strand does not match expected."


def test_extract_orf_start_context_negative_strand():
    """
    Test the extraction of ORF start context on the negative strand.
    """
    orf_object = orf_classes.OrfBase(
        parent_sequence_id='chrom1', start_pos=8, strand='-', codons=['ATG', 'ATG', 'ATG'])

    parent_transcript_seq = Bio.Seq.reverse_complement(
        'GAGAGAATGCATCATCATATATATAT')

    upstream_context, downstream_context = sequence_tools.extract_orf_start_context(orf_object=orf_object,
                                                                                    parent_transcript_seq=parent_transcript_seq)
    assert str(
        upstream_context) == 'ATAT', "The extracted ORF upstream start context on the negative strand does not match expected."
    assert str(
        downstream_context) == 'A', "The extracted ORF downstream start context on the negative strand does not match expected."
