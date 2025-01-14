import pytest

from orf_finding import coordinates, orf_classes


# sys.path.insert(0, os.path.abspath(
#     os.path.join(os.path.dirname(__file__), '..')))


def test_GffInterval_properties():
    gff_interval = coordinates.GffInterval(5, 10)
    assert gff_interval.python_start == 4
    assert gff_interval.python_end == 10


def test_VeliaInterval_properties():
    velia_interval = coordinates.VeliaInterval(5, 10)
    assert velia_interval.python_start == 4
    assert velia_interval.python_end == 10


def test_PythonInterval_properties():
    python_interval = coordinates.PythonInterval(5, 10)
    assert python_interval.python_start == 5
    assert python_interval.python_end == 10


def test_chrom_starts_and_block_sizes_to_interval_objects():
    chrom_starts = [1, 5]
    block_sizes = [2, 3]
    interval_type = 'velia'
    intervals = coordinates.convert_chrom_starts_and_block_sizes_to_interval_objects(
        chrom_starts, block_sizes, interval_type)
    assert len(intervals) == 2
    assert intervals[0].python_start == 0
    assert intervals[0].python_end == 3
    assert intervals[1].python_start == 4
    assert intervals[1].python_end == 8


def test_chrom_starts_and_block_sizes_to_interval_objects_invalid_interval_type():
    chrom_starts = [1, 5]
    block_sizes = [2, 3]
    interval_type = 'invalid'
    with pytest.raises(ValueError):
        coordinates.convert_chrom_starts_and_block_sizes_to_interval_objects(
            chrom_starts, block_sizes, interval_type)


def test_IntervalFactory_make_interval():
    interval_type = 'gff'
    interval_class = coordinates.IntervalFactory.make_interval(interval_type)
    assert interval_class == coordinates.GffInterval

    interval_type = 'velia'
    interval_class = coordinates.IntervalFactory.make_interval(interval_type)
    assert interval_class == coordinates.VeliaInterval

    interval_type = 'python'
    interval_class = coordinates.IntervalFactory.make_interval(interval_type)
    assert interval_class == coordinates.PythonInterval


def test_IntervalFactory_make_interval_invalid_type():
    interval_type = 'invalid'
    with pytest.raises(ValueError):
        coordinates.IntervalFactory.make_interval(interval_type)


def test_convert_chrom_starts_and_block_sizes_to_intervals():
    chrom_starts = [1, 5]
    block_sizes = [2, 3]
    source_type = 'velia'
    destination_type = 'python'
    intervals = coordinates.convert_chrom_starts_and_block_sizes_to_intervals(
        chrom_starts, block_sizes, source_type, destination_type)
    assert len(intervals) == 2
    assert intervals[0][0] == 0
    assert intervals[0][1] == 3
    assert intervals[1][0] == 4
    assert intervals[1][1] == 8


def test_convert_chrom_starts_and_block_sizes_to_intervals_invalid_source_type():
    chrom_starts = [1, 5]
    block_sizes = [2, 3]
    source_type = 'invalid'
    destination_type = 'velia'
    with pytest.raises(ValueError):
        coordinates.convert_chrom_starts_and_block_sizes_to_intervals(
            chrom_starts, block_sizes, source_type, destination_type)


def test_convert_chrom_starts_and_block_sizes_to_intervals_invalid_destination_type():
    chrom_starts = [1, 5]
    block_sizes = [2, 3]
    source_type = 'gff'
    destination_type = 'invalid'
    with pytest.raises(ValueError):
        coordinates.convert_chrom_starts_and_block_sizes_to_intervals(
            chrom_starts, block_sizes, source_type, destination_type)


def test_convert_gff_chrom_starts_and_block_sizes():
    chrom_starts = [1, 5]
    block_sizes = [2, 3]
    source_type = 'gff'
    destination_type = 'python'

    converted_chrom_starts, converted_block_sizes = coordinates.convert_chrom_starts_and_block_sizes(chrom_starts=chrom_starts,
                                                                                                     block_sizes=block_sizes,
                                                                                                     source_type=source_type,
                                                                                                     destination_type=destination_type)
    assert converted_chrom_starts == [0, 4]
    assert converted_block_sizes == [2, 3]


def test_convert_velia_chrom_starts_and_block_sizes():
    chrom_starts = [1, 5]
    block_sizes = [2, 3]
    source_type = 'velia'
    destination_type = 'python'

    converted_chrom_starts, converted_block_sizes = coordinates.convert_chrom_starts_and_block_sizes(chrom_starts=chrom_starts,
                                                                                                     block_sizes=block_sizes,
                                                                                                     source_type=source_type,
                                                                                                     destination_type=destination_type)
    assert converted_chrom_starts == [0, 4]
    assert converted_block_sizes == [3, 4]


def test_compute_orf_genomic_start_real():
    mock_orf = orf_classes.OrfBase(
        1, 418, '+', ['CTG'] + ['AAA'] * (203 // 3 - 2) + ['TAA'])

    mock_transcript = {'assembly.genbank_accession': 'CM000669.2', 'assembly.genome_accession': 'GCA_000001405.28_GRCh38.p13_assembly_report', 'assembly.id': 7, 'assembly.sequence_length': 159345973, 'assembly.ucsc_style_name': 'chr7', 'transcript.block_sizes': '155;81;110;72;126;488', 'transcript.chrom_starts':
                       '127588411;127589083;127589485;127590066;127590963;127591213', 'transcript.end': 127591700, 'transcript.ensembl_id': 'ENST00000000233.10', 'transcript.id': 1, 'transcript.refseq_id': 'rna-NM_001662.4', 'transcript.start': 127588411, 'transcript.strand': '+', 'transcript.transcript_type': 'protein_coding'}

    computed_genomic_start = coordinates.compute_orf_genomic_start(
        mock_orf, mock_transcript)
    assert computed_genomic_start == 127590962


def test_compute_orf_genomic_start_contrived():
    mock_orf = orf_classes.OrfBase(
        'test+', 6, '+', ['ATG', 'AAA', 'CCC', 'GGG', 'TTT'])

    mock_transcript = {'transcript.chrom_starts': '1;19;114;120',
                       'transcript.block_sizes': '6;4;4;100'}

    computed_genomic_start = coordinates.compute_orf_genomic_start(
        mock_orf, mock_transcript)
    assert computed_genomic_start == 18


def test_compute_phases_pos_strand():
    """
    Validating against phases listed for ENST00000246314.10
    """
    py_chrom_starts = [36004384, 36008689, 36008896, 36009474, 36013629, 36013914,
                       36027113, 36034173, 36036176, 36039789, 36040306, 36043446, 36054945, 36055636]
    py_block_sizes = [91, 88, 148, 120, 123, 134,
                      185, 160, 91, 195, 135, 102, 200, 109]
    strand = '+'
    seq_length = 248956422
    computed_ensembl_phases, frames = coordinates.compute_cds_phases_and_frames(py_chrom_starts=py_chrom_starts,
                                                                                py_block_sizes=py_block_sizes,
                                                                                strand=strand,
                                                                                chrom_length=seq_length,
                                                                                phase_style='ensembl')
    assert computed_ensembl_phases == [
        0, 1, 2, 0, 0, 0, 2, 1, 2, 0, 0, 0, 0, 2]
    computed_gencode_phases, frames = coordinates.compute_cds_phases_and_frames(py_chrom_starts=py_chrom_starts,
                                                                                py_block_sizes=py_block_sizes,
                                                                                strand=strand,
                                                                                chrom_length=seq_length,
                                                                                phase_style='gencode')

    assert computed_gencode_phases == [
        0, 2, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 1]


def test_compute_phases_neg_strand():
    """
    Validating against phases listed for ENST00000237305.11
    """
    py_chrom_starts = [134170267, 134170825, 134171022, 134171636, 134172192,
                       134172661, 134173022, 134173273, 134173461, 134174004, 134174510, 134174731]
    py_block_sizes = [168, 90, 156, 96, 124, 113, 132, 84, 105, 76, 76, 76]
    strand = '-'
    seq_length = 170805979
    ensembl_phases, frames = coordinates.compute_cds_phases_and_frames(py_chrom_starts=py_chrom_starts,
                                                                       py_block_sizes=py_block_sizes,
                                                                       strand=strand,
                                                                       chrom_length=seq_length,
                                                                       phase_style='ensembl')

    gencode_phases, frames = coordinates.compute_cds_phases_and_frames(py_chrom_starts=py_chrom_starts,
                                                                       py_block_sizes=py_block_sizes,
                                                                       strand=strand,
                                                                       chrom_length=seq_length,
                                                                       phase_style='gencode')
    assert ensembl_phases == [0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 1, 0]
    assert gencode_phases == [0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 2, 0]


def test_compute_phases_pos_strand2():
    """
    Uses data from https://useast.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000185946;r=1:103525691-103553343;t=ENST00000524631
    """
    py_chrom_starts = [103526070,
                       103527694,
                       103533738,
                       103534773,
                       103535329,
                       103536125,
                       103537341,
                       103541349,
                       103543298,
                       103544940,
                       103546247,
                       103546976,
                       103550940,
                       103551720]
    py_block_sizes = [192, 48, 119, 84, 112,
                      69, 143, 126, 149, 162, 95, 59, 133, 60]
    ensembl_phases = [0, 0, 0, 2, 2, 0, 0, 2, 2, 1, 1, 0, 2, 0]
    gencode_phases = [0, 0, 0, 1, 1, 0, 0, 1, 1, 2, 2, 0, 1, 0]
    computed_ensembl_phases, computed_frames = coordinates.compute_cds_phases_and_frames(py_chrom_starts=py_chrom_starts,
                                                                                         py_block_sizes=py_block_sizes,
                                                                                         strand='+',
                                                                                         chrom_length=248956422,
                                                                                         phase_style='ensembl')
    assert computed_ensembl_phases == ensembl_phases, computed_ensembl_phases

    computed_gencode_phases, computed_frames = coordinates.compute_cds_phases_and_frames(py_chrom_starts=py_chrom_starts,
                                                                                         py_block_sizes=py_block_sizes,
                                                                                         strand='+',
                                                                                         chrom_length=248956422,
                                                                                         phase_style='gencode')
    assert computed_gencode_phases == gencode_phases, computed_gencode_phases


def test_compute_phases_neg_strand2():
    """
    Uses data from https://useast.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000174876;r=1:103687415-103696210;t=ENST00000330330
    """
    py_chrom_starts = [
        103687448,
        103688957,
        103689193,
        103691293,
        103691487,
        103692426,
        103693299,
        103693978,
        103694979,
        103695471
    ]
    py_block_sizes = [190, 126, 119, 100, 123, 134, 231, 198, 147, 168]
    ensembl_phases = [2, 2, 0, 2, 2, 0, 0, 0, 0, 0]

    computed_phases, computed_frames = coordinates.compute_cds_phases_and_frames(py_chrom_starts=py_chrom_starts,
                                                                                 py_block_sizes=py_block_sizes,
                                                                                 strand='-',
                                                                                 chrom_length=248956422,
                                                                                 phase_style='ensembl')
    assert computed_phases == ensembl_phases
