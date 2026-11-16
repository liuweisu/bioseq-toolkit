import pytest

from bioseq import align, assembly, bacteria, proteomics, seqio, simulate, te, variants
from bioseq.benchmark import score_intervals

# --- seqio -----------------------------------------------------------------

def test_fasta_roundtrip():
    recs = [seqio.Record("a", "ACGT" * 30, "first"), seqio.Record("b", "GGCC")]
    parsed = list(seqio.parse_fasta(seqio.format_fasta(recs).splitlines()))
    assert parsed == recs


def test_fastq_parse_and_quality():
    text = "@r1 x\nACGT\n+\nIIII\n@r2\nGG\n+\n++\n"
    recs = list(seqio.parse_fastq(text.splitlines()))
    assert [r.id for r in recs] == ["r1", "r2"]
    assert seqio.mean_quality(recs[0].quality) == 40


def test_reverse_complement_and_gc():
    assert seqio.reverse_complement("AACGTN") == "NACGTT"
    assert seqio.gc_content("GGCCAATT") == 0.5
    assert seqio.gc_content("NNNN") == 0.0


def test_translate():
    assert seqio.translate("ATGGCCTAAGGG") == "MA*G"
    assert seqio.translate("ATGGCCTAAGGG", to_stop=True) == "MA"


# --- alignment -------------------------------------------------------------

def test_global_alignment_identity():
    aln = align.needleman_wunsch("GATTACA", "GATTACA")
    assert aln.identity == 1.0 and aln.score == 7


def test_local_alignment_finds_embedded_match():
    aln = align.smith_waterman("TTTTTGATTACATTTTT", "CCCGATTACACCC")
    assert aln.a == "GATTACA" and aln.start_a == 5


# --- assembly --------------------------------------------------------------

def test_n50():
    stats = assembly.assembly_stats(["A" * 80, "A" * 70, "A" * 50, "A" * 40, "A" * 30, "A" * 20, "A" * 10])
    assert stats.total_length == 300 and stats.n50 == 70 and stats.l50 == 2


def test_ng50_uses_genome_size():
    stats = assembly.assembly_stats(["A" * 100, "A" * 50], genome_size=1000)
    assert stats.ng50 == 0  # the assembly covers less than half the genome


def test_assembler_reconstructs_linear_genome():
    genome = assembly.random_genome(5000, seed=1)
    reads = assembly.simulate_reads(genome, read_length=100, coverage=30, error_rate=0.005, seed=1)
    contigs = assembly.assemble([r.seq for r in reads], k=31, min_count=3, min_length=100)
    best = contigs[0].seq
    assert len(best) > 4800
    assert best in genome or seqio.reverse_complement(best) in genome


def test_assembler_closes_circular_genome():
    genome = assembly.random_genome(4000, seed=2)
    reads = assembly.simulate_reads(genome, coverage=30, error_rate=0.0, seed=2, circular=True)
    contigs = assembly.assemble([r.seq for r in reads], k=31, min_count=2)
    assert len(contigs) == 1 and "circular" in contigs[0].description
    assert len(contigs[0]) == len(genome)


def test_kmer_spectrum_separates_errors():
    genome = assembly.random_genome(3000, seed=3)
    reads = assembly.simulate_reads(genome, coverage=40, error_rate=0.01, seed=3)
    spectrum = assembly.kmer_spectrum(assembly.count_kmers([r.seq for r in reads], 21))
    assert spectrum[1] > spectrum.get(30, 0)  # many singleton error k-mers


# --- transposable elements -------------------------------------------------

@pytest.fixture(scope="module")
def plant():
    return simulate.plant_genome_with_tes(seed=7)


def test_ltr_finder_precision_and_recall(plant):
    pred = [(c.start, c.end) for c in te.find_ltr_retrotransposons(plant.sequence)]
    truth = [(t.start, t.end) for t in plant.truth if t.te_class.startswith("LTR")]
    score = score_intervals(pred, truth)
    assert score.precision >= 0.9 and score.recall >= 0.75


def test_ltr_age_increases_with_divergence():
    young = te.LTRCandidate(0, 5000, 400, 0.999, "AAAAA", True)
    old = te.LTRCandidate(0, 5000, 400, 0.95, "AAAAA", True)
    assert old.insertion_age() > young.insertion_age() > 0


def test_tir_finder_finds_some_elements(plant):
    pred = [(c.start, c.end) for c in te.find_tir_transposons(plant.sequence)]
    truth = [(t.start, t.end) for t in plant.truth if t.te_class.startswith("DNA")]
    assert score_intervals(pred, truth).true_positives >= 3


def test_repeat_mask_close_to_truth(plant):
    mask = te.kmer_repeat_mask(plant.sequence)
    true_fraction = sum(t.end - t.start for t in plant.truth) / len(plant.sequence)
    assert abs(mask.repeat_fraction - true_fraction) < 0.1


def test_classify_ltr_element(plant):
    t = next(t for t in plant.truth if t.family == "RLG_Retand")
    hit = te.classify(plant.sequence[t.start : t.end], plant.library)
    assert hit.family == "RLG_Retand" and hit.te_class == "LTR/Gypsy"


def test_kimura_distance():
    assert te.kimura_2p("ACGTACGTAC", "ACGTACGTAC") == 0
    assert te.kimura_2p("AAAAAAAAAA", "GAAAAAAAAA") > 0


# --- medical: TE insertion in a gene ---------------------------------------

def test_detects_alu_like_insertion_in_exon():
    case = simulate.human_gene_with_te_insertion()
    ins = variants.find_insertions(case.reference, case.patient)
    assert len(ins) == 1
    assert abs(ins[0].ref_position - case.true_position) <= 3
    assert any(ins[0].in_interval(s, e) for s, e in case.exons)
    hit = te.classify(ins[0].sequence, case.library)
    assert hit.family == "AluLike"
    assert te.poly_a_tail(ins[0].sequence)


def test_no_insertions_in_identical_sequences():
    case = simulate.human_gene_with_te_insertion()
    assert variants.find_insertions(case.reference, case.reference) == []


# --- bacteria --------------------------------------------------------------

def test_ori_prediction_and_dnaa_enrichment():
    bac = simulate.bacterial_genome(length=60_000, ori=20_000)
    pred = bacteria.predict_ori(bac.sequence)
    assert abs(pred.ori - bac.ori) < 500
    assert abs(pred.ter - bac.ter) < 3000
    assert pred.enrichment > 5


def test_orf_finder_recovers_planted_genes():
    bac = simulate.bacterial_genome(length=60_000, ori=20_000)
    orfs = bacteria.find_orfs(bac.sequence)
    score = score_intervals([(o.start, o.end) for o in orfs], [(g[0], g[1]) for g in bac.genes], 0.8)
    assert score.recall > 0.85


def test_orf_coordinates_on_minus_strand():
    cds = "ATG" + "GCT" * 120 + "TAA"
    seq = "CC" + seqio.reverse_complement(cds) + "GG"
    orf = bacteria.find_orfs(seq, min_protein_length=100)[0]
    assert orf.strand == "-" and seq[orf.start : orf.end] == seqio.reverse_complement(cds)


# --- proteomics ------------------------------------------------------------

def test_peptide_mass_matches_known_value():
    # Angiotensin II (DRVYIHPF): monoisotopic [M+H]+ 1046.5418
    assert proteomics.mz("DRVYIHPF") == pytest.approx(1046.5418, abs=0.002)


def test_trypsin_rules():
    peps = [p.sequence for p in proteomics.digest("AAAKPAAARAAAAKAAAAAA", missed_cleavages=0, min_length=1)]
    assert peps == ["AAAKPAAAR", "AAAAK", "AAAAAA"]  # no cleavage before proline


def test_isoelectric_point_ordering():
    assert proteomics.isoelectric_point("KKKKRRRR") > 10
    assert proteomics.isoelectric_point("DDDDEEEE") < 4


def test_pmf_identifies_protein():
    data = simulate.pmf_dataset(n_proteins=100)
    hits = proteomics.pmf_search(data.observed_mh, data.database)
    assert hits[0].protein_id == data.true_protein
    assert hits[0].score > 3 * max(hits[1].score, 1)


def test_unknown_amino_acid_raises():
    with pytest.raises(ValueError):
        proteomics.peptide_mass("PEPTIDEB")
