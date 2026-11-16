"""Five short case studies, one per field. Run: ``python examples/case_studies.py``.

All data are simulated with a known truth (see ``bioseq.simulate``), so each
study ends by checking the answer instead of just printing it.
"""

from __future__ import annotations

import time

from bioseq import assembly, bacteria, proteomics, simulate, te, variants
from bioseq.benchmark import score_intervals


def header(title: str) -> None:
    print(f"\n## {title}\n")


def study_assembly() -> None:
    header("1. Bacterial genome assembly from short reads")
    bac = simulate.bacterial_genome(length=60_000, ori=20_000, seed=3)
    reads = assembly.simulate_reads(bac.sequence, coverage=30, error_rate=0.005, seed=3, circular=True)
    print(f"{len(reads)} reads of 100 bp, 0.5% substitution errors, 30x coverage")
    spectrum = assembly.kmer_spectrum(assembly.count_kmers([r.seq for r in reads], 31))
    print(f"k-mer spectrum: {spectrum.get(1, 0)} k-mers seen once (errors), peak near "
          f"{max((m for m in spectrum if m > 5), key=spectrum.get)}x (genomic)")
    for k in (21, 31, 51):
        t = time.time()
        contigs = assembly.assemble([r.seq for r in reads], k=k, min_count=3, min_length=200)
        s = assembly.assembly_stats([c.seq for c in contigs], genome_size=len(bac.sequence))
        print(f"k={k}: {s.n_contigs} contigs, N50 {s.n50:,} bp, NG50 {s.ng50:,} bp, "
              f"total {s.total_length:,} bp ({time.time() - t:.1f}s)")


def study_plant_te() -> None:
    header("2. Plant genome: LTR retrotransposons and DNA transposons")
    rows = {"LTR": [], "TIR": []}
    for seed in range(7, 17):
        g = simulate.plant_genome_with_tes(seed=seed)
        for name, fn, prefix in (
            ("LTR", te.find_ltr_retrotransposons, "LTR"),
            ("TIR", te.find_tir_transposons, "DNA"),
        ):
            pred = [(c.start, c.end) for c in fn(g.sequence)]
            truth = [(t.start, t.end) for t in g.truth if t.te_class.startswith(prefix)]
            rows[name].append(score_intervals(pred, truth))
    print("| finder | precision | recall | exact boundaries |")
    print("| --- | --- | --- | --- |")
    for name, scores in rows.items():
        tp = sum(s.true_positives for s in scores)
        pred = sum(s.predicted for s in scores)
        truth = sum(s.truth for s in scores)
        exact = sum(s.exact for s in scores)
        print(f"| {name} | {tp / pred:.2f} | {tp / truth:.2f} | {exact}/{truth} |")

    g = simulate.plant_genome_with_tes(seed=7)
    print("\nLTR insertion ages (seed 7), molecular clock r = 1.3e-8 subst/site/year:")
    for c in te.find_ltr_retrotransposons(g.sequence):
        hit = te.classify(g.sequence[c.start : c.end], g.library)
        print(f"  {c.start:>7}-{c.end:<7} {hit.family:<11} {hit.te_class:<10} LTR identity {c.identity:.3f}"
              f"  ~{c.insertion_age() / 1e6:.2f} Myr")


def study_medical_insertion() -> None:
    header("3. Human DNA: a new Alu-like insertion in a patient's gene")
    case = simulate.human_gene_with_te_insertion()
    ins = variants.find_insertions(case.reference, case.patient)[0]
    exon = next((i + 1 for i, (s, e) in enumerate(case.exons) if ins.in_interval(s, e)), None)
    hit = te.classify(ins.sequence, case.library)
    print(f"Insertion of {ins.length} bp at reference position {ins.ref_position} "
          f"(truth {case.true_position}), TSD {ins.tsd} ({len(ins.tsd)} bp)")
    print(f"Location: {'exon ' + str(exon) if exon else 'intron'}; "
          f"frameshift: {variants.frameshift(ins.length)}")
    print(f"Best library match: {hit.family} ({hit.te_class}), identity {hit.identity:.2f}, "
          f"poly(A) tail: {te.poly_a_tail(ins.sequence)}")
    print("Interpretation: an exonic SINE insertion with a TSD and poly(A) tail is the signature of "
          "L1-mediated retrotransposition; in a real patient this would be confirmed by PCR across "
          "both junctions.")


def study_bacterial_ori() -> None:
    header("4. Bacterial DNA: origin of replication and genes")
    bac = simulate.bacterial_genome(seed=3)
    pred = bacteria.predict_ori(bac.sequence)
    print(f"Cumulative GC-skew minimum at {pred.ori:,} (true oriC {bac.ori:,}); "
          f"maximum at {pred.ter:,} (true ter {bac.ter:,})")
    print(f"DnaA boxes (<=1 mismatch) within 500 bp of the predicted origin: {pred.dnaa_boxes_near_ori} "
          f"vs {pred.dnaa_boxes_expected:.2f} expected ({pred.enrichment:.0f}x enrichment)")
    orfs = bacteria.find_orfs(bac.sequence)
    s = score_intervals([(o.start, o.end) for o in orfs], [(g[0], g[1]) for g in bac.genes], 0.8)
    print(f"ORFs >= 100 aa: {len(orfs)}; planted genes recovered: {s.true_positives}/{s.truth} "
          f"(precision {s.precision:.2f}: random sequence also contains long ORFs, which is why "
          "real gene finders model codon usage and ribosome binding sites)")


def study_proteomics() -> None:
    header("5. Proteomics: peptide mass fingerprinting")
    data = simulate.pmf_dataset()
    print(f"{len(data.observed_mh)} MALDI peaks searched against {len(data.database)} proteins, 20 ppm")
    hits = proteomics.pmf_search(data.observed_mh, data.database)
    for h in hits[:3]:
        print(f"  {h.protein_id}: {h.matched} peptides, coverage {h.coverage:.0%}, score {h.score}")
    top = data.database[hits[0].protein_id]
    print(f"Top hit is the true protein: {hits[0].protein_id == data.true_protein}. "
          f"MW {proteomics.molecular_weight(top) / 1000:.1f} kDa, pI {proteomics.isoelectric_point(top)}, "
          f"GRAVY {proteomics.gravy(top)}")


if __name__ == "__main__":
    print("# bioseq-toolkit case studies")
    study_assembly()
    study_plant_te()
    study_medical_insertion()
    study_bacterial_ori()
    study_proteomics()
