"""Command line interface: ``bioseq <command> ...``."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

from . import __version__, assembly, bacteria, proteomics, simulate, te
from .seqio import Record, format_fasta, format_fastq, read_fasta, read_fastq


def _records(path: str) -> list[Record]:
    p = Path(path)
    if any(s in p.suffixes for s in (".fq", ".fastq")):
        return read_fastq(p)
    return read_fasta(p)


def cmd_stats(args) -> None:
    recs = _records(args.fasta)
    print(json.dumps(assembly.assembly_stats([r.seq for r in recs], args.genome_size).as_dict(), indent=2))


def cmd_assemble(args) -> None:
    reads = [r.seq for r in _records(args.reads)]
    contigs = assembly.assemble(reads, k=args.k, min_count=args.min_count, min_length=args.min_length)
    sys.stdout.write(format_fasta(contigs))
    stats = assembly.assembly_stats([c.seq for c in contigs])
    print(json.dumps(stats.as_dict()), file=sys.stderr)


def cmd_te(args) -> None:
    rows = []
    library = read_fasta(args.library) if args.library else []
    for rec in read_fasta(args.fasta):
        for kind, cands in (
            ("LTR", te.find_ltr_retrotransposons(rec.seq)),
            ("TIR", te.find_tir_transposons(rec.seq)),
        ):
            for c in cands:
                row = {"seq": rec.id, "type": kind, **asdict(c)}
                if kind == "LTR":
                    row["age_years"] = round(c.insertion_age())
                if library:
                    hit = te.classify(rec.seq[c.start : c.end], library)
                    row["family"] = hit.family if hit else None
                    row["class"] = hit.te_class if hit else None
                rows.append(row)
        mask = te.kmer_repeat_mask(rec.seq)
        print(f"# {rec.id}: {mask.repeat_fraction:.1%} in high-copy k-mer repeats", file=sys.stderr)
    _print_table(rows)


def cmd_ori(args) -> None:
    for rec in read_fasta(args.fasta):
        pred = bacteria.predict_ori(rec.seq)
        orfs = bacteria.find_orfs(rec.seq, min_protein_length=args.min_protein)
        print(json.dumps({
            "seq": rec.id,
            "length": len(rec),
            **asdict(pred),
            "dnaa_enrichment": round(pred.enrichment, 2),
            "orfs": len(orfs),
            "coding_density": round(bacteria.coding_density(orfs, len(rec)), 3),
        }, indent=2))


def cmd_protein(args) -> None:
    seq = args.sequence
    print(json.dumps({
        "length": len(seq),
        "molecular_weight": round(proteomics.molecular_weight(seq), 2),
        "isoelectric_point": proteomics.isoelectric_point(seq),
        "gravy": proteomics.gravy(seq),
        "tryptic_peptides": [
            {"sequence": p.sequence, "mh": round(p.mh, 4), "missed": p.missed_cleavages}
            for p in proteomics.digest(seq, missed_cleavages=args.missed)
        ],
    }, indent=2))


def cmd_simulate(args) -> None:
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    plant = simulate.plant_genome_with_tes(seed=args.seed)
    (out / "plant_fragment.fa").write_text(format_fasta([Record("plant_fragment", plant.sequence)]))
    (out / "te_library.fa").write_text(format_fasta(plant.library))
    _write_tsv(out / "plant_fragment_truth.tsv", [asdict(t) for t in plant.truth])
    bac = simulate.bacterial_genome(seed=args.seed)
    (out / "bacterium.fa").write_text(format_fasta([Record("bacterium", bac.sequence, "circular")]))
    reads = assembly.simulate_reads(bac.sequence[:30000], coverage=30, seed=args.seed)
    (out / "bacterium_reads.fq").write_text(format_fastq(reads))
    print(f"Wrote simulated data to {out}/", file=sys.stderr)


def _write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    cols = list(rows[0])
    lines = ["\t".join(cols)] + ["\t".join(str(r[c]) for c in cols) for r in rows]
    path.write_text("\n".join(lines) + "\n")


def _print_table(rows: list[dict]) -> None:
    if not rows:
        print("no candidates found")
        return
    cols = list(dict.fromkeys(c for r in rows for c in r))
    print("\t".join(cols))
    for r in rows:
        print("\t".join(str(r.get(c, "")) for c in cols))


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="bioseq", description=__doc__)
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(required=True)

    p = sub.add_parser("stats", help="N50/L50 and friends for a FASTA file")
    p.add_argument("fasta")
    p.add_argument("--genome-size", type=int)
    p.set_defaults(func=cmd_stats)

    p = sub.add_parser("assemble", help="de Bruijn assembly of reads (FASTQ/FASTA) to contigs on stdout")
    p.add_argument("reads")
    p.add_argument("-k", type=int, default=31)
    p.add_argument("--min-count", type=int, default=3)
    p.add_argument("--min-length", type=int, default=200)
    p.set_defaults(func=cmd_assemble)

    p = sub.add_parser("te", help="find LTR retrotransposons and TIR transposons")
    p.add_argument("fasta")
    p.add_argument("--library", help="TE consensus library (RepeatMasker-style headers family#Class/Superfamily)")
    p.set_defaults(func=cmd_te)

    p = sub.add_parser("ori", help="predict oriC/ter and ORFs of a bacterial chromosome")
    p.add_argument("fasta")
    p.add_argument("--min-protein", type=int, default=100)
    p.set_defaults(func=cmd_ori)

    p = sub.add_parser("protein", help="mass, pI, hydropathy and tryptic digest of a protein")
    p.add_argument("sequence")
    p.add_argument("--missed", type=int, default=0)
    p.set_defaults(func=cmd_protein)

    p = sub.add_parser("simulate", help="write the simulated example datasets")
    p.add_argument("--outdir", default="data")
    p.add_argument("--seed", type=int, default=7)
    p.set_defaults(func=cmd_simulate)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
