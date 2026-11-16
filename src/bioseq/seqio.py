"""Reading and writing sequence files, plus basic sequence utilities.

Only the standard library is used, so the whole package also runs in the
browser through Pyodide (see ``docs/index.html``).
"""

from __future__ import annotations

import gzip
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

_COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")

CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass(frozen=True)
class Record:
    """One sequence with its identifier. ``quality`` is set for FASTQ only."""

    id: str
    seq: str
    description: str = ""
    quality: str | None = None

    def __len__(self) -> int:
        return len(self.seq)


def _open(path: str | Path) -> TextIO:
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def parse_fasta(lines: Iterable[str]) -> Iterator[Record]:
    """Parse FASTA from any iterable of lines (a file handle, a list, ...)."""
    header: str | None = None
    chunks: list[str] = []
    for line in lines:
        line = line.rstrip("\r\n")
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if header is not None:
                yield _make_record(header, chunks)
            header, chunks = line[1:], []
        else:
            if header is None:
                raise ValueError("FASTA data must start with a '>' header line")
            chunks.append(line.strip())
    if header is not None:
        yield _make_record(header, chunks)


def _make_record(header: str, chunks: list[str]) -> Record:
    rid, _, desc = header.partition(" ")
    return Record(rid, "".join(chunks).upper(), desc)


def read_fasta(path: str | Path) -> list[Record]:
    with _open(path) as fh:
        return list(parse_fasta(fh))


def parse_fastq(lines: Iterable[str]) -> Iterator[Record]:
    it = iter(lines)
    for header in it:
        header = header.rstrip()
        if not header:
            continue
        if not header.startswith("@"):
            raise ValueError(f"Expected FASTQ header starting with '@', got {header[:20]!r}")
        seq = next(it).rstrip()
        plus = next(it).rstrip()
        qual = next(it).rstrip()
        if not plus.startswith("+") or len(qual) != len(seq):
            raise ValueError(f"Malformed FASTQ record {header!r}")
        rid, _, desc = header[1:].partition(" ")
        yield Record(rid, seq.upper(), desc, qual)


def read_fastq(path: str | Path) -> list[Record]:
    with _open(path) as fh:
        return list(parse_fastq(fh))


def format_fasta(records: Iterable[Record], width: int = 60) -> str:
    out: list[str] = []
    for r in records:
        out.append(f">{r.id}" + (f" {r.description}" if r.description else ""))
        out.extend(r.seq[i : i + width] for i in range(0, len(r.seq), width))
    return "\n".join(out) + "\n"


def write_fasta(records: Iterable[Record], path: str | Path, width: int = 60) -> None:
    Path(path).write_text(format_fasta(records, width))


def format_fastq(records: Iterable[Record]) -> str:
    return "".join(f"@{r.id}\n{r.seq}\n+\n{r.quality}\n" for r in records)


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def gc_content(seq: str) -> float:
    """Fraction of G+C among unambiguous bases (N is ignored)."""
    seq = seq.upper()
    acgt = sum(seq.count(b) for b in "ACGT")
    if acgt == 0:
        return 0.0
    return (seq.count("G") + seq.count("C")) / acgt


def translate(seq: str, frame: int = 0, to_stop: bool = False) -> str:
    """Translate DNA to protein with the standard genetic code (unknown codons -> X)."""
    seq = seq.upper()
    protein: list[str] = []
    for i in range(frame, len(seq) - 2, 3):
        aa = CODON_TABLE.get(seq[i : i + 3], "X")
        if aa == "*" and to_stop:
            break
        protein.append(aa)
    return "".join(protein)


def mean_quality(quality: str, offset: int = 33) -> float:
    if not quality:
        return 0.0
    return sum(ord(c) - offset for c in quality) / len(quality)
