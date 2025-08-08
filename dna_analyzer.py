#!/usr/bin/env python3

import argparse
import gzip
import io
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict
from typing import Dict, Generator, Iterable, List, Optional, Tuple


# -----------------------------
# File and sequence IO utilities
# -----------------------------

def _open_maybe_gzip(path: str) -> io.TextIOBase:
    if path == "-":
        return io.TextIOWrapper(sys.stdin.buffer, encoding="utf-8")
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def detect_format(handle: Iterable[str]) -> str:
    """Detect input format from the first non-empty line.

    Returns one of: "fasta", "fastq", or "raw".
    The handle must be a seekable text stream.
    """
    pos = handle.tell()
    try:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith(">"):
                handle.seek(pos)
                return "fasta"
            if line.startswith("@"):
                # Could be FASTQ. We will assume fastq if the next 3 lines are present.
                handle.seek(pos)
                return "fastq"
            # Otherwise likely raw sequence
            handle.seek(pos)
            return "raw"
        handle.seek(pos)
        return "raw"
    except (OSError, io.UnsupportedOperation):
        # Non-seekable (e.g., stdin). Fallback to raw.
        return "raw"


def read_fasta(handle: Iterable[str]) -> Generator[Tuple[str, str], None, None]:
    header = None
    seq_chunks: List[str] = []
    for line in handle:
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield (header, "".join(seq_chunks).upper())
            header = line[1:].strip()
            seq_chunks = []
        else:
            seq_chunks.append(line.strip())
    if header is not None:
        yield (header, "".join(seq_chunks).upper())


def read_fastq(handle: Iterable[str]) -> Generator[Tuple[str, str, str], None, None]:
    while True:
        header = handle.readline()
        if not header:
            break
        if not header.startswith("@"):
            # Skip until a header line is found
            header = header.strip()
            if not header:
                continue
            raise ValueError("Invalid FASTQ: expected '@' at header line, got: " + header[:20])
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not seq or not plus or not qual:
            raise ValueError("Invalid FASTQ: truncated record")
        if not plus.startswith("+"):
            raise ValueError("Invalid FASTQ: '+' line missing")
        yield (header[1:].strip(), seq.strip().upper(), qual.strip())


def read_raw_sequences(handle: Iterable[str]) -> Generator[Tuple[str, str], None, None]:
    """Treat each non-empty line as a separate sequence. Header is a 1-based index."""
    for idx, line in enumerate(handle, start=1):
        seq = line.strip()
        if not seq:
            continue
        yield (f"seq_{idx}", seq.upper())


# -----------------------------
# Sequence analysis primitives
# -----------------------------

DNA_COMPLEMENT: Dict[str, str] = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "R": "Y",  # A or G -> T or C
    "Y": "R",  # C or T -> G or A
    "S": "S",  # G or C
    "W": "W",  # A or T
    "K": "M",  # G or T -> C or A
    "M": "K",  # A or C -> T or G
    "B": "V",  # C/G/T -> G/C/A
    "D": "H",  # A/G/T -> T/C/A
    "H": "D",  # A/C/T -> T/G/A
    "V": "B",  # A/C/G -> T/G/C
    "N": "N",
}

RNA_TRANSCRIPTION = str.maketrans({"T": "U"})

STANDARD_GENETIC_CODE: Dict[str, str] = {
    # Phenylalanine, Leucine
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    # Isoleucine, Methionine
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    # Valine
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    # Serine, Proline
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Threonine, Alanine
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyrosine, Histidine, Glutamine, Asparagine, Lysine, Aspartate, Glutamate
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    # Cysteine, Tryptophan, Arginine, Glycine, Serine, Stop
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODONS = {"ATG"}


def sanitize_dna(seq: str) -> str:
    return re.sub(r"[^ACGTRYSWKMBDHVN]", "N", seq.upper())


def reverse_complement(seq: str) -> str:
    seq = sanitize_dna(seq)
    return "".join(DNA_COMPLEMENT.get(base, "N") for base in reversed(seq))


def transcribe_to_rna(seq: str) -> str:
    return sanitize_dna(seq).translate(RNA_TRANSCRIPTION)


def translate_dna(seq: str, frame: int = 0, to_stop: bool = False) -> str:
    seq = sanitize_dna(seq)
    if frame not in (0, 1, 2):
        raise ValueError("frame must be 0, 1, or 2")
    protein_chars: List[str] = []
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        aa = STANDARD_GENETIC_CODE.get(codon, "X")
        if to_stop and aa == "*":
            break
        protein_chars.append(aa)
    return "".join(protein_chars)


def sequence_stats(seq: str, qualities: Optional[str] = None) -> Dict[str, float]:
    seq = sanitize_dna(seq)
    length = len(seq)
    counts = Counter(seq)
    g = counts.get("G", 0)
    c = counts.get("C", 0)
    a = counts.get("A", 0)
    t = counts.get("T", 0)
    n = counts.get("N", 0)
    gc = g + c
    at = a + t
    gc_pct = (gc / length * 100) if length else 0.0
    at_pct = (at / length * 100) if length else 0.0
    n_pct = (n / length * 100) if length else 0.0
    gc_skew = ((g - c) / (g + c)) if (g + c) else 0.0
    # Shannon entropy over A,C,G,T,N; ignore zeros
    probs = [counts[b] / length for b in ("A", "C", "G", "T", "N") if counts[b] > 0]
    entropy = -sum(p * math.log2(p) for p in probs) if length else 0.0

    mean_q = None
    if qualities is not None:
        # Assume Phred+33
        q_scores = [max(0, ord(ch) - 33) for ch in qualities]
        mean_q = sum(q_scores) / len(q_scores) if q_scores else None

    return {
        "length": length,
        "count_A": a,
        "count_C": c,
        "count_G": g,
        "count_T": t,
        "count_N": n,
        "gc_percent": gc_pct,
        "at_percent": at_pct,
        "n_percent": n_pct,
        "gc_skew": gc_skew,
        "shannon_entropy": entropy,
        "mean_quality": mean_q if mean_q is not None else float("nan"),
    }


def kmer_frequencies(seq: str, k: int, skip_ambiguous: bool = True) -> Dict[str, int]:
    seq = sanitize_dna(seq)
    counts: Dict[str, int] = defaultdict(int)
    if k <= 0:
        return {}
    for i in range(0, len(seq) - k + 1):
        kmer = seq[i : i + k]
        if skip_ambiguous and "N" in kmer:
            continue
        counts[kmer] += 1
    return dict(counts)


# -----------------------------
# ORF finding and motif search
# -----------------------------

def find_orfs(seq: str, min_aa_len: int = 30, both_strands: bool = True) -> List[Dict[str, object]]:
    seq = sanitize_dna(seq)

    def _find_on_strand(s: str, strand: str) -> List[Dict[str, object]]:
        orfs: List[Dict[str, object]] = []
        for frame in (0, 1, 2):
            i = frame
            while i <= len(s) - 3:
                codon = s[i : i + 3]
                if codon in START_CODONS:
                    j = i + 3
                    while j <= len(s) - 3:
                        stop_codon = s[j : j + 3]
                        if stop_codon in STOP_CODONS:
                            aa_len = (j - i) // 3
                            if aa_len >= min_aa_len:
                                dna_seq = s[i : j + 3]
                                prot = translate_dna(s[i:], frame=0, to_stop=True)
                                orfs.append(
                                    {
                                        "strand": strand,
                                        "frame": frame,
                                        "start": i,
                                        "end": j + 3,
                                        "length_nt": j + 3 - i,
                                        "length_aa": aa_len,
                                        "dna": dna_seq,
                                        "protein": prot,
                                    }
                                )
                            i = j + 3
                            break
                        j += 3
                    else:
                        # No stop codon found; move to next start position
                        i += 3
                else:
                    i += 3
        return orfs

    results = _find_on_strand(seq, "+")
    if both_strands:
        rc = reverse_complement(seq)
        # Positions for reverse complement are reported w.r.t. the RC string to avoid confusion.
        results.extend(_find_on_strand(rc, "-"))
    # Sort by decreasing length_aa
    results.sort(key=lambda r: (-r["length_aa"], r["strand"], r["frame"], r["start"]))
    return results


IUPAC_TO_REGEX = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


def iupac_to_regex(pattern: str) -> str:
    regex_chars: List[str] = []
    for ch in pattern.upper():
        regex_chars.append(IUPAC_TO_REGEX.get(ch, re.escape(ch)))
    return "".join(regex_chars)


def find_motifs(seq: str, pattern: str, use_iupac: bool = True, overlapping: bool = True) -> List[Tuple[int, int, str]]:
    seq = sanitize_dna(seq)
    regex = iupac_to_regex(pattern) if use_iupac else pattern
    flags = re.IGNORECASE
    compiled = re.compile(regex, flags)
    matches: List[Tuple[int, int, str]] = []
    if overlapping:
        # Lookahead to capture overlapping matches
        overlap_pattern = re.compile(f"(?=({compiled.pattern}))", flags)
        for m in overlap_pattern.finditer(seq):
            s = m.start(1)
            e = m.end(1)
            matches.append((s, e, seq[s:e]))
    else:
        for m in compiled.finditer(seq):
            matches.append((m.start(), m.end(), m.group(0)))
    return matches


# -----------------------------
# CLI utilities
# -----------------------------

def load_sequences(path: Optional[str], raw_sequence: Optional[str]) -> List[Tuple[str, str, Optional[str]]]:
    """Load sequences from a file or a provided raw sequence.

    Returns list of tuples: (header, sequence, qualities or None)
    """
    if raw_sequence:
        return [("input", sanitize_dna(raw_sequence), None)]
    if not path:
        raise ValueError("Either --input or --sequence must be provided")

    with _open_maybe_gzip(path) as handle:
        fmt = detect_format(handle)
        if fmt == "fasta":
            return [(h, s, None) for h, s in read_fasta(handle)]
        if fmt == "fastq":
            return list(read_fastq(handle))
        return [(h, s, None) for h, s in read_raw_sequences(handle)]


def cmd_stats(args: argparse.Namespace) -> int:
    seqs = load_sequences(args.input, args.sequence)
    results = []
    for header, seq, qual in seqs:
        res = sequence_stats(seq, qualities=qual if args.include_quality else None)
        res["id"] = header
        results.append(res)

    if args.format == "json":
        print(json.dumps(results if args.per_seq else aggregate_stats(results), indent=2))
    elif args.format == "tsv":
        rows = results if args.per_seq else [aggregate_stats(results)]
        headers = list(rows[0].keys()) if rows else []
        print("\t".join(headers))
        for row in rows:
            print("\t".join(str(row.get(h, "")) for h in headers))
    else:
        if args.per_seq:
            for row in results:
                print(format_stats_text(row))
                print("-")
        else:
            print(format_stats_text(aggregate_stats(results)))
    return 0


def aggregate_stats(per_seq: List[Dict[str, float]]) -> Dict[str, float]:
    if not per_seq:
        return {}
    out: Dict[str, float] = {}
    total_len = sum(int(x["length"]) for x in per_seq)
    out["sequences"] = len(per_seq)
    out["total_bases"] = total_len
    for key in ("count_A", "count_C", "count_G", "count_T", "count_N"):
        out[key] = sum(int(x[key]) for x in per_seq)
    out["gc_percent"] = (100.0 * (out["count_G"] + out["count_C"]) / total_len) if total_len else 0.0
    out["at_percent"] = (100.0 * (out["count_A"] + out["count_T"]) / total_len) if total_len else 0.0
    out["n_percent"] = (100.0 * out["count_N"] / total_len) if total_len else 0.0
    # Average of per-seq skew and entropy weighted by length
    def weighted_avg(key: str) -> float:
        num = sum(float(x[key]) * int(x["length"]) for x in per_seq)
        den = total_len if total_len else 1
        return num / den

    out["gc_skew"] = weighted_avg("gc_skew")
    out["shannon_entropy"] = weighted_avg("shannon_entropy")
    # Mean quality averaged over sequences (if present)
    if not math.isnan(per_seq[0].get("mean_quality", float("nan"))):
        out["mean_quality"] = sum(float(x.get("mean_quality", 0.0)) for x in per_seq) / len(per_seq)
    return out


def format_stats_text(stats: Dict[str, float]) -> str:
    lines = []
    order = [
        "id", "length", "count_A", "count_C", "count_G", "count_T", "count_N",
        "gc_percent", "at_percent", "n_percent", "gc_skew", "shannon_entropy", "mean_quality",
    ]
    for key in order:
        if key in stats:
            lines.append(f"{key}: {stats[key]}")
    for key in stats:
        if key not in order:
            lines.append(f"{key}: {stats[key]}")
    return "\n".join(lines)


def cmd_kmer(args: argparse.Namespace) -> int:
    seqs = load_sequences(args.input, args.sequence)
    combined = Counter()
    for _, seq, _ in seqs:
        combined.update(kmer_frequencies(seq, args.k, skip_ambiguous=not args.include_ambiguous))
    if args.top:
        items = combined.most_common(args.top)
    else:
        items = sorted(combined.items(), key=lambda kv: (-kv[1], kv[0]))

    if args.format == "json":
        print(json.dumps({k: v for k, v in items}, indent=2))
    else:
        print("kmer\tcount")
        for k, v in items:
            print(f"{k}\t{v}")
    return 0


def cmd_revcomp(args: argparse.Namespace) -> int:
    seqs = load_sequences(args.input, args.sequence)
    for header, seq, _ in seqs:
        rc = reverse_complement(seq)
        if args.fasta:
            print(f">{header}|revcomp")
            print(rc)
        else:
            print(rc)
    return 0


def cmd_translate(args: argparse.Namespace) -> int:
    seqs = load_sequences(args.input, args.sequence)
    frames: List[int]
    if args.frame == "all":
        frames = [0, 1, 2]
    else:
        frames = [int(args.frame) - 1]

    for header, seq, _ in seqs:
        targets = [("+", seq)]
        if args.both_strands:
            targets.append(("-", reverse_complement(seq)))
        for strand, s in targets:
            for fr in frames:
                prot = translate_dna(s, frame=fr, to_stop=args.to_stop)
                prefix = f">{header}|strand={strand}|frame={fr+1}"
                if args.fasta:
                    print(prefix)
                    print(prot)
                else:
                    print(prefix, prot, sep="\t")
    return 0


def cmd_orfs(args: argparse.Namespace) -> int:
    seqs = load_sequences(args.input, args.sequence)
    for header, seq, _ in seqs:
        orfs = find_orfs(seq, min_aa_len=args.min_aa, both_strands=args.both_strands)
        if args.format == "json":
            payload = {"id": header, "orfs": orfs}
            print(json.dumps(payload, indent=2))
        else:
            print(f"# {header}")
            print("strand\tframe\tstart\tend\tlen_nt\tlen_aa\tprotein")
            for o in orfs:
                print(
                    f"{o['strand']}\t{o['frame']+1}\t{o['start']}\t{o['end']}\t{o['length_nt']}\t{o['length_aa']}\t{o['protein']}"
                )
    return 0


def cmd_motifs(args: argparse.Namespace) -> int:
    seqs = load_sequences(args.input, args.sequence)
    for header, seq, _ in seqs:
        matches = find_motifs(seq, pattern=args.pattern, use_iupac=not args.regex, overlapping=not args.non_overlapping)
        if args.format == "json":
            payload = {"id": header, "matches": [
                {"start": s, "end": e, "match": m} for (s, e, m) in matches
            ]}
            print(json.dumps(payload, indent=2))
        else:
            print(f"# {header}")
            print("start\tend\tmatch")
            for s, e, m in matches:
                print(f"{s}\t{e}\t{m}")
    return 0


# -----------------------------
# Argument parser setup
# -----------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="DNA Sequence Analyzer: stats, k-mers, reverse complement, translation, ORFs, motif search",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    sub = p.add_subparsers(dest="command", required=True)

    # stats
    sp = sub.add_parser("stats", help="Compute sequence statistics")
    sp.add_argument("--input", "-i", type=str, help="FASTA/FASTQ/raw file path or '-' for stdin")
    sp.add_argument("--sequence", "-s", type=str, help="Raw sequence string")
    sp.add_argument("--per-seq", action="store_true", help="Report per-sequence stats instead of aggregate")
    sp.add_argument("--include-quality", action="store_true", help="Compute mean quality if FASTQ")
    sp.add_argument("--format", choices=["text", "json", "tsv"], default="text")
    sp.set_defaults(func=cmd_stats)

    # kmer
    sp = sub.add_parser("kmer", help="Count k-mer frequencies")
    sp.add_argument("--input", "-i", type=str, help="FASTA/FASTQ/raw file path or '-' for stdin")
    sp.add_argument("--sequence", "-s", type=str, help="Raw sequence string")
    sp.add_argument("--k", type=int, default=3, help="k-mer size")
    sp.add_argument("--include-ambiguous", action="store_true", help="Include kmers containing N")
    sp.add_argument("--top", type=int, help="Show only the top N k-mers")
    sp.add_argument("--format", choices=["tsv", "json"], default="tsv")
    sp.set_defaults(func=cmd_kmer)

    # revcomp
    sp = sub.add_parser("revcomp", help="Reverse-complement sequences")
    sp.add_argument("--input", "-i", type=str, help="FASTA/FASTQ/raw file path or '-' for stdin")
    sp.add_argument("--sequence", "-s", type=str, help="Raw sequence string")
    sp.add_argument("--fasta", action="store_true", help="FASTA-style output")
    sp.set_defaults(func=cmd_revcomp)

    # translate
    sp = sub.add_parser("translate", help="Translate DNA to protein")
    sp.add_argument("--input", "-i", type=str, help="FASTA/FASTQ/raw file path or '-' for stdin")
    sp.add_argument("--sequence", "-s", type=str, help="Raw sequence string")
    sp.add_argument("--frame", choices=["1", "2", "3", "all"], default="1", help="Frame to translate")
    sp.add_argument("--both-strands", action="store_true", help="Translate both strands")
    sp.add_argument("--to-stop", action="store_true", help="Stop translation at first stop codon")
    sp.add_argument("--fasta", action="store_true", help="FASTA-style output")
    sp.set_defaults(func=cmd_translate)

    # orfs
    sp = sub.add_parser("orfs", help="Find open reading frames (ATG..stop)")
    sp.add_argument("--input", "-i", type=str, help="FASTA/FASTQ/raw file path or '-' for stdin")
    sp.add_argument("--sequence", "-s", type=str, help="Raw sequence string")
    sp.add_argument("--min-aa", type=int, default=30, help="Minimum ORF length in amino acids")
    sp.add_argument("--both-strands", action="store_true", help="Search both strands")
    sp.add_argument("--format", choices=["text", "json"], default="text")
    sp.set_defaults(func=cmd_orfs)

    # motifs
    sp = sub.add_parser("motifs", help="Find motifs (IUPAC or regex)")
    sp.add_argument("--input", "-i", type=str, help="FASTA/FASTQ/raw file path or '-' for stdin")
    sp.add_argument("--sequence", "-s", type=str, help="Raw sequence string")
    sp.add_argument("--pattern", "-p", type=str, required=True, help="Motif pattern (IUPAC by default)")
    sp.add_argument("--regex", action="store_true", help="Interpret pattern as raw regex instead of IUPAC")
    sp.add_argument("--non-overlapping", action="store_true", help="Do not report overlapping matches")
    sp.add_argument("--format", choices=["text", "json"], default="text")
    sp.set_defaults(func=cmd_motifs)

    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except BrokenPipeError:
        # Allow piping to head, etc.
        return 0
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())