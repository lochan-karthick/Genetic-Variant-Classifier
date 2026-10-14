#!/usr/bin/env python3
"""Annotate high-quality SNVs using spliced, strand-aware CDS sequences.

Reports one representative transcript per VCF record, not exhaustive transcript effects.
"""
import argparse
import csv
import logging
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq

LOG = logging.getLogger("genvar")


def read_annotations(gff_path):
    transcripts = defaultdict(list)
    bins = defaultdict(set)
    with open(gff_path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "CDS" or cols[6] not in ("+", "-"):
                continue
            attributes = dict(part.split("=", 1) for part in cols[8].split(";") if "=" in part)
            parents = attributes.get("Parent", "").split(",")
            if not parents or not parents[0] or cols[7] not in ("0", "1", "2"):
                continue
            chrom, start, end, strand, phase = cols[0], int(cols[3]), int(cols[4]), cols[6], int(cols[7])
            for parent in parents:
                feat = (chrom, start, end, strand, phase)
                transcripts[parent].append(feat)
                for window in range((start - 1) // 10000, (end - 1) // 10000 + 1):
                    bins[(chrom, window)].add(parent)
    return transcripts, bins


def annotate_variant(chrom, pos, ref, alt, transcripts, bins, fasta, sequence_for):
    candidates = sorted(bins.get((chrom, (pos - 1) // 10000), ()))
    overlapping = False
    for parent in candidates:
        segments = transcripts[parent]
        if not any(c == chrom and a <= pos <= b for c, a, b, _, _ in segments):
            continue
        overlapping = True
        try:
            cds, locations, initial_phase, strand = sequence_for(parent)
            offset = next(
                (start_index + (pos - a if strand == "+" else b - pos) - initial_phase)
                for a, b, start_index in locations if a <= pos <= b
            )
            if offset < 0:
                raise ValueError("Variant lies in incomplete 5' codon")
            codon_start = (offset // 3) * 3
            codon = cds[codon_start:codon_start + 3]
            expected_ref = ref if strand == "+" else str(Seq(ref).complement())
            alt_on_cds = alt if strand == "+" else str(Seq(alt).complement())
            if len(codon) != 3 or any(b not in "ACGT" for b in codon):
                raise ValueError("Incomplete or ambiguous reference codon")
            if cds[offset] != expected_ref:
                raise ValueError(f"VCF REF does not match transcript: {ref}")
            alternate_codon = codon[:offset % 3] + alt_on_cds + codon[offset % 3 + 1:]
            original_aa = str(Seq(codon).translate())
            alternate_aa = str(Seq(alternate_codon).translate())
            consequence = "Synonymous" if original_aa == alternate_aa else "Non-Synonymous"
            return [consequence, parent, codon_start // 3 + 1, original_aa, alternate_aa]
        except (ValueError, KeyError, IndexError) as exc:
            LOG.warning("Could not annotate %s:%s for %s: %s", chrom, pos, parent, exc)
            continue
    if overlapping:
        return ["Unclassified", "NA", "NA", "NA", "NA"]
    return ["Non-coding", "NA", "NA", "NA", "NA"]


def main():
    parser = argparse.ArgumentParser(description="Annotate high-quality SNVs against spliced CDS models")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--fasta", required=True)
    args = parser.parse_args()

    logging.basicConfig(filename="error_log_v2.log", level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(message)s")
    transcripts, bins = read_annotations(args.gff)
    fasta = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

    @lru_cache(maxsize=2048)
    def sequence_for(parent):
        features = transcripts[parent]
        chrom, _, _, strand, initial_phase = features[0]
        if any(c != chrom or direction != strand for c, _, _, direction, _ in features):
            raise ValueError("Inconsistent CDS chromosome or strand")
        features = sorted(features, key=lambda feature: feature[1], reverse=(strand == "-"))
        initial_phase = features[0][4]
        pieces, locations = [], []
        cumulative = -initial_phase
        offset = 0
        for i, (c, start, end, direction, phase) in enumerate(features):
            if i and phase != (-cumulative) % 3:
                raise ValueError("Inconsistent internal CDS phase")
            piece = str(fasta[c].seq[start - 1:end]).upper()
            if strand == "-":
                piece = str(Seq(piece).reverse_complement())
            pieces.append(piece)
            locations.append((start, end, offset))
            offset += len(piece)
            cumulative += len(piece)
        return "".join(pieces)[initial_phase:], locations, initial_phase, strand

    counts = Counter()
    discarded_quality = 0
    with open(args.vcf, encoding="utf-8") as vcf, open("Output_Variants_v2.tsv", "w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(["Chrom", "Pos", "Ref", "Alt", "Type", "Transcript", "Protein Location", "Ref AA", "Alt AA"])
        for line in vcf:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                LOG.warning("Skipping malformed VCF record")
                continue
            chrom, pos_str, ref, alt, qual_str = fields[0], fields[1], fields[3].upper(), fields[4].upper(), fields[5]
            try:
                qual = float(qual_str)
                pos = int(pos_str)
            except ValueError:
                LOG.warning("Skipping record with missing or invalid quality/position")
                continue
            if qual <= 20:
                discarded_quality += 1
                continue
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
                LOG.warning("Skipping non-SNV record %s:%s", chrom, pos)
                continue
            try:
                genome_base = str(fasta[chrom].seq[pos - 1:pos]).upper()
            except (KeyError, IndexError) as exc:
                LOG.warning("Missing reference position %s:%s: %s", chrom, pos, exc)
                annotation = ["Unclassified", "NA", "NA", "NA", "NA"]
            else:
                if genome_base != ref:
                    LOG.warning("Reference mismatch at %s:%s, VCF=%s, FASTA=%s", chrom, pos, ref, genome_base)
                    annotation = ["Unclassified", "NA", "NA", "NA", "NA"]
                else:
                    annotation = annotate_variant(chrom, pos, ref, alt, transcripts, bins, fasta, sequence_for)
            writer.writerow([chrom, pos, ref, alt, *annotation])
            counts[annotation[0]] += 1

    LOG.info("Low quality variants excluded: %s", discarded_quality)
    LOG.info("Classifications: %s", dict(counts))
    plt.bar(list(counts), list(counts.values()))
    plt.xlabel("Variant classification")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig("Variant_Numbers_v2.png")
    plt.close()
    print("Results written to Output_Variants_v2.tsv")
    print("Classification counts:", dict(counts))


if __name__ == "__main__":
    main()
