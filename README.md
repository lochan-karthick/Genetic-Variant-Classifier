
# Genetic Variant Classifier

This repository contains a Python-based tool for classifying genetic variants as non-coding, synonymous, or non-synonymous using VCF variant data, GFF gene annotations, and a reference genome in FASTA format.

The project was developed using *Plasmodium falciparum* 3D7 genomic data.

## Overview

The classifier identifies variants overlapping coding regions and evaluates their effects on amino acid sequences. It reconstructs coding sequences using GFF annotations, accounting for exon structure, strand orientation, and reading frames.

Variants are categorised based on whether they fall outside coding regions or alter the translated amino acid sequence. The analysis uses one representative transcript per variant.

## Files

| File | Description |
|---|---|
| `GenVar_classifier.py` | Main Python script for variant classification. |
| `Output_Variants_v2.tsv` | Corrected variant classifications and predicted amino acid changes. |
| `Variant_Numbers_v2.png` | Visualisation of variant classification counts. |
| `error_log_v2.log` | Processing information and warnings. |
| `code_documenation.pdf` | Original coursework documentation. |
| `legacy/` | Original implementation and results preserved for reference. |

## Requirements

- Python 3
- Biopython
- Matplotlib

## Usage

Run the classifier using:

```bash
python3 GenVar_classifier.py --vcf variants.vcf --gff annotations.gff --fasta genome.fasta
```

Replace the example filenames with the paths to your input files.

The script generates a TSV file containing variant classifications, transcript IDs, protein positions, and reference and alternate amino acids. It also produces a classification plot and log file.

## Results

After quality filtering, 8,370 variants were classified:

| Classification | Count |
|---|---:|
| Non-coding | 2,380 |
| Synonymous | 1,992 |
| Non-synonymous | 3,998 |

## Original Implementation

The `legacy/` folder contains the original coursework implementation and outputs. The updated classifier improves coding sequence reconstruction and accounts for negative-strand genes and multi-exon transcripts.

The original VCF, GFF, and reference FASTA input files are not included in this repository.
