#!/bin/env python3
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Align

from pathlib import Path

import pandas as pd
import sys

upstream = 'cgccatactgtgcgttataccgccagtaatgctgctcgttttgccggattatgggaaagaaataatctcataaacgaaaaattaaaaagagaagaggtttgatttaacttattgataataaagttaaaaaaacaaataaatacaagacaa'
fims = 'ttggggccaaactgtccatatcataaataagttacgtattttttctcaagcataaaaatattaaaaaacgacaaaaagcatctaactgtttgatatgtaaattatttctattgtaaattaatttcacatcacctccgctatatgtaaagctaacgtttctgtggctcgacgcatcttcctcattcttctctccaaaaaccacctcatgcaatataaacatctataaataaagataacaatagaatattaagccaacaaataaactgaaaaagtttgtccgcgatgctttcctctatgagtcaaaatggccccaa'
downstream = 'atgtttcatcttttgggggaaaactgtgcagtgttggcagtcaaactcgtttacaaaacaaagtgtacagaacgactgcccatgtcgatttagaaatagttttttgaaaggaaagcagcatgaaaattaaaactctggcaatcgttgttc'

on_seq = (upstream + fims + downstream).upper()
off_seq = (upstream + str(Seq(fims).reverse_complement()) + downstream).upper()

aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.gap_score = -3
aligner.mismatch_score = -1
aligner.match_score = 1

def check_alignment(read_seq, ref, score_thr = 130):
    '''Returns the read sequence, alignment score, whether it matches the top strand, and the coordinates in the ref'''

    alignments = aligner.align(read_seq, ref)

    if alignments.score > score_thr:

        best_alignment = alignments[0]  # Gets highest scoring alignment
        read_start = int(best_alignment.coordinates[0][0])  # Start position in read
        read_end = int(best_alignment.coordinates[0][-1])   # End position in read
        ref_start = int(best_alignment.coordinates[1][0])  # Start position in reference
        ref_end = int(best_alignment.coordinates[1][-1])   # End position in reference
        mutations = get_mutations(alignments[0])

        return [read_seq, alignments.score, True, read_start, read_end, ref_start, ref_end, mutations]

    # If there's no match forward, try reverse-complement
    rev_read_seq = str(Seq(read_seq).reverse_complement())
    alignments = aligner.align(rev_read_seq, ref)

    if alignments.score > score_thr:

        best_alignment = alignments[0]  # Gets highest scoring alignment
        read_start = int(best_alignment.coordinates[0][0])  # Start position in read
        read_end = int(best_alignment.coordinates[0][-1])   # End position in read
        ref_start = int(best_alignment.coordinates[1][0])  # Start position in reference
        ref_end = int(best_alignment.coordinates[1][-1])   # End position in reference
        mutations = get_mutations(alignments[0])

        return [rev_read_seq, alignments.score, False, read_start, read_end, ref_start, ref_end, mutations]

    return []

def get_mutations(alignment):
    '''Extract the mutations detected in the alignment.'''
    mutations = []

    ref_idx = alignment.indices[1]
    ref_seq = str(alignment[0])
    alt_seq = str(alignment[1])

    for pos, ref, alt in zip(ref_idx, ref_seq, alt_seq):
        if ref != alt:
            mutations.append(ref + str(pos) + alt)

    return ', '.join(mutations)

input_folder = Path(sys.argv[1])
output_folder = Path(sys.argv[2])

# Get all fastq files recursively
files = input_folder.rglob('*.fastq')

for file in files:

    print(f"Processing {file}")

    hits = []

    for read_id, record in enumerate(SeqIO.parse(file, "fastq")):

        read_seq = str(record.seq)
        quality = record.letter_annotations["phred_quality"]

        # Align to ON orientation
        on_hit = check_alignment(read_seq, on_seq)
        if on_hit:
            print(f'On hit for {read_id}')
            hits.append([read_id, True] + on_hit)

        # Align to OFF orientation
        off_hit = check_alignment(read_seq, off_seq)
        if off_hit:
            print(f'Off hit for {read_id}')
            hits.append([read_id, False] + off_hit)

        break

    df = pd.DataFrame(hits, columns=['read_id', 'on', 'read_seq', 'score', 'top_strand', 'read_start', 'read_end', 'ref_start', 'ref_end', 'mutations'])

    name = Path(file).stem
    df['name'] = name
    df.to_csv(output_folder / f'fims_{name}.tsv', sep = '\t', index = False)
