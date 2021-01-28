__author__ = 'Amber Biology'

import pytest, os
from amber_biotools.sequence_alignments import SequenceAlignment
from biotite.sequence import ProteinSequence, NucleotideSequence, GeneralSequence, Alphabet
from biotite.sequence.io import fasta


sequences = {
    'adalimumab':'/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/adalimumab_FASTA.txt',
    'ab_full':'/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/antibody_full_FASTA.txt',
    'emd72000':'/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/emd72000_FASTA.txt',
    'inx021':'/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/inx021_FASTA.txt',
    'rituximab':'/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/rituximab_FASTA.txt'}

chains = {}
select_heavy = 0
select_light = 1

for sequence in sequences:
    chains[sequence] = []
    seq_data = fasta.FastaFile.read(sequences[sequence])
    for header, string in seq_data.items():
        chains[sequence].append(ProteinSequence(string))

sequence_stats = {'adalimumab':[4, 3498], 'ab_full':[6, 7126], 'emd72000':[21, 14573],
    'inx021':[4, 537], 'rituximab':[9, 14323]}


def test_banner():
    print('\nStarting Sequence Alignment Tests ...')

def test_sequence_alignment_init():
    print('Testing sequence alignment init ...')
    seq1 = chains['adalimumab'][select_heavy]
    seq2 = chains['rituximab'][select_heavy]
    limit = None
    sa = SequenceAlignment(seq1, seq2,  residue_limit=limit, local=False)
    assert sa.get_alignment_count() == 1
    assert sa.get_score(0) == 873
    gs = sa.get_gapped_sequences(0)
    gs0 = 'DIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQRYNRAPYTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGE-'
    gs1 = 'QIVLSQSPAILSASPGEKVTMTCRASSSV-SYIHWFQQKPGSSPKPWIYATSNLASGVPVRFSGSGSGTSYSLTISRVEAEDAATYYCQQWTSNPPTFGGGTKLEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC'
    assert gs[0] == gs0
    assert gs[1] == gs1
    trace = sa.get_trace_map(0)
    assert trace[0][0] == 0
    assert trace[1][0] == 0
    assert trace[0][-1] == 212
    assert trace[1][-1] == 211