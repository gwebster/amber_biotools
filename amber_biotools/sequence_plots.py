__author__ = 'Amber Biology'

from biotite.sequence import ProteinSequence, NucleotideSequence, GeneralSequence
import matplotlib.pyplot as plt
import numpy as np
import math

class SequencePlotter:

    def __init__(self, sequence):
        self.sequence = sequence
        self.figwidth = 10.0
        self.figheight = 5.0
        self.padfactor = 0.2
        self.figdpi = 200
        self.fontsize = 12
        self.rowsize = 50
        self.fontcolor = 'black'
        self.setup_plot()

    def setup_plot(self):
        self.nrows = math.ceil(len(self.sequence) / self.rowsize)
        self.figure = plt.figure(figsize=(self.figwidth, self.figheight), dpi=self.figdpi)
        self.xpadding = self.rowsize * self.padfactor
        self.ypadding = self.nrows * self.padfactor
        self.ax = plt.axis([0, self.rowsize+self.xpadding, 0, self.nrows+self.ypadding])
        self.xmin = self.ax[0]
        self.xmax = self.ax[1]
        self.ymin = self.ax[2]
        self.ymax = self.ax[3]
        plt.axis('off')

    def show_plot(self):
        plt.show()

    def plot_residue(self, nresidue):
        pos = self.get_position(nresidue)
        plt.text(pos[0], pos[1], self.sequence[nresidue-1], color=self.fontcolor, fontsize=self.fontsize, ha='center')

    def get_position(self, nresidue):
        x = self.xmin + self.xpadding/2.0 + (nresidue % self.rowsize)
        row_offset = int(math.ceil(nresidue / self.rowsize)) - 1
        y = self.ymax - self.ypadding/2.0 - row_offset
        return (x, y)


if __name__ == '__main__':
    from biotite.sequence.io import fasta
    sequences = {
        'adalimumab': '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/adalimumab_FASTA.txt',
        'ab_full': '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/antibody_full_FASTA.txt',
        'emd72000': '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/emd72000_FASTA.txt',
        'inx021': '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/inx021_FASTA.txt',
        'rituximab': '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/rituximab_FASTA.txt'}

    chains = {}

    for sequence in sequences:
        chains[sequence] = []
        seq_data = fasta.FastaFile.read(sequences[sequence])
        for header, string in seq_data.items():
            chains[sequence].append(ProteinSequence(string))

    use = 'adalimumab'
    select_chain = 0
    sequence = chains[use][select_chain]
    sp = SequencePlotter(sequence)
    for n in range(1, 106):
        sp.plot_residue(n)
    sp.show_plot()