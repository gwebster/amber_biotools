__author__ = 'Amber Biology'

import amber_biotools as amber
from amber_biotools.human_mhc2_data import human_mhc2_alleles
from biotite.sequence.io.fasta import FastaFile
from biotite.sequence import ProteinSequence
import matplotlib.pyplot as plt
import os, math, pickle


class MHCLibrary:

    default_epitope_threshold = 0.2

    @staticmethod
    def generate_human_mhc_library():
        mhclib = MHCLibrary('Human MHCII Alleles', human_mhc2_alleles)
        return mhclib

    @staticmethod
    def profile_antibody_db_immunogenicity(fasta_file_path, save_folder):
        mhclib = MHCLibrary('Human MHCII Alleles', human_mhc2_alleles)
        fasta_file = FastaFile.read(fasta_file_path)
        antibodies = {}
        antibodies[AntibodyChain.LIGHT] = {}
        antibodies[AntibodyChain.HEAVY] = {}
        for header, string in fasta_file.items():
            if header.find('VL') > -1:
                antibodies[AntibodyChain.LIGHT][header] = ProteinSequence(string)
            elif header.find('VH') > -1:
                antibodies[AntibodyChain.HEAVY][header] = ProteinSequence(string)
            else:
                continue
        stats = {AntibodyChain.LIGHT:{'min':9.9e9, 'max':-9.9e9, 'min_ab':None, 'max_ab':None},
                 AntibodyChain.HEAVY:{'min':9.9e9, 'max':-9.9e9, 'min_ab':None, 'max_ab':None}}
        scores = {AntibodyChain.LIGHT:[], AntibodyChain.HEAVY:[]}
        for chain in antibodies:
            for ab_chain in antibodies[chain]:
                profile = mhclib.profile_immunogenicity(antibodies[chain][ab_chain])
                summary = mhclib.summarize_immunogenicity(profile)
                scores[chain].append(summary['mean_score'])
                if summary['mean_score'] < stats[chain]['min']:
                    stats[chain]['min'] = summary['mean_score']
                    stats[chain]['min_ab'] = ab_chain
                if summary['mean_score'] > stats[chain]['max']:
                    stats[chain]['max'] = summary['mean_score']
                    stats[chain]['max_ab'] = ab_chain
            stats[chain]['mean'] = sum(scores[chain]) / len(scores[chain])
        for chain in scores:
            dev2 = 0.0
            for score in scores[chain]:
                dev2 += (score - stats[chain]['mean'])**2
            stats[chain]['stdev'] = (dev2 / len(scores[chain]))**0.5
        save_file = os.path.join(save_folder, 'antibody_db_immunogenicity.dat')
        with open(save_file, 'wb') as pickle_file:
            pickle.dump(stats, pickle_file)
        return stats

    def __init__(self, name, mhc_data):
        self.data = mhc_data
        self.name = name
        allele = list(self.data.keys())[0]
        residue = list(self.data[allele].keys())[0]
        self.epitope_length = len(self.data[allele][residue])
        self.max_scores = {}
        self.compile_allele_max_scores()

    def compile_allele_max_scores(self):
        for allele in self.alleles():
            max_pos = [-1.0e8] * self.epitope_length
            for residue in self.data[allele].keys():
                for i in range(0, self.epitope_length):
                    if self.data[allele][residue][i] > max_pos[i]:
                        max_pos[i] = self.data[allele][residue][i]
            self.max_scores[allele] = sum(max_pos)

    def alleles(self):
        return sorted(self.data.keys())

    def get(self,allele,residue,position):
        return self.data[allele][residue][position]

    def score(self, epitope):
        result = []
        for allele in self.alleles():
            score = 0.0
            for i in range(0, len(epitope)):
                score += self.get(allele, epitope[i], i)
            if score < 0.0: score = 0.0
            frac = score/self.max_scores[allele]
            result.append((frac, allele))
        return sorted(result, reverse=True)

    def profile_immunogenicity(self, protein_sequence, threshold = default_epitope_threshold):
        result = {'epitope_threshold':threshold, 'scores':[], 'alleles':[], 'sequence':[]}
        sequence_length = len(protein_sequence)
        for n in range(0, sequence_length):
            result['sequence'].append(protein_sequence[n])
            if (n + 9) < sequence_length:
                peptide = str(protein_sequence[n:n + 9])
            else:
                peptide = 'G' * self.epitope_length
            epitopes = self.score(peptide)
            alleles = []
            total = 0.0
            n_alleles_above_threshold = 0
            for epitope in epitopes:
                if epitope[0] > threshold:
                    total += epitope[0]
                    n_alleles_above_threshold += 1
                    alleles.append((epitope[1], epitope[0]))
            if n_alleles_above_threshold > 0:
                score = (total / n_alleles_above_threshold) * 100.0
            else:
                score = 0.0
            result['scores'].append(score)
            result['alleles'].append(alleles)
        return result

    def summarize_immunogenicity(self, profile_results):
        alleles = []
        sum_score = 0.0
        length = len(profile_results['scores'])
        for n in range(0, length):
            for allele in profile_results['alleles'][n]:
                sum_score += allele[1]
                if not allele[0] in alleles:
                    alleles.append(allele[0])
        mean = sum_score / length
        result = {'sequence_length':length}
        result['score'] = sum_score
        result['mean_score'] = mean
        result['n_alleles'] = len(alleles)
        result['alleles'] = list(alleles)
        return result

    def generate_immunogenicity_plot(self, title, profile_results, save_folder=None, row_size=50, figwidth=10.0,
                                     aspect=0.25, dpi=100, text_y_offset=0.8, y_title_offset=1.0, manual_file_name=None):
        sequence = ['X'] + profile_results['sequence']
        x = list(range(0, len(sequence)))
        y = [0.0] + profile_results['scores']
        nalleles = [0]
        for na in range(0, len(sequence)-1):
            nalleles.append(len(profile_results['alleles'][na]))
        maxscore = max(y)
        ymax = maxscore + (maxscore * 0.05)
        nplots = math.ceil(len(x) / row_size)
        height = nplots * (figwidth * aspect)
        fig, ax = plt.subplots(nplots, 1, figsize=(figwidth, height), dpi=dpi)
        fig.suptitle(title, y=y_title_offset)
        xpos = 1
        nplot = 0
        nallele_offset = ymax + text_y_offset
        while xpos < len(x) + 1:
            xmin = xpos - 0.5
            if xpos + row_size > len(x) + 1:
                xend = -1
                xmax = xpos + row_size - 0.5
            else:
                xend = xpos + row_size
                xmax = xpos + row_size - 0.5
            ax[nplot].bar(x[xpos:xend], y[xpos:xend], color='salmon',edgecolor='white',width=1.0,align='center')
            ax[nplot].set_xlim(xmin, xmax)
            ax[nplot].set_ylim(0, ymax)
            dot_line_end = xpos + row_size - 9
            if xend == -1:
                seq_tick_end = len(sequence)
                ax[nplot].set_xticks(x[xpos + 9::10])
            else:
                seq_tick_end = xend
                ax[nplot].set_xticks(x[xpos + 9:xend:10])
            for pos in range(xpos+9, dot_line_end, 10):
                ax[nplot].axvline(x=pos, linestyle='dashed', color='lightgray')
            for pos in range(xpos, seq_tick_end):
                ax[nplot].text(pos, text_y_offset, sequence[pos], color='royalblue',
                               fontweight='bold', horizontalalignment='center')
                if nalleles[pos] > 0:
                    ax[nplot].text(pos, nallele_offset, str(nalleles[pos]), color='salmon',
                                   fontweight='normal', fontsize=8, horizontalalignment='center')
            xpos += row_size
            if xpos >= len(x):
                break
            nplot += 1
        plt.tight_layout()
        if save_folder == None:
            plt.show()
        else:
            if manual_file_name == None:
                file_stem = amber.standardize_file_stem(title)
                save_file = '{}.png'.format(file_stem)
            else:
                save_file = manual_file_name
            save_path = os.path.join(save_folder, save_file)
            plt.savefig(save_path, bbox_inches='tight')
        return

    def generate_epitope_report(self, profile):
        result = []
        for i in range(0,len(profile['scores'])):
            score = profile['scores'][i]
            if score == 0.0:
                continue
            alleles = profile['alleles'][i]
            position = i + 1
            ninemer = ''.join(profile['sequence'][i:i+9])
            nhits = len(alleles)
            txt = '\nPosition: {:5d}, Epitope = {}, Score = {:6.2f} for {:4d} allele(s)'.format(position, ninemer, score/100.0, nhits)
            txt += '\nAllele          Score    Allele          Score    Allele          Score'
            result.append(txt)
            nitem = 0
            txt = ''
            for nh in range(0, nhits):
                nitem += 1
                if nitem > 3:
                    nitem = 1
                    result.append(txt)
                    txt = ''
                txt += '{:12}  {:6.2f}    '.format(alleles[nh][0], alleles[nh][1])
            if len(txt) > 0:
                result.append(txt)
        summary = self.summarize_immunogenicity(profile)
        txt = '\nRaw Immunogenicity Score = {:.2f}, Normalized Immunogenicity Score = {:.2f}'.format(summary['score'], summary['mean_score'])
        result.append(txt)
        return '\n'.join(result)

    def __str__(self):
        txt = ["{\n"]
        for allele in sorted(self.data.keys()):
            txt.append("   '%s': {\n" % allele)
            for residue in sorted(self.data[allele].keys()):
                txt.append("      '%s': %s,\n" % (residue, self.data[allele][residue]))
            txt.append("   },\n")
        txt.append("}\n")
        return "".join(txt)


if __name__ == "__main__":
    mhclib = MHCLibrary('Human MHCII Alleles', human_mhc2_alleles)
    #print(mhclib.max_scores)
    #score = mhclib.score('VQLVESGAE')
    #for entry in score:
    #    assert entry[1] in mhclib.alleles()
    #assert len(score) == 50
    #assert score[0] == (0.40425531914893614, 'HLA-DRB1_0410')
    #assert score[1] == (0.372093023255814, 'HLA-DRB1_0806')
    #assert score[49] == (0.0, 'HLA-DRB1_0101')
    #print(score)
    from biotite.sequence.io import fasta
    from amber_biotools.antibody import Antibody, AntibodyChain
    fapath = '/Users/gordon/Google Drive/Python/amber_tools/data/inx021_FASTA.txt'
    seq_data = fasta.FastaFile.read(fapath)
    chains = []
    for header, string in seq_data.items():
        chains.append(string)
    light_chain = chains[0]
    heavy_chain = chains[1]
    ab = Antibody('My antibody', light_chain, heavy_chain)
    vl = ab.get_variable_region(AntibodyChain.LIGHT)
    vh = ab.get_variable_region(AntibodyChain.HEAVY)
    profile = mhclib.profile_immunogenicity(vh[2])
    print(profile['alleles'])
    print(profile['scores'])
    print(profile.keys())
    si = mhclib.summarize_immunogenicity(profile)
    print(si)
    sandbox = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox'
    #approved_abs = '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/approved_ab_sequences.txt'
    #abs = MHCLibrary.profile_antibody_db_immunogenicity(approved_abs, save_folder=sandbox)
    #for chain in abs:
    #    print(abs[chain])

    mhclib.generate_immunogenicity_plot('Test nalleles', profile, dpi=200, save_folder=sandbox)

    '''a_file = '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/human_albumin_FASTA.txt'
    e_file = '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/human_epo_FASTA.txt'
    a_data = fasta.FastaFile.read(a_file)
    e_data = fasta.FastaFile.read(e_file)
    for header, string in a_data.items():
        a_str = string
    for header, string in e_data.items():
        e_str = string
    albumin = ProteinSequence(a_str)
    epo = ProteinSequence(e_str)
    a_profile = mhclib.profile_immunogenicity(albumin)
    e_profile = mhclib.profile_immunogenicity(epo)
    print(mhclib.summarize_immunogenicity(a_profile))
    print(mhclib.summarize_immunogenicity(e_profile))'''
    #ab_db_datafile = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox/antibody_db_immunogenicity.dat'
    #with open(ab_db_datafile, 'rb') as pfile:
    #    ab_db_data = pickle.load(pfile)
    print()
    abl_profile = mhclib.profile_immunogenicity(ab.get_variable_region(AntibodyChain.LIGHT)[2])
    abh_profile = mhclib.profile_immunogenicity(ab.get_variable_region(AntibodyChain.HEAVY)[2])
    print(mhclib.summarize_immunogenicity(abl_profile))
    print(mhclib.summarize_immunogenicity(abh_profile))
    er = mhclib.generate_epitope_report(abh_profile)
    print(er)