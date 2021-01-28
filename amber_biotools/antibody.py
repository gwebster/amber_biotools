__author__ = 'Amber Biology'

import amber_biotools as amber
from amber_biotools import AntibodyChain
from amber_biotools.antibody_germ_line import AntibodyGermLine
from biotite.sequence import ProteinSequence
from biotite.sequence.io import fasta
import biotite.sequence.align as bioalign
from amber_biotools.sequence_alignments import MultipleSequenceAlignment, BlastAlignment
from amber_biotools.immunogenicity import MHCLibrary
import re, time

default_subsitution_matrix = bioalign.SubstitutionMatrix.std_protein_matrix()
human_antibody_germ_line = AntibodyGermLine()


class Antibody:

    modifications = {'n-linked-glycosylation': 'N[ACDEFGHIKLMNQRSTVWY][ST]',
                          'o-linked-glycosylation': '(P[S/T])|([S/T]..P)',
                          'deamidation': '[NQ][AGST]'}

    human_mhc_lib = MHCLibrary.generate_human_mhc_library()

    def __init__(self, name, light_chain_sequence, heavy_chain_sequence):
        self.name = name
        self.light_sequence = light_chain_sequence
        self.heavy_sequence = heavy_chain_sequence
        self.sequences = {AntibodyChain.LIGHT:self.light_sequence, AntibodyChain.HEAVY:self.heavy_sequence}
        self.structure = None
        self.cdr = {AntibodyChain.LIGHT:{}, AntibodyChain.HEAVY:{}}
        self.cdr_scan()

    def cdr_scan(self):
        self.scan_cdr_light()
        self.scan_cdr_heavy()

    def scan_cdr_light(self):
        cys_pos = 22
        slack = 5
        cdr_scan = {1:(), 2:(), 3:()}
        sequence = str(self.light_sequence)
        cdr1 = re.compile('C.{8,20}W[YLF][QL]')
        match = cdr1.finditer(sequence, cys_pos - slack)
        cdr1_start = -1
        for m in match:
            cdr1_start = m.start() + 1
            cdr1_end = m.end() - 3
            cdr1 = sequence[cdr1_start:cdr1_end]
            if cys_pos - slack <= cdr1_start <= cys_pos + slack:
                break
        # found cdr1?
        if cdr1_start > -1:
            cdr_scan[1] = (cdr1_start + 1, cdr1_end, cdr1, (cdr1_start, cdr1_end))
        else:
            return cdr_scan
        cdr2_start = cdr1_end + 15
        cdr2_end = cdr2_start + 7
        cdr2 = sequence[cdr2_start:cdr2_end]
        # found cdr2
        cdr_scan[2] = (cdr2_start + 1, cdr2_end, cdr2, (cdr2_start, cdr2_end))
        cdr3_start = cdr2_end + 32
        cdr3_end = -1
        cdr3_tail = re.compile('FG.G')
        match = cdr3_tail.finditer(sequence, cdr3_start)
        for m in match:
            cdr3_end = m.start()
        # found cdr3?
        if cdr3_end > -1:
            cdr3 = sequence[cdr3_start:cdr3_end]
            cdr_scan[3] = (cdr3_start + 1, cdr3_end, cdr3, (cdr3_start, cdr3_end))
        self.cdr[AntibodyChain.LIGHT] = cdr_scan
        return

    def scan_cdr_heavy(self):
        cys_pos = 22
        slack = 5
        cdr_scan = {1:(), 2:(), 3:()}
        sequence = str(self.heavy_sequence)
        cdr1 = re.compile('C.{10,20}W[VIA]')
        match = cdr1.finditer(sequence, cys_pos - slack)
        cdr1_start = -1
        for m in match:
            cdr1_start = m.start() + 4
            cdr1_end = m.end() - 2
            cdr1 = sequence[cdr1_start:cdr1_end]
            if cys_pos - slack <= cdr1_start <= cys_pos + slack:
                break
        # found cdr1?
        if cdr1_start > -1:
            cdr_scan[1] = (cdr1_start+1, cdr1_end, cdr1, (cdr1_start, cdr1_end))
        else:
            return cdr_scan
        cdr2_start = cdr1_end + 15
        cdr2_end = -1
        cdr2_tail = re.compile('[KR][LIVFTA][TSIA]')
        match = cdr2_tail.finditer(sequence, cdr2_start)
        for m in match:
            cdr2_end = m.start()
            break
        # found cdr2?
        cdr2 = sequence[cdr2_start:cdr2_end]
        if cdr2_end > -1:
            cdr_scan[2] = (cdr2_start + 1, cdr2_end, cdr2, (cdr2_start, cdr2_end))
        else:
            return cdr_scan
        cdr3_pattern = re.compile('C.{5,30}WG.G')
        cdr3_start = -1
        match = cdr3_pattern.finditer(sequence, cdr2_end)
        for m in match:
            cdr3_start = m.start() + 3
            cdr3_end = m.end() - 4
        # found cdr3?
        if cdr3_end > -1:
            cdr3 = sequence[cdr3_start:cdr3_end]
            cdr_scan[3] = (cdr3_start+1, cdr3_end, cdr3, (cdr3_start, cdr3_end))
        self.cdr[AntibodyChain.HEAVY] = cdr_scan
        return

    def get_cdr(self, chain, cdr_number):
        return self.cdr[chain][cdr_number]

    def get_variable_region(self, chain):
        if chain == AntibodyChain.LIGHT:
            after_cdr3 = 10
        elif chain == AntibodyChain.HEAVY:
            after_cdr3 = 11
        else:
            return ''
        variable_region_end = self.get_cdr(chain, 3)[1] + after_cdr3
        variable_region = ProteinSequence(self.sequences[chain][:variable_region_end])
        return (0, variable_region_end, variable_region)

    def load_structure(self, structure, light_chain_id, heavy_chain_id):
        if not isinstance(light_chain_id, tuple):
            light_chain_id = (light_chain_id, 1)
        if not isinstance(heavy_chain_id, tuple):
            heavy_chain_id = (heavy_chain_id, 1)
        self.structure = structure
        self.structure_chain_id = {}
        self.structure_chain_id[AntibodyChain.LIGHT] = light_chain_id
        self.structure_chain_id[AntibodyChain.HEAVY] = heavy_chain_id
        return

    def has_structure(self):
        return not self.structure == None

    def map_structure_to_sequence(self):
        light_chain_id = self.structure_chain_id[AntibodyChain.LIGHT]
        heavy_chain_id = self.structure_chain_id[AntibodyChain.HEAVY]
        alignments = {}
        alignments[AntibodyChain.LIGHT] = self.structure.align_with_sequence(self.light_sequence, default_subsitution_matrix, light_chain_id)
        alignments[AntibodyChain.HEAVY] = self.structure.align_with_sequence(self.heavy_sequence, default_subsitution_matrix, heavy_chain_id)
        result = {}
        nblank = -1
        for chain in [AntibodyChain.LIGHT, AntibodyChain.HEAVY]:
            result[chain] = {}
            chain_id = self.structure_chain_id[chain]
            struc_chain = self.structure.chains[chain_id]['sequence']
            residue_ids = list(struc_chain.keys())
            alignment = alignments[chain].alignments[0]
            for pair_map in alignment.trace:
                seq_pos = pair_map[1]
                struc_pos = pair_map[0]
                if not struc_pos == -1:
                    residue_id = residue_ids[struc_pos]
                    residue_name = struc_chain[residue_id][0]
                else:
                    residue_id = -1
                    residue_name = ' '
                if seq_pos == -1:
                    seq_pos = nblank
                    nblank -= 1
                result[chain][seq_pos] = (struc_pos, residue_id, residue_name)
        return result

    def get_sequence_string(self, chain):
        return str(self.sequences[chain])

    def scan_modification_sites(self, chain, modifications=modifications):
        chain_string = self.get_sequence_string(chain)
        result = {}
        for mod_type in modifications:
            result[mod_type] = []
        for mod_type in modifications:
            matches = re.finditer(modifications[mod_type], chain_string)
            for match in matches:
                result[mod_type].append((match.start(), match.group()))
        return result

    def compute_solvent_accessibility(self):
        if self.structure == None:
            return
        chain_ids = list(self.structure_chain_id.values())
        result = self.structure.map_solvent_accessibility(chain_ids)
        return result

    def plot_solvent_accessibility(self, sa_results, save_folder=None, figwidth = 10.0, aspect = 0.2, y_title_offset=1.0,
                                   dpi=100, row_size=25, text_y_offset=10.0, title_dict=None, manual_file_name_dict=None):
        if self.structure == None:
            return
        self.structure.plot_solvent_accessibility(sa_results, save_folder=save_folder, figwidth=figwidth, aspect=aspect,
                    y_title_offset=y_title_offset, dpi=dpi, row_size=row_size, text_y_offset=text_y_offset,
                    title_dict=title_dict, manual_file_name_dict=manual_file_name_dict)
        return

    def get_homologous_human_frameworks(self, chain, nresults=10):
        scan_sequence = self.get_variable_region(chain)[2]
        scan = human_antibody_germ_line.scan_human_frameworks(scan_sequence, chain)
        return scan[:nresults]

    def align_with_human_frameworks(self, chain, nframeworks=10):
        scans = self.get_homologous_human_frameworks(chain, nframeworks)
        if chain == AntibodyChain.LIGHT:
            this_label = '{} VL'.format(self.name)
        elif chain == AntibodyChain.HEAVY:
            this_label = '{} VH'.format(self.name)
        else:
            return
        this_sequence = self.get_variable_region(chain)[2]
        labels = [this_label]
        sequences = [this_sequence]
        for scan in scans:
            framework = scan[2]
            labels.append(framework)
            sequences.append(human_antibody_germ_line.get_sequence(chain, framework))
        alignment = MultipleSequenceAlignment(labels, sequences)
        return alignment

    def generate_human_alignment_plot(self, chain, nframeworks=10, figsize=(8.0, 11.0), dpi=100, save_path=None,
                                      show_numbers=True, show_line_position=False,type_based=False):
        alignment = self.align_with_human_frameworks(chain, nframeworks=nframeworks)
        alignment.plot_alignment(figsize=figsize, dpi=dpi, save_path=save_path, show_numbers=show_numbers,
                                 show_line_position=show_line_position,type_based=type_based)
        return

    def profile_immunogenicity(self):
        result = {}
        light_sequence = self.get_variable_region(AntibodyChain.LIGHT)[2]
        heavy_sequence = self.get_variable_region(AntibodyChain.HEAVY)[2]
        result[AntibodyChain.LIGHT] = Antibody.human_mhc_lib.profile_immunogenicity(light_sequence)
        result[AntibodyChain.HEAVY] = Antibody.human_mhc_lib.profile_immunogenicity(heavy_sequence)
        return result

    def generate_immunogenicity_plots(self, file_name_dict, save_folder=None, row_size=50, figwidth=10.0,
                                     aspect=0.25, dpi=100, text_y_offset=0.8):
        subtitle = {AntibodyChain.LIGHT:'(light chain)', AntibodyChain.HEAVY:'(heavy chain)'}
        profile = self.profile_immunogenicity()
        for chain in profile:
            title = 'Epitopes: {} {}'.format(self.name, subtitle[chain])
            self.human_mhc_lib.generate_immunogenicity_plot(title, profile[chain], save_folder=save_folder, row_size=row_size,
                         figwidth=figwidth, aspect=aspect, dpi=dpi, text_y_offset=text_y_offset, manual_file_name=file_name_dict[chain])
        return

    def generate_epitope_reports(self):
        result = {}
        profile = self.profile_immunogenicity()
        for chain in profile:
            result[chain] = Antibody.human_mhc_lib.generate_epitope_report(profile[chain])
        return result

    def generate_blast_alignments(self, max_results=100, save_file_stem=None, pause=60):
        suffix = {AntibodyChain.LIGHT: '_vl.dat', AntibodyChain.HEAVY: '_vh.dat'}
        blast_alignments = {}
        for chain in self.sequences:
            sequence = self.get_variable_region(chain)[2]
            if save_file_stem == None:
                save_file = None
            else:
                save_file = save_file_stem + suffix[chain]
            blast_alignments[chain] = BlastAlignment.blast_sequence(sequence, max_results=max_results, save_file=save_file)
            # pause to avoid violating BLAST server rules for too frequent requests
            time.sleep(pause)
        return blast_alignments

    def generate_canonical_sequence_plots(self, blast_alignments_dict, file_name_dict, save_folder=None, row_size=50, figwidth=10.0,
                                     aspect=0.25, dpi=100, text_y_offset=2.5, y_title_offset=1.0):
        chains = {AntibodyChain.LIGHT:'(light chain)', AntibodyChain.HEAVY:'(heavy chain)'}
        for chain in self.sequences:
            blast = blast_alignments_dict[chain]
            title = 'Canonical Sequence Analysis: {} {}'.format(self.name, chains[chain])
            file_name = file_name_dict[chain]
            blast.generate_canonical_sequence_plot(title, save_folder=save_folder, row_size=row_size, figwidth=figwidth,
                aspect=aspect, dpi=dpi, text_y_offset=text_y_offset, y_title_offset=y_title_offset, manual_file_name=file_name)
        return

    def generate_organism_pie_charts(self, blast_alignments_dict, file_name_dict, explode_organism='Homo sapiens', save_folder=None,
                                     figsize=(10.0, 10.0), dpi=100, nlegend=5, explode_level=0.1):
        chains = {AntibodyChain.LIGHT:'(light chain)', AntibodyChain.HEAVY:'(heavy chain)'}
        for chain in self.sequences:
            blast = blast_alignments_dict[chain]
            title = 'Organism Frequency in Alignments: {} {}'.format(self.name, chains[chain])
            file_name = file_name_dict[chain]
            blast.generate_organism_pie_chart(title, explode_organism=explode_organism, save_folder=save_folder, figsize=figsize,
                                dpi=dpi, nlegend=nlegend, explode_level=explode_level, manual_file_name=file_name)
        return

    def format_sequence(self, chain, variable_region=True, highlight='cdr', add_numbers=True):
        if variable_region:
            sequence = self.get_variable_region(chain)[2]
        else:
            sequence = self.sequences[chain]
        if highlight == 'cdr':
            highlight = []
            for cdr in self.cdr[chain]:
                start = self.cdr[chain][cdr][0]
                end = self.cdr[chain][cdr][1] + 1
                for n in range(start, end):
                    highlight.append(n)
        return amber.format_sequence(sequence, highlight=highlight, add_numbers=add_numbers)

    def scan_cdr_residues(self, chain, residues=['W', 'M']):
        result = {1:{}, 2:{}, 3:{}}
        cdr = self.cdr[chain]
        sequence = self.sequences[chain]
        for ncdr in cdr:
            result[ncdr]['cdr'] = '{:d}-{}-{:d}'.format(cdr[ncdr][0], cdr[ncdr][2], cdr[ncdr][1])
            start = cdr[ncdr][3][0]
            end = cdr[ncdr][3][1]
            for residue in residues:
                result[ncdr][residue] = []
                for n in range(start, end):
                    if sequence[n] == residue:
                        result[ncdr][residue].append(n)
        return result





if __name__ == '__main__':
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
    ab = Antibody(use, chains[use][int(AntibodyChain.LIGHT)], chains[use][int(AntibodyChain.HEAVY)])
    print(ab.cdr)
    print(ab.get_variable_region(AntibodyChain.LIGHT)[2])
    print(ab.get_variable_region(AntibodyChain.HEAVY)[2])
    #s = Structure.load_from_pdb('adalimumab', '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/4nyl.pdb')
    #ab.load_structure(s, ('L',1), ('H',1))
    print(ab.scan_modification_sites(AntibodyChain.LIGHT))
    print(ab.scan_modification_sites(AntibodyChain.HEAVY))
    print(ab.cdr[AntibodyChain.LIGHT])
    print(ab.scan_cdr_residues(AntibodyChain.LIGHT))
    print(ab.cdr[AntibodyChain.HEAVY])
    print(ab.scan_cdr_residues(AntibodyChain.HEAVY))
    #abmap = ab.map_structure_to_sequence()
    #for chain in abmap:
    #    for pos in abmap[chain]:
    #        print(ab.sequences[chain][pos], end='')
    #    print()
    #    for pos in abmap[chain]:
    #        print(abmap[chain][pos][2], end='')
    #    print()
    #scans = ab.get_homologous_human_frameworks(AntibodyChain.HEAVY)
    #for scan in scans:
    #    print(scan)
    #align = ab.align_with_human_frameworks(AntibodyChain.HEAVY)
    #print(align)
    #plot_path = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox/ab_full_heavy_align.png'
    #ab.generate_human_alignment_plot(AntibodyChain.HEAVY, save_path=plot_path, dpi=200)
    #profile = ab.profile_immunogenicity()
    #print(len(profile[AntibodyChain.HEAVY]['scores']))
    #print(len(profile[AntibodyChain.HEAVY]['sequence']))
    sandbox = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox/'
    #ab.generate_immunogenicity_plots(dpi=200, save_folder=sandbox)
    #align = ab.generate_blast_alignments(max_results=500, save_file_stem=sandbox + 'blast_adalimumab')
    #ab.generate_canonical_sequence_plots(align, save_folder=sandbox)
    #print()
    #print(ab.format_sequence(AntibodyChain.HEAVY))
    #print()
    #print(ab.format_sequence(AntibodyChain.HEAVY, highlight=[49, 50, 51, 52, 53, 54, 55]))
    #ab.generate_human_alignment_plot(AntibodyChain.LIGHT, save_path=sandbox+'ab_human_align_light.png')
    epitopes = ab.generate_epitope_reports()
    for chain in epitopes:
        print(epitopes[chain])
        print('\n----------------------------\n')