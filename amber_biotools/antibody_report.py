__author__ = 'Amber Biology'

import amber_biotools as amber
import amber_biotools.antibody as amber_antibody
import amber_biotools.structure as amber_structure
import amber_biotools.sequence_alignments as amber_align
import os, pickle

introduction_text = ['This report contains a {}-based analysis of the antibody "{}".',
                     'The supplied light and heavy chain sequence variable regions (shown below)',
                     'have approximate molecular weights of {:.1f} and {:.1f} respectively.',
                     'Most of the analyses featured in this report are computed using the',
                     'data from the antibody variable regions. If a 3D structure of the',
                     'antibody is also available, the solvent-accessible surface analysis',
                     'is performed for the entire structure.']

class AntibodyReport:

    def __init__(self, name, light_sequence, heavy_sequence, folder_path):
        self.name = name
        self.antibody = amber_antibody.Antibody(self.name, light_sequence, heavy_sequence)
        report_folder_stem = amber.standardize_file_stem('{}_antibody_report'.format(self.name))
        self.report_folder = os.path.join(folder_path, report_folder_stem)
        if not os.path.exists(self.report_folder):
            os.mkdir(self.report_folder)

    def load_structure(self, pdb_file_path, light_chain_id, heavy_chain_id):
        structure = amber_structure.Structure.load_from_pdb(self.name, pdb_file_path)
        self.antibody.load_structure(structure, light_chain_id, heavy_chain_id)

    def generate_introduction(self):
        analysis = 'sequence'
        if self.antibody.has_structure():
            analysis = 'structure and sequence'
        light_mwt = self.antibody.get_variable_region(amber_antibody.AntibodyChain.LIGHT)[2].get_molecular_weight()
        heavy_mwt = self.antibody.get_variable_region(amber_antibody.AntibodyChain.HEAVY)[2].get_molecular_weight()
        intro_string = ' '.join(introduction_text)
        introduction = intro_string.format(analysis, self.name, light_mwt, heavy_mwt)
        intro_file_path = os.path.join(self.report_folder, '01_introduction.txt')
        with open(intro_file_path, 'w') as intro_file:
            intro_file.write(introduction)
        return

    def generate_formatted_sequences(self):
        formatted_sequence = self.antibody.format_sequence(amber_antibody.AntibodyChain.LIGHT)
        sequence_file_path = os.path.join(self.report_folder, '02_light_chain_sequence.txt')
        with open(sequence_file_path, 'w') as sequence_file:
            sequence_file.write(formatted_sequence)
        formatted_sequence = self.antibody.format_sequence(amber_antibody.AntibodyChain.HEAVY)
        sequence_file_path = os.path.join(self.report_folder, '03_heavy_chain_sequence.txt')
        with open(sequence_file_path, 'w') as sequence_file:
            sequence_file.write(formatted_sequence)
        return

    def generate_blast_alignments(self):
        save_file_stem = os.path.join(self.report_folder, 'blast')
        blast = self.antibody.generate_blast_alignments(max_results=500, save_file_stem=save_file_stem)
        return blast

    def reload_blast_alignments(self):
        blast_file_names = {amber_antibody.AntibodyChain.LIGHT: 'blast_vl.dat', amber_antibody.AntibodyChain.HEAVY: 'blast_vh.dat'}
        blast = {}
        for chain in blast_file_names:
            data_file_path = os.path.join(self.report_folder, blast_file_names[chain])
            with open(data_file_path, 'rb') as pfile:
                pdata = pickle.load(pfile)
            blast[chain] = amber_align.BlastAlignment(pdata[0], pdata[1])
        return blast

    def generate_canonical_sequence_plots(self, blast_alignments_dict, row_size=50, figwidth=10.0, aspect=0.25, dpi=100,
                                          text_y_offset=2.5, y_title_offset=1.0):
        chain_ids = {amber_antibody.AntibodyChain.LIGHT: ('(light chain)', '04', 'light'), amber_antibody.AntibodyChain.HEAVY: ('(heavy chain)', '06', 'heavy')}
        use_name = amber.standardize_file_stem(self.name)
        file_names = {}
        for chain in chain_ids:
            file_names[chain] = '{}_csa_{}_{}.png'.format(chain_ids[chain][1], use_name, chain_ids[chain][2])
        self.antibody.generate_canonical_sequence_plots(blast_alignments_dict, file_names, save_folder=self.report_folder, row_size=row_size,
                    figwidth=figwidth, aspect=aspect, dpi=dpi, text_y_offset=text_y_offset, y_title_offset=y_title_offset)
        return

    def generate_organism_pie_charts(self, blast_alignments_dict, explode_organism='Homo sapiens', figsize=(10.0, 10.0),
                                     dpi=100, nlegend=5, explode_level=0.1):
        chain_ids = {amber_antibody.AntibodyChain.LIGHT:('05', 'light'), amber_antibody.AntibodyChain.HEAVY:('07', 'heavy')}
        use_name = amber.standardize_file_stem(self.name)
        file_names = {}
        for chain in chain_ids:
            file_names[chain] = '{}_organisms_{}_{}.png'.format(chain_ids[chain][0], use_name, chain_ids[chain][1])
        self.antibody.generate_organism_pie_charts(blast_alignments_dict, file_names, explode_organism=explode_organism, save_folder=self.report_folder,
                                    figsize=figsize, dpi=dpi, nlegend=nlegend, explode_level=explode_level)
        return

    def generate_immunogenicity_plots(self, row_size=50, figwidth=10.0, aspect=0.25, dpi=100, text_y_offset=0.8):
        chain_ids = {amber_antibody.AntibodyChain.LIGHT: ('14', 'light'), amber_antibody.AntibodyChain.HEAVY: ('16', 'heavy')}
        use_name = amber.standardize_file_stem(self.name)
        file_names = {}
        for chain in chain_ids:
            file_names[chain] = '{}_epitopes_{}_{}.png'.format(chain_ids[chain][0], use_name, chain_ids[chain][1])
        self.antibody.generate_immunogenicity_plots(file_names, save_folder=self.report_folder, row_size=row_size,
                                        figwidth=figwidth, aspect=aspect, dpi=dpi, text_y_offset=text_y_offset)
        return

    def generate_epitope_reports(self):
        chain_ids = {amber_antibody.AntibodyChain.LIGHT: ('15', 'light'), amber_antibody.AntibodyChain.HEAVY: ('17', 'heavy')}
        use_name = amber.standardize_file_stem(self.name)
        result = self.antibody.generate_epitope_reports()
        for chain in result:
            file_name = '{}_epitope_report_{}_{}'.format(chain_ids[chain][0], use_name, chain_ids[chain][1])
            file_path = os.path.join(self.report_folder, file_name)
            with open(file_path, 'w') as rfile:
                rfile.write(result[chain])
        return

    def generate_solvent_accessibility_plots(self, figwidth=10.0, aspect=0.2, dpi=100, row_size=25, text_y_offset=10.0, y_title_offset=1.0):
        light_id = self.antibody.structure_chain_id[amber_antibody.AntibodyChain.LIGHT][0]
        heavy_id = self.antibody.structure_chain_id[amber_antibody.AntibodyChain.HEAVY][0]
        use_name = amber.standardize_file_stem(self.name)
        chain_ids = {light_id: ('12', 'light'), heavy_id: ('13', 'heavy')}
        file_names = {}
        titles = {}
        sa_results = self.antibody.compute_solvent_accessibility()
        for chain in sa_results:
            file_names[chain] = '{}_sa_{}_{}.png'.format(chain_ids[chain][0], use_name, chain_ids[chain][1])
            titles[chain] = 'Solvent Accesibility Plot: {} (chain: {})'.format(self.name, chain)
        self.antibody.plot_solvent_accessibility(sa_results, save_folder=self.report_folder, figwidth=figwidth, aspect=aspect,
                                        y_title_offset=y_title_offset, dpi=dpi, row_size=row_size, text_y_offset=text_y_offset,
                                        title_dict=titles, manual_file_name_dict=file_names)
        return sa_results

    def generate_cdr_hotspot_reports(self, field_format='{:10}'):
        chains = {}
        chains['Light Chain'] = self.antibody.scan_cdr_residues(amber_antibody.AntibodyChain.LIGHT)
        chains['Heavy Chain'] = self.antibody.scan_cdr_residues(amber_antibody.AntibodyChain.HEAVY)
        use_name = amber.standardize_file_stem(self.name)
        result_list = []
        for chain in chains:
            data_dict = chains[chain]
            result_list.append('{}: '.format(chain))
            for ncdr in data_dict:
                line =  '      CDR {:d}:  {:30}'.format(ncdr, data_dict[ncdr]['cdr'])
                for residue in data_dict[ncdr]:
                    if residue == 'cdr':
                        continue
                    line += '    {}: '.format(residue)
                    poslist = []
                    for pos in data_dict[ncdr][residue]:
                        poslist.append(str(pos+1))
                    if len(poslist) == 0:
                        line += field_format.format('None')
                    else:
                        line += field_format.format(', '.join(poslist))
                result_list.append(line)
            result_list.append(' ')
        result_text = '\n'.join(result_list)
        file_path = os.path.join(self.report_folder, '08_cdr_residues_{}.txt'.format(use_name))
        with open(file_path, 'w') as rfile:
            rfile.write(result_text)
        return

    def generate_deamidation_reports(self, field_format='{:>15}', line_items=4):
        chains = {}
        chains['Light Chain'] = self.antibody.scan_modification_sites(amber_antibody.AntibodyChain.LIGHT)
        chains['Heavy Chain'] = self.antibody.scan_modification_sites(amber_antibody.AntibodyChain.HEAVY)
        use_name = amber.standardize_file_stem(self.name)
        result_list = []
        for chain in chains:
            data_list = chains[chain]['deamidation']
            line_list = ['{}: '.format(chain)]
            nitems = 0
            for item in data_list:
                field = '{:d}-{}-{:d}'.format(item[0]+1, item[1], item[0]+len(item[1]))
                line_list.append(field_format.format(field))
                nitems += 1
                if nitems == line_items:
                    line = ''.join(line_list)
                    result_list.append(line)
                    line_list = ['             ']
                    nitems = 0
            if nitems > 0:
                line = ''.join(line_list)
                result_list.append(line)
            result_list.append(' ')
        result_text = '\n'.join(result_list)
        file_path = os.path.join(self.report_folder, '09_deamidation_{}.txt'.format(use_name))
        with open(file_path, 'w') as rfile:
            rfile.write(result_text)
        return

    def generate_o_glycosylation_reports(self, field_format='{:>15}', line_items=4):
        chains = {}
        chains['Light Chain'] = self.antibody.scan_modification_sites(amber_antibody.AntibodyChain.LIGHT)
        chains['Heavy Chain'] = self.antibody.scan_modification_sites(amber_antibody.AntibodyChain.HEAVY)
        use_name = amber.standardize_file_stem(self.name)
        result_list = []
        for chain in chains:
            data_list = chains[chain]['o-linked-glycosylation']
            line_list = ['{}: '.format(chain)]
            nitems = 0
            for item in data_list:
                field = '{:d}-{}-{:d}'.format(item[0]+1, item[1], item[0]+len(item[1]))
                line_list.append(field_format.format(field))
                nitems += 1
                if nitems == line_items:
                    line = ''.join(line_list)
                    result_list.append(line)
                    line_list = ['             ']
                    nitems = 0
            if nitems > 0:
                line = ''.join(line_list)
                result_list.append(line)
            result_list.append(' ')
        result_text = '\n'.join(result_list)
        file_path = os.path.join(self.report_folder, '10_o_glycosylation_{}.txt'.format(use_name))
        with open(file_path, 'w') as rfile:
            rfile.write(result_text)
        return

    def generate_n_glycosylation_reports(self, field_format='{:>15}', line_items=4):
        chains = {}
        chains['Light Chain'] = self.antibody.scan_modification_sites(amber_antibody.AntibodyChain.LIGHT)
        chains['Heavy Chain'] = self.antibody.scan_modification_sites(amber_antibody.AntibodyChain.HEAVY)
        use_name = amber.standardize_file_stem(self.name)
        result_list = []
        for chain in chains:
            data_list = chains[chain]['n-linked-glycosylation']
            line_list = ['{}: '.format(chain)]
            nitems = 0
            if len(data_list) == 0:
                line_list.append('None')
            for item in data_list:
                field = '{:d}-{}-{:d}'.format(item[0]+1, item[1], item[0]+len(item[1]))
                line_list.append(field_format.format(field))
                nitems += 1
                if nitems == line_items:
                    line = ''.join(line_list)
                    result_list.append(line)
                    line_list = ['             ']
                    nitems = 0
            if nitems > 0 or len(data_list) == 0:
                line = ''.join(line_list)
                result_list.append(line)
            result_list.append(' ')
        result_text = '\n'.join(result_list)
        file_path = os.path.join(self.report_folder, '11_n_glycosylation_{}.txt'.format(use_name))
        with open(file_path, 'w') as rfile:
            rfile.write(result_text)
        return

    def generate_human_alignment_plots(self):
        chain_ids = {amber_antibody.AntibodyChain.LIGHT: ('18', 'light'), amber_antibody.AntibodyChain.HEAVY: ('19', 'heavy')}
        use_name = amber.standardize_file_stem(self.name)
        for chain in chain_ids:
            file_name = '{}_human_alignment_{}_{}'.format(chain_ids[chain][0], use_name, chain_ids[chain][1])
            file_path = os.path.join(self.report_folder, file_name)
            self.antibody.generate_human_alignment_plot(chain, save_path=file_path)
        return



if __name__ == '__main__':
    from biotite.sequence.io.fasta import FastaFile
    from biotite.sequence import ProteinSequence
    sequence_file = '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/adalimumab_FASTA.txt'
    structure_file = '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/4nyl.pdb'
    sandbox = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox'
    sequence_data = FastaFile.read(sequence_file)
    sequences = []
    for header, sequence in sequence_data.items():
        sequences.append(ProteinSequence(sequence))
    report = AntibodyReport('Adalimumab', sequences[0], sequences[1], sandbox)
    report.load_structure(structure_file, 'L', 'H')
    report.generate_introduction()
    report.generate_formatted_sequences()
    #report.generate_blast_alignments()
    blast = report.reload_blast_alignments()
    report.generate_canonical_sequence_plots(blast)
    report.generate_organism_pie_charts(blast)
    report.generate_immunogenicity_plots()
    report.generate_solvent_accessibility_plots()
    report.generate_cdr_hotspot_reports()
    report.generate_deamidation_reports()
    report.generate_o_glycosylation_reports()
    report.generate_n_glycosylation_reports()
    report.generate_epitope_reports()
    report.generate_human_alignment_plots()