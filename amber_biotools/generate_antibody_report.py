__author__ = 'Amber Biology'

from biotite.sequence.io.fasta import FastaFile
from biotite.sequence import ProteinSequence
from amber_biotools.antibody_report import AntibodyReport
import os

# --- fill out this section ------------------------

# enter the antibody name/identifier
antibody_name = 'CAMPATH-1G IGG2B'

# enter existing folder path into which antibody report folder will be created
antibody_report_folder = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox'

# paste the light chain sequence here
light_chain_sequence = 'DIKMTQSPSFLSASVGDRVTLNCKASQNIDKYLNWYQQKLGESPKLLIYNTNNLQTGIPSRFSGSGSGTDFTLTISSLQPEDVATYFCLQHISRPRTFGTGTKLELKRANAAPTVSIFPPSTEQLATGGASVVCLMNKFYPRDISVKWK'


# paste the heavy chain sequence here
heavy_chain_sequence = 'EVKLLESGGGLVQPGGSMRLSCAGSGFTFTDFYMNWIRQPAGKAPEWLGFIRDKAKGYTTEYNPSVKGRFTISRDNTQNMLYLQMNTLRAEDTATYYCAREGHTAAPFDYWGQGVMVTVSSAQTTAPSVYPLAPGCGDTTSSTVTLGCLVK'


# paste the path to the PDB structure file here (if one is available)
# and specify the light and heavy chain ids
structure_file_path = '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/1bfo.pdb'
structure_light_chain_id = 'A'
structure_heavy_chain_id = 'B'

if len(structure_file_path) > 0 and os.path.exists(structure_file_path):
    structure_available = True
else:
    structure_available = False



# --- run this section ------------------------------------

if __name__ == '__main__':
    print('\nAssembling antibody ...')
    light_chain = ProteinSequence(light_chain_sequence)
    heavy_chain = ProteinSequence(heavy_chain_sequence)
    print('Initializing report ...')
    report = AntibodyReport(antibody_name, light_chain, heavy_chain, antibody_report_folder)
    if structure_available:
        report.load_structure(structure_file_path, structure_light_chain_id, structure_heavy_chain_id)
    # comment out this next step if loading pre-saved BLAST from file
    print('Running BLAST ...')
    blast = report.generate_blast_alignments()
    # uncomment this next step if re-running report and loading pre-saved BLAST from file
    #blast = report.reload_blast_alignments()
    print('Generating introduction ...')
    report.generate_introduction()
    print('Generating formatted sequences ...')
    report.generate_formatted_sequences()
    print('Generating canonical sequence plots ...')
    report.generate_canonical_sequence_plots(blast)
    print('Generating organism pie charts ...')
    report.generate_organism_pie_charts(blast)
    print('Generating CDR hotspot summaries ...')
    report.generate_cdr_hotspot_reports()
    print('Generating deamidation summaries ...')
    report.generate_deamidation_reports()
    print('Generating glycosylation summaries ...')
    report.generate_o_glycosylation_reports()
    report.generate_n_glycosylation_reports()
    if structure_available:
        print('Generating solvent-accessible surface plots ...')
        report.generate_solvent_accessibility_plots()
    print('Generating immunogenicity plots ...')
    report.generate_immunogenicity_plots()
    print('Generating epitope summaries ...')
    report.generate_epitope_reports()
    print('Generating human alignment plots ...')
    report.generate_human_alignment_plots()
    print('Report successfully generated.\n')