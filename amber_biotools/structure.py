__author__ = 'Amber Biology'

import amber_biotools as amber
from biotite.structure.io.pdb import PDBFile
from biotite.structure import AtomArray
import biotite.structure as biostruc
import biotite.sequence.align as bioalign
from biotite.sequence import ProteinSequence, NucleotideSequence, GeneralSequence, Alphabet
import biotite.database.rcsb as biorcsb
from amber_biotools.sequence_map import SequenceMap, SequenceType, backbone_atoms, modified_nucleotide_sequence
from amber_biotools.sequence_alignments import SequenceAlignment
import matplotlib.pyplot as plt
import numpy as np
import os, math

sequence_map = SequenceMap()
default_protein_subsitution_matrix = bioalign.SubstitutionMatrix.std_protein_matrix()
default_nucleotide_subsitution_matrix = bioalign.SubstitutionMatrix.std_nucleotide_matrix()
atom_field_map = {'chain_id':'chain_id', 'residue_id':'res_id', 'insert_code':'ins_code', 'residue_name':'res_name',
                  'hetero':'hetero', 'atom_name':'atom_name', 'element':'element', 'atom_id':'atom_id',
                  'b_factor':'b_factor', 'occupancy':'occupancy', 'charge':'charge', 'xyz':'coord'}

class Structure:

    # static methods

    @staticmethod
    def load_from_pdb(name, pdb_file_path, alias_map=None):
        pdb_file = PDBFile.read(pdb_file_path)
        atom_array = pdb_file.get_structure(model=1)
        return Structure(name, atom_array, alias_map)

    @staticmethod
    def write_to_pdb(pdb_file_path, atom_array):
        pdb_file = PDBFile()
        pdb_file.set_structure(atom_array)
        pdb_file.write(pdb_file_path)
        return

    @staticmethod
    def combine_atom_arrays(atom_array_list):
        result = AtomArray(0)
        for atom_array in atom_array_list:
            result += atom_array
        return result

    @staticmethod
    def compile_chains(atom_array):
        result = {}
        for atom in atom_array:
            chain_id = atom.chain_id
            residue_id = atom.res_id
            if not chain_id in result:
                result[chain_id] = []
            if not residue_id in result[chain_id]:
                result[chain_id].append(residue_id)
        return result

    @staticmethod
    def rcsb_sequence_query(sequence, search_scope, min_fractional_identity=0.0):
        if not search_scope in ['protein', 'dna', 'rna']:
            return
        else:
            query = biorcsb.SequenceQuery(sequence, search_scope, min_identity=min_fractional_identity)
            return sorted(biorcsb.search(query))

    @staticmethod
    def rcsb_structure_query(pdb_id, chain=None):
        query = biorcsb.StructureQuery(pdb_id, chain)
        return sorted(biorcsb.search(query))

    @staticmethod
    def fetch_from_rcsb(pdb_id, file_format, download_folder, ok_to_overwrite=False):
        if not file_format in ['pdb', 'pdbx', 'cif', 'mmcif', 'mmtf', 'fasta']:
            return
        if not os.path.exists(download_folder):
            return
        fetch_file = biorcsb.fetch(pdb_id, file_format, target_path=download_folder, overwrite=ok_to_overwrite)
        return os.path.join(download_folder, os.path.basename(fetch_file))

    # instance methods

    def __init__(self, name, atom_array, alias_map=None):
        self.name = name
        self.atoms = atom_array
        self.clear_selection()
        if alias_map == None:
            self.alias_map = None
            self.analyze()
        else:
            self.alias_map = alias_map
            self.run_alias_map(alias_map)
        return

    def parse_chain_id(self, chain_name, chain_index):
        if isinstance(chain_name, tuple):
            return chain_name
        else:
            return (chain_name, chain_index)

    def clear_selection(self):
        self.selection = []
        self.selection_map = []
        self.selection_chains = {}

    def get_selection_size(self):
        return len(self.selection)

    def add_to_selection(self, atom_list, atom_map, overwrite=True):
        if len(atom_list) == 0:
            return
        if overwrite:
            self.selection = biostruc.array(atom_list)
            self.selection_map = atom_map
        else:
            self.selection += biostruc.array(atom_list)
            self.selection_map += atom_map
        self.analyze_selection()

    def analyze_selection(self):
        self.analyze(use_selection=True)

    def run_alias_map_on_selection(self, alias_map):
        self.run_alias_map(alias_map, use_selection=True)

    def convert_selection_to_structure(self, name, alias_map=None):
        return Structure(name, self.selection, alias_map)

    def get_valid_selection_keys(self):
        return ['sphere', 'box', 'atom_name', 'atom_number', 'residue_name', 'residue_id', 'element', 'hetero']

    def select(self, select_dict, overwrite=True):
        select_list = []
        select_map = []
        # setup sphere/box to save time if they are included in selection criteria
        if 'sphere' in select_dict:
            x = select_dict['sphere'][0]
            y = select_dict['sphere'][1]
            z = select_dict['sphere'][2]
            sphere_centroid = biostruc.Atom([x, y, z])
            sphere_radius = select_dict['sphere'][3]
        if 'box' in select_dict:
            box_xmin = select_dict['box'][0] - select_dict['box'][3]
            box_xmax = select_dict['box'][0] + select_dict['box'][3]
            box_ymin = select_dict['box'][1] - select_dict['box'][4]
            box_ymax = select_dict['box'][1] + select_dict['box'][4]
            box_zmin = select_dict['box'][2] - select_dict['box'][5]
            box_zmax = select_dict['box'][2] + select_dict['box'][5]
        # chain_id, atom_name, atom_id, residue_name, residue_id, element, hetero, sphere, box
        for natom in range(0, len(self.atoms)):
            checklist = {}
            if 'chain_id' in select_dict:
                checklist['chain_id'] = False
                if not isinstance(select_dict['chain_id'], list):
                    select_dict['chain_id'] = [select_dict['chain_id']]
                for chain_id in select_dict['chain_id']:
                    chain_id = self.parse_chain_id(chain_id, 1)
                    srch_chain = chain_id
                    if natom >= self.chains[srch_chain]['start'] and natom <= self.chains[srch_chain]['end']:
                        checklist['chain_id'] = True
                        break
            if 'atom_name' in select_dict:
                checklist['atom_name'] = False
                if not isinstance(select_dict['atom_name'], list):
                    select_dict['atom_name'] = [select_dict['atom_name']]
                for atom_name in select_dict['atom_name']:
                    if self.atoms[natom].atom_name == atom_name:
                        checklist['atom_name'] = True
                        break
            if 'atom_number' in select_dict:
                checklist['atom_number'] = False
                start = select_dict['atom_number'][0]
                end = select_dict['atom_number'][1]
                if natom >= start and natom <= end:
                    checklist['atom_number'] = True
            if 'residue_name' in select_dict:
                checklist['residue_name'] = False
                if not isinstance(select_dict['residue_name'], list):
                    select_dict['residue_name'] = [select_dict['residue_name']]
                for residue_name in select_dict['residue_name']:
                    if self.atoms[natom].res_name == residue_name:
                        checklist['residue_name'] = True
                        break
            if 'residue_id' in select_dict:
                checklist['residue_id'] = False
                if not isinstance(select_dict['residue_id'], list):
                    select_dict['residue_id'] = [select_dict['residue_id']]
                if len(select_dict['residue_id']) == 1:
                    select_dict['residue_id'].append(select_dict['residue_id'][0])
                start = select_dict['residue_id'][0]
                end = select_dict['residue_id'][1]
                if self.atoms[natom].res_id >= start and self.atoms[natom].res_id <= end:
                    checklist['residue_id'] = True
            if 'element' in select_dict:
                checklist['element'] = False
                if not isinstance(select_dict['element'], list):
                    select_dict['element'] = [select_dict['element']]
                for element in select_dict['element']:
                    if self.atoms[natom].element == element:
                        checklist['element'] = True
                        break
            if 'hetero' in select_dict:
                checklist['hetero'] = False
                if self.atoms[natom].hetero == select_dict['hetero']:
                    checklist['hetero'] = True
            if 'sphere' in select_dict:
                checklist['sphere'] = False
                if biostruc.distance(sphere_centroid, self.atoms[natom]) <= sphere_radius:
                    checklist['sphere'] = True
            if 'box' in select_dict:
                checklist['box'] = False
                if self.atoms[natom].coord[0] >= box_xmin and self.atoms[natom].coord[0] <= box_xmax and \
                self.atoms[natom].coord[1] >= box_ymin and self.atoms[natom].coord[1] <= box_ymax and \
                self.atoms[natom].coord[2] >= box_zmin and self.atoms[natom].coord[2] <= box_zmax:
                    checklist['box'] = True
            # selected?
            if not False in checklist.values():
                #print('selected:', natom, self.atoms[natom])
                #new_atom = self.atoms[natom].copy()
                select_list.append(self.atoms[natom])
                select_map.append(natom)
        self.add_to_selection(select_list, select_map, overwrite)
        return

    def run_alias_map(self, alias_map, use_selection=False):
        if alias_map == None:
            return
        atom_set = self.atoms
        if use_selection:
            if self.selection == None or len(self.selection) == 0:
                return
            atom_set = self.selection
        for n in range(0, len(atom_set)):
            if atom_set[n].res_name in alias_map:
                new_name = alias_map[atom_set[n].res_name]
                new_atom = atom_set[n].copy()
                new_atom.res_name = new_name
                atom_set[n] = new_atom
        self.analyze(use_selection)
        return

    def analyze(self, use_selection=False):
        if use_selection:
            if self.selection == None or len(self.selection) == 0:
                return
            atom_set = self.selection
            self.selection_chains = {}
            chain_set = self.selection_chains
        else:
            atom_set = self.atoms
            self.chains = {}
            chain_set = self.chains
        if len(atom_set) == 0:
            return
        chain_ids = []
        chain_starts = biostruc.get_chain_starts(atom_set)
        residue_starts = biostruc.get_residue_starts(atom_set)
        # initialize chain_set and first chain key for edge case where len(residue_starts) = 1
        # i.e. the first residue in the current chain is also the last one
        current_chain_key = (atom_set[residue_starts[0]].chain_id, 1)
        chain_set[current_chain_key] = {}
        current_chain = -1
        for n in range(0, len(residue_starts)):
            if n == len(residue_starts) - 1:
                chain_set[current_chain_key]['end'] = len(atom_set) - 1
            position = residue_starts[n]
            if position in chain_starts:
                if not current_chain == -1:
                    chain_set[current_chain_key]['end'] = position - 1
                current_chain = atom_set[position].chain_id
                chain_ids.append(current_chain)
                current_chain_key = (current_chain, chain_ids.count(current_chain))
                chain_set[current_chain_key] = {}
                chain_set[current_chain_key]['type'] = position
                chain_set[current_chain_key]['start'] = position
                chain_set[current_chain_key]['end'] = -1
                chain_set[current_chain_key]['sequence'] = {}
            residue_id = atom_set[position].res_id
            residue_name = atom_set[position].res_name
            residue_marker = residue_name
            if len(residue_name) == 1:
                residue_symbol = residue_name
                residue_name = 'XXX'
            else:
                residue_symbol = 'X'
            if n+1 == len(residue_starts):
                residue_end = len(residue_starts) - 1
            else:
                residue_end = residue_starts[n+1] - 1
            chain_set[current_chain_key]['sequence'][residue_id] = [residue_symbol, residue_name, position, residue_end, residue_marker]
        # fill out residue names/symbols
        for chain in chain_set:
            try:
                sequence = []
                for position in chain_set[chain]['sequence']:
                    sequence.append(chain_set[chain]['sequence'][position][4])
                sequence_type = sequence_map.detect_sequence_type(sequence)
                chain_set[chain]['type'] = sequence_type
                if sequence_type != SequenceType.UNKNOWN:
                    for position in chain_set[chain]['sequence']:
                        map = sequence_map.get_map(sequence_type, chain_set[chain]['sequence'][position][4])
                        chain_set[chain]['sequence'][position][0] = map[0]
                        chain_set[chain]['sequence'][position][1] = map[1]
                        # remove redundant sequence marker
                        del chain_set[chain]['sequence'][position][4]
            except:
                chain_set[chain]['type'] = SequenceType.UNKNOWN
        return

    def get_chain(self, chain_name, chain_index=1, include_heteroatoms=False, use_selection=False):
        chain_id = self.parse_chain_id(chain_name, chain_index)
        atom_set = self.atoms
        chain_set = self.chains
        if use_selection:
            if self.selection == None or len(self.selection) == 0:
                return
            atom_set = self.selection
            chain_set = self.selection_chains
        atom_list = []
        for n in range(chain_set[chain_id]['start'], chain_set[chain_id]['end'] + 1):
            if include_heteroatoms:
                atom_list.append(atom_set[n])
            else:
                if not atom_set[n].hetero:
                    atom_list.append(atom_set[n])
        if len(atom_list) == 0:
            return
        else:
            return biostruc.array(atom_list)

    def get_sequence(self, chain_name, chain_index=1, skip_unknown_residues=False, convert_u=False, use_selection=False):
        chain_id = self.parse_chain_id(chain_name, chain_index)
        chain_set = self.chains
        if use_selection:
            if self.selection == None or len(self.selection) == 0:
                return
            chain_set = self.selection_chains
        chain_data = chain_set[chain_id]
        if chain_set[chain_id]['type'] == SequenceType.UNKNOWN:
            return
        sequence_list = []
        residue_id_list = []
        for residue in chain_data['sequence']:
            if skip_unknown_residues and chain_data['sequence'][residue][0] =='X':
                continue
            else:
                sequence_list.append(chain_data['sequence'][residue][0])
                residue_id_list.append(residue)
        sequence = ''.join(sequence_list)
        if convert_u:
            sequence = sequence.replace('U', 'T')
        if len(sequence) == 0:
            return
        alphabet = Alphabet(set(sequence))
        if chain_data['type'] == SequenceType.PROTEIN:
            try:
                return (ProteinSequence(sequence), SequenceType.PROTEIN, residue_id_list)
            except:
                return (GeneralSequence(alphabet, sequence), SequenceType.UNKNOWN, residue_id_list)
        elif chain_data['type'] == SequenceType.NUCLEOTIDE:
            try:
                return (NucleotideSequence(sequence), SequenceType.NUCLEOTIDE, residue_id_list)
            except:
                return (GeneralSequence(alphabet, sequence), SequenceType.UNKNOWN, residue_id_list)
        else:
            return (GeneralSequence(alphabet, sequence), SequenceType.UNKNOWN, residue_id_list)

    def list_chains_of_type(self, sequence_type, use_selection=False):
        chain_set = self.chains
        if use_selection:
            if self.selection == None or len(self.selection) == 0:
                return
            chain_set = self.selection_chains
        result = []
        for chain in chain_set:
            if chain_set[chain]['type'] == sequence_type:
                result.append(chain)
        return result

    def get_chain_length(self, chain_name, chain_index=1, use_selection=False):
        chain_id = self.parse_chain_id(chain_name, chain_index)
        chain_set = self.chains
        if use_selection:
            if self.selection == None or len(self.selection) == 0:
                return
            chain_set = self.selection_chains
        return 1 + (chain_set[chain_id]['end'] - chain_set[chain_id]['start'])

    def __str__(self):
        return self.atoms.__str__()

    def save_as_pdb(self, pdb_file_path, use_selection=False):
        atom_set = self.atoms
        if use_selection:
            if self.selection == None or len(self.selection) == 0:
                return
            atom_set = self.selection
        pdb_file = PDBFile()
        pdb_file.set_structure(atom_set)
        pdb_file.write(pdb_file_path)
        return

    def get_backbone(self, chain_name, chain_index=1, use_selection=False):
        chain_id = self.parse_chain_id(chain_name, chain_index)
        atom_set = self.atoms
        chain_set = self.chains
        if use_selection:
            if self.selection == None or len(self.selection) == 0:
                return
            atom_set = self.selection
            chain_set = self.selection_chains
        chain_data = chain_set[chain_id]
        chain_type = chain_data['type']
        if not chain_type in [SequenceType.PROTEIN, SequenceType.NUCLEOTIDE]:
            return
        n_start = chain_data['start']
        n_end = chain_data['end']
        valid_atoms = backbone_atoms[chain_type]
        atom_list = []
        for n in range(n_start, n_end):
            if not atom_set[n].hetero:
                if atom_set[n].atom_name in valid_atoms:
                    atom_list.append(atom_set[n])
        return biostruc.array(atom_list)

    def get_contacts(self, chain1_id, chain2_id, radius=4.0, use_selection=False, add_to_selection=False, overwrite_selection=False):
        chain1_id = self.parse_chain_id(chain1_id, 1)
        chain2_id = self.parse_chain_id(chain2_id, 1)
        atom_set = self.atoms
        chain_set = self.chains
        if use_selection:
            if self.selection == None or len(self.selection) == 0:
                return
            atom_set = self.selection
            chain_set = self.selection_chains
        contacts = []
        atom_list1 = []
        atom_list2 = []
        start1 = chain_set[chain1_id]['start']
        end1 = chain_set[chain1_id]['end'] + 1
        start2 = chain_set[chain2_id]['start']
        end2 = chain_set[chain2_id]['end'] + 1
        for natom1 in range(start1, end1):
            for natom2 in range(start2, end2):
                distance = biostruc.distance(atom_set[natom1], atom_set[natom2])
                if distance <= radius:
                    atom_list1.append(atom_set[natom1])
                    atom_list2.append(atom_set[natom2])
                    result = {}
                    result['distance'] = distance
                    result['data'] = {}
                    result['data'][chain1_id] = {}
                    result['data'][chain1_id]['residue_id'] = atom_set[natom1].res_id
                    result['data'][chain1_id]['residue_name'] = atom_set[natom1].res_name
                    result['data'][chain1_id]['atom_number'] = natom1
                    result['data'][chain1_id]['atom_name'] = atom_set[natom1].atom_name
                    result['data'][chain2_id] = {}
                    result['data'][chain2_id]['residue_id'] = atom_set[natom2].res_id
                    result['data'][chain2_id]['residue_name'] = atom_set[natom2].res_name
                    result['data'][chain2_id]['atom_number'] = natom2
                    result['data'][chain2_id]['atom_name'] = atom_set[natom2].atom_name
                    contacts.append(result)
        if add_to_selection:
            atom_list = atom_list1 + atom_list2
            self.add_to_selection(atom_list, overwrite_selection)
        return contacts

    def superimpose(self, self_chain_id, other, other_chain_id, use_self_selection=False,
                    substitution_matrix=default_protein_subsitution_matrix, convert_u=False,
                    use_other_selection=False, self_residue_range=None, use_backbone=True):
        self_chain_id = self.parse_chain_id(self_chain_id, 1)
        other_chain_id = self.parse_chain_id(other_chain_id, 1)
        self_sequence = self.get_sequence(self_chain_id, use_selection=use_self_selection, convert_u=convert_u)
        other_sequence = other.get_sequence(other_chain_id, use_selection=use_other_selection, convert_u=convert_u)
        if self_residue_range == None:
            self_residue_range = [self_sequence[2][0], self_sequence[2][-1]]
            #print('self_residue_range ---', self_residue_range)
        sequence_alignment = SequenceAlignment(self_sequence[0], other_sequence[0], substitution_matrix=substitution_matrix)
        trace = sequence_alignment.alignments[0].trace
        mapped_sequences = self.map_aligned_sequences(self_sequence, other_sequence, trace, self_residue_range)
        #print('mapped sequences ---', len(mapped_sequences[0]), len(mapped_sequences[1]))
        #print(mapped_sequences[0])
        #print(mapped_sequences[1])
        # assemble atom arrays
        if use_backbone:
            self_chain_atoms = self.get_backbone(self_chain_id, use_selection=use_self_selection)
            other_chain_atoms = other.get_backbone(other_chain_id, use_selection=use_other_selection)
        else:
            self_chain_atoms = self.get_chain(self_chain_id, use_selection=use_self_selection)
            other_chain_atoms = other.get_chain(other_chain_id, use_selection=use_other_selection)
        self_atoms = {}
        other_atoms = {}
        for atom in self_chain_atoms:
            if not atom.res_id in self_atoms:
                self_atoms[atom.res_id] = []
            if len(atom.ins_code) == 0:
                self_atoms[atom.res_id].append(atom)
        for atom in other_chain_atoms:
            if not atom.res_id in other_atoms:
                other_atoms[atom.res_id] = []
            if len(atom.ins_code) == 0:
                other_atoms[atom.res_id].append(atom)
        self_select = []
        other_select = []
        for n in range(0, len(mapped_sequences[0])):
            self_residue_id = mapped_sequences[0][n]
            other_residue_id = mapped_sequences[1][n]
            # skip if residue id not in atoms (e.g. a skipped residue of heteroatoms)
            if not self_residue_id in self_atoms:
                continue
            if not other_residue_id in other_atoms:
                continue
            #print('Adding:', self_residue_id, len(self_atoms[self_residue_id]), '   ', other_residue_id, len(other_atoms[other_residue_id]))
            for atom in self_atoms[self_residue_id]:
                self_select.append(atom)
            for atom in other_atoms[other_residue_id]:
                other_select.append(atom)
        #print('lengths ---', len(self_select), len(other_select))
        self_array = biostruc.array(self_select)
        other_array = biostruc.array(other_select)
        transform = biostruc.superimpose(self_array, other_array)
        # prepare results
        result = {}
        rmsd_before = biostruc.rmsd(self_array, other_array)
        rmsd_after = biostruc.rmsd(self_array, transform[0])
        result['rmsd_before'] = rmsd_before
        result['rmsd_after'] = rmsd_after
        #print('rmsd ---', rmsd_before, rmsd_after)
        # apply transform to all atoms in other chain
        unshifted_other = other.get_chain(other_chain_id, include_heteroatoms=True, use_selection=use_other_selection)
        shifted_other = biostruc.superimpose_apply(unshifted_other, transform[1])
        result['transformed_structure'] = shifted_other
        result['transform'] = transform[1]
        return result

    def map_aligned_sequences(self, get_sequence1, get_sequence2, alignment_trace, residue_range=None):
        result = [[], []]
        for map_pair in alignment_trace:
            if -1 in map_pair:
                continue
            else:
                pos1 = map_pair[0]
                pos2 = map_pair[1]
                residue1_id = get_sequence1[2][pos1]
                residue2_id = get_sequence2[2][pos2]
                if not residue_range == None:
                    if not (residue1_id >= residue_range[0] and residue1_id <= residue_range[1]):
                        continue
                    else:
                        result[0].append(residue1_id)
                        result[1].append(residue2_id)
        return result

    def edit_selection(self, edit_dict, edit_selection_only=False):
        for natom in range(0, len(self.selection)):
            new_atom = self.selection[natom].copy()
            for field_name in edit_dict:
                attr_name = atom_field_map[field_name]
                setattr(new_atom, attr_name, edit_dict[field_name])
            self.selection[natom] = new_atom
            if not edit_selection_only:
                nmap = self.selection_map[natom]
                self.atoms[nmap] = new_atom
        self.analyze()
        self.analyze(use_selection=True)
        return

    def map_solvent_accessibility(self, chain_id_list):
        for nchain in range(0, len(chain_id_list)):
            chain_id_list[nchain] = self.parse_chain_id(chain_id_list[nchain], 1)
        atom_list = []
        for chain_id in chain_id_list:
            atom_list += self.get_chain(chain_id, include_heteroatoms=False)
        atom_set = biostruc.array(atom_list)
        chains = Structure.compile_chains(atom_set)
        atom_sasa = biostruc.sasa(atom_set)
        res_sasa = biostruc.apply_residue_wise(atom_set, atom_sasa, np.sum)
        result = {}
        nscore = 0
        for chain_id in chains:
            if not chain_id in result:
                result[chain_id] = {}
                result[chain_id]['residue_id'] = []
                result[chain_id]['score'] = []
            for residue_id in chains[chain_id]:
                result[chain_id]['residue_id'].append(residue_id)
                result[chain_id]['score'].append(res_sasa[nscore])
                nscore += 1
        return result

    def plot_solvent_accessibility(self, sa_results, save_folder=None, figwidth = 10.0, aspect = 0.2, y_title_offset=1.0,
                                   dpi=100, row_size=25, text_y_offset=10.0, title_dict=None, manual_file_name_dict=None):
        for chain in sa_results:
            chain_id = (chain, 1)
            x = sa_results[chain]['residue_id']
            y = sa_results[chain]['score']
            ymax = max(y)
            nplots = math.ceil(len(x) / row_size)
            height = nplots * (figwidth * aspect)
            fig, ax = plt.subplots(nplots, 1, figsize=(figwidth, height), dpi=dpi)
            if not title_dict == None:
                fig.suptitle(title_dict[chain], y=y_title_offset)
            xpos = 0
            nplot = 0
            while xpos < len(x):
                xmin = x[xpos]
                if xpos + row_size > len(x):
                    xend = -1
                    xmax = x[xpos] + row_size - 1
                else:
                    xend = xpos + row_size
                    xmax = x[xpos + row_size] - 1
                ax[nplot].plot(x[xpos:xend], y[xpos:xend])
                ax[nplot].set_xlim(xmin, xmax)
                ax[nplot].set_ylim(0, ymax)
                ax[nplot].set_xticks(x[xpos:xend])
                ax[nplot].set_yticks(np.arange(0.0, ymax, 50.0))
                seq_ticks = []
                if xend == -1:
                    xtick_end = len(x)
                else:
                    xtick_end = xend
                for xp in range(xpos, xtick_end):
                    seq_ticks.append(self.chains[chain_id]['sequence'][x[xp]][0])
                ntick = 0
                for pos in x[xpos:xend]:
                    ax[nplot].axvline(x=pos, linestyle='dashed', color='lightgray')
                    ax[nplot].text(pos, ymax+text_y_offset, seq_ticks[ntick], color='blue',
                                   fontweight='bold', horizontalalignment='center')
                    ntick += 1
                #plt.xlabel("Residue")
                #plt.ylabel("SASA")
                xpos += row_size
                if xpos >= len(x):
                    break
                nplot += 1
            plt.tight_layout()
            if save_folder == None:
                plt.show()
            else:
                file_stem = amber.standardize_file_stem(self.name)
                if manual_file_name_dict == None:
                    save_file = 'sa_{}_{}.png'.format(file_stem, chain)
                else:
                    save_file = manual_file_name_dict[chain]
                save_path = os.path.join(save_folder, save_file)
                plt.savefig(save_path, bbox_inches='tight')
        return

    def list_residue_types(self, chain_name, chain_index=1, use_selection=False):
        chain_id = self.parse_chain_id(chain_name, chain_index)
        chain_set = self.chains
        if use_selection:
            if self.selection == None or len(self.selection) == 0:
                return []
            chain_set = self.selection_chains
        chain_data = chain_set[chain_id]
        result = []
        for position in chain_data['sequence']:
            if not chain_data['sequence'][position][1] in result:
                result.append(chain_data['sequence'][position][1])
        return sorted(result)

    def format_modified_nucleotides(self, chain_name, chain_index=1, use_selection=False):
        chain_id = self.parse_chain_id(chain_name, chain_index)
        atom_set = self.atoms
        chain_set = self.chains
        if use_selection:
            if self.selection == None or len(self.selection) == 0:
                return []
            atom_set = self.selection
            chain_set = self.selection_chains
        chain_data = chain_set[chain_id]
        for position in chain_data['sequence']:
            residue_name = chain_data['sequence'][position][1]
            if residue_name in modified_nucleotide_sequence:
                chain_data['sequence'][position][0] = modified_nucleotide_sequence[residue_name]
                atom_start = chain_data['sequence'][position][2]
                atom_end = chain_data['sequence'][position][3]
                for natom in range(atom_start, atom_end+1):
                    new_atom = atom_set[natom].copy()
                    new_atom.hetero = False
                    atom_set[natom] = new_atom
        return

    def align_with_sequence(self, other_sequence, substitution_matrix, chain_name, chain_index=1, use_selection=False,
                            skip_unknown_residues=False, convert_u=False):
        chain_id = self.parse_chain_id(chain_name, chain_index)
        self_sequence = self.get_sequence(chain_id, use_selection=use_selection, skip_unknown_residues=skip_unknown_residues, convert_u=convert_u)
        sequence_alignment = SequenceAlignment(self_sequence[0], other_sequence, substitution_matrix=substitution_matrix)
        return sequence_alignment

if __name__ == '__main__':
    pdb_path = '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/rituximab.pdb'
    #sequence_map.add_modified_nucleotides()
    #dnm_alias_map = {'DC':'CYT', 'DG':'GUA', 'DA':'ADE', 'DT':'THY'}
    s = Structure.load_from_pdb('My Structure', pdb_path, alias_map=None)
    #s.run_alias_map(alias_map)
    #s.analyze()
    #print(s.atoms)
    for chain in s.chains:
        print(chain, s.chains[chain]['type'], s.chains[chain]['start'], s.chains[chain]['end'], s.chains[chain]['sequence'])
    #print(s.get_sequence('B'))
    print(s.get_sequence('H'))
    #print(s.get_sequence('B', 2))
    #print(s.get_sequence('A', 2))
    #print(s.list_chains_of_type(SequenceType.NUCLEOTIDE))
    lc = s.get_chain('L')
    hc = s.get_chain('H')
    ab_atoms = Structure.combine_atom_arrays([lc, hc])
    ab_path = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox/ab_atoms.pdb'
    Structure.write_to_pdb(ab_path, ab_atoms)
    pdb_path1 = '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/1f7u.pdb'
    pdb_path2 = '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/1ttt.pdb'
    s1 = Structure.load_from_pdb('My Structure', pdb_path1, alias_map=None)
    s2 = Structure.load_from_pdb('My Structure', pdb_path2, alias_map=None)
    print(s1.chains)
    print(s2.chains['A', 2])
    #print(s2.get_chain('A', 4, include_heteroatoms=True))
    a4 = s2.get_chain('A', 4, include_heteroatoms=True)
    print(len(a4))
    #for atom in a4:
    #    print(atom.hetero, atom)
    #print(s2.atoms[90:100])
    #rna1 = s1.get_backbone('B')
    #print(s2.get_backbone('D'))
    #print(s2.get_chain('A'))
    #print()
    #select_dict = {'atom_number': [90, 100]}
    #select_dict = {'chain_id':['A'], 'residue_id':[390, 405], 'atom_name':['CA', 'N', 'O', 'C'], 'element':['N']}
    #select_dict = {'chain_id':[('A', 2)], 'hetero':True}
    #select_dict = {'chain_id': [('A', 1)]}
    #select_dict = {'chain_id': ['A'], 'box':[85.427, -1.143, 22.857, 5.0, 15.0, 10.0]}
    #s2.select(select_dict)
    #print(s2.selection)
    #print(s2.selection[0].coord[1])
    #sphere_path = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox/sphere.pdb'
    #s2.save_as_pdb(sphere_path, use_selection=True)
    #s2.analyze(use_selection=True)
    #print(s.selection_chains)
    #pdb_path = '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/rituximab.pdb'
    #s = Structure.load_from_pdb('My Structure', pdb_path, alias_map=None)
    #print(len(s.atoms))
    #contacts = s.get_contacts('H', 'L', 4.0, add_to_selection=True, overwrite_selection=True)
    #print(s.get_selection_size())
    #contacts_path = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox/rituximab_contacts.pdb'
    #s.save_as_pdb(contacts_path, use_selection=True)
    #ab_light = 'DIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQRYNRAPY'
    #srch = Structure.rcsb_sequence_query(ab_light, search_scope='protein', min_fractional_identity=0.92)
    #print(srch)
    #srch = Structure.rcsb_structure_query('1ttt', chain='A')
    #print(srch)
    #sandbox_path = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox'
    #pdb_file = Structure.fetch_from_rcsb('1dnm', 'pdb', sandbox_path)
    #print(pdb_file)
    #fasta_file = Structure.fetch_from_rcsb('1ttt', 'fasta', sandbox_path)
    #print(fasta_file)
    '''print('\n\n')
    pdb_path1 = '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/1f7u.pdb'
    pdb_path2 = '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/1ttt.pdb'
    s1 = Structure.load_from_pdb('My Structure', pdb_path1, alias_map=None)
    s2 = Structure.load_from_pdb('My Structure', pdb_path2, alias_map=None)
    s1.format_modified_nucleotides('B')
    s2.format_modified_nucleotides('D')
    matrix = default_nucleotide_subsitution_matrix
    s2_transform = s1.superimpose('B', s2, 'D', self_residue_range=None, use_backbone=True, substitution_matrix=matrix, convert_u=True)
    Structure.write_to_pdb('/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox/shifted_1ttt.pdb', s2_transform['transformed_structure'])
    s1.select({'chain_id':[('B', 1)], 'residue_id':[901]})
    print('\nSelection (pre-edit) ...')
    for n in range(0, len(s1.selection)):
        print(s1.selection[n], s1.selection[n].hetero)
    print('\nAtoms (pre-edit) ...')
    for n in range(0, 21):
        print(s1.atoms[n], s1.atoms[n].hetero)
    print('\n\n\nSelection (post-edit) ...')
    s1.edit_selection({'hetero':False})
    for n in range(0, len(s1.selection)):
        print(s1.selection[n], s1.selection[n].hetero)
    print('\nAtoms (post-edit) ...')
    for n in range(0, 21):
        print(s1.atoms[n], s1.atoms[n].hetero)
    pdb_path1 = '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/1f7u.pdb'
    pdb_path2 = '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/1ttt.pdb'
    s1 = Structure.load_from_pdb('My Structure', pdb_path1, alias_map=None)
    s1.format_modified_nucleotides('B')
    arg_trna = NucleotideSequence('GGCTCCGTGGCGCAATGGATAGCGCATTGGACTTCAAATTCAAAGGTTCCGGGTTCGAGTCCCGGCGGAGTCG')
    trna = NucleotideSequence('GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'.replace('U', 'T'))
    align = s1.align_with_sequence(trna, default_nucleotide_subsitution_matrix, 'B', convert_u=True)
    gs = align.get_gapped_sequences(0)
    print(gs[0])
    print(gs[1])
    pdb_path = '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/rituximab.pdb'
    s = Structure.load_from_pdb('My Structure', pdb_path, alias_map=None)
    adalimumab_h = 'EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLDMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKI'
    protseq = ProteinSequence(adalimumab_h)
    align = s.align_with_sequence(protseq, default_protein_subsitution_matrix, 'H')
    gs = align.get_gapped_sequences(0)
    print(gs[0])
    print(gs[1])
    plot_file = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox/rituximab_align.png'
    align.plot_alignment(0, 'Rituximab', 'Adalimumab', save_path=None)'''
    pdb_path = '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/rituximab.pdb'
    s = Structure.load_from_pdb('My Structure', pdb_path, alias_map=None)
    sa = s.map_solvent_accessibility(['H', 'L'])
    for chain in sa:
        for key in sa[chain]:
            print(chain, key, sa[chain][key])
    save_folder = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox/'
    s.plot_solvent_accessibility(sa, save_folder=save_folder, dpi=200)

