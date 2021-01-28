__author__ = 'Amber Biology'

import pytest, os
from amber_biotools.structure import Structure, sequence_map
import biotite.structure as biostruc
from biotite.structure.io.pdb import PDBFile
from biotite.structure import AtomArray
from biotite.sequence import ProteinSequence, NucleotideSequence, GeneralSequence, Alphabet
from biotite.sequence.align import align_optimal
from amber_biotools.sequence_map import SequenceMap, SequenceType, backbone_atoms, modified_nucleotide_sequence

# structure file paths
structures = {
    'rituximab':'/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/rituximab.pdb',
    'pdb_1f7u':'/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/1f7u.pdb',
    'pdb_1ttt':'/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/1ttt.pdb',
    'pdb_1dnm':'/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/1dnm.pdb',
    'pdb_1eqr':'/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/structures/1eqr.pdb'}

structure_stats = {'rituximab':[4, 3498], 'pdb_1f7u':[6, 7126], 'pdb_1ttt':[21, 14573],
    'pdb_1dnm':[4, 537], 'pdb_1eqr':[9, 14323]}

# build atom array for testing
pdb_file = PDBFile.read(structures['rituximab'])
atom_array = pdb_file.get_structure()[0]

def test_banner():
    print('\nStarting Structure Tests ...')

def test_structure_init():
    print('Testing structure init ...')
    structure = Structure('My structure', atom_array)
    assert len(structure.atoms) == 3498
    assert structure.atoms[0].__str__() == '    L       3  VAL N      N        -3.788  -21.408   34.882'
    assert structure.atoms[-1].__str__() == 'HET H     411  HOH O      O        -6.262  -15.345  -10.257'
    assert list(structure.chains.keys()) == [('L', 1), ('H', 1), ('L', 2), ('H', 2)]
    assert structure.chains[('L', 1)]['start'] ==  0
    assert structure.chains[('L', 1)]['end'] == 1605
    assert len(structure.chains[('L', 1)]['sequence']) == 211
    assert structure.chains[('H', 2)]['start'] ==  3387
    assert structure.chains[('H', 2)]['end'] == 3497
    assert len(structure.chains[('H', 2)]['sequence']) == 111

def test_read_from_pdb():
    print('Testing read structure from PDB ...')
    for structure in structures:
        s = Structure.load_from_pdb('My Structure', structures[structure], alias_map=None)
        assert len(s.chains) == structure_stats[structure][0]
        assert len(s.atoms) == structure_stats[structure][1]

def test_get_sequence():
    print('Testing get sequence ...')
    s = Structure.load_from_pdb('My Structure', structures['rituximab'], alias_map=None)
    assert len(s.get_sequence('H')[0]) == 221
    assert len(s.get_sequence('H')[0]) == len(s.get_sequence('H')[2])
    assert str(s.get_sequence('H')[0][:50]) == 'XVQLQQPGAELVKPGASVKMSCKASGYTFTSYNMHWVKQTPGRGLEWIGA'
    assert s.get_sequence('H')[1] == SequenceType.PROTEIN
    assert s.get_sequence('H')[2][:20] == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    assert len(s.get_sequence('L')[0]) == 211
    assert len(s.get_sequence('L')[0]) == len(s.get_sequence('L')[2])
    assert str(s.get_sequence('L', 1)[0][:50]) == 'VLSQSPAILSASPGEKVTMTCRASSSVSYIHWFQQKPGSSPKPWIYATSN'
    assert s.get_sequence('L', 1)[1] == SequenceType.PROTEIN
    assert s.get_sequence('L', 1)[2][:20] == [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
    # test skip unknown residues option
    assert str(s.get_sequence('H', skip_unknown_residues=True)[0][:50]) == 'VQLQQPGAELVKPGASVKMSCKASGYTFTSYNMHWVKQTPGRGLEWIGAI'
    assert s.get_sequence('H', skip_unknown_residues=True)[2][:20] == [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
    # sequence with multiple chain types
    s = Structure.load_from_pdb('My Structure', structures['pdb_1ttt'], alias_map=None)
    assert str(s.get_sequence('D')[0]) == 'GCGGAUUUAXCUCAGXXGGGAGAGCXCCAGAXUXAAXAXXUGGAGXUCXUGUGXXCGXUCCACAGAAUUCGCACCA'
    assert s.chains[('D', 1)]['type'] == SequenceType.NUCLEOTIDE
    assert s.get_sequence('D', 2) == None
    assert s.chains[('D', 2)]['type'] == SequenceType.UNKNOWN

def test_sequence_map():
    print('Testing sequence map ...')
    s = Structure.load_from_pdb('My Structure', structures['pdb_1f7u'], alias_map=None)
    assert str(s.get_sequence('B', 1)[0]) == 'XUCCUCGUXXCCCAAXGGXCACGGCXXCUGGCUICGAACCAGAAGAXUXCAGGXXCAXGUCCUGGCGGGGAAGCCA'
    sequence_map.add_modified_nucleotides()
    s = Structure.load_from_pdb('My Structure', structures['pdb_1f7u'], alias_map=None)
    assert str(s.get_sequence('B', 1)[0]) == 'PUCCUCGUK#CCCAADGGDCACGGCLPCUGGCUICGAACCAGAAGADU?CAGGTPCA"GUCCUGGCGGGGAAGCCA'
    sequence_map.reset()
    s = Structure.load_from_pdb('My Structure', structures['pdb_1f7u'], alias_map=None)
    assert str(s.get_sequence('B', 1)[0]) == 'XUCCUCGUXXCCCAAXGGXCACGGCXXCUGGCUICGAACCAGAAGAXUXCAGGXXCAXGUCCUGGCGGGGAAGCCA'

def test_alias_map():
    print('Testing alias map ...')
    dnm_alias_map = {'DC': 'CYT', 'DG': 'GUA', 'DA': 'ADE', 'DT': 'THY'}
    s = Structure.load_from_pdb('My Structure', structures['pdb_1dnm'], alias_map=None)
    assert s.get_sequence('A') == None
    assert s.get_sequence('B') == None
    s = Structure.load_from_pdb('My Structure', structures['pdb_1dnm'], alias_map=dnm_alias_map)
    assert str(s.get_sequence('A')[0]) == 'CGCAAGCTGGCG'
    assert str(s.get_sequence('B')[0]) == 'CGCAAGCTGGCG'

def test_get_chain():
    print('Testing get chain ...')
    s = Structure.load_from_pdb('My Structure', structures['rituximab'], alias_map=None)
    assert len(s.get_chain('H')) == 1632
    assert len(s.get_chain('L')) == 1606
    atoms = s.get_chain('L')
    assert atoms[0].res_id == 3
    assert atoms[0].res_name == 'VAL'
    assert atoms[0].atom_name == 'N'
    assert atoms[-1].res_id == 213
    assert atoms[-1].res_name == 'CYS'
    assert atoms[-1].atom_name == 'OXT'

def test_list_chains_of_type():
    print('Testing list chains of type ...')
    s = Structure.load_from_pdb('My Structure', structures['pdb_1ttt'], alias_map=None)
    sequence_types = [SequenceType.PROTEIN, SequenceType.NUCLEOTIDE, SequenceType.UNKNOWN]
    chain_type_lists = {}
    for sequence_type in sequence_types:
        chain_type_lists[sequence_type] = s.list_chains_of_type(sequence_type)
    for chain_id in s.list_chains_of_type(sequence_type):
        assert chain_id in chain_type_lists[sequence_type]

def test_get_chain_length():
    print('Testing get chain length ...')
    s = Structure.load_from_pdb('My Structure', structures['pdb_1f7u'], alias_map=None)
    for chain in s.chains:
        start = s.chains[chain]['start']
        end = s.chains[chain]['end']
        chain_length = 1 + (end - start)
        assert s.get_chain_length(chain) == chain_length

def test_save_as_pdb():
    print('Testing save as PDB ...')
    s = Structure.load_from_pdb('My Structure', structures['rituximab'], alias_map=None)
    ab_path = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox/ab_atoms.pdb'
    if os.path.exists(ab_path):
        os.remove(ab_path)
    assert not os.path.exists(ab_path)
    s.save_as_pdb(ab_path)
    assert os.path.exists(ab_path)
    ab = Structure.load_from_pdb('My Structure', ab_path, alias_map=None)
    for chain in ab.chains:
        assert len(ab.chains[chain]) == len(s.chains[chain])

def test_get_backbone():
    print('Testing get backbone ...')
    # protein
    s = Structure.load_from_pdb('My Structure', structures['rituximab'], alias_map=None)
    backbone = s.get_backbone('L')
    for atom in backbone:
        assert atom.atom_name in backbone_atoms[SequenceType.PROTEIN]
    assert len(backbone) == 844
    assert str(backbone[0]) == '    L       3  VAL N      N        -3.788  -21.408   34.882'
    assert str(backbone[-1]) == '    L     213  CYS O      O       -26.711  -26.544  -25.697'
    # DNA
    dnm_alias_map = {'DC': 'CYT', 'DG': 'GUA', 'DA': 'ADE', 'DT': 'THY'}
    s = Structure.load_from_pdb('My Structure', structures['pdb_1dnm'], alias_map=dnm_alias_map)
    backbone = s.get_backbone('A')
    for atom in backbone:
        assert atom.atom_name in backbone_atoms[SequenceType.NUCLEOTIDE]
    assert str(backbone[0]) == "    A       1  CYT O5'    O        21.140   36.574   24.570"
    assert str(backbone[-1]) == "    A      12  GUA C1'    C        18.044   20.326  -11.644"
    # RNA
    s = Structure.load_from_pdb('My Structure', structures['pdb_1ttt'], alias_map=None)
    backbone = s.get_backbone('D')
    for atom in backbone:
        assert atom.atom_name in backbone_atoms[SequenceType.NUCLEOTIDE]
    assert str(backbone[0]) == "    D       1  G   P      P        98.728   -7.754   26.443"
    assert str(backbone[-1]) == "    D      76  A   C1'    C        98.404   -3.535    8.665"

def test_selections():
    print('Testing selections ...')
    s = Structure.load_from_pdb('My Structure', structures['rituximab'], alias_map=None)
    s.select({'chain_id':['L']})
    assert len(s.selection) == 1606
    assert len(s.selection_chains[('L', 1)]['sequence']) == 211
    s.select({'chain_id':[('H', 1)]})
    assert len(s.selection) == 1640
    assert len(s.selection_chains[('H', 1)]['sequence']) == 221
    s.clear_selection()
    assert len(s.selection) == 0
    s.select({'chain_id': ['L']})
    s.select({'chain_id': [('H', 1)]}, overwrite=False)
    assert len(s.selection) == (1640 + 1606)
    assert len(s.selection_chains[('L', 1)]['sequence']) == 211
    assert len(s.selection_chains[('H', 1)]['sequence']) == 221
    s.clear_selection()
    s.select({'chain_id': [('H', 1)], 'hetero':False})
    assert len(s.selection) == 1632
    assert len(s.selection_chains[('H', 1)]['sequence']) == 220
    s.select({'chain_id': [('H', 1)], 'residue_id': [100, 199]})
    assert len(s.selection) == 716
    assert len(s.selection_chains[('H', 1)]['sequence']) == 100
    s.select({'chain_id': [('H', 1)], 'residue_name': ['SER', 'THR', 'GLY']})
    assert len(s.selection) == 430
    assert len(s.selection_chains[('H', 1)]['sequence']) == 75
    s.select({'chain_id': [('H', 1)], 'residue_id': [100, 199], 'atom_name':['CA', 'N', 'C', 'O']})
    assert len(s.selection) == 400
    assert len(s.selection_chains[('H', 1)]['sequence']) == 100
    s.select({'chain_id': [('H', 1)], 'residue_id': [100, 199], 'element':['S']})
    assert len(s.selection) == 1
    assert len(s.selection_chains[('H', 1)]['sequence']) == 1
    s.select({'chain_id': [('H', 1)], 'element':['S']})
    assert len(s.selection) == 7
    assert len(s.selection_chains[('H', 1)]['sequence']) == 7

def test_use_selection():
    print('Testing using selections ...')
    s = Structure.load_from_pdb('My Structure', structures['rituximab'], alias_map=None)
    s.select({'chain_id': ['L']})
    selection = s.get_sequence('L', 1, use_selection=True)
    assert str(selection[0][:50]) == 'VLSQSPAILSASPGEKVTMTCRASSSVSYIHWFQQKPGSSPKPWIYATSN'
    assert selection[1] == SequenceType.PROTEIN
    assert selection[2][:20] == [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
    assert len(s.get_chain('L', use_selection=True)) == 1606
    assert s.list_chains_of_type(SequenceType.PROTEIN) == [('L', 1), ('H', 1)]
    assert s.list_chains_of_type(SequenceType.PROTEIN, use_selection=True) == [('L', 1)]
    s.select({'chain_id': ['H']}, overwrite=False)
    assert s.list_chains_of_type(SequenceType.PROTEIN) == [('L', 1), ('H', 1)]
    for chain in s.selection_chains:
        start = s.selection_chains[chain]['start']
        end = s.selection_chains[chain]['end']
        chain_length = 1 + (end - start)
        assert s.get_chain_length(chain, use_selection=True) == chain_length
    dnm_alias_map = {'DC': 'CYT', 'DG': 'GUA', 'DA': 'ADE', 'DT': 'THY'}
    s = Structure.load_from_pdb('My Structure', structures['pdb_1dnm'], alias_map=None)
    assert s.get_sequence('A') == None
    s.select({'chain_id': ['A']})
    assert s.get_sequence('A', use_selection=True) == None
    s.run_alias_map_on_selection(dnm_alias_map)
    assert str(s.get_sequence('A', use_selection=True)[0]) == 'CGCAAGCTGGCG'

def test_edit():
    print('Testing edit ...')
    s = Structure.load_from_pdb('My Structure', structures['pdb_1ttt'], alias_map=None)
    s.select({'chain_id':('D', 1), 'residue_id':10})
    for n in range(0, len(s.selection)):
        nmap = s.selection_map[n]
        assert s.selection[n].hetero == True
        assert s.selection[n].res_name == '2MG'
        assert s.atoms[nmap].hetero == True
        assert s.atoms[nmap].res_name == '2MG'
    s.edit_selection({'hetero':False, 'residue_name':'GUA'})
    for n in range(0, len(s.selection)):
        nmap = s.selection_map[n]
        assert s.selection[n].hetero == False
        assert s.selection[n].res_name == 'GUA'
        assert s.atoms[nmap].hetero == False
        assert s.atoms[nmap].res_name == 'GUA'
    # test edit selection only
    s.select({'chain_id': ('D', 1), 'residue_id': 16})
    for n in range(0, len(s.selection)):
        nmap = s.selection_map[n]
        assert s.selection[n].hetero == True
        assert s.selection[n].res_name == 'H2U'
        assert s.atoms[nmap].hetero == True
        assert s.atoms[nmap].res_name == 'H2U'
    s.edit_selection({'hetero':False, 'residue_name':'URA'}, edit_selection_only=True)
    for n in range(0, len(s.selection)):
        nmap = s.selection_map[n]
        assert s.selection[n].hetero == False
        assert s.selection[n].res_name == 'URA'
        assert s.atoms[nmap].hetero == True
        assert s.atoms[nmap].res_name == 'H2U'

#@pytest.mark.skip()
def test_solvent_accessibility():
    print('Testing solvent_accessibility ...')
    s = Structure.load_from_pdb('My Structure', structures['rituximab'], alias_map=None)
    sa = s.map_solvent_accessibility(['H', 'L'])
    assert len(sa) == 2
    assert sa['H']['residue_id'][0] == 2
    assert sa['H']['score'][0] == pytest.approx(79.66533)
    assert sa['L']['residue_id'][0] == 3
    assert sa['L']['score'][0] == pytest.approx(120.10853)
    s = Structure.load_from_pdb('My Structure', structures['pdb_1f7u'], alias_map=None)
    sa = s.map_solvent_accessibility(['A', 'B'])
    assert len(sa) == 2
    assert sa['A']['residue_id'][0] == 2
    assert sa['A']['score'][0] == pytest.approx(103.47849)
    assert sa['B']['residue_id'][0] == 902
    assert sa['B']['score'][0] == pytest.approx(262.43134)


def test_list_residue_types():
    print('Testing list residue types ...')
    s = Structure.load_from_pdb('My Structure', structures['rituximab'], alias_map=None)
    rt= ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PCA', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    assert s.list_residue_types('H') == rt
    s = Structure.load_from_pdb('My Structure', structures['pdb_1ttt'], alias_map=None)
    rt = ['1MA', '2MG', '5MC', '5MU', '7MG', 'ADE', 'CYT', 'GUA', 'H2U', 'M2G', 'OMC', 'OMG', 'PSU', 'URA', 'YYG']
    assert s.list_residue_types('D') == rt

def test_format_modified_nucleotides():
    print('Testing format modified nucleotides ...')
    for test_rna in [('pdb_1f7u', ('B',1)), ('pdb_1ttt', ('D',1))]:
        test_pdb = test_rna[0]
        test_chain = test_rna[1]
        for use_selection in [False, True]:
            s1 = Structure.load_from_pdb('My Structure', structures[test_pdb], alias_map=None)
            if use_selection:
                s1.select({'chain_id':test_chain})
                chain_set = s1.selection_chains[test_chain]['sequence']
                atom_set = s1.selection
            else:
                chain_set = s1.chains[test_chain]['sequence']
                atom_set = s1.atoms
            for position in chain_set:
                atom_start = chain_set[position][2]
                atom_end = chain_set[position][3]
                res_name = chain_set[position][1]
                if res_name in modified_nucleotide_sequence:
                    if not res_name == 'INO':
                        assert chain_set[position][0] == 'X'
                        for natom in range(atom_start, atom_end+1):
                            assert atom_set[natom].hetero == True
            if use_selection:
                s1.format_modified_nucleotides(test_chain, use_selection=True)
            else:
                s1.format_modified_nucleotides(test_chain)
            for position in chain_set:
                atom_start = chain_set[position][2]
                atom_end = chain_set[position][3]
                res_name = chain_set[position][1]
                if res_name in modified_nucleotide_sequence:
                    assert chain_set[position][0] == modified_nucleotide_sequence[res_name]
                    for natom in range(atom_start, atom_end+1):
                        assert atom_set[natom].hetero == False

@pytest.mark.skip()
def test_contacts():
    print('Testing contacts ...')
    s = Structure.load_from_pdb('My Structure', structures['rituximab'], alias_map=None)
    contacts = s.get_contacts('H', 'L', 4.0, add_to_selection=True, overwrite_selection=True)
    assert len(s.atoms) == 3498
    assert len(s.selection) == 536
    assert contacts[0]['data'][('H', 1)]['atom_number'] == 1864
    assert contacts[0]['data'][('L', 1)]['atom_number'] == 663
    assert contacts[0]['distance'] == pytest.approx(3.6832654)
    assert contacts[-1]['data'][('H', 1)]['atom_number'] == 3213
    assert contacts[-1]['data'][('L', 1)]['atom_number'] == 900
    assert contacts[-1]['distance'] == pytest.approx(3.9771187)
    contacts_path = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox/rituximab_contacts.pdb'
    if os.path.exists(contacts_path):
        os.remove(contacts_path)
    s.save_as_pdb(contacts_path, use_selection=True)
    assert os.path.exists(contacts_path)
    sel = Structure.load_from_pdb('Contacts', contacts_path)
    assert len(sel.atoms) == 536










if __name__ == '__main__':
    pass