__author__ = 'Amber Biology'

from enum import Enum
import copy
from biotite.sequence import Alphabet

class SequenceType(Enum):
    NUCLEOTIDE = 1
    PROTEIN = 2
    CUSTOM = 3
    UNKNOWN = 4

sequence_type_names = {SequenceType.NUCLEOTIDE:'Nucleotide', SequenceType.PROTEIN:'Protein',
                       SequenceType.CUSTOM:'Custom', SequenceType.UNKNOWN:'Unknown'}

sequence_map = {SequenceType.NUCLEOTIDE: {1:{'A': 'ADE', 'C': 'CYT', 'T': 'THY', 'G': 'GUA', 'U': 'URA', 'I':'INO'}, 3:{} },
                SequenceType.PROTEIN: {1:{'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY', 'H': 'HIS',
                                    'I': 'ILE', 'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN',
                                    'R': 'ARG', 'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR', }, 3:{}},
                SequenceType.CUSTOM:{1:{}, 3:{}}, SequenceType.UNKNOWN:{1:{}, 3:{}}}

modified_nucleotides = [('P', 'PSU'), ('K', '1MG'), ('#', '2MG'), ('D', 'H2U'), ('L', 'M2G'), ('?', '5MC'), ('T', '5MU'),
                        ('"', '1MA'), ("'", '3MC'), ('7', '7MC')]

modified_nucleotide_sequence = {'PSU':'U', '1MG':'G', '2MG':'G', 'H2U':'U', 'M2G':'G', '5MC':'C', '5MU':'U', '1MA':'A', '3MC':'C', '7MC':'C',
                                'INO':'A', 'OMC':'C', 'OMG':'G', 'YYG':'G', '7MG':'G'}

backbone_atoms = {SequenceType.NUCLEOTIDE:["P", "O1P", "O2P", "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", "C3'", "O3'", "O2'"],
                  SequenceType.PROTEIN:["N", "CA", "C", "O"]}




class SequenceMap:

    def __init__(self):
        self.build_sequence_map()

    def reset(self):
        self.build_sequence_map()

    def build_sequence_map(self):
        self.map = {}
        self.map = copy.deepcopy(sequence_map)
        for seq_type in self.map:
            for one in self.map[seq_type][1]:
                three = self.map[seq_type][1][one]
                self.map[seq_type][3][three] = one

    def add(self, sequence_type, one_letter, three_letter):
        if sequence_type == SequenceType.UNKNOWN:
            return
        self.map[sequence_type][1][one_letter] = three_letter
        self.map[sequence_type][3][three_letter] = one_letter

    def add_modified_nucleotides(self):
        for nuc in modified_nucleotides:
            self.add(SequenceType.NUCLEOTIDE, nuc[0], nuc[1])

    def get_map(self, sequence_type, symbol):
        if symbol in self.map[sequence_type][1]:
            one_letter = symbol
            three_letter = self.map[sequence_type][1][symbol]
        elif symbol in self.map[sequence_type][3]:
            one_letter = self.map[sequence_type][3][symbol]
            three_letter = symbol
        else:
            if len(symbol) == 1:
                return (symbol, 'XXX')
            else:
                return ('X', symbol)
        return (one_letter, three_letter)

    def detect_sequence_type(self, sequence):
        scores = {SequenceType.PROTEIN:0, SequenceType.NUCLEOTIDE:0, SequenceType.CUSTOM:0, SequenceType.UNKNOWN:0}
        match = {SequenceType.PROTEIN:0.0, SequenceType.NUCLEOTIDE:0.0, SequenceType.CUSTOM:0.0, SequenceType.UNKNOWN:0.0}
        length_key = len(sequence[0])
        if length_key != 1 and length_key != 3:
            return
        for seq_type in SequenceType:
            valid_residues = self.get_valid_residues(seq_type)
            residues = list(valid_residues[1]) + list(valid_residues[3])
            unique = set([])
            for residue in sequence:
                if residue in residues:
                    unique.add(residue)
                    scores[seq_type] += 2
            if len(residues) == 0:
                match[seq_type] = 0.0
            else:
                match[seq_type] = len(unique) / len(residues)
        if scores[SequenceType.PROTEIN] == 0 and scores[SequenceType.NUCLEOTIDE] == 0 and scores[SequenceType.CUSTOM] == 0:
            return SequenceType.UNKNOWN
        if scores[SequenceType.PROTEIN] == scores[SequenceType.NUCLEOTIDE] and scores[SequenceType.PROTEIN] > scores[SequenceType.CUSTOM]:
            if scores[SequenceType.PROTEIN] == scores[SequenceType.NUCLEOTIDE]:
                if match[SequenceType.PROTEIN] > match[SequenceType.NUCLEOTIDE]:
                    return SequenceType.PROTEIN
                else:
                    return SequenceType.NUCLEOTIDE
            else:
                return max(scores, key=scores.get)
        else:
            return max(scores, key=scores.get)


    def is_valid_sequence(self, sequence_type, sequence, threshold=100):
        residues = self.map[sequence_type][1].keys()
        sequence_length = len(sequence)
        count = 0
        for n in range(0, sequence_length):
            if sequence[n] in residues:
                count += 1
        if 100.0 * (count/sequence_length) >= threshold:
            return True
        else:
            return False

    def get_valid_residues(self, sequence_type):
        result = {}
        result[1] = list(self.map[sequence_type][1].keys())
        result[3] = list(self.map[sequence_type][3].keys())
        return result

    def generate_alphabet(self, sequence_type):
        letters = list(self.map[sequence_type][1].keys())
        return Alphabet(letters)



if __name__ == '__main__':
    from biotite.sequence import ProteinSequence
    sm = SequenceMap()
    print(sm.get_map(SequenceType.NUCLEOTIDE, 'A'))
    print(sm.get_map(SequenceType.NUCLEOTIDE, 'GUA'))
    print(sm.get_map(SequenceType.PROTEIN, 'S'))
    print(sm.get_map(SequenceType.PROTEIN, 'TYR'))
    print(sm.get_map(SequenceType.PROTEIN, 'B'))
    print(sm.get_map(SequenceType.PROTEIN, 'PCA'))
    seq = ProteinSequence('QIVLSQSPAILSASPGEKVTMTCRASSSVSYIHWFQQKXPGSSPKPWIYA')
    print(sm.is_valid_sequence(SequenceType.PROTEIN, seq, threshold=90))
    print(sm.is_valid_sequence(SequenceType.NUCLEOTIDE, seq))
    print(sm.get_valid_residues(SequenceType.PROTEIN))
    print(sm.get_valid_residues(SequenceType.NUCLEOTIDE))
    sm.add_modified_nucleotides()
    print(sm.get_valid_residues(SequenceType.NUCLEOTIDE))
    print(sm.detect_sequence_type(seq))
    seq = ['URA', 'ADE', 'THY', 'CYT', 'THY', 'GUA', 'ADE', 'CYT', 'GUA']
    print(sm.detect_sequence_type(seq))
    seq = 'TACGTACGTACGTGGCTAGCATCGATCGATC'
    print(sm.detect_sequence_type(seq))
    print(sm.generate_alphabet(SequenceType.NUCLEOTIDE))


