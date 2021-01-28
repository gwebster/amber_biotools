__author__ = 'Amber Biology'

import amber_biotools as amber
from biotite.sequence import ProteinSequence, NucleotideSequence, GeneralSequence, Alphabet
import biotite.sequence.align as bioalign
import biotite.sequence.graphics as graphics
import matplotlib.pyplot as plt
from biotite.application.blast import BlastWebApp
import biotite.database.entrez as entrez
from amber_biotools.sequence_map import SequenceType, sequence_type_names
import os, pickle, math

default_protein_subsitution_matrix = bioalign.SubstitutionMatrix.std_protein_matrix()
default_nucleotide_subsitution_matrix = bioalign.SubstitutionMatrix.std_nucleotide_matrix()
clustal_path = '/Applications/clustal/clustalo'

# entrez database access functions

def entrez_fetch_by_id(entrez_id_list, fasta_file_folder=None, suffix='fasta', ret_type='fasta', sequence_type=SequenceType.PROTEIN):
    db_name = entrez.get_database_name(sequence_type_names[sequence_type])
    fasta_files = entrez.fetch(entrez_id_list, fasta_file_folder, suffix=suffix, ret_type=ret_type, db_name=db_name)
    sequences = {}
    for file in fasta_files:
        try:
            fasta_file = fasta.FastaFile.read(file)
        except:
            continue
        for header, sequence in fasta_file.items():
            blank = header.index(' ')
            seq_id = header[:blank]
            try:
                title = header[blank+1:header.rindex('[')].strip()
            except:
                title = header[blank+1:]
            try:
                organism = header[header.rindex('[') + 1: header.rindex(']')]
            except:
                organism = 'unidentified'
            sequences[seq_id] = {}
            sequences[seq_id]['title'] = title
            sequences[seq_id]['organism'] = organism
            sequences[seq_id]['header'] = header
            if sequence_type == SequenceType.PROTEIN:
                sequences[seq_id]['sequence'] = ProteinSequence(sequence)
            elif sequence_type == SequenceType.NUCLEOTIDE:
                sequences[seq_id]['sequence'] = NucleotideSequence(sequence)
            else:
                sequences[seq_id]['sequence'] = GeneralSequence(sequence)
    return sequences

def entrez_search(query_list, max_hits=20, sequence_type=SequenceType.PROTEIN):
    db_name = entrez.get_database_name(sequence_type_names[sequence_type])
    search_query = None
    operator = None
    nitem = 0
    while nitem < len(query_list):
        query = entrez.SimpleQuery(query_list[nitem], query_list[nitem+1])
        if search_query == None:
            search_query = query
        if operator == 'AND':
            search_query = entrez.CompositeQuery(operator, search_query, query)
        elif operator == 'OR':
            search_query = search_query | query
        elif operator == 'NOT':
            search_query = search_query ^ query
        if nitem+2 < len(query_list):
            operator = query_list[nitem+2]
        nitem += 3
    ids = entrez.search(search_query, db_name, number=max_hits)
    items = entrez_fetch_by_id(ids, sequence_type=sequence_type)
    return (str(search_query), items)

# simple sequence alignment

class SequenceAlignment:

    def __init__(self, sequence1, sequence2, residue_limit=None, substitution_matrix=default_protein_subsitution_matrix, max_align=100, local=False):
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.substitution_matrix = substitution_matrix
        self.run_alignments(residue_limit, max_align, local)

    def run_alignments(self, residue_limit=None, max_align=100, local=False):
        self.local = local
        self.max_align = max_align
        if residue_limit == None:
            seq1 = self.sequence1
            seq2 = self.sequence2
        else:
            seq1 = self.sequence1[:residue_limit]
            seq2 = self.sequence2[:residue_limit]
        self.alignments = bioalign.align_optimal(seq1, seq2, self.substitution_matrix, local=local, max_number=max_align)

    def get_alignment_count(self):
        return len(self.alignments)

    def get(self, n_alignment):
        return self.alignments[n_alignment]

    def get_gapped_sequences(self, n_alignment):
        return self.alignments[n_alignment].get_gapped_sequences()

    def get_score(self, n_alignment):
        return self.alignments[n_alignment].score

    def get_trace(self, n_alignment):
        return self.alignments[n_alignment].trace

    def get_trace_map(self, n_alignment):
        trace_map = [[], []]
        trace = self.alignments[n_alignment].trace
        for pair_map in trace:
            query_pos = pair_map[0]
            hit_pos = pair_map[1]
            if query_pos == -1 or hit_pos == -1:
                continue
            else:
                trace_map[0].append(query_pos)
                trace_map[1].append(hit_pos)
        return trace_map

    def plot_alignment(self, n_alignment, label1, label2, figsize=(8.0, 2.5), dpi=200, save_path=None, show_numbers=True, show_line_position=False):
        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax = fig.add_subplot(111)
        graphics.plot_alignment_similarity_based(ax, self.alignments[n_alignment], matrix=self.substitution_matrix, labels=[label1, label2],
            show_numbers=show_numbers, show_line_position=show_line_position)
        fig.tight_layout()
        if save_path == None:
            plt.show()
        else:
            plt.savefig(save_path)
        return

# BLAST search and alignment

class BlastAlignment:

    @staticmethod
    def blast_sequence(sequence, save_file=None, max_results=100):
        blast_sequence = str(sequence)
        blaster = BlastWebApp('blastp', sequence)
        blaster.set_max_results(max_results)
        blaster.start()
        blaster.join()
        blaster.clean_up()
        alignments = blaster.get_alignments()
        result = (sequence, alignments)
        if not save_file == None:
            with open(save_file, 'wb') as pickle_file:
                pickle.dump(result, pickle_file)
        return BlastAlignment(sequence, alignments)

    def __init__(self, sequence, alignments):
        self.sequence = sequence
        self.alignments = alignments
        self.residue_frequency = {}
        self.organism_frequency = {}
        self.compile_blast_alignments()

    def compile_blast_alignments(self):
        chain_length = len(self.sequence)
        for n in range(0, chain_length):
            self.residue_frequency[n] = {}
            self.residue_frequency[n][self.sequence[n]] = 0
        for alignment in self.alignments:
            # compile organism
            hit = alignment.hit_definition
            try:
                organism = hit[hit.rindex('[')+1: hit.rindex(']')]
            except:
                organism = 'unidentified'
            if not organism in self.organism_frequency:
                self.organism_frequency[organism] = 0
            self.organism_frequency[organism] += 1
            # compile frequency
            trace = alignment.trace
            gapped_sequences = alignment.get_gapped_sequences()
            query_sequence = gapped_sequences[0]
            hit_sequence = gapped_sequences[1]
            for pair_map in trace:
                query_pos = pair_map[0]
                hit_pos = pair_map[1]
                if query_pos == -1 or query_sequence[query_pos] == '-':
                    continue
                if hit_pos == -1 or hit_sequence[hit_pos] == '-':
                    continue
                if not hit_sequence[hit_pos] in self.residue_frequency[query_pos]:
                    self.residue_frequency[query_pos][hit_sequence[hit_pos]] = 0
                self.residue_frequency[query_pos][hit_sequence[hit_pos]] += 1
        return

    def filter_alignments_by_score(self, score_threshold=500):
        result = []
        for nalign in range(0, len(self.alignments)):
            if self.alignments[nalign].score >= score_threshold:
                result.append((nalign, self.alignments[nalign].score, self.alignments[nalign].hit_id))
        return result

    def filter_alignments_by_organism(self, search_organism):
        search_organism = search_organism.lower()
        result = []
        for nalign in range(0, len(self.alignments)):
            hit = self.alignments[nalign].hit_definition
            try:
                organism = hit[hit.rindex('[')+1: hit.rindex(']')]
                organism_field = organism.lower()
            except:
                organism = 'unidentified'
                organism_field = 'unidentified'
            if organism_field.find(search_organism) > -1:
                result.append((nalign, organism, self.alignments[nalign].hit_id))
        return result

    def alignment_count(self):
        return len(self.alignments)

    def generate_canonical_sequence_plot(self, title, save_folder=None, row_size=50, figwidth=10.0, aspect=0.25, dpi=100,
                                         text_y_offset=2.5, y_title_offset=1.0, manual_file_name=None):
        sequence = ['X'] + list(self.sequence)
        x = list(range(0, len(self.sequence) + 1))
        y = [0.0]
        nalign = self.alignment_count()
        for n in self.residue_frequency:
            residue = self.sequence[n]
            y.append((self.residue_frequency[n][residue] / nalign) * 100.0)
        maxscore = max(y)
        ymax = maxscore + (maxscore * 0.05)
        nplots = math.ceil(len(x) / row_size)
        height = nplots * (figwidth * aspect)
        fig, ax = plt.subplots(nplots, 1, figsize=(figwidth, height), dpi=dpi)
        fig.suptitle(title, y=y_title_offset)
        xpos = 1
        nplot = 0
        #print(len(x), x)
        #print(len(y), y)
        while xpos < len(x) + 1:
            xmin = xpos - 0.5
            if xpos + row_size > len(x) + 1:
                xend = len(x)
                xmax = xpos + row_size - 0.5
            else:
                xend = xpos + row_size
                xmax = xpos + row_size - 0.5
            #print('xpos, xend =', xpos, xend)
            ax[nplot].bar(x[xpos:xend], y[xpos:xend], color='green',edgecolor='white',width=1.0,align='center')
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
                #print('dle-------', pos, dot_line_end)
                ax[nplot].axvline(x=pos, linestyle='dashed', color='lightgray')
            for pos in range(xpos, seq_tick_end):
                ax[nplot].text(pos, text_y_offset, sequence[pos], color='lightsteelblue',
                               fontweight='bold', horizontalalignment='center')
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

    def generate_organism_pie_chart(self, title, explode_organism='Homo sapiens', save_folder=None, figsize=(10.0, 10.0),
                                    dpi=100, nlegend=5, explode_level=0.2, print_title=True, manual_file_name=None):
        total = 0
        for organism in self.organism_frequency:
            total += self.organism_frequency[organism]
        percentages = []
        for organism in self.organism_frequency:
            count = self.organism_frequency[organism]
            pct = (count / total) * 100.0
            percentages.append((count, organism, pct))
        percentages = sorted(percentages, reverse=True)
        fmt = '{} ({:.2f}%)'
        counts = []
        organisms = []
        legend_labels = []
        if len(percentages) < nlegend:
            nlegend = len(percentages)
        for n in range(0, nlegend):
            data = percentages[n]
            counts.append(data[0])
            organisms.append(data[1])
            legend_labels.append(fmt.format(percentages[n][1], percentages[n][2]))
        explode = [0] * nlegend
        #explode_index = organisms.index(explode_organism)
        #if explode_index > -1:
        explode[0] = explode_level
        fig1, ax1 = plt.subplots(figsize=figsize)
        if print_title:
            fig1.suptitle(title)
        wedges, texts = ax1.pie(counts, explode=explode, labels=organisms, autopct=None, labeldistance=None,
                shadow=False, startangle=90)
        ax1.axis('equal')
        ax1.legend(wedges, legend_labels,
                  title="Organisms",
                  loc="best")
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


# multiple sequence alignment


class MultipleSequenceAlignment:

    def __init__(self, label_list, sequence_list, matrix=default_protein_subsitution_matrix, gap_penalty=-10,
                 terminal_penalty=True, distances=None, guide_tree=None):
        self.labels = label_list
        self.substitution_matrix = matrix
        self.sequences = sequence_list
        self.residue_frequency = {}
        multialign = bioalign.align_multiple(self.sequences, matrix=matrix, gap_penalty=gap_penalty,
                     terminal_penalty=terminal_penalty, distances=distances, guide_tree=guide_tree)
        self.alignments = multialign[0]

    def plot_alignment(self, figsize=(8.0, 11.0), dpi=100, save_path=None, show_numbers=True, show_line_position=False, type_based=False):
        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax = fig.add_subplot(111)
        if type_based:
            graphics.plot_alignment_type_based(ax, self.alignments, labels=self.labels,
                show_numbers=show_numbers, show_line_position=show_line_position)
        else:
            graphics.plot_alignment_similarity_based(ax, self.alignments, matrix=self.substitution_matrix, labels=self.labels,
                show_numbers=show_numbers, show_line_position=show_line_position)
        fig.tight_layout()
        if save_path == None:
            plt.show()
        else:
            plt.savefig(save_path, bbox_inches='tight')
        return

    def __str__(self, line_length=50):
        maxlen = 0
        for label in self.labels:
            if len(label) > maxlen:
                maxlen = len(label)
        label_fmt = '{{:{:d}}}'.format(maxlen)
        result = []
        gap_seqs = self.alignments.get_gapped_sequences()
        seq_len = len(gap_seqs[0])
        pos_fmt = '{{:{:d}}}'.format(int(math.log10(seq_len)+1))
        fmt = label_fmt + ' ' + pos_fmt + ' ' + ' {} ' + pos_fmt
        npos = 0
        ndash = [0] * len(self.labels)
        while npos < seq_len:
            for nseq in range(0, len(self.labels)):
                if (npos + line_length) >= seq_len:
                    end = seq_len
                else:
                    end = npos + line_length
                ndash[nseq] += gap_seqs[nseq][npos:end].count('-')
                line = fmt.format(self.labels[nseq], npos+1, gap_seqs[nseq][npos:end], end-ndash[nseq])
                result.append(line)
            result.append(' ')
            npos += line_length
        return '\n'.join(result)

    def compile_alignments(self):
        chain_length = 0
        for seq in self.sequences:
            if len(seq) > chain_length:
                chain_length = len(seq)
        for n in range(0, chain_length):
            self.residue_frequency[n] = {}
        for alignment in self.alignments:
            # compile frequency
            trace = alignment.trace
            gapped_sequences = alignment.get_gapped_sequences()
            query_sequence = gapped_sequences[0]
            hit_sequence = gapped_sequences[1]
            for pair_map in trace:
                query_pos = pair_map[0]
                hit_pos = pair_map[1]
                if query_pos == -1 or query_sequence[query_pos] == '-':
                    continue
                if hit_pos == -1 or hit_sequence[hit_pos] == '-':
                    continue
                if not hit_sequence[hit_pos] in self.residue_frequency[query_pos]:
                    self.residue_frequency[query_pos][hit_sequence[hit_pos]] = 0
                self.residue_frequency[query_pos][hit_sequence[hit_pos]] += 1
        return




if __name__ == '__main__':
    from biotite.sequence.io import fasta
    sequences = {
        'adalimumab': '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/adalimumab_FASTA.txt',
        'ab_full': '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/antibody_full_FASTA.txt',
        'emd72000': '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/emd72000_FASTA.txt',
        'inx021': '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/inx021_FASTA.txt',
        'rituximab': '/Users/gordon/Google Drive/amber/Python/amber_biotools/test_data/sequences/rituximab_FASTA.txt'}
    # load sequences
    chains = {}
    select_heavy = 0
    select_light = 1
    for sequence in sequences:
        chains[sequence] = []
        seq_data = fasta.FastaFile.read(sequences[sequence])
        for header, string in seq_data.items():
            chains[sequence].append(ProteinSequence(string))
    # align
    seq1 = chains['adalimumab'][select_heavy]
    seq2 = chains['rituximab'][select_heavy]
    limit = None
    sa = SequenceAlignment(seq1, seq2, residue_limit=limit, local=False)
    print(sa.get_alignment_count())
    for n in range(0, sa.get_alignment_count()):
    #for n in range(0, 1):
        print(sa.get_score(n))
        gs = sa.get_gapped_sequences(n)
        print(gs[0])
        print(gs[1])
        #trace = sa.get_trace_map(n)
        #print(trace[0])
        #print(trace[1])
        print()


    ab_light = 'DIVMTQSPSFLSASVGDRVTITCKASSNLGHAVAWYQQKPGKSPKLLIYSASNRYTGVPDRFSGSGSGTDFTLTISSLQPEDFADYFCQQYDDYPYTFGGGTKLEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC'
    save_file_path = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox/ab_light_blast.dat'
    #ba = BlastAlignment.blast_sequence(ab_light, save_file=save_file_path, max_results=1000)
    #with open(save_file_path, 'rb') as pfile:
    #    blast = pickle.load(pfile)
    #print(blast)

    with open(save_file_path, 'rb') as pfile:
        blast = pickle.load(pfile)
    ba = BlastAlignment(blast[0], blast[1])
    #for align in ba.alignments:
    #    gs = align.get_gapped_sequences()
    #    print()
    #    print(gs[0])
    #    print(gs[1])
    print('\n')
    print(ba.residue_frequency)
    print(ba.filter_alignments_by_score(1030))
    print(ba.filter_alignments_by_organism('mulatta'))
    print(ba.organism_frequency)
    sandbox = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox'
    ba.generate_canonical_sequence_plot('Test Canonical Sequence Plot', dpi=200, save_folder=sandbox)
    #ba.generate_organism_pie_chart('Light Chain Alignments', dpi=200, save_folder=sandbox)
    print()

    #save_folder = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox'
    #sequences = fetch_entrez_by_id(['ABO30531.1','ELL39179.1'])
    #for seq_id in sequences:
    #    print(seq_id, sequences[seq_id])
    #seq1 = sequences['ABO30531.1']['sequence']
    #seq2 = sequences['ELL39179.1']['sequence']
    #sa = SequenceAlignment(seq1, seq2, residue_limit=None, local=False)
    #fig_path = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox/ef_alignments.png'
    #sa.plot_alignment(0, 'ABO30531.1', 'ELL39179.1', figsize=(8.0, 5.0), save_path=fig_path)

    #items = entrez_search(['homo', 'Organism', "AND", "EF1A", "Gene Name", 'AND', '1-alpha', 'Title'], max_hits=20)
    #items = entrez_search(['pyrococcus', 'Organism', "AND", 'trna synthetase', 'Title', "AND", 'leucyl', 'Title'], max_hits=20)
    #items = entrez_search(['trna synthetase', 'Title', "AND", 'leucyl', 'Title'], max_hits=200)
    #print(items[0])
    #hits = items[1]
    #for item in hits:
    #    print(item, hits[item]['organism'], hits[item]['title'])

    '''synthetase_ids = ['BAA95667.1', 'GHM89476.1', 'RCH87981.1']
    synthetases = entrez_fetch_by_id(synthetase_ids)
    labels = []
    sequences = []
    for item in synthetases:
        labels.append(item)
        sequences.append(synthetases[item]['sequence'])
    ma = MultipleSequenceAlignment(labels, sequences)
    print(len(ma))
    print(ma[0])
    print(ma[0].trace)
    save_path = '/Users/gordon/Google Drive/amber/Python/amber_biotools/sandbox/synth_align.png'
    ma.plot_alignment(figsize=(8.0, 15.0), dpi=200, save_path=save_path, type_based=False)'''



