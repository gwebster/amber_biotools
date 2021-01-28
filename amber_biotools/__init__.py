__author__ = 'Amber Biology'

from enum import IntEnum
import re, math

class AntibodyChain(IntEnum):
    LIGHT = 0
    HEAVY = 1

def standardize_file_stem(raw_file_stem):
    file_name = raw_file_stem.lower()
    try:
        last_dot = file_name.rindex('.')
        file_name = file_name[:last_dot]
    except:
        pass
    file_name = re.sub(r'\s+', '_', file_name)
    file_name = re.sub(r'[^a-zA-Z0-9_.]', '_', file_name)
    file_name = re.sub(r'_+', '_', file_name)
    if file_name[-1] == '_':
        file_name = file_name[:-1]
    return file_name

def format_sequence(sequence, row_size=50, add_numbers=True, highlight=[]):
    result = []
    seq_len = len(sequence)
    seg_fmt = ' {{:<{:d}}} '.format(row_size)
    numlen = int(math.log10(seq_len) + 1)
    pos_fmt = '{{:{:d}d}}'.format(numlen)
    blank = ' ' * numlen
    if add_numbers:
        fmt = pos_fmt + seg_fmt + pos_fmt
        ulfmt = blank + seg_fmt + blank
    else:
        fmt = seg_fmt
        ulfmt = seg_fmt
    pos = 0
    while pos < seq_len:
        end = pos + row_size
        if end >= seq_len:
            segment = str(sequence[pos:])
            end = seq_len
        else:
            segment = str(sequence[pos:end])
        ul_list = []
        for npos in range(pos, end):
            if npos+1 in highlight:
                ul_list.append('_')
            else:
                ul_list.append(' ')
        #ul_list.append('\n')
        ul_segment = ''.join(ul_list)
        if len(highlight) > 0:
            result.append(ulfmt.format(ul_segment))
        if add_numbers:
            result.append(fmt.format(pos+1, segment, end))
        else:
            result.append(fmt.format(segment))
        pos += row_size
    return '\n'.join(result)




if __name__ == '__main__':
    from biotite.sequence import ProteinSequence
    print(standardize_file_stem('My structure   -001 _____'))
    ab_light = ProteinSequence('DIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQRYNRAPYTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGE')
    print(format_sequence(ab_light, highlight=[1, 2, 3, 4, 49, 50, 51, 52], add_numbers=True))
    #print(format_sequence(ab_light, underline=[]))







