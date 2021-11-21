def MakeSeqDict(FilePath, no_enter=False):
    seq_dict = {}
    f = open(FilePath, 'r')
    F = f.read()
    f.close()
    seq_list = F.split('>')
    seq_list.remove('')
    for unit in seq_list:
        header = unit.split('\n')[0]
        if no_enter:
            sequence = ''.join(unit.split('\n')[1:]).replace(' ', '')
        else:
            sequence = '\n'.join(unit.split('\n')[1:]).replace(' ', '')
        seq_dict[header] = sequence
    return seq_dict
