import os
import filelist
from tqdm import tqdm


def make_empty_species_sequence_dict(gene_list, alignment_directory):
    species_list = []

    f = open('{0}/{1}'.format(alignment_directory, gene_list[0]), 'r')
    F = f.read()
    f.close()

    species_sequence_list = F.split('>')
    species_sequence_list.remove('')

    for i in species_sequence_list:
        species_name = i.split('\n')[0]  # Edit according to the format of header
        species_list.append(species_name)

    species_sequence_dict = {}

    for species in species_list:
        species_sequence_dict[species] = ''

    return species_sequence_dict


def concatenate(empty_species_sequence_dict, gene_list, alignment_directory):
    for file_name in tqdm(gene_list):
        f = open('{0}/{1}'.format(alignment_directory, file_name), 'r')
        F = f.read()
        f.close()

        species_sequence_list2 = F.split('>')
        species_sequence_list2.remove('')

        for i in species_sequence_list2:
            species_name2 = i.split('\n')[0]  # Edit according to the format of header
            sequence2 = '\n'.join(i.split('\n')[1:])
            empty_species_sequence_dict[species_name2] += sequence2

    U_list = []  # numbers of columns that have unusual amino acid
    for key in empty_species_sequence_dict:
        for i in range(len(empty_species_sequence_dict[key])):
            if empty_species_sequence_dict[key][i] == 'U' and i not in U_list:
                U_list.append(i)
    U_list.sort()

    # delete the columns with unusual amino acid
    for key in empty_species_sequence_dict:
        no_U_sequence = ''
        if U_list:
            for i in range(len(U_list)):
                if i == 0:
                    no_U_sequence += empty_species_sequence_dict[key][0:U_list[i]] + empty_species_sequence_dict[key][U_list[i]+1:U_list[i+1]]
                elif i == len(U_list) - 1:
                    no_U_sequence += empty_species_sequence_dict[key][U_list[i]+1:]
                else:
                    no_U_sequence += empty_species_sequence_dict[key][U_list[i]+1:U_list[i+1]]
            empty_species_sequence_dict[key] = no_U_sequence

    return empty_species_sequence_dict


def write(made_species_sequence_dict):
    g = open('./Alignment/Protein/Concatenated/concatenated.fas', 'w')
    for key in made_species_sequence_dict:
        header = '>' + key + '\n'
        sequence = made_species_sequence_dict[key]
        g.write(header)
        g.write(sequence)
    g.close()


if __name__ == '__main__':
    os.system('mkdir alignment_concatenated')
    gene_list = filelist.mklist_full(extension='.fas', directory='./alignment_automated')
    species_sequence_dict = make_empty_species_sequence_dict(gene_list)
    made_species_sequence_dict = concatenate(species_sequence_dict, gene_list)
    write(made_species_sequence_dict)
