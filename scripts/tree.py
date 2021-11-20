import re
import os


def delete_branch_length(tree_with_branch_length):
    f = open(tree_with_branch_length, 'r')
    tree = f.readline()
    f.close()

    no_branch_length_tree = ''
    for i in tree:
        if re.match('[a-zA-Z,_();]', i):
            no_branch_length_tree += i

    g = open('./SpeciesTree/no_branch_length.nwk', 'w')
    g.write(no_branch_length_tree)
    g.close()


def root(Rscript, outgroup):
    f = open(Rscript, 'r')
    F = f.read()
    f.close()

    F1 = F.replace('>>>>', outgroup)

    g = open('/home/ldy9381/programs/ComparativeGenomics/scripts/RootTree.R', 'w')
    g.write(F1)
    g.close()

    os.system('Rscript /home/ldy9381/programs/ComparativeGenomics/scripts/RootTree.R')


def label(unrooted_tree_file, target_list):
    f = open(unrooted_tree_file, 'r')
    F = f.read()
    f.close()

    for species in target_list:
        F = F.replace(species, species + ' #1')

    g = open('./SpeciesTree/unrooted_labeled.nwk', 'w')
    g.write(F)
    g.close()

