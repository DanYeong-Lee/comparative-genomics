import sys
import os
sys.path.append(os.path.abspath("./"))

sys.path.append("/data/project/")
from scripts.__main__ import main
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--species', help='file name of species list', default='input/species_list.txt')
    parser.add_argument('-w', '--work', help='work you wanna run',
                        choices=['all', 'download', 'ortholog', 'reformat', 'align', 'codonalign', 'taas', 'model', 'tree', 'PS', 'ultra','cafe'])
    parser.add_argument('-a', '--alignment', help='alignment tool', choices=['PRANK', 'MAFFT', 'MUSCLE'],
                        default='PRANK')
    parser.add_argument('-t', '--threads', help='number of threads', type=int, default=4)
    parser.add_argument('-T', '--tree', help='tree building program', choices=['iqtree', 'raxml', 'timetree'],
                        default='timetree')
    parser.add_argument('-m', '--modeltest', help='modeltest tool', choices=['modelfinder', 'prottest', 'modeltest-ng'],
                        default='modelfinder')
    args = parser.parse_args()

    main(args)