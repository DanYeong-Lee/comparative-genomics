import os
import filelist
import multiprocessing


def make_ctl(gene_name, align_tool):
    f = open("ctl/codeml.ctl", 'r')
    F = f.readlines()
    f.close()

    # H0 control file
    for line in F:
        if line.startswith('seqfile ='):
            F[F.index(line)] = 'seqfile = ../Alignment/Codon/Trimmed/{0}_pal2nal_gblocks_paml/{1}.fas\n'.format(align_tool, gene_name)  # directory of sequence file
        elif line.startswith('      outfile = '):
            F[F.index(line)] = '      outfile = ./Results/codeml_results/{0}_H0\n'.format(gene_name)  # output directory
        elif line.startswith('     treefile ='):
            F[F.index(line)] = '     treefile = ../SpeciesTree/unrooted_labeled.nwk\n'  # directory of labeled_unrooted_tree file
        elif line.startswith('        model ='):
            F[F.index(line)] = '        model = 2\n'
        elif line.startswith('      NSsites ='):
            F[F.index(line)] = '      NSsites = 2  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;\n'
        elif line.startswith('    fix_omega ='):
            F[F.index(line)] = '    fix_omega = 1  * 1: omega or omega_1 fixed, 0: estimate \n'
        elif line.startswith('        omega ='):
            F[F.index(line)] = '        omega = 1. * initial or fixed omega, for codons or codon-based AAs\n'

    g = open('Results/PositiveSelection/Control_files/{0}_H0.ctl'.format(gene_name), 'w')
    for line in F:
        g.write(line)

    # H1 control file
    for line in F:
        if line.startswith('      outfile = '):
            F[F.index(line)] = '      outfile = ./Results/codeml_results/{0}_H1\n'.format(gene_name)  # output directory
        elif line.startswith('    fix_omega ='):
            F[F.index(line)] = '    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate \n'

    g = open('Results/PositiveSelection/Control_files/{0}_H1.ctl'.format(gene_name), 'w')
    for line in F:
        g.write(line)


def make_done_list(AllGeneList):
    raw_done_list = filelist.mklist(extension='_H0', directory='Results/PositiveSelection/Results/codeml_results')
    complete_list = []
    incomplete_list = []
    for gene in raw_done_list:
        if gene + '_H1' not in os.listdir('Results/PositiveSelection/Results/codeml_results'):
            incomplete_list.append(gene)
        else:
            if os.path.getsize('Results/PositiveSelection/Results/codeml_results/' + gene + '_H0') >= os.path.getsize('Results/PositiveSelection/Results/codeml_results/' + gene + '_H1'):
                incomplete_list.append(gene)
            else:
                complete_list.append(gene)
    for gene in complete_list:
        AllGeneList.remove(gene)
    for gene in incomplete_list:
        os.system('rm Results/PositiveSelection/Results/codeml_results/' + gene + '_*')

    return AllGeneList


def run_codeml(gene_name):
    os.system('codeml ./Control_files/{0}_H0.ctl'.format(gene_name))
    os.system('codeml ./Control_files/{0}_H1.ctl'.format(gene_name))

