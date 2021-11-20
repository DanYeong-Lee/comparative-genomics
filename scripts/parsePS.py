import os
import filelist
import csv
from scipy import stats
from Bio.Seq import Seq


def make_csv():
    header_site = ['\t', 'lnLH0', 'lnLH1', 'LRT', 'p_value', 'foreground w', '\t', 'PS site', 'Amino acid', 'P']
    header = ['\t', 'lnLH0', 'lnLH1', 'LRT', 'p_value', 'foreground w']
    h = open('Results/PositiveSelection/Results/summary/positive_selection_result_site.csv', 'w')
    wr = csv.writer(h)
    wr.writerow(header_site)
    h.close()
    j = open('Results/PositiveSelection/Results/summary/positive_selection_result.csv', 'w')
    wr1 = csv.writer(j)
    wr1.writerow(header)
    j.close()


def parse(gene_name):
    result_directory = 'Results/PositiveSelection/Results/codeml_results/'
    H0_file = gene_name + '_H0'
    H1_file = gene_name + '_H1'
    if H0_file in os.listdir(result_directory) and H1_file in os.listdir(result_directory) and os.path.getsize(
            result_directory + H1_file) > os.path.getsize(result_directory + H0_file):
        f = open(result_directory + H0_file, 'r')
        F = f.readlines()
        f.close()

        for line in F:
            if line.startswith('lnL(ntime:'):
                H0_likelihood = float(line.split('):')[1].strip().split(' ')[0])
                H0_np = int(line.split('np: ')[1].split(')')[0])

        g = open(result_directory + H1_file, 'r')
        G = g.readlines()
        g.close()

        for line in G:
            if line.startswith('lnL(ntime:'):
                H1_likelihood = float(line.split('):')[1].strip().split(' ')[0])
                H1_np = int(line.split('np: ')[1].split(')')[0])

            elif line.startswith('foreground w'):
                w_candidates = line.split('foreground w')[1].strip()
                w_candidate_real = w_candidates.replace('  ', ' ')
                w = w_candidate_real.split(' ')[3]

            elif line.startswith('Bayes Empirical Bayes (BEB) analysis'):
                BEB_index = G.index(line)

            elif line.startswith('The grid'):
                BEB_next_index = G.index(line)

        try:
            site_list = G[BEB_index + 2:BEB_next_index]
            while '\n' in site_list:
                site_list.remove('\n')
            for i in site_list:
                site_list[site_list.index(i)] = i.strip()
            LRT = 2 * (H1_likelihood - H0_likelihood)
            df = H1_np - H0_np

            p_value = (1 - stats.chi2.cdf(LRT, df)) / 2

            h = open('Results/PositiveSelection/Results/summary/positive_selection_result_site.csv', 'a')
            wr = csv.writer(h)

            for i in range(len(site_list)):
                if i == 0:
                    line = [gene_name, H0_likelihood, H1_likelihood, LRT, p_value, w, '\t'] + site_list[i].split(' ')
                    wr.writerow(line)
                else:
                    line = ['\t', '\t', '\t', '\t', '\t', '\t', '\t'] + site_list[i].split(' ')
                    wr.writerow(line)
            wr.writerow('\t')
            h.close()

            j = open('Results/PositiveSelection/Results/summary/positive_selection_result.csv', 'a')
            wr1 = csv.writer(j)
            line1 = [gene_name, H0_likelihood, H1_likelihood, LRT, p_value, w]
            wr1.writerow(line1)
            j.close()

        except:
            pass


def choose():
    os.system('Rscript scripts/FDR.R')
    f = open('Results/PositiveSelection/Results/summary/positive_selection_adjusted.csv', 'r')
    header = f.readline()
    F = f.readlines()
    f.close()

    PS_lines = []

    for i in F:
        if float(i.split(',')[7]) < 0.1 and float(i.split(',')[6]) > 1.0:
            PS_lines.append(i)

    g = open('Results/PositiveSelection/Results/summary/positively_selected_genes.csv', 'w')
    g.write(header)
    for i in PS_lines:
        g.write(i)
    g.close()

    PS_gene_list = []
    for line in PS_lines:
        gene_name = line.split(',')[1].strip('"')
        PS_gene_list.append(gene_name)

    return PS_gene_list


def TAAS_PS(result_directory, psList):
    f = open(result_directory + 'TAAS_genes.txt', 'r')
    F = f.read()
    f.close()

    TAAS_gene_list = F.split('\n')
    TAAS_gene_list.remove('')

    TAAS_PS_list = []
    for gene in TAAS_gene_list:
        if gene in psList:
            TAAS_PS_list.append(gene)

    h = open(result_directory + 'TAAS&PS_genes.txt', 'w')
    header = '>Genes that selected positively & have TAAS sites ' + '(Total: ' + str(len(TAAS_PS_list)) + ')\n\n'
    h.write(header)
    for gene in TAAS_PS_list:
        line = gene + '\n'
        h.write(line)
    h.close()

    return TAAS_PS_list


def overlap(TAAS_PS_list, trim):
    site_overlap_dict = {}
    if trim == 'Gblocks':
        TAAS_directory = 'Results/TAAS/Gblocks/'
    elif trim == 'NotTrimmed':
        TAAS_directory = 'Results/TAAS/NotTrimmed/'

    for gene in TAAS_PS_list:
        g = open('Results/PositiveSelection/Results/codeml_results/{}_H1'.format(gene), 'r')
        G = g.readlines()
        g.close()

        h = open('Results/Alignment/Codon/Trimmed/PRANK_pal2nal_gblocks/{}.fas'.format(gene), 'r')
        H = h.read()
        h.close()

        codon_seq_list = H.split('>')
        codon_seq_list.remove('')
        species_seq_dict = {}
        for unit in codon_seq_list:
            species_name = unit.split('\n')[0]
            cod_seq = Seq(''.join(unit.split('\n')[1:]).replace(' ', ''))
            pep_seq = cod_seq.translate()
            species_seq_dict[species_name] = pep_seq

        for line in G:
            if line.startswith('Bayes Empirical Bayes (BEB) analysis'):
                BEB_index = G.index(line)
            elif line.startswith('The grid'):
                BEB_next_index = G.index(line)

        site_list = G[BEB_index + 2: BEB_next_index]
        while '\n' in site_list:
            site_list.remove('\n')
        for i in site_list:
            site_list[site_list.index(i)] = i.strip()

        ps_site_list = []
        for i in site_list:
            ps_site_list.append(i)

        species_list = list(species_seq_dict.keys())

        f = open(TAAS_directory + 'TAAS_sites/{}.csv'.format(gene), 'r')
        species = f.readline().strip().split(',')[1:]
        F = f.readlines()
        f.close()

        TAAS_site_dict = {}
        for k in F:
            line = k.strip()
            aa_species_dict = {}
            site_number = line.split(',')[0]
            aa_list = line.split(',')[1:]
            for i in species:
                aa_species_dict[i] = aa_list[species.index(i)]
            TAAS_site_dict[site_number] = aa_species_dict

        overlap_ps_site_list = []
        for i in ps_site_list:
            sn = i.split(' ')[0]
            rep_aa = i.split(' ')[1]
            p_value = i.split(' ')[2]

            for TAAS_site_number in TAAS_site_dict:
                w = 0
                for species in species_list:
                    if TAAS_site_dict[TAAS_site_number][species] == species_seq_dict[species][int(sn)-1]:
                        w += 1
                if w == len(species_list) and rep_aa == TAAS_site_dict[TAAS_site_number][species_list[0]] and i not in overlap_ps_site_list:
                    overlap_ps_site_list.append(i)
        if overlap_ps_site_list:
            site_overlap_dict[gene] = overlap_ps_site_list

    return site_overlap_dict


def write(site_overlap_dict, trim):
    if trim == 'Gblocks':
        TAAS_directory = 'Results/TAAS/Gblocks/'
    elif trim == 'NotTrimmed':
        TAAS_directory = 'Results/TAAS/NotTrimmed/'

    f = open(TAAS_directory + 'TAAS&PS_overlapping.txt', 'w')
    header = 'Genes that have positively selected TAAS sites ' + '(Total: ' + str(len(site_overlap_dict)) + ')\n\n'
    f.write(header)
    for gene in site_overlap_dict:
        header = '>' + gene + '' + '\n'
        f.write(header)
        for site in site_overlap_dict[gene]:
            line = site + '\n'
            f.write(line)
        f.write('\n')
    f.close()


def parsePS_main():
    os.makedirs('Results/PositiveSelection/Results/summary', exist_ok=True)
    gene_list = filelist.mklist(directory='Results/PositiveSelection/Results/codeml_results', extension='_H0')
    make_csv()
    for i in gene_list:
        parse(i)
    PS_list = choose()
    TAAS_PS_Gblocks_list = TAAS_PS('Results/TAAS/Gblocks/', PS_list)
    TAAS_PS_NotTrimmed_list = TAAS_PS('Results/TAAS/NotTrimmed/', PS_list)
    GOOD_dict_gblocks = overlap(TAAS_PS_Gblocks_list, 'Gblocks')
    Good_dict_nottrimmed = overlap(TAAS_PS_NotTrimmed_list, 'NotTrimmed')
    write(GOOD_dict_gblocks, 'Gblocks')
    write(Good_dict_nottrimmed, 'NotTrimmed')


if __name__ == '__main__':
    parsePS_main()
