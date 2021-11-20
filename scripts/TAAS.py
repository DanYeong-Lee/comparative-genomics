import os
import csv
import argparse
import filelist
import parsefas
import multiprocessing


def TAAS(gene_name, target_list, control_list, alignment_directory, format):
    if format == 'Gblocks':
        file_format = 'fas'
        result_directory = 'Results/TAAS/Gblocks/'
    elif format == 'NotTrimmed':
        file_format = 'best.fas'
        result_directory = 'Results/TAAS/NotTrimmed/'

    file_path = '{0}/{1}.{2}'.format(alignment_directory, gene_name, file_format)  # Edit(Target files extension)
    Species_dic = parsefas.MakeSeqDict(file_path, no_enter=True)

    Index = []

    for site_number in range(len(Species_dic[target_list[0]])):
        target_seq = [Species_dic[species][site_number] for species in target_list]
        control_seq = [Species_dic[species][site_number] for species in control_list]

        n = 0
        for base in target_seq:
            if base in control_seq:
                n += 1

        if n == 0:
            Index.append(site_number)

    if len(Index) != 0:
        K = ['\t'] + target_list + control_list
        f = open(result_directory + 'TAAS_sites/%s.csv' % gene_name, 'w')
        wr = csv.writer(f)
        wr.writerow(K)
        for i in Index:
            Seq = [str(i + 1)]
            for j in target_list:
                Seq.append(Species_dic[j][i])
            for j in control_list:
                Seq.append(Species_dic[j][i])
            wr.writerow(Seq)
        f.close()

        g = open(result_directory + 'TAAS_genes.txt', 'a')
        line = gene_name + '\n'
        g.write(line)
        g.close()


def analyse(target_list, format):
    if format == 'Gblocks':
        result_directory = 'Results/TAAS/Gblocks/'
    elif format == 'NotTrimmed':
        result_directory = 'Results/TAAS/NotTrimmed/'

    t1_gene = []
    t1_site = 0
    t2_gene = []
    t2_site = 0
    t3_gene = []
    t3_site = 0
    t4_gene = []
    t4_site = 0
    TAAS_list = filelist.mklist('.csv', directory=result_directory + 'TAAS_sites')
    for gene in TAAS_list:
        f = open(result_directory + 'TAAS_sites/%s.csv' % gene, 'r')
        F_all = f.read()
        F = F_all.split('\n')[1:-1]
        f.close()
        for line in F:
            all_base_list = line.split(',')[1:]
            target_base_list = all_base_list[0:len(target_list)]
            control_base_list = all_base_list[len(target_list):]
            if target_base_list.count(target_base_list[0]) == len(target_base_list) and control_base_list.count(
                    control_base_list[0]) == len(control_base_list):
                t1_site += 1
                if gene not in t1_gene:
                    t1_gene.append(gene)
            elif target_base_list.count(target_base_list[0]) == len(target_base_list):
                t2_site += 1
                if gene not in t2_gene:
                    t2_gene.append(gene)
            elif control_base_list.count(control_base_list[0]) == len(control_base_list):
                t3_site += 1
                if gene not in t3_gene:
                    t3_gene.append(gene)
            else:
                t4_site += 1
                if gene not in t4_gene:
                    t4_gene.append(gene)

    all_sites = t1_site + t2_site + t3_site + t4_site

    g = open(result_directory + 'TAAS_summary.txt', 'w')
    header = 'Total (TAAS genes: ' + str(len(TAAS_list)) + ', TAAS sites: ' + str(all_sites) + ')\n\n'
    line1 = '>Type 1\n' + 'Genes: ' + str(len(t1_gene)) + '\n' + 'Sites: ' + str(t1_site) + ' (' + str(
        round(t1_site / all_sites, 4) * 100) + '%)\n\n'
    line2 = '>Type 2\n' + 'Genes: ' + str(len(t2_gene)) + '\n' + 'Sites: ' + str(t2_site) + ' (' + str(
        round(t2_site / all_sites, 4) * 100) + '%)\n\n'
    line3 = '>Type 3\n' + 'Genes: ' + str(len(t3_gene)) + '\n' + 'Sites: ' + str(t3_site) + ' (' + str(
        round(t3_site / all_sites, 4) * 100) + '%)\n\n'
    line4 = '>Type 4\n' + 'Genes: ' + str(len(t4_gene)) + '\n' + 'Sites: ' + str(t4_site) + ' (' + str(
        round(t4_site / all_sites, 4) * 100) + '%)\n\n'
    g.write(header)
    g.write(line1)
    g.write(line2)
    g.write(line3)
    g.write(line4)
    g.close()

    h = open(result_directory + 'TAAS_genes_type.txt', 'w')
    h.write('>Type 1\n')
    for i in t1_gene:
        line = i + '\n'
        h.write(line)
    h.write('\n')

    h.write('>Type 2\n')
    for i in t2_gene:
        line = i + '\n'
        h.write(line)
    h.write('\n')

    h.write('>Type 3\n')
    for i in t3_gene:
        line = i + '\n'
        h.write(line)
    h.write('\n')

    h.write('>Type 4\n')
    for i in t4_gene:
        line = i + '\n'
        h.write(line)
    h.write('\n')

    h.close()


def RunTAAS(align_tool, TargetList, ControlList, threads):
    os.makedirs('Results/TAAS/Gblocks/TAAS_sites', exist_ok=True)
    os.makedirs('Results/TAAS/NotTrimmed/TAAS_sites', exist_ok=True)
    g = open('Results/TAAS/Gblocks/TAAS_genes.txt', 'w')
    g.close()
    h = open('Results/TAAS/NotTrimmed/TAAS_genes.txt', 'w')
    h.close()

    Gene = filelist.mklist('.fas', directory='Results/Alignment/Protein/Trimmed/{}_gblocks'.format(align_tool))

    global TAAS_wrap_gblocks, TAAS_wrap_nottrimmed

    def TAAS_wrap_gblocks(gene_name):
        TAAS(gene_name, TargetList, ControlList, 'Results/Alignment/Protein/Trimmed/{}_gblocks'.format(align_tool), 'Gblocks')

    def TAAS_wrap_nottrimmed(gene_name):
        TAAS(gene_name, TargetList, ControlList, 'Results/Alignment/Protein/NotTrimmed/{}'.format(align_tool), 'NotTrimmed')

    pool = multiprocessing.Pool(processes=threads)
    pool.map(TAAS_wrap_gblocks, Gene)
    pool.map(TAAS_wrap_nottrimmed, Gene)
    pool.close()
    pool.join()

    analyse(TargetList, 'Gblocks')
    analyse(TargetList, 'NotTrimmed')

