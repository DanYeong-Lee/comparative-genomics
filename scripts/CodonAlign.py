from scripts import filelist
import os
import multiprocessing
import argparse


def pal2nal_normal(gene_name, align_tool):
    os.system(
        'pal2nal.pl Results/Alignment/Protein/NotTrimmed/{1}/{0}.best.fas Results/All_sequences/CDS_sequence/{0}.fasta -output fasta > Results/Alignment/Codon/NotTrimmed/{1}_pal2nal/{0}.fas'.format(
            gene_name, align_tool))


def pal2nal_error(gene_name, align_tool):
    os.system(
        'pal2nal.pl Results/Alignment/Protein/NotTrimmed/{1}/{0}.best.fas Results/Alignment/Codon/NotTrimmed/{1}_pal2nal_errorfix/{0}.fasta -output fasta > Results/Alignment/Codon/NotTrimmed/{1}_pal2nal/{0}.fas'.format(
            gene_name, align_tool))


def errorlist(align_tool):
    pal2nal_list = filelist.mklist(extension='.fas',
                                   directory='Results/Alignment/Codon/NotTrimmed/{}_pal2nal'.format(align_tool))
    error_list = []
    for file in pal2nal_list:
        if os.path.getsize('Results/Alignment/Codon/NotTrimmed/{0}_pal2nal/{1}.fas'.format(align_tool, file)) == 0:
            error_list.append(file)
    return error_list


def errordict(errorlist):
    all_dict = {}
    for gene in errorlist:
        nuc_dict = {}

        f = open('Results/All_sequences/CDS_sequence/{0}.fasta'.format(gene))
        F = f.read()
        f.close()

        unit_list = F.split('>')
        unit_list.remove('')

        for unit in unit_list:
            species = unit.split('\n')[0]
            nucseq = '\n'.join(unit.split('\n')[1:])
            nuc_dict[species] = nucseq

        all_dict[gene] = nuc_dict

    return all_dict


def fixinputnuc(errordict, align_tool):
    for gene in errordict:
        for species in errordict[gene]:
            nuc_len = len(errordict[gene][species].replace('\n', ''))
            if nuc_len % 3 == 2:
                errordict[gene][species] = errordict[gene][species].strip() + 'N\n'
            elif nuc_len % 3 == 1:
                errordict[gene][species] = errordict[gene][species] + 'NN\n'

        f = open('Results/Alignment/Codon/NotTrimmed/{0}_pal2nal_errorfix/{1}.fasta'.format(align_tool, gene), 'w')
        for species in errordict[gene]:
            unit = '>' + species + '\n' + errordict[gene][species]
            f.write(unit)
        f.close()

        pal2nal_error(gene, align_tool)


def convert(gene, align_tool):
    f = open('Results/Alignment/Codon/Trimmed/{0}_pal2nal_gblocks/{1}.fas'.format(align_tool, gene))
    F = f.read()
    f.close()

    g = open('Results/Alignment/Codon/Trimmed/{0}_pal2nal_gblocks_paml/{1}.fas'.format(align_tool, gene), 'w')
    seq_list = F.split('>')
    seq_list.remove('')

    seq_dict = {}
    for unit in seq_list:
        species = unit.split('\n')[0]
        sequence = '\n'.join(unit.split('\n')[1:]).replace(' ', '')
        seq_dict[species] = sequence

    species_num = len(seq_dict)
    seq_len = len(seq_dict[list(seq_dict.keys())[0]].replace('\n', ''))
    g.write('  {0}   {1}\n'.format(species_num, seq_len))
    for i in seq_dict:
        unit = i + '\n' + seq_dict[i]
        g.write(unit)
    g.close()


def codonalign(align_tool, threads):
    os.makedirs('Results/Alignment/Codon/NotTrimmed/{}_pal2nal'.format(align_tool), exist_ok=True)
    os.makedirs('Results/Alignment/Codon/NotTrimmed/{}_pal2nal_errorfix'.format(align_tool), exist_ok=True)
    os.makedirs('Results/Alignment/Codon/Trimmed/{}_pal2nal_gblocks'.format(align_tool), exist_ok=True)
    os.makedirs('Results/Alignment/Codon/Trimmed/{}_pal2nal_gblocks_paml'.format(align_tool), exist_ok=True)

    gene_list = filelist.mklist(extension='.best.fas', directory='Results/Alignment/Protein/NotTrimmed/{}'.format(align_tool))

    global wrap1, wrap2

    def wrap1(gene_name):
        pal2nal_normal(gene_name, align_tool)

    pool1 = multiprocessing.Pool(processes=threads)
    pool1.map(wrap1, gene_list)
    pool1.close()
    pool1.join()

    error_list = errorlist(align_tool)
    error_dict = errordict(error_list)
    fixinputnuc(error_dict, align_tool)

    for gene in gene_list:
        os.system('Gblocks Results/Alignment/Codon/NotTrimmed/{1}_pal2nal/{0}.fas -t=c -p=n'.format(gene, align_tool))
        os.system(
            'cp Results/Alignment/Codon/NotTrimmed/{1}_pal2nal/{0}.fas-gb Results/Alignment/Codon/Trimmed/{1}_pal2nal_gblocks/{0}.fas'.format(
                gene, align_tool))
        os.system('rm Results/Alignment/Codon/NotTrimmed/{1}_pal2nal/{0}.fas-gb'.format(gene, align_tool))

    gblocked_list = filelist.mklist(extension='.fas', directory='Results/Alignment/Codon/Trimmed/{}_pal2nal_gblocks'.format(align_tool))

    blank_gene_list = []
    for gene in gblocked_list:
        with open('Results/Alignment/Codon/Trimmed/{1}_pal2nal_gblocks/{0}.fas'.format(gene, align_tool), 'r') as f:
            f.readline()
            if f.readline() == '\n':
                blank_gene_list.append(gene)

    with open('Results/Alignment/Codon/Trimmed/Blank_genes.txt', 'w') as g:
        for gene in blank_gene_list:
            os.system('rm Results/Alignment/Codon/Trimmed/{1}_pal2nal_gblocks/{0}.fas'.format(gene, align_tool))
            line = gene + '\n'
            g.write(line)

    new_gblocked_list = filelist.mklist(extension='.fas', directory='Results/Alignment/Codon/Trimmed/{}_pal2nal_gblocks'.format(align_tool))
    for gene in new_gblocked_list:
        convert(gene, align_tool)
