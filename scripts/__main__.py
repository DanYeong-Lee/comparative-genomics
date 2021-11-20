import os
import argparse
import multiprocessing
import re
import parsefas
import filelist
from scripts import download_genomes, primary_transcript, reformat, TAAS, CodonAlign, codeml, parsePS, tree, align


def MakeList(species_file):
    with open(species_file) as f:
        raw_species_list = f.readlines()
        raw_species_list = [x.strip().replace(" ", "_") for x in raw_species_list]
    species_list = [x.replace("*", "").replace("#", "") for x in raw_species_list]
    
    target_list = [x.replace("#", "") for x in raw_species_list if x.endswith("#")]
    control_list = [x.replace("*", "") for x in raw_species_list if not x.endswith("#")]
    outgroup = [x.replace("*", "") for x in raw_species_list if x.endswith("*")][0]

    return species_list, target_list, control_list, outgroup

def DownloadGenome():
    download_genomes.main()
    
    pool = multiprocessing.Pool(processes=len(species_list))
    pool.map(primary_transcript.primary_transcript, species_list)
    pool.close()
    pool.join()

def FindOrthologs():
    os.makedirs('Results/Program_results', exist_ok=True)

    os.system(
        'orthofinder -t {0} -a {0} -og -o Results/Program_results/OrthoFinder -f Results/Genome_data/Protein/primary_transcripts'.format(
            args.threads))  # TODO

def Reformat():
    os.makedirs('Results/All_sequences/Protein_sequence', exist_ok=True)
    os.makedirs('Results/All_sequences/CDS_sequence', exist_ok=True)
    OF_dir = os.listdir('Results/Program_results/OrthoFinder')[0]
    sin_ort_list = filelist.mklist_full('.fa', f"Results/Program_results/OrthoFinder/{OF_dir}/Single_Copy_Orthologue_Sequences", absolute=True)
    print(sin_ort_list[0])
    reformat.main(species_list, sin_ort_list)


def Align():
    gene_list = filelist.mklist_full('.fasta', 'Results/All_sequences/Protein_sequence', absolute=True)
    protein_aligner = align.ProteinAligner(args.alignment)
    pool = multiprocessing.Pool(processes=args.threads)
    pool.map(protein_aligner.align, gene_list)
    pool.close()
    pool.join()

    if args.alignment == 'prank':
        os.system('rm -r tmpdirprankmsa*')

    trimmed_list = filelist.mklist_full('fas', f"Results/Alignment/Protein/Trimmed/{args.alignment.upper()}_gblocks", absolute=True)
    blank_gene_list = []
    for gene in trimmed_list:
        with open(gene) as f:
            f.readline()
            if f.readline() == '\n':
                blank_gene_list.append(gene.split('/')[-1].split('.')[0])

    with open('Results/Alignment/Protein/Trimmed/Blank_genes.txt', 'w') as g:
        for gene in blank_gene_list:
            os.system(f"rm Results/Alignment/Protein/Trimmed/{args.alignment.upper()}_gblocks/{gene}.fas")
            line = gene + '\n'
            g.write(line)

    
def codon_align():
    CodonAlign.codonalign(args.alignment, args.threads)


def TAAS_analysis():
    TAAS.RunTAAS(args.alignment, target_list, control_list, args.threads)


def BuildTree():
    os.makedirs('Results/SpeciesTree', exist_ok=True)

    if args.tree == 'timetree':
        tree.delete_branch_length('input/timetree.nwk')
        tree.root('scripts/RootTree_source.R', outgroup)
        tree.label('Results/SpeciesTree/unrooted.nwk', target_list)
        f = open('input/timetree.nwk')
        F = f.read()
        f.close()

        ultra = re.sub("'[0-9]+'", '', F)

        g = open('Results/SpeciesTree/Ultrametric.nwk', 'w')
        g.write(ultra)
        g.close()


def PositiveSelection():
    os.makedirs('Results/PositiveSelection/Results/codeml_results', exist_ok=True)
    os.makedirs('Results/PositiveSelection/Control_files', exist_ok=True)
    os.makedirs('Results/PositiveSelection/Results/summary', exist_ok=True)

    gene_list = filelist.mklist(extension='.fas',
                                directory='Results/Alignment/Codon/Trimmed/{}_pal2nal_gblocks_paml'.format(args.alignment))
    if 'Genes_with_PSC.txt' not in os.listdir('Results/PositiveSelection'):
        psc_gene_dict = {}
        sc_list = ['TGA', 'TAA', 'TAG']

        for file in gene_list:
            file_path = 'Results/Alignment/Codon/Trimmed/{0}_pal2nal_gblocks/{1}.fas'.format(args.alignment, file)  # TODO
            seq_dict = parsefas.MakeSeqDict(file_path, no_enter=True)

            for species in seq_dict:
                seq = seq_dict[species].replace(' ', '')
                for i in range(0, int(len(seq) / 3)):
                    codon = seq[3 * i:3 * (i + 1)]
                    if codon in sc_list:
                        psc_gene_dict[file] = ': ' + species + ', ' + str(3 * i + 1) + ' ~ ' + str(
                            3 * (i + 1)) + ', ' + codon
                        break
                if file in psc_gene_dict:
                    break

        h = open('Results/PositiveSelection/Genes_with_PSC.txt', 'w')

        for psc_gene in psc_gene_dict:
            gene_list.remove(psc_gene)
            line = psc_gene + psc_gene_dict[psc_gene]
            h.write(line)
            h.write('\n')

        h.close()
    else:
        f = open('Results/PositiveSelection/Genes_with_PSC.txt')
        F = f.readlines()
        for psc_gene in F:
            gene_list.remove(psc_gene.split(':')[0])

    if not os.listdir('Results/PositiveSelection/Control_files'):
        for i in gene_list:
            codeml.make_ctl(i, args.alignment)

    new_gene_list = codeml.make_done_list(gene_list)

    os.chdir('Results/PositiveSelection')

    pool2 = multiprocessing.Pool(processes=args.threads)
    pool2.map(codeml.run_codeml, new_gene_list)
    pool2.close()
    pool2.join()

    os.chdir('../../')

    parsePS.parsePS_main()

    
def CAFE():
    cafe.cafe_main(args.threads, target_list)
    parseCAFE.parseCAFE_main(target_list)


def main(MainArgs):
    global args
    args = MainArgs
    global species_list, target_list, control_list, outgroup
    species_list, target_list, control_list, outgroup = MakeList(args.species)

    if args.work == 'download':
        DownloadGenome()
    elif args.work == 'ortholog':
        FindOrthologs()
    elif args.work == 'reformat':
        Reformat()
    elif args.work == 'align':
        Align()
    elif args.work == 'codonalign':
        codon_align()
    elif args.work == 'taas':
        TAAS_analysis()
    elif args.work == 'tree':
        BuildTree()
    elif args.work == 'PS':
        PositiveSelection()
    elif args.work == 'cafe':
        CAFE()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--species', help='file name of species list', default='input/species_list.txt')
    parser.add_argument('-w', '--work', help='work you wanna run',
                        choices=['all', 'download', 'ortholog', 'reformat', 'align', 'codonalign', 'taas', 'tree', 'PS', 'cafe'])
    parser.add_argument('-a', '--alignment', help='alignment tool', choices=['prank', 'mafft', 'muscle'],
                        default='prank')
    parser.add_argument('-t', '--threads', help='number of threads', type=int, default=4)
    parser.add_argument('-T', '--tree', help='tree building program', choices=['iqtree', 'raxml', 'timetree'],
                        default='timetree')
    parser.add_argument('-d', '--database', help='Genome database[NCBI or Ensembl]', choices=['ncbi', 'ensembl'], default='ncbi')
    args = parser.parse_args()

    main()
