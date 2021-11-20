import os
import argparse
import multiprocessing
import re
import parsefas
import filelist
from scripts import download_genomes, primary_transcript, reformat, TAAS, concatenate, CodonAlign, codeml, parsePS, \
    tree, align


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
    os.makedirs('./Program_results', exist_ok=True)

    os.system(
        'orthofinder -t {0} -a {0} -og -o ./Program_results/OrthoFinder -f ./Genome_data/Protein/primary_transcripts'.format(
            args.threads))  # TODO

def Reformat():
    os.makedirs('./All_sequences/Protein_sequence', exist_ok=True)
    os.makedirs('./All_sequences/CDS_sequence', exist_ok=True)
    OF_dir = os.listdir('./Program_results/OrthoFinder')[0]
    sin_ort_list = filelist.mklist_full('.fa', f"./Program_results/OrthoFinder/{OF_dir}/Single_Copy_Orthologue_Sequences", absolute=True)
    print(sin_ort_list[0])
    reformat.main(species_list, sin_ort_list)


def Align():
    gene_list = filelist.mklist_full('.fasta', './All_sequences/Protein_sequence', absolute=True)
    protein_aligner = align.ProteinAligner(args.alignment)
    pool = multiprocessing.Pool(processes=args.threads)
    pool.map(protein_aligner.align, gene_list)
    pool.close()
    pool.join()

    if args.alignment == 'prank':
        os.system('rm -r tmpdirprankmsa*')

    trimmed_list = filelist.mklist_full('fas', f"./Alignment/Protein/Trimmed/{args.alignment.upper()}_gblocks", absolute=True)
    blank_gene_list = []
    for gene in trimmed_list:
        with open(gene) as f:
            f.readline()
            if f.readline() == '\n':
                blank_gene_list.append(gene.split('/')[-1].split('.')[0])

    with open('./Alignment/Protein/Trimmed/Blank_genes.txt', 'w') as g:
        for gene in blank_gene_list:
            os.system(f"rm ./Alignment/Protein/Trimmed/{args.alignment.upper()}_gblocks/{gene}.fas")
            line = gene + '\n'
            g.write(line)

    
def codon_align():
    CodonAlign.codonalign(args.alignment, args.threads)


def TAAS_analysis():
    TAAS.RunTAAS(args.alignment, target_list, control_list, args.threads)


def FindModel():
    os.makedirs('./Alignment/Protein/Concatenated', exist_ok=True)
    gene_list = filelist.mklist_full(extension='.fas',
                                     directory='./Alignment/Protein/Trimmed/{}_gblocks'.format(args.alignment))
    species_sequence_dict = concatenate.make_empty_species_sequence_dict(gene_list,
                                                                         './Alignment/Protein/Trimmed/{}_gblocks'.format(
                                                                             args.alignment))
    made_species_sequence_dict = concatenate.concatenate(species_sequence_dict, gene_list,
                                                         './Alignment/Protein/Trimmed/{}_gblocks'.format(
                                                             args.alignment))
    concatenate.write(made_species_sequence_dict)

    os.makedirs('./EvolutionaryModel', exist_ok=True)
    if args.modeltest == 'prottest':
        os.makedirs('./Program_results/ProtTest3', exist_ok=True)
        os.system(
            'java -jar /home/ldy9381/programs/prottest-3.4.2/prottest-3.4.2.jar -i ./Alignment/Protein/Concatenated/concatenated.fas -o ./Program_results/ProtTest3/protest.txt -all-distributions -BIC -tc 0.5 -threads %d' % args.threads)
        f = open('./Program_results/ProtTest3/prottest.txt', 'r')
        F = f.readlines()
        f.close()
        for line in F:
            if line.startswith('Best model according to BIC:'):
                BestModel = line.split('Best model according to BIC:')[1].strip()

    elif args.modeltest == 'modelfinder':
        os.makedirs('./Program_results/ModelFinder', exist_ok=True)
        os.system(
            'iqtree -s ./Alignment/Protein/Concatenated/concatenated.fas -o {0} -st AA -pre ./Program_results/ModelFinder/modelfinder -T {1} -m TESTONLY'.format(
                outgroup, args.threads))
        os.system('gzip -d ./Program_results/ModelFinder/*.gz')
        f = open('./Program_results/ModelFinder/modelfinder.iqtree', 'r')
        F = f.readlines()
        f.close()
        for line in F:
            if line.startswith('Best-fit model according to BIC:'):
                BestModel = line.split('Best-fit model according to BIC:')[1].strip()

    elif args.modeltest == 'modeltest-ng':
        os.makedirs('./Program_results/ModelTest-NG', exist_ok=True)
        os.system(
            'modeltest-ng -i ./alignment_concatenated/concatenated.fas -d aa -o ./Program_results/ModelTest-NG/modeltest -t ml -T raxml -p {}'.format(
                args.threads))
        f = open('./Program_results/ModelTest-NG/modeltest.out', 'r')
        F = f.readlines()
        f.close()

        for line in F:
            if line.startswith('Best model according to BIC'):
                BestModel = F[F.index(line) + 2].split('Model:')[1].strip()

    h = open('./EvolutionaryModel/Model.txt', 'w')
    line = 'Best model according to BIC: ' + BestModel
    h.write(line)
    h.close()


def BuildTree():
    os.makedirs('./SpeciesTree', exist_ok=True)

    if args.tree == 'iqtree':
        f = open('./EvolutionaryModel/Model.txt', 'r')
        F = f.readline()
        f.close()
        BestModel = F.split('Best model according to BIC:')[1].strip()
        os.makedirs('./Program_results/IQ-TREE', exist_ok=True)
        os.system(
            'iqtree -s ./Alignment/Protein/Concatenated/concatenated.fas -st AA -pre ./Program_results/IQ-TREE/iqtree -o {0} -T {1} -m {2} -B 1000'.format(
                outgroup, args.threads, BestModel))
        # TODO

    elif args.tree == 'raxml':
        f = open('./EvolutionaryModel/Model.txt', 'r')
        F = f.readline()
        f.close()
        BestModel = F.split('Best model according to BIC:')[1].strip()
        os.makedirs('./Program_results/RAxML-NG', exist_ok=True)
        os.system(
            'raxml-ng --all --msa ./Alignment/Protein/Concatenated/concatenated.fas --msa-format FASTA --data-type AA --prefix ./Program_results/RAxML-NG/raxml --outgroup {0} --bs-trees autoMRE --model {1} --threads {2}'.format(
                outgroup, BestModel, args.threads))
        # TODO

    elif args.tree == 'timetree':
        tree.delete_branch_length('./input/timetree.nwk')
        tree.root('/home/ldy9381/programs/ComparativeGenomics/scripts/RootTree_source.R', outgroup)
        tree.label('./SpeciesTree/unrooted.nwk', target_list)
        f = open('./input/timetree.nwk')
        F = f.read()
        f.close()

        ultra = re.sub("'[0-9]+'", '', F)

        g = open('./SpeciesTree/Ultrametric.nwk', 'w')
        g.write(ultra)
        g.close()


def PositiveSelection():
    os.makedirs('./PositiveSelection/Results/codeml_results', exist_ok=True)
    os.makedirs('./PositiveSelection/Control_files', exist_ok=True)
    os.makedirs('./PositiveSelection/Results/summary', exist_ok=True)

    gene_list = filelist.mklist(extension='.fas',
                                directory='./Alignment/Codon/Trimmed/{}_pal2nal_gblocks_paml'.format(args.alignment))
    if 'Genes_with_PSC.txt' not in os.listdir('./PositiveSelection'):
        psc_gene_dict = {}
        sc_list = ['TGA', 'TAA', 'TAG']

        for file in gene_list:
            file_path = './Alignment/Codon/Trimmed/{0}_pal2nal_gblocks/{1}.fas'.format(args.alignment, file)  # TODO
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

        h = open('./PositiveSelection/Genes_with_PSC.txt', 'w')

        for psc_gene in psc_gene_dict:
            gene_list.remove(psc_gene)
            line = psc_gene + psc_gene_dict[psc_gene]
            h.write(line)
            h.write('\n')

        h.close()
    else:
        f = open('./PositiveSelection/Genes_with_PSC.txt')
        F = f.readlines()
        for psc_gene in F:
            gene_list.remove(psc_gene.split(':')[0])

    if not os.listdir('./PositiveSelection/Control_files'):
        for i in gene_list:
            codeml.make_ctl(i, args.alignment)

    new_gene_list = codeml.make_done_list(gene_list)

    os.chdir('./PositiveSelection')

    pool2 = multiprocessing.Pool(processes=args.threads)
    pool2.map(codeml.run_codeml, new_gene_list)
    pool2.close()
    pool2.join()

    os.chdir('../')

    parsePS.parsePS_main()


def UltrametricTree():
    os.makedirs('./Program_results/MCMCTree/codeml', exist_ok=True)
    os.makedirs('./Program_results/MCMCTree/result1', exist_ok=True)
    mcmctree.strict_codeml()
    mcmctree.prepare()
    os.chdir('./Program_results/MCMCTree/result1')
    mcmctree.mcmc()
    os.chdir('../result2')
    mcmctree.mcmc()


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
    elif args.work == 'model':
        FindModel()
    elif args.work == 'tree':
        BuildTree()
    elif args.work == 'PS':
        PositiveSelection()
    elif args.work == 'ultra':
        UltrametricTree()
    elif args.work == 'cafe':
        CAFE()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--species', help='file name of species list', default='input/species_list.txt')
    parser.add_argument('-w', '--work', help='work you wanna run',
                        choices=['all', 'download', 'ortholog', 'reformat', 'align', 'codonalign', 'taas', 'model', 'tree', 'PS', 'ultra',
                                 'cafe'])
    parser.add_argument('-a', '--alignment', help='alignment tool', choices=['prank', 'mafft', 'muscle'],
                        default='prank')
    parser.add_argument('-t', '--threads', help='number of threads', type=int, default=4)
    parser.add_argument('-T', '--tree', help='tree building program', choices=['iqtree', 'raxml', 'timetree'],
                        default='timetree')
    parser.add_argument('-m', '--modeltest', help='modeltest tool', choices=['modelfinder', 'prottest', 'modeltest-ng'],
                        default='modelfinder')
    parser.add_argument('-d', '--database', help='Genome database[NCBI or Ensembl]', choices=['ncbi', 'ensembl'], default='ncbi')
    args = parser.parse_args()

    main()
