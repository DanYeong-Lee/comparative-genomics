from scripts import filelist
import multiprocessing
import os


def makedict(species_list):
    all_species_CDS_dict = {}
    all_species_TCDS_dict = {}
    for species_name in species_list:
        q_CDS = open('./Genome_data/CDS/%s.fna' % species_name, 'r')
        Q_CDS = q_CDS.read()
        q_CDS.close()
        CDS_species_sequence_list = Q_CDS.split('>lcl')
        CDS_species_sequence_list.remove('')

        cds_dict = {}  # {cds id: sequence}

        for unit in CDS_species_sequence_list:
            cds_id = unit.split('_cds_')[1].split(' [')[0]
            cds_seq_lines = unit.split('\n')[1:]
            while cds_seq_lines.count('') > 1:
                cds_seq_lines.remove('')
            cds_sequence = '\n'.join(cds_seq_lines)
            cds_dict[cds_id] = cds_sequence

        all_species_CDS_dict[species_name] = cds_dict

        q_TCDS = open('./Genome_data/Translated_CDS/%s.faa' % species_name, 'r')
        Q_TCDS = q_TCDS.read()
        q_TCDS.close()
        TCDS_species_sequence_list = Q_TCDS.split('>lcl')
        TCDS_species_sequence_list.remove('')

        tcds_dict = {}

        for unit in TCDS_species_sequence_list:
            tcds_id = unit.split('_prot_')[1].split(' [')[0]
            tcds_seq_lines = unit.split('\n')[1:]
            while tcds_seq_lines.count('') > 1:
                tcds_seq_lines.remove('')
            tcds_sequence = '\n'.join(tcds_seq_lines)
            tcds_dict[tcds_id] = tcds_sequence

        all_species_TCDS_dict[species_name] = tcds_dict

    return all_species_CDS_dict, all_species_TCDS_dict


def match(protein_file_name, all_species_CDS_dict, all_species_TCDS_dict):
    f = open('./Single_Copy_Orthologue_Sequences/Protein_sequence/%s' % protein_file_name, 'r')
    F = f.read()
    f.close()
    gene_symbol = protein_file_name.split('.')[0]
    species_sequence_list = F.split('>')
    species_sequence_list.remove('')

    species_pt_sequence_dict = {}  # {species_name : protein_sequence}

    for i in species_sequence_list:
        species_name = i.split('\n')[0]
        seq_lines = i.split('\n')[1:]
        while seq_lines.count('') > 1:
            seq_lines.remove('')
        sequence = '\n'.join(seq_lines)
        species_pt_sequence_dict[species_name] = sequence

    for species_name in species_pt_sequence_dict:
        for tcds_id in all_species_TCDS_dict[species_name]:
            if all_species_TCDS_dict[species_name][tcds_id] == species_pt_sequence_dict[species_name]:
                ID = tcds_id  # protein_ID
                break

        coding_sequence = all_species_CDS_dict[species_name][ID]

        g = open('./PositiveSelection/Materials/matched_CDS/%s_CDS.faa' % gene_symbol, 'a')
        header = '>' + species_name + '\n'
        g.write(header)
        g.write(coding_sequence)
        g.close()


def wrap(protein_file_name):
    match(protein_file_name, all_species_CDS_dict, all_species_TCDS_dict)


if __name__ == '__main__':
    protein_file_list = filelist.mklist_full(extension='.fasta', directory='./Single_Copy_Orthologue_Sequences/Protein_sequence/')
    species_list = filelist.mklist(extension='.fna', directory='./Genome_data/CDS')
    os.system('mkdir matched_CDS')

    all_species_CDS_dict, all_species_TCDS_dict = makedict(species_list)

    pool = multiprocessing.Pool(processes=15)
    pool.map(wrap, protein_file_list)
    pool.close()
    pool.join()
