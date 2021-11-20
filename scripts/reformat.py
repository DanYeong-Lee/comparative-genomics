# Modify the format of Single_Copy_Orthologue_Sequences from OrthoFinder output
# File name: Temporary value -> Gene symbol(of Homo sapiens)
# Headers: Protein ID -> Species name

# We should know which species have this protein ID & Gene symbol of Homo sapiens

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import multiprocessing


def get_cds_dict(species_list):
    cds_dict = {}
    for species in species_list:
        per_species_dict = {}
        for record in SeqIO.parse(f"Results/Genome_data/CDS/primary_transcripts/{species}.fa", 
                                 "fasta"):
            per_species_dict[record.id] = record.seq
        cds_dict[species] = per_species_dict
    
    return cds_dict


def get_gene_symbol_dict():
    gene_symbol_dict = {}
    for record in SeqIO.parse("Results/Genome_data/CDS/Homo_sapiens.fna", "fasta"):
        if "protein_id" not in record.description:
            continue
        gene_symbol = record.description.split("gene=")[-1].split("]")[0]
        protein_id = record.description.split("protein_id=")[-1].split("]")[0]
        
        gene_symbol_dict[protein_id] = gene_symbol

    return gene_symbol_dict


def reformat_one_gene(input_file, cds_dict, gene_symbol_dict):
    cds_records = []
    protein_records = []
    
    for record in SeqIO.parse(input_file, "fasta"):
        for species, cdss in cds_dict.items():
            if record.id in cdss.keys():
                the_species = species
                cds = cdss[record.id]
                break
        if the_species == "Homo_sapiens":
            gene_symbol = gene_symbol_dict[record.id]
        
        protein_records.append(SeqRecord(record.seq, id=the_species, description=""))
        cds_records.append(SeqRecord(cds, id=the_species, description=""))
    SeqIO.write(protein_records, f"Results/All_sequences/Protein_sequence/{gene_symbol}.fasta", "fasta")
    SeqIO.write(cds_records, f"Results/All_sequences/CDS_sequence/{gene_symbol}.fasta", "fasta")
    

def main(species_list, file_list):
    cds_dict = get_cds_dict(species_list)
    gene_symbol_dict = get_gene_symbol_dict()
    
    global wrapper
    def wrapper(input_file):
        reformat_one_gene(input_file, cds_dict, gene_symbol_dict)
    
    pool = multiprocessing.Pool(processes=16)
    pool.map(wrapper, file_list)
    pool.close()
    pool.join()
    