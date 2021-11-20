import pandas as pd
import os
import multiprocessing
import argparse
import itertools

def get_species_list(species_file):
    with open(species_file) as f:
        species_list = f.readlines()
        species_list = [species.strip().replace("*", "").replace("#", "") for species in species_list]
        
    return species_list

def get_ftp_urls(summary_file, species_list):
    df = pd.read_csv(summary_file, sep="\t", skiprows=2, header=None)
    df.index = df[7]
    df = df[19]
    
    return list(df[species_list])

def download(species, ftp_url):
    assembly_name = ftp_url.split("/")[-1]
    ftp_url += "/"
    ftp_url += assembly_name
    
    os.system(f"wget {ftp_url + '_cds_from_genomic.fna.gz'} "
              f"-O Genome_data/CDS/{species}.fna.gz")
    os.system(f"wget {ftp_url + '_protein.faa.gz'} "
              f"-O Genome_data/Protein/{species}.faa.gz")
    
    #os.system(f"wget {ftp_url + '_genomic.gtf.gz'} "
    #         f"-O {os.path.join(genome_path, 'GTF.gtf.gz')}")

    
def main():
    os.makedirs("Genome_data/CDS/primary_transcripts", exist_ok=True)
    os.makedirs("Genome_data/Protein/primary_transcripts", exist_ok=True)
    
    assembly_summary = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"
    os.system(f"wget {assembly_summary} -O Genome_data/assembly_summary_refseq.txt")
    
    species_list = get_species_list("input/species_list.txt")
    urls_list = get_ftp_urls("Genome_data/assembly_summary_refseq.txt", species_list)

    
    pool = multiprocessing.Pool(processes=len(species_list))
    pool.starmap(download, zip([x.replace(" ", "_") for x in species_list], urls_list))
    pool.close()
    pool.join()
    
    os.system('gzip -d Genome_data/CDS/*.gz')
    os.system('gzip -d Genome_data/Protein/*.gz')
    
    
if __name__ == "__main__":
    main()