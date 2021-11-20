from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def primary_transcript(species):
    cds_primary_dict = {}

    for record in SeqIO.parse(f"Results/Genome_data/CDS/{species}.fna", "fasta"):
        if "protein_id" not in record.description:
            continue
        gene_name = record.description.split("gene=")[-1].split("]")[0]
        protein_id = record.description.split("protein_id=")[-1].split("]")[0]
        value = [protein_id, record.seq]

        try:
            if len(cds_primary_dict[gene_name][-1]) < len(record.seq):
                cds_primary_dict[gene_name] = value
        except:
            cds_primary_dict[gene_name] = value

    protein_dict = SeqIO.to_dict(SeqIO.parse(f"Results/Genome_data/Protein/{species}.faa", "fasta"))

    cds_records = []
    protein_records = []

    for gene_name, values in cds_primary_dict.items():
        cds_records.append(SeqRecord(values[1], id=values[0], description=""))
        protein_record = protein_dict[values[0]]
        protein_record.description = ""
        protein_records.append(protein_record)

    SeqIO.write(cds_records, f"Results/Genome_data/CDS/primary_transcripts/{species}.fa", "fasta")
    SeqIO.write(protein_records, f"Results/Genome_data/Protein/primary_transcripts/{species}.fa", "fasta")