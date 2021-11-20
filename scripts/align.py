import os
import filelist
import parsefas
import numpy as np


class ProteinAligner:
    def __init__(self, align_tool):
        self.align_tool = align_tool.upper()
        self.alignment_path = f"./Alignment/Protein/NotTrimmed/{self.align_tool}/"  # Alignment output 파일이 저장되는 경로
        self.trimmed_path = f"./Alignment/Protein/Trimmed/{self.align_tool}_gblocks/"
        os.makedirs(self.alignment_path, exist_ok=True)
        os.makedirs(self.trimmed_path, exist_ok=True)
        self.tool_dict = {'PRANK': self.prank, 'MAFFT': self.mafft, 'MUSCLE': self.muscle}

    def prank(self, gene_path):
        gene_name = gene_path.split('/')[-1].split('.fasta')[0]
        os.system(f"prank -protein -F -quiet -d={gene_path} -o={os.path.join(self.alignment_path, gene_name)}")

    def mafft(self, gene_path):
        gene_name = gene_path.split('/')[-1].split('.fasta')[0]
        os.system(f"mafft-linsi --anysymbol --quiet {gene_path} > {os.path.join(self.alignment_path, gene_name)}.best.fas")

    def muscle(self, gene_path):
        gene_name = gene_path.split('/')[-1].split('.fasta')[0]
        os.system(f"muscle -in {gene_path} -out {os.path.join(self.alignment_path, gene_name)}.best.fas -seqtype protein -quiet")

    def replace_seleno(self, alignment_path):
        with open(alignment_path) as f:
            F = f.readlines()
        newfile = [line.replace('U', 'X') if not line.startswith('>') else line for line in F]

        with open(alignment_path, 'w') as g:
            for line in newfile:
                g.write(line)

    def reorder(self, alignment_path):
        with open(alignment_path) as f:
            F = f.read().split('>')[1:]
        F.sort()

        with open(alignment_path, 'w') as g:
            for unit in F:
                new_unit = '>' + unit
                g.write(new_unit)

    def gblocks(self, alignment_path):
        gene_name = alignment_path.split('/')[-1].split('.best.fas')[0]
        os.system(f"Gblocks {alignment_path} -t=p -p=n")
        os.system(f"cp {alignment_path}-gb ./Alignment/Protein/Trimmed/{self.align_tool}_gblocks/{gene_name}.fas")
        os.system(f"rm {alignment_path}-gb")

    def align(self, gene_path):
        gene_name = gene_path.split('/')[-1].split('.fasta')[0] + '.best.fas'
        self.tool_dict[self.align_tool](gene_path)
        aligned_path = os.path.join(self.alignment_path, gene_name)
        self.replace_seleno(aligned_path)
        self.reorder(aligned_path)
        self.gblocks(aligned_path)


class CodonAligner:
    def __init__(self, align_tool):
        self.align_tool = align_tool.upper()
        self.prot_align_path = f"./Alignment/Protein/NotTrimmed/{self.align_tool}"  # Protein alignment 파일이 있는 경로
        self.codon_align_path = f"./Alignment/Codon/NotTrimmed/{self.align_tool}_pal2nal"  # Codon alignment 파일이 저장될 경로
        self.codon_trimmed_path = f"./Alignment/Codon/Trimmed/{self.align_tool}_pal2nal_gblocks"
        os.makedirs(self.codon_align_path, exist_ok=True)
        os.makedirs(self.codon_trimmed_path, exist_ok=True)
        os.makedirs(self.codon_trimmed_path + '_paml', exist_ok=True)

        self.gene_list = filelist.mklist('.best.fas', self.prot_align_path)  # 정렬된 단백질 서열의 유전자 이름 리스트

        with open('./Alignment/Codon/filtered_cds.txt', 'w') as f:
            f.write('Filtered gene list\n\n')

    def align(self, gene):
        seq_dict = parsefas.MakeSeqDict(f"./All_sequences/CDS_sequence/{gene}.fasta", no_enter=True)
        remainder = np.array([len(seq) % 3 for seq in seq_dict.values()])
        if np.all(remainder == 0):
            self.pal2nal(gene)
            self.gblocks(gene)
            self.convert_to_paml(gene)

        else:
            with open('./Alignment/Codon/filtered_cds.txt', 'a') as f:
                f.write(f"{gene} (Not multiple of 3)\n")

    def pal2nal(self, gene):
        os.system(f"pal2nal.pl {self.prot_align_path}/{gene}.best.fas ./All_sequences/CDS_sequence/{gene}.fasta -output fasta > {self.codon_align_path}/{gene}.fas -nogap -nomismatch")

    def gblocks(self, gene):
        os.system(f"Gblocks {self.codon_align_path}/{gene}.fas -t=c -p=n")
        os.system(f"cp {self.codon_align_path}/{gene}.fas-gb {self.codon_trimmed_path}/{gene}.fas")
        os.system(f"rm {self.codon_align_path}/{gene}.fas-gb")

    def convert_to_paml(self, gene):
        if os.path.exists(f"{self.codon_trimmed_path}/{gene}.fas"):
            seq_dict = parsefas.MakeSeqDict(f"{self.codon_trimmed_path}/{gene}.fas")
            for key, seq in seq_dict.items():
                seq_dict[key] = seq.replace(' ', '')
            species_num = len(seq_dict)
            seq_len = len(list(seq_dict.values())[0].replace('\n', ''))
            header = f"  {species_num}   {seq_len}\n"
            with open(f"{self.codon_trimmed_path}_paml/{gene}.fas", 'w') as f:
                f.write(header)
                for key, seq in seq_dict.items():
                    unit = f"{key}\n{seq}"
                    f.write(unit)
        else:
            with open('./Alignment/Codon/filtered_cds.txt', 'a') as f:
                f.write(f"{gene} (Pal2nal error)\n")


if __name__ == '__main__':
    al = CodonAligner('prank')
    for gene in al.gene_list:
        al.align(gene)



