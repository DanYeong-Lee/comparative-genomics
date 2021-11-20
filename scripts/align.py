import os
import filelist
import parsefas
import numpy as np


class ProteinAligner:
    def __init__(self, align_tool):
        self.align_tool = align_tool.upper()
        self.alignment_path = f"Results/Alignment/Protein/NotTrimmed/{self.align_tool}/"  # Alignment output 파일이 저장되는 경로
        self.trimmed_path = f"Results/Alignment/Protein/Trimmed/{self.align_tool}_gblocks/"
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
        os.system(f"cp {alignment_path}-gb Results/Alignment/Protein/Trimmed/{self.align_tool}_gblocks/{gene_name}.fas")
        os.system(f"rm {alignment_path}-gb")

    def align(self, gene_path):
        gene_name = gene_path.split('/')[-1].split('.fasta')[0] + '.best.fas'
        self.tool_dict[self.align_tool](gene_path)
        aligned_path = os.path.join(self.alignment_path, gene_name)
        self.replace_seleno(aligned_path)
        self.reorder(aligned_path)
        self.gblocks(aligned_path)




