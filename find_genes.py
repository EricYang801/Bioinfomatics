from pyfaidx import Fasta

genome = Fasta("Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa")

seq = genome["1"][1471765 - 1:1497848]

print(">ATAD3B chr1:1471765-1497848 strand=+")
print(str(seq).upper())