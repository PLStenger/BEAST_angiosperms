from Bio.Data import CodonTable
from Bio import SeqIO
import statistics

standard_table = CodonTable.unambiguous_dna_by_id[11]
codon2aa = standard_table.forward_table

def count_syn_non_syn_sites(seq):
    seq = seq.upper()
    n_sites = 0
    s_sites = 0
    bases = ['A', 'T', 'C', 'G']
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) != 3 or "N" in codon: continue
        aa_ref = codon2aa.get(codon, "X")
        for pos in range(3):
            for base in bases:
                if base == codon[pos]: continue
                mut_codon = codon[:pos] + base + codon[pos+1:]
                aa_mut = codon2aa.get(mut_codon, "X")
                if aa_mut == "X": continue
                if aa_mut == aa_ref:
                    s_sites += 1/3
                else:
                    n_sites += 1/3
    return s_sites, n_sites

s_all = []
n_all = []

for record in SeqIO.parse("/Users/stengerpierre-louis/Downloads/Vvinifera_145_cds.fa", "fasta"):
    s, n = count_syn_non_syn_sites(str(record.seq))
    s_all.append(s)
    n_all.append(n)

# Statistiques pour S
print("S (sites synonymes)")
print(f"Moyenne : {statistics.mean(s_all):.2f}")
print(f"Médiane : {statistics.median(s_all):.2f}")
print(f"Écart-type : {statistics.stdev(s_all):.2f}")

# Statistiques pour N
print("\nN (sites non synonymes)")
print(f"Moyenne : {statistics.mean(n_all):.2f}")
print(f"Médiane : {statistics.median(n_all):.2f}")
print(f"Écart-type : {statistics.stdev(n_all):.2f}")
