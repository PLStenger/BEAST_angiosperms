import re

keep_taxa = {
    "Arabidopsis_thaliana", "Arabidopsis_lyrata", "Arabidopsis_halleri",
    "Capsella_rubella", "Thellungiella_parvula", "Brassica_oleracea",
    "Brassica_napus_C", "Brassica_rapa", "Brassica_juncea", "Brassica_nigra"
}

input_file = "/Users/stengerpierre-louis/Documents/INRAE_PaleoLab/03_angiosperms_origins/04_beauti/SpeciesTreeAlignment_filtered.fa.nex"
output_file = "/Users/stengerpierre-louis/Documents/INRAE_PaleoLab/03_angiosperms_origins/08_Beast_by_family/SpeciesTreeAlignment_filtered_brassicaceae.nex"

with open(input_file) as infile:
    lines = infile.readlines()

# Variables pour parser le NEXUS en interleave
in_matrix = False
taxa_seqs = {}
taxa_order = []
nchar = None
header_lines = []
footer_lines = []
i = 0
while i < len(lines):
    line = lines[i]
    # Detect début matrix
    if re.match(r'^\s*matrix\s*$', line, re.IGNORECASE):
        in_matrix = True
        i += 1
        # Parse blocs interleaved
        while i < len(lines):
            line = lines[i]
            if re.match(r'^\s*;\s*$', line):
                in_matrix = False
                break
            # Séquence potentielle dans format interleave :
            # - lignes de séquences avec taxon seq_part
            # - ou lignes vides répétées séparant blocs
            if line.strip() == "":
                i += 1
                continue
            parts = line.rstrip().split(None, 1)
            if len(parts) == 2:
                taxon, seqpart = parts
                if taxon not in taxa_seqs:
                    taxa_seqs[taxon] = seqpart.strip()
                    taxa_order.append(taxon)
                else:
                    taxa_seqs[taxon] += seqpart.strip()
            else:
                # ligne qui peut contenir uniquement des séquences (interleave)
                # à appliquer à toutes séquences précédemment enregistrées
                seqpart = line.strip()
                if seqpart != "":
                    for t in taxa_order:
                        taxa_seqs[t] += seqpart
            i += 1
        break
    else:
        header_lines.append(line)
    i += 1

# Récupérer le footer après ; du matrix
footer_lines = lines[i+1:]

# Récupérer nchar actuel pour recalcul dimensions plus tard
for hl in header_lines:
    match = re.match(r'^\s*dimensions\s+ntax\s*=\s*\d+\s+nchar\s*=\s*(\d+)\s*;?', hl, re.IGNORECASE)
    if match:
        nchar = int(match.group(1))
        break

# Filtrer taxa selon liste keep_taxa et séquences valides (non vide)
filtered_taxa_seqs = {taxon: seq for taxon, seq in taxa_seqs.items() if taxon in keep_taxa and len(seq.replace("-", "").replace("X","").replace("?","")) > 0}

# Mettre à jour nchar selon la longueur séquence (on prend la 1ère séquence car elles doivent être toutes même longueur)
if filtered_taxa_seqs:
    nchar_new = len(next(iter(filtered_taxa_seqs.values())))
else:
    nchar_new = 0
ntax_new = len(filtered_taxa_seqs)

# Ecrire fichier NEXUS en sequential (une séquence par ligne après MATRIX)
with open(output_file, "w") as outfile:
    # écrire header avec mise à jour dimensions
    for line in header_lines:
        if re.match(r'^\s*dimensions\s+ntax\s*=\s*\d+\s+nchar\s*=\s*\d+\s*;?', line, re.IGNORECASE):
            outfile.write(f"dimensions ntax={ntax_new} nchar={nchar_new};\n")
        else:
            outfile.write(line)
    # écrire matrix
    outfile.write("matrix\n")
    for taxon in filtered_taxa_seqs:
        outfile.write(f"{taxon:<25} {filtered_taxa_seqs[taxon]}\n")
    outfile.write(";\n")
    # écrire footer (reste du fichier)
    for line in footer_lines:
        outfile.write(line)

print(f"Fichier filtré écrit dans '{output_file}' avec {ntax_new} taxa et {nchar_new} caractères par séquence.")


