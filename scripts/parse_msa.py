# Script que toma como primer argumento un alineamiento y extrae el subalineamiento
# donde esta nuestra muestra de interes

import sys
import re
from Bio import SeqIO, AlignIO

infile = sys.argv[1]
name = re.sub("_aln.*", "", infile)
file = open(infile, "r")
output_name = str(name) + "_subaln.fasta"
output = open(output_name, "w")


my_dict = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
align = AlignIO.read(infile, "fasta")
# Mi muestra de referencia es la ultima del alineamiento
id=align[len(align)-1].id

value = my_dict.get(id)

longitud = []

sec = ""


for character in value.seq:
    sec += str(character)


for character in range(0, len(sec)):
    if sec[character] != "-":
        longitud.append(character)


start = longitud[0]
end = longitud[len(longitud)-1]

# Porque el alineamiento no empieza desde la pos 0
subalign = align[:,start:end+1]
AlignIO.write(subalign, output, "fasta")
