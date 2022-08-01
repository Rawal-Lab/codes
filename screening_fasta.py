from Bio import SeqIO
import sys
import os

fasta_path = os.path.abspath(sys.argv[1])
file_name = os.path.basename(fasta_path).split(".")[0]

with open(fasta_path) as f:
    seq_records = [s for s in SeqIO.parse(f, "fasta")]
    passed_sequence_records = [seq_record for seq_record in seq_records if len(seq_record.seq) > 9]

visited = []
final_seq_records = []
for i, seqrd in enumerate(passed_sequence_records):
    replace_with = f"SEQ_{i}_{seqrd.id}"
    if not seqrd.id in visited:
        visited.append(seqrd.id)
        seqrd.id = replace_with
        final_seq_records.append(seqrd)
    else:
        print(f"WARNING: Repeated ID --> {seqrd.id}")

with open(f"{file_name}_formatted.fasta", "w") as f:
    SeqIO.write(final_seq_records, f, format="fasta")

invalid_amino_acids = ["B", "b", "X", "x", "Z", "z"]
res = ""
current_id = ""
with open(f"{file_name}_formatted.fasta") as f:
    for i, line in enumerate(f):
        if not line.strip().startswith(">"):
            for char in invalid_amino_acids:
                if char in line:
                    print(f"WARNING: Removing invalid character {char} in ID: {current_id}")
                    line = line.replace(char, "")
        else:
            current_id = line.strip(">").split(" ")[0]
        res += line

with open(f"{file_name}_formatted.fasta", "w") as f:
    f.write(res)