from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess, os, tempfile

def run_muscle(input_sequences, muscle_exe="muscle", need_log=False):
    if type(input_sequences) is str:
        input_file_name = input_sequences
    elif type(input_sequences) is list:
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as input_file:
            records = [SeqRecord(Seq(seq), id=f"Seq{i+1}") for i, seq in enumerate(input_sequences)]
            SeqIO.write(records, input_file, "fasta")
            input_file_name = input_file.name
    else:
        raise ValueError("input_sequences must be either a path or a list of strings")

    output_file = tempfile.NamedTemporaryFile(delete=False)
    output_file_name = output_file.name
    output_file.close()

    if need_log:
        subprocess.run([muscle_exe, "-align", input_file_name, "-output", output_file_name])
    else:
        with open(os.devnull, 'w') as devnull:
            subprocess.run([muscle_exe, "-align", input_file_name, "-output", output_file_name],
                        stdout=devnull, stderr=devnull)

    alignment = list(SeqIO.parse(output_file_name, "fasta"))

    if input_file.__class__ == tempfile._TemporaryFileWrapper:
        os.remove(input_file_name)
    os.remove(output_file_name)

    return alignment

if __name__ == "__main__":

    sequences = [
        "ATGCTACGATCG",
        "ATGCTACGTGC",
        "ATGCTAACGATCG"
    ]

    muscle = "muscle"
    
    alignment = run_muscle(sequences, muscle_exe=muscle, need_log=False)

    for record in alignment:
        print(f">{record.id}\n{record.seq}")