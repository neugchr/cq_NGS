import argparse

parser = argparse.ArgumentParser(
    description="parse a fasta file"
)

# # positional argument mandatory
parser.add_argument("infile", type=argparse.FileType("r"), help=("turns multi-fasta files into single sequence fastas, the sequence is in one line"))
args = parser.parse_args()



class Sequence:
    def __init__(self, description, sequence):
        self.description = description
        self.sequence = sequence

    def __str__(self):
        ostr = f"\n--------------------------\n{self.description}\n{self.sequence}\n"
        return(ostr)


def parse_file(fh):
    sequences = []
    description = ''
    sequence = ''
    for line in fh:
        if line[0] == '>':
            #write the previous sequence
            if description != '':
                sequences.append( Sequence(description, sequence) )
            description = line
            sequence = ''
        else:
            sequence += line
    # now write the last sequence as well
    sequences.append( Sequence(description, sequence) )
    return(sequences)



with args.infile as file:
    seqs = parse_file(file)


for seq in seqs:
    fname = seq.description[1:15].replace(" ", "_")
    with open(f"{fname}.fasta", "w+") as fh:
        outstring = seq.description.replace("\n", "")
        outstring += '\n'+seq.sequence.replace("\n", "")
        fh.write(outstring)
