import argparse

parser = argparse.ArgumentParser(
    description="parse a fasta file"
)

# # positional argument mandatory
parser.add_argument("infile", type=argparse.FileType("r"), help=("input fasta file"))
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


with open('results.txt', 'w+') as fh:
    for i in range(len(seqs)-1,-1,-1):
        print(str(seqs[i]))
        fh.write(str(seqs[i]))
        seqs.pop(i)
