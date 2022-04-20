
class Sequence:
    def __init__(self, description, sequence):
        self.description = description
        self.sequence = sequence
        self.suffixes = []

    def __str__(self):
        ostr = f"\n--------------------------\n{self.description}\n{self.sequence}\n"
        for i in self.suffixes:
            ostr += i+'\n'
        return(ostr)

    def create_suffixes(self):
        for i in range(len(self.sequence)):
            self.suffixes.append(self.sequence[i:])
        self.suffixes.sort()


def parse_file(file):
    sequences = []
    with open(file, "r") as fh:
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

seqs = parse_file('unknown.fasta')
#seqs = parse_file('shortened.fasta')
with open('results.txt', 'w+') as fh:
    for i in range(len(seqs)-1,-1,-1):
        seqs[i].create_suffixes()
        print(str(seqs[i]))
        fh.write(str(seqs[i]))
        seqs.pop(i)
