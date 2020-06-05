def translateThat(dna):

    dna1to3 = dna[0:3]
    dna3to6 = dna[3:6]
    dnarest = dna[-3:]

# table gotten from https://www.geeksforgeeks.org/dna-protein-python-3/
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '', 'TAG': '',
        'TGC': 'C', 'TGT': 'C', 'TGA': '', 'TGG': 'W',
    }

    protein = ""
    if len(dna1to3) % 3 == 0:
        for i in range(0, len(dna1to3), 3):
            codon = dna1to3[i:i + 3]
            protein += table[codon]

    if len(dna3to6) % 3 == 0:
        for i in range(0, len(dna3to6), 3):
            codon = dna3to6[i:i + 3]
            protein += table[codon]

    if len(dnarest) % 3 == 0:
        for i in range(0, len(dnarest), 3):
            codon = dnarest[i:i + 3]
            protein += table[codon]

    return protein


def mutate():
    DnaFile = open("DNAFile.txt", "r+")

    normalDNA = DnaFile.read()
    normalDNA = normalDNA.replace("\n", "")

    normalDNA = normalDNA.replace("a", "A")

    DnaFile.close()

    normal = open("normalDNA.txt", "w")

    normal.write(normalDNA)

    normal.close()

    DnaFile = open("DNAFile.txt", "r+")

    normalDNA = DnaFile.read()
    normalDNA = normalDNA.replace("\n", "")

    normalDNA = normalDNA.replace("a", "T")

    DnaFile.close()

    mutated = open("mutatedDNA.txt", "w")

    mutated.write(normalDNA)

    mutated.close()


def txtTranslate():
    normalDNA = open("normalDNA.txt", "r")
    mutated = open("mutatedDNA.txt", "r")

    return normalDNA.read() + "\n" + mutated.read()


print(txtTranslate())
