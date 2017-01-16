def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

def CountDict(Text, k):
    Count = {}
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        Count[i] = PatternCount(Pattern, Text)
    return Count

def remove_duplicates(lista):
    listaUnicos = []
    for elem in lista:
        if elem not in listaUnicos:
            listaUnicos.append(elem)
    return listaUnicos

def FrequentWords(Text, k):
    FrequentPatterns = []
    Count = CountDict(Text, k)
    m = max(Count.values())
    for i in Count:
        if Count[i] == m:
            FrequentPatterns.append(Text[i:i+k])
    FrequentPatternsNoDuplicates = remove_duplicates(FrequentPatterns)
    return FrequentPatternsNoDuplicates

def reverse(text):
    reverse_text_list = []
    for i in range(len(text)):
        reverse_text_list.append(text[len(text)-1-i])
    return "".join(reverse_text_list)

def complement(Nucleotide):
    comp = '' # output variable
    # your code here
    if Nucleotide == "A":
        comp = "T"
    elif Nucleotide == "T":
        comp = "A"
    elif Nucleotide == "G":
        comp = "C"
    elif Nucleotide == "C":
        comp = "G"
    else:
        comp = Nucleotide

    return comp

def ReverseComplement(Pattern):
    revComp = '' # output variable
    # your code here
    rev = reverse(Pattern)
    for nucleotide in rev:
        revComp = revComp + complement(nucleotide)

    return revComp

def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    # your code here
    for i in range(len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions

def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(symbol, Genome[0:n//2])
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

def Skew(Genome):
    skew = {} #initializing the dictionary
    # your code here
    n = len(Genome)
    skew[0] = 0
    for i in range(1,n+1):
        if Genome[i-1] == "G":
            skew[i] = skew[i-1]+1
        elif Genome[i-1] == "C":
            skew[i] = skew[i-1]-1
        else:
            skew[i] = skew[i-1]

    return skew

def MinimumSkew(Genome):
    positions = [] # output variable
    # your code here
    skew = Skew(Genome)
    minValue = skew[0]
    positions.append(minValue)
    for i in range(1,len(Genome)):
        if skew[i] < minValue:
            minValue = skew[i]
            positions = [i]
        if skew[i] == minValue and i not in positions:
            positions.append(i)
    return positions

def MaximumSkew(Genome):
    positions = [] # output variable
    # your code here
    skew = Skew(Genome)
    maxValue = skew[0]
    positions.append(maxValue)
    for i in range(1,len(Genome)):
        if skew[i] > maxValue:
            maxValue = skew[i]
            positions = [i]
        if skew[i] == maxValue and i not in positions:
            positions.append(i)
    return positions

def HammingDistance(s1, s2):
    distancia = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            distancia = distancia + 1
    return distancia

# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Pattern, Text, d):
    positions = [] # output variable
    # your code here
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions

def ApproximatePatternCount(Pattern, Text, d):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)],Pattern) <= d :
            count = count+1
    return count

# print(PatternCount("TGT", "ACTGTACGATGATGTGTGTCAAAG"))
# Text = "TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT"
# print(FrequentWords(Text, 3))
# print( ReverseComplement("GATTACA") )

#Genome = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
#print(MinimumSkew(Genome))

s1 = "CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA"
s2 = "CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG"
print(HammingDistance(s1,s2))

#Genome = "CATTCCAGTACTTCGATGATGGCGTGAAGA"
#print(MaximumSkew(Genome))
