
from math import *
import random

def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return  count

def CountWithPseudocounts(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return  count

def Profile(Motifs):
    profile = {}
    count = Count(Motifs)
    k = len(Motifs[0])
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            total = count["A"][j] + count["C"][j] + count["G"][j] + count["T"][j]
            profile[symbol].append(count[symbol][j]/total)

    return profile

def ProfileWithPseudocounts(Motifs):
    profile = {}
    count = CountWithPseudocounts(Motifs)
    k = len(Motifs[0])
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            total = count["A"][j] + count["C"][j] + count["G"][j] + count["T"][j]
            profile[symbol].append(count[symbol][j]/total)

    return profile

def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol

    return consensus


def Score(Motifs):
    score = 0
    k = len(Motifs[0])
    t = len(Motifs)
    consensus = Consensus(Motifs)

    for j in range(k):
        cont = 0
        for i in range(t):
            if Motifs[i][j] != consensus[j]:
                cont += 1
        score += cont

    return score

def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return p

def ProfileMostProbablePattern(Text, k, Profile):
    """ Explicación de porque poner la pMax a -1
    (Sacado del comentario de Cris Lawrence ) a la tarea:
    I too did this in PyCharm and got the expected answer but failed the test until
    I also set p to -1 initially as someone else did.
    I guess the issue is that by setting p = 0 for starters,
    you could wind up with just an empty string if all possibilities
    occur with 0 probability.  The algorithm spec says return the first best.
    If all have 0 probability, that would be the first 5-mer in the sequence.
    Setting the probability to less than 0 guarantees you return something.
    I guess that's why an initial value for p of 0 doesn't work. """
    pMax = -1
    mostProbPattern = ""
    for i in range(0, len(Text)-k + 1):
        pattern = Text[i:i+k]
        p = Pr(pattern, Profile )
        if p > pMax:
            mostProbPattern = pattern
            pMax = p

    return mostProbPattern

def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs

def Entropy(Profile):
    totalEntropy = 0
    k = len(Profile["A"])

    for j in range(k):
        colEntropy = 0
        for i in ["A","C","G","T"]:
            p = Profile[i][j]
            if p == 0:
                colEntropy += 0
            else:
                colEntropy +=  p * log(p,2)

        totalEntropy += -1 * colEntropy

    return totalEntropy

# Input:  A profile matrix Profile and a list of strings Dna
# Output: Motifs(Profile, Dna)
def Motifs(Profile, Dna,k):
    # insert your code here
    motifs = []
    for i in range(len(Dna)):
        kmer = ProfileMostProbablePattern(Dna[i], k, Profile)
        motifs.append(kmer)

    return motifs

def RandomMotifs(Dna, k, t):
    motifs = []
    for i in range(t):
        init = random.randint(0,len(Dna[0])-k-1)
        motifs.append(Dna[i][init:init+k])
    return motifs

def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna,k)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs

def Normalize(Probabilities):
    P = 0
    for key in Probabilities:
        P += Probabilities[key]
    for key in Probabilities:
        Probabilities[key] = Probabilities[key] / P
    return Probabilities


def WeightedDie(Probabilities):
    """
        >>> lineaProb = {"A":0.1, "C":0.1, "G":0, "T":0.8}
        >>> lineaProb
            {'A': 0.1, 'C': 0.1, 'T': 0.8, 'G': 0}
        >>> sorted(lineaProb.items(), key = lambda x: x[1])
            [('G', 0), ('A', 0.1), ('C', 0.1), ('T', 0.8)]
    """

    lineaProb = {}
    probAcumulada = 0
    for key in Probabilities.keys():
        lineaProb[key] = probAcumulada + Probabilities[key]
        probAcumulada = probAcumulada + Probabilities[key]

    p = random.uniform(0, 1)

    lineaProbOrd = sorted(lineaProb.items(), key = lambda x: x[1])

    for tupla in lineaProbOrd:
        if p <= tupla[1]:
            return tupla[0]


def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)

    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)


def GibbsSampler(Dna, k , t, N):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    for j in range(N):
        i = random.randint(0,t-1)
        M.pop(i)
        Profile =  ProfileWithPseudocounts(M)
        newMotif =  ProfileGeneratedString(Dna[i], Profile, k)
        M.insert(i,newMotif)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M

    return BestMotifs


probabilities = {1:0.22, 2:0.54, 3:0.58, 4:0.36, 5:0.3}
print(Normalize(probabilities))

# Dna = ["AAGCCAAA","AATCCTGG","GCTACTTG","ATGTTTTG"]
#
# M = ["CCA","CCT","CTT","TTG"]
#
#
# Profile = ProfileWithPseudocounts(M)
# newM = Motifs(Profile, Dna, 3)
#
# print(newM)







# k = 3
# t = 5
# Dna = ["TTACCTTAAC","GATGTCTGTC","ACGGCGTTAG","CCCTAACGAG","CGTCAGAGGT"]
#
#
# print(random.randint(0,len(Dna[0])-k))
# print(RandomMotifs(Dna, k, t))

# Profile = {"A":[0.4,0.3,0.0,0.1,0.0,0.9],
#            "C":[0.2,0.3,0.0,0.4,0.0,0.1],
#            "G":[0.1,0.3,1.0,0.1,0.5,0.0],
#            "T":[0.3,0.1,0.0,0.4,0.5,0.0]}
#
# Text = "CAGTGA"
# print(Pr(Text,Profile))
# Motifs = ["AACGTA",
#            "CCCGTT",
#            "CACCTT",
#            "GGATTA",
#            "TTCCGG"]
#
# print(Score(Motifs))
