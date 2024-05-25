import argparse
import hashlib
from typing import List

# inputs
parser = argparse.ArgumentParser(prog="Vanheusden", 
                                 description="kmer distribution of reads\ncol 1: count, col 2: number of kmers that appeared count times")
parser.add_argument("-k", "--kmerLength")
#parser.add_argument('-i', '--ignoreDirectionality', action='store_true', required=False)
parser.add_argument("reads")
args = parser.parse_args()

k = int(args.kmerLength)
#ignoreDirectionality = args.ignoreDirectionality
output = args.reads.split(".")[0] + ".histo.tsv"
input = args.reads
bases = ["A", "T", "C", "G"]

# histogram
distribution = {}
for i in range(1, 500):
    distribution[i] = 0

# other parameters
width = 1000
initMin = 9999
numFunctions = 10

# for ignoreDirectionality
def ReverseComplement(pattern: str) -> str:
    reverse_complement = ""
    start_idx = len(pattern) - 1
    for i in range(start_idx, -1, -1):
        if pattern[i] == "A":
            reverse_complement += "T"
        if pattern[i] == "T":
            reverse_complement += "A"
        if pattern[i] == "C":
            reverse_complement += "G"
        if pattern[i] == "G":
            reverse_complement += "C"
    return reverse_complement

# Count-min sketch
countMinSketch = [[0 for j in range(1000)] for i in range(numFunctions)]

# choose between strand and reverse complement
def ChooseStrand(kmer: str) -> str:
    kmerRev = ReverseComplement(kmer)
    if hash(kmer) < hash(kmerRev):
        return kmer
    return kmerRev

# return a list of hash values
def Hashes(kmer: str) -> List[int]:
    hashes = []

    # initialize hash functions
    functions = []
    functions.append(hashlib.sha1())
    functions.append(hashlib.sha224())
    functions.append(hashlib.sha256())
    functions.append(hashlib.sha384())
    functions.append(hashlib.sha512())
    functions.append(hashlib.sha3_224())
    functions.append(hashlib.sha3_256())
    functions.append(hashlib.sha3_384())
    functions.append(hashlib.sha3_512())
    functions.append(hashlib.md5())

    # update hash functions
    code = kmer.encode("utf-8")
    for i in range(numFunctions):
        functions[i].update(code)

    # compute hash values
    for i in range(numFunctions):
        hashes.append(int(functions[i].hexdigest(), 16) % width)
    return hashes

# update countMinSketch
def updateCount(kmer: str):
    i = 0
    for j in Hashes(kmer):
        countMinSketch[i][j] += 1
        i += 1

# return counts for a kmer
def getCount(kmer: str) -> int:
    min = initMin
    i = 0
    for j in Hashes(kmer):
        if countMinSketch[i][j] < min:
            min = countMinSketch[i][j]
        i += 1
    return min

# fill out countMinSketch
with open(input, "r") as file:
    for line in file:
        if line[0] in bases:
            line = line.strip()
            for i in range(len(line) - k + 1):
                kmer = ChooseStrand(line[i : i + k])
                updateCount(kmer)

# read values from countMinSketch
with open(input, "r") as file:
    for line in file:
        if line[0] in bases:
            line = line.strip()
            for i in range(len(line) - k + 1):
                kmer = ChooseStrand(line[i : i + k])
                distribution[getCount(kmer)] += 1

# write to tsv
with open(output, "w") as file:
    counter = 0
    for count in distribution:
        if distribution[count] == 0:
            counter += 1
        if counter >= 10:
            break
        file.write(str(count) + "\t" + str(distribution[count]) + "\n")
