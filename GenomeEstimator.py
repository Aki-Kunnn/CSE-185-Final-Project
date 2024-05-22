
import hashlib

k = 10
ignoreDirectionality = True
output = "out.tsv"
bases = ["A", "T", "C", "G"]
AorC = ["A", "C"]

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

frequencies = {}
for i in range(1, 200):
    frequencies[i] = 0

sha1 = hashlib.sha1()
sha224 = hashlib.sha224()
sha256 = hashlib.sha256()

list1 = [0] * 1000
list224 = [0] * 1000
list256 = [0] * 1000

def Hash(kmer: str):
    idx1 = sha1.update(kmer).hexdigest() % 1000
    idx224 = sha224.update(kmer).hexdigest() % 1000
    idx256 = sha256.update(kmer).hexdigest() % 1000
    list1[idx1] += 1
    list224[idx224] += 1
    list256[idx256] += 1

def HashCount(kmer: str) -> int:
    idx1 = sha1.update(kmer).hexdigest() % 1000
    idx224 = sha224.update(kmer).hexdigest() % 1000
    idx256 = sha256.update(kmer).hexdigest() % 1000
    min = list1[idx1]
    if list224[idx224] < min:
        min = list224[idx224]
    if list256[idx256] < min:
        min = list256[idx256]
    return min

with open("testing.txt", "r") as file:
    for line in file:
        if line[0] in bases:
            line = line[:-1]
            for i in range(len(line) - k + 1):
                kmer = line[i : i + k]
                if kmer[0] not in AorC:
                    kmer = ReverseComplement(kmer)
                Hash(kmer)
    
with open("testing.txt", "r") as file:
    for line in file:
        if line[0] in bases:
            line = line[:-1]
            for i in range(len(line) - k + 1):
                kmer = line[i : i + k]
            if kmer[0] not in AorC:
                    kmer = ReverseComplement(kmer)
            frequencies[HashCount(kmer)] += 1

