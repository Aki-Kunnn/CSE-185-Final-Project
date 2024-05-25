import argparse
import hashlib
import random
from typing import List

bases = ["A", "T", "C", "G"]

def main():
    parser = argparse.ArgumentParser(prog="Vanheusden Genome Size Estimator")
    
    subparsers = parser.add_subparsers(dest = "function")
    parserGenerate = subparsers.add_parser("generate")
    parserEstimate = subparsers.add_parser("estimate")

    parserGenerate.add_argument("-g", "--genomeLength", required=True, type = int)
    parserGenerate.add_argument("-n", "--numberOfReads", required=True, type = int)
    parserGenerate.add_argument("-l", "--lengthOfReads", required=True, type = int)
    parserGenerate.add_argument("-e", "--errorRate", required=True, type=float)

    parserEstimate.add_argument("-k", "--kmerLength", required=True, type=int)
    parserEstimate.add_argument("reads")

    args = parser.parse_args()


    if args.function == "generate":
        g = args.genomeLength
        n = args.numberOfReads
        l = args.lengthOfReads
        e = args.errorRate 
        output = "g" + str(g) + "_n" + str(n) + "_l" + str(l) + "_e" + str(e) + ".fasta"

        # used for generating errors in reads
        def GenerateError(nucleotide: str) -> str:
            bases = ["A", "T", "C", "G"]
            bases.remove(nucleotide)
            return bases[random.randint(0, 2)]

        # generate genome
        genome = ""
        for i in range(g):
            genome += bases[random.randint(0, 3)]

        # generate reads
        reads = []
        for i in range(n):
            start = random.randint(0, len(genome) - l + 1)
            read = genome[start : start + l]
            for i in range(len(read)):
                rng = random.randint(0, 1 / e)
                if rng == 0:
                    read = read[:i] + GenerateError(read[i]) + read[i+1:]
            reads.append(read)

        # write to 
        with open(output, "w") as file:
            file.write(">header\n")
            for read in reads:
                file.write(read + "\n")

    elif args.function == "estimate":
        k = args.kmerLength
        input = args.reads
        output = "k" + str(k) + "_histogram.tsv"

        distribution = {i: 0 for i in range(1, 500)}
        width = 10000
        numFunctions = 10

        # used for ignoreDirectionality
        def ReverseComplement(pattern: str) -> str:
            reverse_complement = ""
            start_idx = len(pattern) - 1
            for i in range(start_idx, -1, -1):
                if pattern[i] == "A":
                    reverse_complement += "T"
                elif pattern[i] == "T":
                    reverse_complement += "A"
                elif pattern[i] == "C":
                    reverse_complement += "G"
                elif pattern[i] == "G":
                    reverse_complement += "C"
            return reverse_complement
        
        # Count-min sketch
        countMinSketch = [[0 for j in range(width)] for i in range(numFunctions)]

        # choose which strand will represent the sequence
        def ChooseStrand(kmer: str) -> str:
            kmerRev = ReverseComplement(kmer)
            codeA = kmer.encode("utf-8")
            codeB = kmerRev.encode("utf-8")
            
            sha1 = hashlib.sha1()
            sha1.update(codeA)
            a = int(sha1.hexdigest(), 16) % width

            sha1 = hashlib.sha1()
            sha1.update(codeB)
            b = int(sha1.hexdigest(), 16) % width
            if a < b:
                return kmer
            return kmerRev

        # return a list of hash values
        def Hashes(kmer: str) -> List[int]:
            hashes = []
            functions = [hashlib.sha1(), hashlib.sha224(), hashlib.sha256(), 
                        hashlib.sha384(), hashlib.sha512(), hashlib.sha3_224(), 
                        hashlib.sha3_256(), hashlib.sha3_384(), hashlib.sha3_512(), 
                        hashlib.md5()]
            code = kmer.encode("utf-8")
            for i in range(numFunctions):
                functions[i].update(code)
                hashes.append(int(functions[i].hexdigest(), 16) % width)
            return hashes
        
        # update countMinSketch
        def updateCount(kmer: str):
            for i, j in enumerate(Hashes(kmer)):
                countMinSketch[i][j] += 1

        # return counts for a kmer
        def getCount(kmer: str) -> int:
            return min(countMinSketch[i][j] for i, j in enumerate(Hashes(kmer)))

        # fill out countMinSketch
        total = 0
        with open(input, "r") as file:
            for line in file:
                if line[0] in bases:
                    line = line.strip()
                    for i in range(len(line) - k + 1):
                        kmer = ChooseStrand(line[i : i + k])
                        updateCount(kmer)
                        total += 1

        # read values from countMinSketch
        with open(input, "r") as file:
            for line in file:
                if line[0] in bases:
                    line = line.strip()
                    for i in range(len(line) - k + 1):
                        kmer = ChooseStrand(line[i : i + k])
                        distribution[getCount(kmer)] += 1

        kmerCoverage = 0
        peak = 0
        for abundance, numKmers in distribution.items():
            if peak < numKmers:
                peak = numKmers
                kmerCoverage = abundance
        genomeSize = int(total / kmerCoverage)
        print("Total number of kmers:\t" + str(total))
        print("Kmer coverage:\t" + str(kmerCoverage))
        print("Estimated genome size:\t" + str(genomeSize))

        # write to tsv
        with open(output, "w") as file:
            counter = 0
            for abundance, numKmers in distribution.items():
                if numKmers == 0:
                    counter += 1
                else:
                    counter = 0
                if counter >= 10:
                    break
                file.write(f"{abundance}\t{numKmers}\n")

    else:
        print("Please use \"vanheusden generate\" or \"vanheusden estimate\"")

if __name__ == "__main__":
    main()