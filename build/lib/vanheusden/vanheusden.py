import argparse
import random
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time
from typing import List
from collections import defaultdict

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
        output = "g" + str(g) + "_n" + str(n) + "_l" + str(l) + "_e" + str(e) + ".fastq"

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
                rng = random.randint(0, int(1 / e))
                if rng == 0:
                    read = read[:i] + GenerateError(read[i]) + read[i+1:]
            reads.append(read)

        # write to 
        with open(output, "w") as file:
            file.write("@header\n")
            for read in reads:
                file.write(read + "\n" + "+" + "\n" + "I" * len(read) + "\n")

    elif args.function == "estimate":
        start = time.time()

        k = args.kmerLength
        input = args.reads
        readLen = 999
        reads = 0

        with open(input) as file:
            for line in file:
                if line[0] in bases:
                    reads += 1
                    if len(line) - 1 < readLen:
                        readLen = len(line) - 1

        print("k:\t" + str(k))
        print("Reads: " + str(reads))
        print("Read Length:\t" + str(readLen))
        #samplingRate = int(lines / 10000)


        samplingRate = 10
        print("Sampling Rate:\t" + str(samplingRate))

        bucketWidth = 10000000
        numHashFunctions = 10

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
        countMinSketch = [[0 for j in range(bucketWidth)] for i in range(numHashFunctions)]

        # choose which strand will represent the sequence
        def ChooseStrand(kmer: str) -> str:
            kmerRev = ReverseComplement(kmer)
            if hash(kmer) < hash(kmerRev):
                return kmer
            return kmerRev

        # return a list of hash values
        def Hashes(kmer: str) -> List[int]:
            hashes = []
            for i in range(numHashFunctions):
                hashes.append(hash(str(i) + kmer) % bucketWidth)
            return hashes
        
        # update countMinSketch
        def updateCount(kmer: str):
            for i, j in enumerate(Hashes(kmer)):
                countMinSketch[i][j] += 1

        # return counts for a kmer
        def getCount(kmer: str) -> int:
            return min(countMinSketch[i][j] for i, j in enumerate(Hashes(kmer)))

        # fill out countMinSketch
        def getRatio(k: int) -> int:
            distribution = defaultdict(int)
            total = 0
            counter = 0
            with open(input, "r") as file:
                for line in file:
                    if line[0] in bases:
                        if counter % samplingRate == 0:
                            line = line.strip()
                            for i in range(len(line) - k + 1):
                                kmer = ChooseStrand(line[i : i + k])
                                updateCount(kmer)
                                total += 1
                        counter += 1

            # read values from countMinSketch
            counter = 0
            with open(input, "r") as file:
                for line in file:
                    if line[0] in bases:
                        if counter % samplingRate == 0:
                            line = line.strip()
                            for i in range(len(line) - k + 1):
                                kmer = ChooseStrand(line[i : i + k])
                                abundance = getCount(kmer)
                                distribution[abundance] += 1
                        counter += 1

            distribution = dict(sorted(distribution.items()))
            valleyY = float("inf")
            valleyX = -1
            abundance = list(distribution.keys())
            numKmers = list(distribution.values())
            errors = 0
            for i in range(len(abundance)):
                x = abundance[i]
                y = numKmers[i]
                if valleyY > y:
                    errors += y
                    valleyY = y
                    valleyX = x
                else:
                    break

            output = "k" + str(k) + "_histogram.tsv"
            with open(output, "w") as file:
                for x, y in distribution.items():
                    file.write(str(x) + "\t" + str(y) + "\n")

            print(numKmers[valleyX - 1])
            print((errors - numKmers[valleyX - 1]))

            if valleyX >= len(abundance):
                return 999

            ratio = numKmers[valleyX - 1] / (errors - numKmers[valleyX - 1])
            return ratio

        """
        k = 10
        ratio = 999
        counter = 0
        for i in range(20, readLen, 20):
            curr = getRatio(i)
            print(str(i) + "\t" + str(curr))
            if curr <= ratio:
                ratio = curr
                k = i
            else:
                counter += 1
                if counter >= 3:
                    break
            # Count-min sketch
            countMinSketch = [[0 for j in range(bucketWidth)] for i in range(numHashFunctions)]
        """

        pct1 = int(0.01 * reads)
        pct10 = int(0.10 * reads)
        pct25 = int(0.25 * reads)
        pct50 = int(0.5 * reads)
        pct75 = int(0.75 * reads)

        distribution = defaultdict(int)
        total = 0
        counter = 0
        with open(input, "r") as file:
            for line in file:
                if line[0] in bases:
                    if counter % samplingRate == 0:
                        line = line.strip()
                        for i in range(len(line) - k + 1):
                            kmer = ChooseStrand(line[i : i + k])
                            updateCount(kmer)
                            total += 1
                    counter += 1
                    if counter == pct1:
                        end = time.time()
                        print("Kmer Hashing:\t1% in " + str(round(end - start, 2)) + "s")
                    if counter == pct10:
                        end = time.time()
                        print("Kmer Hashing:\t10% in " + str(round(end - start, 2)) + "s")
                    if counter == pct25:
                        end = time.time()
                        print("Kmer Hashing:\t25% in " + str(round(end - start, 2)) + "s")
                    if counter == pct50:
                        end = time.time()
                        print("Kmer Hashing:\t50% in " + str(round(end - start, 2)) + "s")
                    if counter == pct75:
                        end = time.time()
                        print("Kmer Hashing:\t75% in " + str(round(end - start, 2)) + "s")

        end = time.time()
        print("Kmer Hashing:\t100% in " + str(round(end - start, 2)) + "s")
        print("Kmer Counting...")
        # read values from countMinSketch
        counter = 0
        with open(input, "r") as file:
            for line in file:
                if line[0] in bases:
                    if counter % samplingRate == 0:
                        line = line.strip()
                        for i in range(len(line) - k + 1):
                            kmer = ChooseStrand(line[i : i + k])
                            abundance = getCount(kmer)
                            distribution[abundance] += 1
                    counter += 1

        distribution = dict(sorted(distribution.items()))

        # write to tsv
        output = "k" + str(k) + "_histogram.tsv"
        with open(output, "w") as file:
            for x, y in distribution.items():
                file.write(str(x) + "\t" + str(y) + "\n")

        valleyY = float("inf")
        valleyX = -1
        abundance = list(distribution.keys())
        numKmers = list(distribution.values())
        errors = 0
        for i in range(len(abundance)):
            x = abundance[i]
            y = numKmers[i]
            if valleyY > y:
                errors += y
                valleyY = y
                valleyX = x
            else:
                break
        peakY = -1
        if valleyX < len(abundance):
            peakY = max(numKmers[valleyX:])
        else:
            peakY = numKmers[0]
        peakX = abundance[numKmers.index(peakY)]

        genomeSize = int(total / peakX)
        print("Total number of kmers:\t" + str(total))
        print("Kmer coverage:\t" + str(peakX))
        #print("Ratio:\t" + str(ratio))
        print("Estimated genome size:\t" + str(genomeSize))

        # generate graph
        df = pd.DataFrame()
        df["Abundance"] = distribution.keys()
        df["Number of Kmers"] = distribution.values()
        sns.set_style("whitegrid")
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        sns.barplot(data=df, x="Abundance", y="Number of Kmers", ax=ax)
        ax.set_title("k = " + str(k))

        numXlabels = 20
        skip = int(1 + len(abundance) / numXlabels)
        #print("skip: " + str(skip))
        for i, label in enumerate(ax.get_xticklabels()):
            if i % skip == 0:
                label.set_visible(True)
            else:
                label.set_visible(False)

        plt.savefig("k" + str(k) + "_histogram.pdf")

        end = time.time()
        print("Complete in " + str(round(end - start, 2)) + "s")

    else:
        print("Please use \"vanheusden generate\" or \"vanheusden estimate\"")

if __name__ == "__main__":
    main()