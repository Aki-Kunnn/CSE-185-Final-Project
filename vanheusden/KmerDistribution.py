# vanheusden/kmer_analysis.py

import argparse
import hashlib
from typing import List

def main():
    """
    Main function to parse arguments, read input file, compute k-mer distribution, 
    and write results to output file.
    """
    parser = argparse.ArgumentParser(prog="Vanheusden", 
                                     description="kmer distribution of reads\ncol 1: count, col 2: number of kmers that appeared count times")
    parser.add_argument("-k", "--kmerLength", required=True, type=int, help="Length of the k-mers")
    parser.add_argument("reads", help="Input file containing DNA sequencing reads")
    args = parser.parse_args()

    k = args.kmerLength
    output = args.reads.split(".")[0] + ".histo.tsv"
    input = args.reads
    bases = ["A", "T", "C", "G"]

    distribution = {i: 0 for i in range(1, 500)}
    width = 1000
    initMin = 9999
    numFunctions = 10

    def ReverseComplement(pattern: str) -> str:
        """
        Compute the reverse complement of a DNA sequence.
        
        Args:
            pattern (str): DNA sequence.
        
        Returns:
            str: Reverse complement of the DNA sequence.
        """
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

    countMinSketch = [[0 for j in range(1000)] for i in range(numFunctions)]

    def ChooseStrand(kmer: str) -> str:
        """
        Choose the lexicographically smaller k-mer between a k-mer and its reverse complement.
        
        Args:
            kmer (str): Original k-mer.
        
        Returns:
            str: Lexicographically smaller k-mer.
        """
        kmerRev = ReverseComplement(kmer)
        if hash(kmer) < hash(kmerRev):
            return kmer
        return kmerRev

    def Hashes(kmer: str) -> List[int]:
        """
        Generate a list of hash values for a k-mer using multiple hash functions.
        
        Args:
            kmer (str): k-mer to be hashed.
        
        Returns:
            List[int]: List of hash values.
        """
        hashes = []
        functions = [
            hashlib.sha1(), hashlib.sha224(), hashlib.sha256(), 
            hashlib.sha384(), hashlib.sha512(), hashlib.sha3_224(), 
            hashlib.sha3_256(), hashlib.sha3_384(), hashlib.sha3_512(), 
            hashlib.md5()
        ]
        code = kmer.encode("utf-8")
        for i in range(numFunctions):
            functions[i].update(code)
            hashes.append(int(functions[i].hexdigest(), 16) % width)
        return hashes

    def updateCount(kmer: str):
        """
        Update the Count-Min Sketch with the given k-mer.
        
        Args:
            kmer (str): k-mer to be added to the Count-Min Sketch.
        """
        for i, j in enumerate(Hashes(kmer)):
            countMinSketch[i][j] += 1

    def getCount(kmer: str) -> int:
        """
        Retrieve the count for a k-mer from the Count-Min Sketch.
        
        Args:
            kmer (str): k-mer whose count is to be retrieved.
        
        Returns:
            int: Count of the k-mer.
        """
        return min(countMinSketch[i][j] for i, j in enumerate(Hashes(kmer)))

    with open(input, "r") as file:
        for line in file:
            if line[0] in bases:
                line = line.strip()
                for i in range(len(line) - k + 1):
                    kmer = ChooseStrand(line[i : i + k])
                    updateCount(kmer)

    with open(input, "r") as file:
        for line in file:
            if line[0] in bases:
                line = line.strip()
                for i in range(len(line) - k + 1):
                    kmer = ChooseStrand(line[i : i + k])
                    distribution[getCount(kmer)] += 1

    with open(output, "w") as file:
        counter = 0
        for count, num_kmers in distribution.items():
            if num_kmers == 0:
                counter += 1
            if counter >= 10:
                break
            file.write(f"{count}\t{num_kmers}\n")
