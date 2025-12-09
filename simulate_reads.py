import argparse
import random
import math

def read_fasta(path): # reads a FASTA file and returns only the sequence
    with open(path) as f:
        f.readline()  # skip FASTA header line
        seq = f.read().replace("\n", "").upper()
    return seq

def main ():
    ap = argparse.ArgumentParser()
    # positional arguements
    # fasta: mutated genome FASTA
    # out_fasta: output FASTQ file path
    ap.add_argument("fasta") 
    ap.add_argument("out_fastq")
    
    # optional parameters:
    # depth: sequencing coverage (30x)
    # readlen: read length (100bp)
    # seed: makes simulation reproducible
    ap.add_argument("--depth", type=float, default=30.0)
    ap.add_argument("--readlen", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    
    random.seed(args.seed) # fix randomness so the same FASTQ is produced every time
    seq = read_fasta(args.fasta) # load mutated genome sequence
    L = len(seq)
    n_reads = math.ceil(args.depth * L / args.readlen) # calculate number of reads required to achieve requested depth:
    # coverage = (total bases sequenced) / genome length
    # total bases = n_reads * readlen
    # so: n_reads = depth * genome_length / read_length
    
    with open(args.out_fastq, "w") as out:
        for i in range(n_reads):
            start = random.randint(0, L - args.readlen) # randomly choose a start position for a read
            read_seq = seq[start:start+args.readlen] # extract a read of length readlen
            # creates dummy high quality score
            qual = "I" * args.readlen
            out.write(f"@read_{i}_pos_{start+1}\n") # write FASTQ entry
            out.write(read_seq + "\n+\n" + qual + "\n")

if __name__ == "__main__":
    main()