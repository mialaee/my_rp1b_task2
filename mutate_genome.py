# input reference FASTA -> output (1) mutated FASTA and (2) "truth" file of all SNPs/indels.

import random
import argparse

BASES = ["A", "C", "G", "T"] # allowed nucleotides for muatations

def read_fasta(path): # Read a FASTA file, returning the chromosome name and sequence as a string
    with open(path) as f:
        header = f.readline().strip()
        seq = f.read().replace("\n", "").upper()
    chrom = header[1:].split()[0]
    return chrom, seq

def write_fasta(path, chrom, seq): # Write a FASTA file, wrapping at 60 characters per line
    with open(path, "w") as f:
        f.write(f">{chrom}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")
             
def mutate_genome(seq, n_snps=300, n_indels=20, max_indel=10, seed=42):
    random.seed(seed) # fix randomness for reproducibility
    seq = list(seq) # convert sequence to list so bases can be modified
    L = len(seq)
        
    snp_positions = set() # step 1: pick random SNP positions
    while len(snp_positions) < n_snps:
        pos = random.randint(0, L-1)
        snp_positions.add(pos)
            
# step 2: pick random INDEL positions (avoid ends of genomes)
    indel_positions = set()
    while len(indel_positions) < n_indels:
        pos = random.randint(10, L-11) # avoid edge issues with deletion
        if pos not in snp_positions: # do not allow overlap 
            indel_positions.add(pos)
                
    variants =[] # list of all mutations
        
# step 3: apply SNP mutations     
    for pos in snp_positions:
        ref = seq[pos] # reference base
        alts = [b for b in BASES if b != ref] # avoid choosing the same base
        alt = random.choice(alts)
        seq[pos] = alt # apply mutation
# VCF is 1-based
        variants.append(("SNP", pos+1, ref, alt))
    

# INDELs (apply after SNPs to keep indexing simple-ish)
# Build new_seq while tracking indels

    new_seq = []
    indel_dict = {p: None for p in indel_positions}
    i = 0
    while i < len(seq):
        if i in indel_dict:
                # decide insertion vs deletion 
            indel_len = random.randint(1, max_indel)
            ins_or_del = random.choice(["INS", "DEL"])
            ref_base = seq[i]
                
            if ins_or_del == "INS":
                    # insertion after current base
                inserted = "".join(random.choice(BASES) for _ in range(indel_len))
                new_seq.append(ref_base + inserted)
                variants.append(("INS", i+1, ref_base, ref_base + inserted))
                i += 1
            else:
                    # deletion starting at this base
                deleted = "".join(seq[i:i+indel_len])
                new_seq.append(ref_base) # keep first base as anchor
                variants.append(("DEL", i+1, ref_base + deleted[1:], ref_base))
                i += indel_len
        else:
            new_seq.append(seq[i])
            i += 1
                
    mutated_seq = "".join(new_seq)
    return mutated_seq, variants

def write_vcf(path, chrom, variants): # write a VCF file describing all introduced variants
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for vtype, pos, ref, alt in variants:
            info = f"TYPE={vtype}"
            f.write(f"{chrom}\t{pos}\t/\t{ref}\t{alt}\t.\tPASS\t{info}\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref_fasta")
    ap.add_argument("out_prefix")
    ap.add_argument("--snps", type=int, default=300)
    ap.add_argument("--indels", type=int, default=20)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    chrom, seq = read_fasta(args.ref_fasta) # load reference genome
    mutated_seq, variants = mutate_genome(seq, args.snps, args.indels, seed=args.seed) # apply random mutations
    write_fasta(args.out_prefix + ".mut.fasta", chrom, mutated_seq) # write mutated genome and truth VCF
    write_vcf(args.out_prefix + ".truth.vcf", chrom, variants)

if __name__ == "__main__":
    main()