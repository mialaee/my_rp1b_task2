import argparse  # allows the script to accept command-line arguments
import gzip  # allows reading .gz compressed VCF files

def load_vcf(path):
    variants = set()   

    # auto-detect gzip vs plain text
    if path.endswith(".gz"):
        opener = gzip.open    # use gzip.open for compressed files
        mode = "rt"      # "rt" = read text mode
    else:
        opener = open    # use normal open() for uncompressed files
        mode = "r"
   
    with opener(path, mode) as f:  # open file with the appropriate method
        for line in f:
            if line.startswith("#"):   # skip header lines that start with '#'
                continue
            chrom, pos, vid, ref, alt, *_ = line.strip().split("\t") # split VCF line into the standard files
            for a in alt.split(","):
                variants.add((chrom, pos, ref, a))

    return variants  # return the set of variants found in VCF


def main(): # create command-line argument parser
    ap = argparse.ArgumentParser(
        description="Compare truth and test VCF files and compute TP/FP/FN, precision, recall."
    )
    ap.add_argument("truth_vcf", help="Ground-truth VCF (from mutate_genome.py)")
    ap.add_argument("test_vcf", help="VCF to evaluate (e.g. bcftools.vcf)")
    args = ap.parse_args()

    truth = load_vcf(args.truth_vcf) #load truth and test variants into sets
    test = load_vcf(args.test_vcf)

    # Set operations
    tp = len(truth & test) # compute True Positives: variants present in both
    fp = len(test - truth) # compute False Positives: variants in test but not in truth
    fn = len(truth - test) # compute False Negatives: variants in TRUTH  in missing from test

    print(f"Loaded {len(truth)} truth variants") # print summary
    print(f"Loaded {len(test)} test variants")
    print(f"TP = {tp}")
    print(f"FP = {fp}")
    print(f"FN = {fn}")

    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    rec  = tp / (tp + fn) if (tp + fn) > 0 else 0.0

    print(f"Precision = {prec:.3f}")
    print(f"Recall    = {rec:.3f}")

# if the script is run directly (not imported), execute main ()
if __name__ == "__main__":
    main()
