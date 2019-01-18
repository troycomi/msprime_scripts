# Converting simulated bedfile of introgressed haplotyes to YMb coordinates.
# File format
#   chr         hap.start   hap.end     len    Ne  ID
# sim_0.task_1  2048116     4132969   2000000
# input from stdin and output to stdout
import sys

max_len = int(sys.argv[1])  # specify the max simulated size

for line in sys.stdin:
    tokens = line.strip().split('\t')
    start = int(tokens[1])
    end = int(tokens[2])
    tokens.append("")  # overwrite below
    for max_chrom_len in range(5000000, max_len+1, 1000000):
        if start < max_chrom_len:
            if end > max_chrom_len:
                tokens[2] = str(max_chrom_len)
            else:
                tokens[2] = str(end)

            # add the window size to the end of the haplotype list
            tokens[-1] = str(max_chrom_len//1000000)
            sys.stdout.write('\t'.join(tokens) + '\n')
