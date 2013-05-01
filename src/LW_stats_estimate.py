from pbcore.io import FastaIO
import numpy as np
from math import exp, log
import sys

fn = sys.argv[1]
genome_size = int(sys.argv[2])


f = FastaIO.FastaReader(fn)
seq_lengths = []
for r in f:
    seq_lengths.append(len(r.sequence))
seq_lengths = np.array(seq_lengths)


total = sum(seq_lengths)
coverage_array = []
print "cutoff\ttotal_base\ttotol/seed\tcov\tcontig_count\tcontig_len/genome_size"
for x in range(3000,10000,500):
    psum = sum(seq_lengths[seq_lengths>x])
    coverage = 0.5 * psum / genome_size # we loss 50% bases after the pre-assembly step
    contig_count = coverage * genome_size / x * exp( -coverage )
    contig_length = (exp(coverage) - 1) * x /coverage
    print "%d\t%d\t%0.2f\t%0.2f\t%0.2f\t%0.2f" % (x, psum, 1.0*total/psum, coverage,  contig_count,  contig_length/genome_size)
    coverage_array.append( [x, psum, 1.0*total/psum, coverage,  contig_count, contig_length/genome_size] )
coverage_array = np.array( coverage_array )


print 
print "recommended cutoff (total/seed > 3, LW contig # <0.25, LW contig length > 0.25x genome)"
print "cutoff\ttotal_base\ttotol/seed\tcov\tcontig_count\tcontig_len/genome_size" 

for l in coverage_array[ (coverage_array[...,2]>3) & (coverage_array[...,4]<0.25) & (coverage_array[...,5]>0.25),...]:
    print "%d\t%d\t%0.2f\t%0.2f\t%0.2f\t%0.2f" % tuple(l)
