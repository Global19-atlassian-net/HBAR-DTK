#!/usr/bin/env python

from pbcore.io.FastqIO import FastqReader, FastqWriter, FastqRecord
import shlex
import sys
import subprocess
import os
import re

usage = "usage: circulization.py initial_contigs.fastq 20000 /tmp circulaized_contigs.fastq"

try:
    fastq_f = FastqReader(sys.argv[1])
    prepostfix_size = int(sys.argv[2])
    tmp_dir = sys.argv[3]
    output_fn = sys.argv[4]
except:
    print usage
    sys.exit(1)

prefix_N = re.compile("^[Nn]+")
postfix_N = re.compile("[Nn]+$")

prefix_fn = os.path.join(tmp_dir,"prefix.fa")
postfix_fn = os.path.join(tmp_dir,"postfix.fa")

with FastqWriter( open(output_fn, "w") ) as output_fh:

    for r in fastq_f:
        r_id = r.name
        r_seq = r.sequence 
        r_qv = r.quality
        m = prefix_N.search(r_seq)
        if m:
            prefix_trim = m.end()
        else:
            prefix_trim = 0
        m = postfix_N.search(r_seq)
        if m:
            postfix_trim = m.start()
        else:
            postfix_trim = len(r_seq)
        
        r_seq = r_seq[prefix_trim: postfix_trim]
        r_qv = r_qv[prefix_trim: postfix_trim]

        if len(r_seq) < prepostfix_size * 2:
            output_fh.writeRecord(FastqRecord(r_id, r_seq, qualityString = r_qv))
            continue

        prefix = r_seq[:prepostfix_size]
        postfix = r_seq[-prepostfix_size:]
        with open(prefix_fn,"w") as pref:
            print >> pref, ">"+r_id+":prefix"
            print >> pref, prefix
        
        with open(postfix_fn,"w") as post:
            print >> post, ">"+r_id+":postfix"
            print >> post, postfix
        
        out = subprocess.check_output(shlex.split("blasr %s %s -m 4 -bestn 1 -noSplitSubreads" % ( prefix_fn, postfix_fn ) ))
        if len(out) == 0:
            output_fh.writeRecord(FastqRecord(r_id, r_seq, qualityString = r_qv))
            continue
        out = out.strip().split()

        identity = float(out[3])
        q_strand = out[4]
        q_start = int(out[5])
        q_end = int(out[6])
        t_strand = out[8]
        t_start = int(out[9])
        t_end = int(out[10])

        if q_strand == "0" and t_strand == "0" and identity > 98 and q_start == 0 and t_end == prepostfix_size:
            start = int( 0.5*(q_start + q_end) ) #avoid using "/" since the semenatic has been changed in Python 3
            end = t_end - t_start - start
            r_seq = r_seq[start:-end]
            r_qv = r_qv[start:-end]
            output_fh.writeRecord(FastqRecord(r_id + ":circularization_trimmed:%d:%d" % (start,end), r_seq, qualityString = r_qv))
        else:
            output_fh.writeRecord(FastqRecord(r_id, r_seq, qualityString = r_qv))

        os.remove( os.path.join(tmp_dir,"prefix.fa") )
        os.remove( os.path.join(tmp_dir,"postfix.fa") )

