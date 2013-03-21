#################################################################################$$
# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#################################################################################$$

import sys
import glob
import pkg_resources
import uuid

from collections import Counter
from multiprocessing import Pool
from pbtools.pbdagcon.q_sense import *
import os

__p4revision__ = "$Revision$"
__p4change__ = "$Change$"
revNum = int(__p4revision__.strip("$").split(" ")[1].strip("#"))
changeNum = int(__p4change__.strip("$").split(":")[-1])
__version__ = "%s-r%d-c%d" % ( pkg_resources.require("pbtools.pbhgap")[0].version, revNum, changeNum )

group_id = int(sys.argv[1])
all_norm_fasta = sys.argv[2]
seed_fasta = sys.argv[3]
bestn = int(sys.argv[4])
tmpdir = sys.argv[5]
num_chunk = int(sys.argv[6])
min_cov = int(sys.argv[7])
max_cov = int(sys.argv[8])
trim_align = int(sys.argv[9])
trim_plr = int(sys.argv[10])
nproc = int(sys.argv[11])

pa_dir = "/%s/pa_%s" % (tmpdir, str(uuid.uuid4()))
os.system("mkdir -p %s" % pa_dir)

def get_ec_reads(input_data):
    ref_data, read_data, cov = input_data
    ref_id = ref_data[0]
    ref_name = "%s/ref_%s.fa" % (pa_dir, ref_id)

    if len(ref_data[1]) < 500:
        return ref_id, "", 0

    with open(ref_name, "w") as f:
        print >>f, ">"+ref_data[0]
        print >>f, ref_data[1]

    read_name = "%s/read_%s.fa" % (pa_dir, ref_id)
    with open(read_name, "w") as f:
        for d in read_data:
            print >>f, ">"+d[0]
            print >>f, d[1]
    
    #os.system("cd %s; q-sense.py r read_%s.fa ref_%s.fa --min_cov %d --max_cov %d --cname %s --max_n_reads 500  -o %s_consensus --nproc %d" % (pa_dir, ref_id, ref_id, min_cov, max_cov, ref_id, ref_id, nproc) )  


    input_fasta_name = os.path.join(pa_dir, "read_%s.fa" % ref_id)
    ref_fasta_name = os.path.join(pa_dir, "ref_%s.fa" % ref_id)
    out_file_name = "%s_consensus" % ref_id
    out_dir_name = pa_dir
    prefix = out_file_name.split(".")
    if len(prefix) > 1:
        prefix = ".".join(prefix[:-1])
    else:
        prefix = ".".join(prefix)
    full_prefix = os.path.join(out_dir_name, prefix)
    hp_corr = False
    mark_lower_case = False
    max_n_reads = 500
    entropy_th = 0.65
    niter = 2
    dump_dag_info = False
    consensus_seq_name = ref_id
    s = generate_consensus(input_fasta_name, ref_fasta_name, full_prefix, consensus_seq_name,
                           hp_corr, 
                           niter, 
                           max_n_reads, 
                           entropy_th,
                           dump_dag_info,
                           min_cov,
                           max_cov,
                           mark_lower_case,
                           nproc)

    try:
        os.system("rm %s/ref_%s.fa %s/read_%s.fa %s/%s_consensus*" % (pa_dir, ref_id, pa_dir, ref_id, pa_dir, ref_id))
    except IOError:
        return ref_id, "", 0

 
    if s != None:
        return ref_id, s, cov
    else:
        return ref_id, "", cov



exe_pool = Pool(4)


"""0x239fb832/0_590 0x722a1e26 -1843 81.6327 0 62 590 590 0 6417 6974 9822 254 11407 -74.5375 -67.9 1"""
query_to_target = {}


for fn in glob.glob("qf_*.dat"):
    with open(fn) as f:
        for l in f:
            d = l.strip().split()
            id1, id2 = d[:2]
            id1 = id1.split("/")[0]
            if id1 == id2:
                continue

            if hash(id2) % num_chunk != group_id:
                continue
            
            #if id1 in query_to_target: continue #not the best hits
            if int(d[2]) > -1000: continue
            query_to_target.setdefault(id1, [])
            query_to_target[id1].append( (int(d[2]),l) )

target_to_query = {}

for id1 in query_to_target:

    query_to_target[id1].sort()
    rank = 0
    for s, ll in query_to_target[id1][:bestn]:
        l = ll.strip()
        d = l.split()
        id1, id2 = d[:2]
        target_to_query.setdefault(id2,[])
        target_to_query[id2].append( ( (rank, int(d[2])), l ) )
        rank += 1


from pbcore.io import FastaIO

f_s = FastaIO.FastaReader(all_norm_fasta)
query_data = {}
for s in f_s:
    id1 = s.name
    if id1 not in query_to_target:
        continue
    query_data[id1]=s.sequence


f_l = FastaIO.FastaReader(seed_fasta)
target_data = {}
for s in f_l:
    id2 = s.name
    if hash(id2) % num_chunk != group_id:
        continue
    target_data[id2]=s.sequence


ec_data = []
base_count = Counter()
for id2 in target_to_query:
    if len(target_to_query[id2])<10:
        continue
    if id2 not in target_data:
        continue
    #if id2 != "0x885e665a": continue
    ref_data = (id2, target_data[id2]) 
    ref_len = len(target_data[id2])
    base_count.clear()
    base_count.update( target_data[id2] )
    if 1.0*base_count.most_common(1)[0][1]/ref_len > 0.8:  # don't do preassmbly if a read is of >80% of the same base
        continue
    read_data = []
    
    query_alignment = target_to_query[id2]
    query_alignment.sort() # get better alignment
    total_bases = 0
    max_cov_bases = max_cov * ref_len * 1.2
    min_cov_bases = min_cov * ref_len * 3
    
    for rank_score, l in query_alignment:
        rank, score = rank_score
        #if rank >= 1 and total_bases > min_cov_bases:
        #    break
        l = l.split()
        id1 = l[0].split("/")[0]
        t_s = int(l[5]) + trim_align
        t_e = int(l[6]) - trim_align
        if t_e - t_s < 400:
            continue
        #if float(l[2]) / (t_e - t_s) > -2.5: #don't use low identity alignment
        #    continue
        total_bases += t_e - t_s
        if total_bases > max_cov_bases:
            break
        q_seq = query_data[id1][t_s:t_e]
        read_data.append( ( "%s/0/%d_%d" % (id1, t_s, t_e), q_seq) )

    if len(read_data) > 5:
        ec_data.append( (ref_data, read_data, 1.0 * total_bases/ref_len) )
print len(ec_data)
#ec_data = ec_data[:20]
with open("pre_assembled_reads_%02d.fa" % group_id, "w") as f:
    #for x in ec_data:
    i = 1
    for res in exe_pool.imap(get_ec_reads, ec_data):
        #res = get_ec_reads(x)
        if len(res[1]) > 500:
            print >>f, ">"+res[0]
            print >>f, res[1][trim_plr:-trim_plr]
            #f.flush()    
            if i % 100 == 0:
                os.fsync(f.fileno()) 
        print "+", i, res[0], len(res[1]), res[2]
        if i % 100 == 0:
            sys.stdout.flush()
        i += 1


