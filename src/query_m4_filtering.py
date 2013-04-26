#!/usr/bin/env python

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


# this script filter out alignment hit that are not the top n hit from a set of m4 file
#
import sys
import heapq

fofn = sys.argv[1]
num_chunk = int(sys.argv[2])
chunk_n = int(sys.argv[3])
best_n = int(sys.argv[4])
length_cutoff = int(sys.argv[5])
output_fn = sys.argv[6]

query_data = {}
with open(fofn) as fnlist:
    for fn in fnlist:
        fn = fn.strip()
        with open(fn) as m4_file:
            for l in m4_file:
                l = l.strip().split()
                q_id, t_id = l[:2]
                if q_id == t_id:
                    continue
                if hash(q_id) % num_chunk != chunk_n:
                    continue
                score = int(l[2])
                if score > -1000: continue
                t_l = int(l[11])
                if t_l < length_cutoff : continue
                query_data.setdefault(q_id, [])
                if len(query_data[q_id]) < best_n:
                    heapq.heappush( query_data[q_id], (-score, l) )
                else:
                    heapq.heappushpop( query_data[q_id], (-score, l) )

with open(output_fn,"w") as out_f:
    for q_id in query_data:
        query_data[q_id].sort(key=lambda x:-x[0])
        for score, l in query_data[q_id]:
            print >>out_f, " ".join(l)
