#!/usr/bin/env python

#################################################################################$$
# Copyright (c) 2011,2012, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this 
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice, 
#   this list of conditions and the following disclaimer in the documentation 
#   and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its contributors 
#   may be used to endorse or promote products derived from this software 
#   without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#################################################################################$$

import os
import sys
import argparse
import logging
import md5
import pbcore.io
import pkg_resources
import tempfile

from pbcore.util.ToolRunner import PBToolRunner

try:
    __p4revision__ = "$Revision: #17 $"
    __p4change__ = "$Change: 121572 $"
    revNum = int(__p4revision__.strip("$").split(" ")[1].strip("#"))
    changeNum = int(__p4change__.strip("$").split(":")[-1])
    __version__ = "%s-r%d-c%d" % ( pkg_resources.require("pbtools.pbdagcon")[0].version, revNum, changeNum )
except:
    __version__ = "pbtools.pbdagcon-github"

class H5FofnToFastaRunner(PBToolRunner):
    
    def __init__(self):
        desc = ['Tool for extracting file from a list of .bas.h5 file',
                'Notes: For all command-line arguments, default values are listed in [].']

        super(H5FofnToFastaRunner, self).__init__('\n'.join(desc))

        self.parser.add_argument('h5_fofn_file', metavar='input.fofn',
                                 help='input .bas.h5/.bax.h5 fofn filename')
        self.parser.add_argument('out_dir', metavar = "out_dir", default='.', 
                                 help='output directory')
        self.parser.add_argument('--min_length', type=int, dest='minlen', default=500,
                                 help='min read length [%(default)s]')
        self.parser.add_argument('--min_seed_length', type=int, dest='min_seed_len', default=4000,
                                 help='min seed read length [%(default)s]')
        self.parser.add_argument('--min_read_score', type=float, dest='minrs', default=0.75,
                                 help='min read score, valid only when used with --readType=Raw [%(default)s]')
        self.parser.add_argument('--tmp_dir', type=str, dest='tmp_dir', default=".",
                                 help='temporary directory used for intermediate fasta file')
        self.parser.add_argument('--replace', type=bool, dest='replace', default=False,
                                 help='replace existing files [%(default)s]')

    def getVersion(self):
        return __version__

    def validateArgs(self):
        pass

    def run(self):

        with open(self.args.h5_fofn_file) as f:
            for l in f:
                l = l.strip()
                assert l.endswith(".bas.h5") or l.endswith(".bax.h5")
                movie_name = ".".join( os.path.basename(l).split(".")[:-2] ) # remove ".bas.h5" or ".bax.h5" 
                movie_name_md5 = md5.md5(movie_name).hexdigest()[:8]
                out_query_name = os.path.join(self.args.out_dir,  movie_name_md5 + "_q.fa")
                out_target_name = os.path.join(self.args.out_dir,  movie_name_md5 + "_t.fa")
                if not self.args.replace and os.path.exists(out_query_name) and os.path.exists(out_target_name):
                    continue
                tmp_prefix = os.path.join(self.args.tmp_dir, "tmp_%s" % movie_name_md5)
                os.system("bash5tools.py %s --readType subreads --outType fasta --minReadScore %f --outFilePref %s" % (l, self.args.minrs, tmp_prefix))
                with open(out_query_name, "w") as o_q:
                    with open(out_target_name, "w") as o_t:
                        fq = pbcore.io.FastaReader("%s.fasta" % tmp_prefix)
                        for r in fq:
                            if len(r.sequence) < self.args.minlen:
                                continue
                            zmw_id = r.name.split("/")[1]
                            subread_start = r.name.split("/")[2].split("_")[0]
                            r_id = movie_name_md5 + "_" + zmw_id + "_" +subread_start
                            print >> o_q, ">" + r_id
                            print >> o_q, r.sequence
                            if len(r.sequence) >= self.args.min_seed_len:
                                print >> o_t, ">" + r_id
                                print >> o_t, r.sequence
                        fq.file.close()
                os.system("rm %s.fasta" % tmp_prefix)


        
if __name__ == '__main__':
    sys.exit(H5FofnToFastaRunner().start())
