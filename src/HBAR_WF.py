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


import sys
import os
import uuid

from contextlib import nested
import logging
import zlib
import glob
import pbcore.io
import pkg_resources


from pypeflow.common import * 
from pypeflow.task import PypeTaskBase
from pypeflow.task import PypeTask, PypeShellTask
from pypeflow.controller import PypeWorkflow
from pypeflow.data import PypeLocalFile, makePypeLocalFile, fn

try:
    __p4revision__ = "$Revision$"
    __p4change__ = "$Change$"
    revNum = int(__p4revision__.strip("$").split(" ")[1].strip("#"))
    changeNum = int(__p4change__.strip("$").split(":")[-1])
    __version__ = "%s-r%d-c%d" % ( pkg_resources.require("pbtools.hbar-dtk")[0].version, revNum, changeNum )
except:
    __version__ = "pbtools.hbar-dtk-github"

def prepare_data(self):

    config = self.config
    RQ_threshold = config["RQ_threshold"]
    
    with nested(open(fn(self.input_fofn)), 
                open(fn(self.normalized_fasta),'w'),
                open("id_map.dat","w")) as (f, nfasta, id_map_f): 
        i = 0
        for l in f:
            l = l.strip()
            file_ext = os.path.basename(l).split(".")[-1]
            if file_ext == "h5":
                os.system("bash5tools.py --readType subreads --outType fasta --minReadScore %f %s --outFilePref tmp" % (RQ_threshold, l) )
            elif file_ext in ["fa", "fasta"]:
                os.system("cp %s tmp.fasta" % l )
            elif file_ext in ["fastq"]:
                fq = pbcore.io.FastqReader(l)
                with open("tmp.fasta","w") as tmpf:
                    for r in fq:
                        print >> tmpf, ">" + r.name
                        print >> tmpf, r.sequence
                fq.file.close()

            tmpfa = pbcore.io.FastaReader("tmp.fasta")
            for r in tmpfa:
                r_id = "%s_%07d" %  (hex(zlib.adler32(r.name+r.sequence) & 0xffffffff), i)
                print >> id_map_f, "%s %s" % (r.name, r_id)
                nfasta.write( ">%s\n" % r_id )
                nfasta.write( "%s\n" % r.sequence.upper() )
                i+=1               
            tmpfa.file.close()
            os.system("rm tmp.fasta")

def prepare_seed_reads(self):

    config = self.config
    length_cutoff = config["length_cutoff"]

    f = pbcore.io.FastaReader(fn(self.normalized_fasta))
    with open(fn(self.seed_fasta),'w') as sfasta:
        for r in f:
            if len(r.sequence) > length_cutoff:
                sfasta.write( ">%s\n" % r.name )
                sfasta.write( "%s\n" % r.sequence.upper() )

def dist_map(self):

    config = self.config
    dist_map_num_chunk = config["dist_map_num_chunk"]
    directory_for_dist_map = config["directory_for_dist_map"]
    sge_option_ck = config["sge_option_ck"]
    sge_option_dm = config["sge_option_dm"]
    install_prefix = config["install_prefix"]
    blasr_opt = config["blasr_opt"]

    #set_up_script = "fastasplit %s %s/ -c %d" % (fn(self.seed_fasta), directory_for_dist_map, dist_map_num_chunk)
    #os.system(set_up_script)
    fasta_file = pbcore.io.FastaReader(fn(self.seed_fasta))
    out_files = []
    for i in range(dist_map_num_chunk):
        out_files.append( open( "%s/%s_chunk_%07d" % (directory_for_dist_map, os.path.basename(fn(self.seed_fasta)), i), "w"))
    for s in fasta_file:
        g = hash(s.name) % dist_map_num_chunk
        out_file = out_files[g]
        out_file.write(">%s\n" % s.name)
        out_file.write("%s\n" % s.sequence)
    for i in range(dist_map_num_chunk):
        out_files[i].close()
    fasta_file.file.close()
    

    align_script_template = """\
. {install_prefix}/bin/activate 
cd %s/%s
blasr {blasr_opt} -m 4 -out m4_%s.dat %s %s 
""".format(install_prefix = install_prefix, blasr_opt=blasr_opt)

    job_name = "dist_map_"+str(uuid.uuid4())
    i = 0
    for chunk_name in glob.glob("%s/%s_chunk_*" % ( directory_for_dist_map, os.path.basename(fn(self.seed_fasta))) ):
        script = align_script_template % (os.getcwd(), directory_for_dist_map, 
                                          os.path.basename(chunk_name), 
                                          fn(self.normalized_fasta), 
                                          os.path.basename(chunk_name))
        with open("scripts/dist_map_%02d.sh" % i,"w") as f:
            print >>f, script
        os.system("qsub -N {jn} {sge_option_dm} -o {cwd}/sge_log -j y\
                -S /bin/bash scripts/dist_map_{jid:02d}.sh".format(jn=job_name+"_%02d" % i, 
                                                                   cwd=os.getcwd(), 
                                                                   sge_option_dm = sge_option_dm, 
                                                                   jid=i))
        i += 1
    
    with open("scripts/mapping_done.sh","w") as f:
        print >>f, "echo done > %s" % fn(self.m4_data_done)
    os.system("""qsub -sync y {sge_option_ck} -hold_jid "{jn}*" -o {cwd}/sge_log -j y\
               -S /bin/bash scripts/mapping_done.sh""".format(jn=job_name, cwd=os.getcwd(), sge_option_ck=sge_option_ck))

def m4_filtering(self):
    config = self.config
    dist_map_num_chunk = config["dist_map_num_chunk"]
    directory_for_dist_map = config["directory_for_dist_map"]
    sge_option_ck = config["sge_option_ck"]
    sge_option_mf = config["sge_option_mf"]
    install_prefix = config["install_prefix"]
    bestn = config["bestn"]
    length_cutoff_pr = config["length_cutoff_pr"]

    os.system('find %s/%s -name "m4*.dat" > %s/%s/m4.fofn' % (os.getcwd(), directory_for_dist_map,
                                                              os.getcwd(), directory_for_dist_map) )

    #set_up_script = "cp filterM4Query.py %s/" % directory_for_dist_map
    #os.system(set_up_script)

    m4filtering_script_template = """\
. {install_prefix}/bin/activate 
cd %s/%s
filterM4Query.py m4.fofn {n_chunk} %d {bestn} {length_cutoff_pr} qf_m4_%02d.dat 
""".format(install_prefix = install_prefix, n_chunk = dist_map_num_chunk, bestn = bestn, length_cutoff_pr = length_cutoff_pr)
    job_name = "m4filtering_"+str(uuid.uuid4())
    for i in range(dist_map_num_chunk):
        script = m4filtering_script_template % (os.getcwd(), directory_for_dist_map, i, i) 
        with open("scripts/m4filtering_%02d.sh" % i,"w") as f:
            print >>f, script
        os.system("qsub -N {jn} {sge_option_mf}  -o {cwd}/sge_log -j y\
                -S /bin/bash scripts/m4filtering_{jid:02d}.sh".format(jn=job_name+"_%02d" % i, 
                                                                   cwd=os.getcwd(), 
                                                                   sge_option_mf = sge_option_mf, 
                                                                   jid=i))
    with open("scripts/m4filtering_done.sh","w") as f:
        print >>f, "echo done > %s" % fn(self.m4filtering_done)
    os.system("""qsub -sync y {sge_option_ck} -hold_jid "{jn}*" -o {cwd}/sge_log -j y\
            -S /bin/bash scripts/m4filtering_done.sh""".format(jn=job_name, cwd=os.getcwd(), sge_option_ck=sge_option_mf))

def get_preassembled_reads(self):

    config = self.config
    directory_for_dist_map = config["directory_for_dist_map"]
    sge_option_ck = config["sge_option_ck"]
    sge_option_pa = config["sge_option_pa"]
    bestn = config["bestn"]
    tmpdir = config["tmpdir"]
    install_prefix = config["install_prefix"]
    num_chunk = config["preassembly_num_chunk"]
    min_cov = config["min_cov"]
    max_cov = config["max_cov"]
    trim_align = config["trim_align"]
    trim_plr = config["trim_plr"]
    q_nproc = config["q_nproc"]

    #set_up_script = "cp generate_preassemble_reads.py %s/" % directory_for_dist_map
    #os.system(set_up_script)
    SGE_script_template = """. %s/bin/activate
cd %s/%s
echo start: `date` > %01d"_job.log"
hostname >> %01d"_job.log"
ls -l m4*.dat >> %01d"_job.log"
%s >> %01d"_job.log"
echo end: `date` >> %01d"_job.log"
"""

    job_name = "preassembly_"+str(uuid.uuid4())
    for j_id in range(0, num_chunk):
        #TODO: use real template lib

        g_plr_str = "generate_preassemble_reads.py %01d %s %s %d %s %d %d %d %d %d %d" % ( j_id, fn(self.normalized_fasta), fn(self.seed_fasta), 
                                                                                      bestn, tmpdir, num_chunk, min_cov, max_cov, trim_align, trim_plr, q_nproc )

        script = SGE_script_template % ( install_prefix, os.getcwd(), directory_for_dist_map,
                                         j_id, j_id, j_id, g_plr_str, j_id, j_id )


        with open("scripts/preassembly_%02d.sh" % j_id,"w") as f:
            print >>f, script

        os.system("qsub -N {jn} {sge_option_pa} -o {cwd}/sge_log -j y\
                -S /bin/bash scripts/preassembly_{jid:02d}.sh".format(jn=job_name+"_%02d" % j_id, 
                                                                      cwd=os.getcwd(), 
                                                                      sge_option_pa = sge_option_pa, 
                                                                      jid=j_id))

    with open("scripts/preassembly_done.sh","w") as f:
        print >>f, "echo done > %s" % fn(self.preassembly_done)
    os.system("""qsub -sync y  {sge_option_ck} -hold_jid "{jn}*" -o {cwd}/sge_log -j y -S /bin/bash scripts/preassembly_done.sh""".format(jn=job_name, cwd=os.getcwd(), sge_option_ck = sge_option_ck))

def run_CA(self):

    config = self.config
    install_prefix = config["install_prefix"]
    sge_option_ca = config["sge_option_ca"]
    sge_option_ck = config["sge_option_ck"]
    directory_for_dist_map = config["directory_for_dist_map"]

    set_up_script = "cat %s/pre_assembled_reads_*.fa > CA/pre_assembled_reads.fa" % (directory_for_dist_map, )
    os.system(set_up_script)

    f = pbcore.io.FastaReader("CA/pre_assembled_reads.fa")
    with open("CA/pre_assembled_reads.fastq","w") as fq:
        for r in f:
            print >>fq, "@"+r.name
            print >>fq, r.sequence
            print >>fq, "+"
            print >>fq, "".join( [chr(33+24)] * len(r.sequence) )

    job_name = "CA_"+str(uuid.uuid4())
    os.system("cd CA/;fastqToCA -technology sanger -type sanger -libraryname preassembled_long -reads pre_assembled_reads.fastq  > pre_assembled_reads.frg")
    #os.system("cp asm.spec CA/")
    script = """. {install_prefix}/bin/activate
cd {cwd}/CA
runCA -d . -p asm -s {install_prefix}/etc/asm.spec  pre_assembled_reads.frg """.format(install_prefix=install_prefix, cwd=os.getcwd())
    with open("scripts/runCA.sh",'w') as f:
        print >>f, script
    os.system("qsub {sge_option_ca} -N {jn} -o {cwd}/sge_log -j y -S /bin/bash scripts/runCA.sh".format(jn=job_name, cwd=os.getcwd(), sge_option_ca=sge_option_ca))
    with open("scripts/CA_done.sh","w") as f:
        print >>f, "echo done > %s" % fn(self.CA_done)
    os.system("""qsub -sync y {sge_option_ck} -hold_jid "{jn}*" -o {cwd}/sge_log -j y -S /bin/bash scripts/CA_done.sh""".format(jn=job_name, cwd=os.getcwd(), sge_option_ck = sge_option_ck))

def quiver_reseq(self):
    config = self.config
    sge_option_ck = config["sge_option_ck"]
    sge_option_qv = config["sge_option_qv"]
    big_tmpdir = config["big_tmpdir"]

    try:
        os.makedirs("quiver_reseq")
    except:
        pass

    SEYMOUR_HOME = config["SEYMOUR_HOME"]
    if SEYMOUR_HOME == None:
        print "SEYMOUR_HOME not set, bypass quiver consensus step"
        return 0

    job_name = "QuiverReq_"+str(uuid.uuid4())
    quiver_script = """#!/bin/bash
export SEYMOUR_HOME=%s
. $SEYMOUR_HOME/etc/setup.sh 
cd %s/quiver_reseq
cp ../CA/9-terminator/asm.ctg.fasta .
referenceUploader -c -p $PWD -n assembly -f asm.ctg.fasta --skipIndexUpdate
compareSequences.py --info --useGuidedAlign --algorithm=blasr --nproc=24 --noXML --h5mode=w --h5fn=out.cmp.h5 --minAccuracy=0.70 --minLength=200 -x -nCandidates 50 -x -minMatch 12 -x -bestn 1 -x -minPctIdentity 70.0 %s assembly/
loadPulses %s out.cmp.h5 -metrics DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag -byread
cmph5tools.py sort out.cmp.h5 --tmp %s
variantCaller.py --algorithm quiver -j 16 --referenceFilename assembly/sequence/assembly.fasta  --parameters best -o output.gff  -o output.fasta -o output.fastq -q 0  -X 80 -x 5 --mapQvThreshold 0 out.cmp.h5
""" % (SYMOURE_HOME, os.getcwd(), fn(self.input_fofn), fn(self.input_fofn), big_tmpdir)
    with open("scripts/quiver_reseq.sh", "w") as f:
        print >>f, quiver_script
    os.system( """qsub -sync y  {sge_option_qv} -N {jn} -o {cwd}/sge_log -j y -S /bin/bash scripts/quiver_reseq.sh """.format(jn=job_name, cwd=os.getcwd(), sge_option_qv = sge_option_qv) )
    with open("scripts/quiver_done.sh","w") as f:
        print >>f, "echo done > %s" % fn(self.Quiver_done)
    os.system("bash scripts/quiver_done.sh")



def run_HGAP(config):

    global prepare_data
    global prepare_seed_reads
    global dist_map
    global generate_preassemble_reads
    global run_CA
    global quiver_reseq
    
    directory_for_dist_map = "dist_map"

    config["install_prefix"] = sys.prefix
    config["directory_for_dist_map"] = directory_for_dist_map

    input_fofn_fn = config["input_fofn_fn"]
    #tmpdir = config["tmpdir"]

    #prepration the distribute mapping directory 
    #try:
        #os.makedirs("%s/ec_data" % directory_for_dist_map)
        #os.makedirs("/%s/ec_data" % tmpdir)
    #except:
        #pass

    try:
        os.makedirs("%s" % directory_for_dist_map)
    except:
        pass
    try:
        os.makedirs("scripts")
    except:
        pass
    try:
        os.makedirs("CA")
    except:
        pass
    try:
        os.makedirs("sge_log")
    except:
        pass

    input_fofn = makePypeLocalFile(input_fofn_fn)
    normalized_fasta = makePypeLocalFile("all_norm.fa")
    seed_fasta = makePypeLocalFile("seeds.fa")

    wf = PypeWorkflow()
    prepare_data_task = PypeTask(inputDataObjs={"input_fofn":input_fofn},
                                 outputDataObjs={"normalized_fasta":normalized_fasta},
                                 config = config ) (prepare_data)

    prepare_seed_reads_task = PypeTask(inputDataObjs = {"normalized_fasta":normalized_fasta},
                                       outputDataObjs = {"seed_fasta":seed_fasta},
                                       config = config)(prepare_seed_reads)
                
    wf.addTasks([prepare_data_task, prepare_seed_reads_task])


    m4_data_done = makePypeLocalFile("%s/m4_data_done" % directory_for_dist_map)
    dist_map_task = PypeTask(inputDataObjs = {"normalized_fasta":normalized_fasta, "seed_fasta":seed_fasta},
                    outputDataObjs = {"m4_data_done":m4_data_done},
                    config = config) (dist_map)

    m4filtering_done = makePypeLocalFile("%s/m4filtering_done" % directory_for_dist_map)
    m4filtering_task = PypeTask(inputDataObjs = {"m4_data_done":m4_data_done},
                       outputDataObjs = {"m4filtering_done":m4filtering_done},
                       config = config) (m4_filtering)

    preassembly_done = makePypeLocalFile("%s/preassembly_done" % directory_for_dist_map)
    get_preassembled_reads_task = PypeTask( inputDataObjs = {"normalized_fasta" : normalized_fasta, 
                                                            "seed_fasta" : seed_fasta, 
                                                            "m4filtering_done" : m4filtering_done},
                                            outputDataObjs = {"preassembly_done" : preassembly_done},
                                            config = config ) (get_preassembled_reads)

    wf.addTasks([dist_map_task, m4filtering_task, get_preassembled_reads_task])


    CA_done = makePypeLocalFile("CA_done")
    run_CA_task = PypeTask( inputDataObjs = {"preassembly_done" : preassembly_done},
                  outputDataObjs = {"CA_done": CA_done},
                  config = config )(run_CA)

    wf.addTasks([run_CA_task])
        
    Quiver_done = makePypeLocalFile("Quiver_done")
    quiver_reseq_task = PypeTask( inputDataObjs = {"CA_done": CA_done, "input_fofn":input_fofn},
                                  outputDataObjs = {"Quiver_done": Quiver_done},
                                  config = config) ( quiver_reseq )

    wf.addTasks([quiver_reseq_task])
    if config["target"] == "all":
        wf.refreshTargets([Quiver_done])
    elif config["target"] == "draft_assembly":
        wf.refreshTargets([CA_done])
    elif config["target"] == "pre_assembly":
        wf.refreshTargets([preassembly_done])



if __name__ == '__main__':
    import ConfigParser
    config = ConfigParser.RawConfigParser()
    if len(sys.argv) < 2:
        print "you need to specify a configuration file"
        print "example: HGAP.py HGAP_run.cfg"
        sys.exit(1)
    config.read(sys.argv[1])

    length_cutoff = config.getint('General', 'length_cutoff')
    input_fofn_fn = config.get('General', 'input_fofn')
    
    length_cutoff_pr = config.getint('General', 'length_cutoff_pr')
    
    RQ_threshold = 0.80
    if config.has_option('General', 'RQ_threshold'):
        RQ_threshold = config.getfloat('General', 'RQ_threshold')

    bestn = 12
    if config.has_option('General', 'bestn'):
        bestn = config.getint('General', 'bestn')

    blasr_opt = """-nCandidates 50 -minMatch 12 -maxLCPLength 15 -bestn 4 -minPctIdentity 70.0 -maxScore -1000 -m 4 -nproc 4"""
    if config.has_option('General', 'blasr_opt'):
        blasr_opt = config.get('General', 'blasr_opt')

    tmpdir = "/tmp"
    if config.has_option('General', 'tmpdir'):
        tmpdir = config.get('General', 'tmpdir')

    big_tmpdir = "/tmp"
    if config.has_option('General', 'big_tmpdir'):
        big_tmpdir = config.get('General', 'big_tmpdir')
    
    if not config.has_option('General', 'SEYMOUR_HOME'):
        print """ SEYMOUR_HOME not found in the configuration file, quiver can not be run """
        SEYMOUR_HOME = None
    else:
        SEYMOUR_HOME = config.get('General', 'SEYMOUR_HOME')

    if config.has_option('General', 'target'):
        target = config.get('General', 'target')
        if target not in ["all", "draft_assembly", "pre_assembly"]:
            print """ Target has to be "all", "draft_assembly" or "pre_assembly". You have an unknown target %s in the configuration file.  """ % target
            sys.exit(1)
    else:
        print """ No target specified, assuming a "pre-assembly" only target """
        target = "pre_assembly"

    preassembly_num_chunk = 12
    if config.has_option('General', 'preassembly_num_chunk'):
        preassembly_num_chunk = config.getint('General', 'preassembly_num_chunk')

    dist_map_num_chunk = 12
    if config.has_option('General', 'dist_map_num_chunk'):
        dist_map_num_chunk = config.getint('General', 'dist_map_num_chunk')
        #if dist_map_num_chunk < bestn / 2.0:
        #    dist_map_num_chunk = int((bestn + 1)/ 2)

    min_cov = 4
    if config.has_option('General', 'min_cov'):
        min_cov = config.getint('General', 'min_cov')

    max_cov = 60 
    if config.has_option('General', 'max_cov'):
        max_cov = config.getint('General', 'max_cov')

    trim_align = 50
    if config.has_option('General', 'trim_align'):
        trim_align = config.getint('General', 'trim_align')

    trim_plr = 50
    if config.has_option('General', 'trim_plr'):
        trim_plr = config.getint('General', 'trim_plr')

    q_nproc = 4
    if config.has_option('General', 'q_proc'):
        q_nproc = config.getint('General', 'q_proc')

    SYMOURE_HOME = config.get("General", "SEYMOUR_HOME")
    hgap_config = {"input_fofn_fn" : input_fofn_fn,
                   "length_cutoff" : length_cutoff,
                   "length_cutoff_pr" : length_cutoff_pr,
                   "sge_option_ck": config.get('General', 'sge_option_ck'),
                   "sge_option_dm": config.get('General', 'sge_option_dm'),
                   "sge_option_mf": config.get('General', 'sge_option_mf'),
                   "sge_option_pa": config.get('General', 'sge_option_pa'),
                   "sge_option_ca": config.get('General', 'sge_option_ca'),
                   "sge_option_qv": config.get('General', 'sge_option_qv'),
                   "bestn" : bestn,
                   "blasr_opt" : blasr_opt,
                   "RQ_threshold" : RQ_threshold,
                   "tmpdir" : tmpdir,
                   "SEYMOUR_HOME" : SEYMOUR_HOME,
                   "target" : target,
                   "dist_map_num_chunk": dist_map_num_chunk,
                   "preassembly_num_chunk": preassembly_num_chunk,
                   "big_tmpdir": big_tmpdir,
                   "min_cov": min_cov,
                   "max_cov": max_cov,
                   "trim_align": trim_align,
                   "trim_plr": trim_plr,
                   "q_nproc": q_nproc}


    run_HGAP(hgap_config)
