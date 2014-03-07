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


from pypeflow.common import * 
from pypeflow.data import PypeLocalFile, makePypeLocalFile, fn
from pypeflow.task import PypeTask, PypeThreadTaskBase, PypeTaskBase
from pypeflow.controller import PypeWorkflow, PypeThreadWorkflow
import glob
import sys
import os
import re
import time
import logging
import uuid

log = 0
if log:
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

def run_script(job_data, job_type = "SGE" ):
    if job_type == "SGE":
        job_name = job_data["job_name"]
        cwd = job_data["cwd"]
        sge_option = job_data["sge_option"]
        script_fn = job_data["script_fn"]
        sge_cmd="qsub -N {job_name} {sge_option} -o {cwd}/sge_log -j y\
                 -S /bin/bash {script}".format(job_name=job_name,  
                                               cwd=os.getcwd(), 
                                               sge_option=sge_option, 
                                               script=script_fn)

        os.system( sge_cmd )
    elif job_type == "local":
        os.system( "bash %s" % job_data["script_fn"] )

def wait_for_file(filename, task = None, job_name = ""):
    while 1:
        time.sleep(30)
        if os.path.exists(filename):
            break

        if task != None:
            if task.shutdown_event != None and task.shutdown_event.is_set(): 
                os.system("qdel %s" % job_name)
                break

def split_fofn( input_fofn, out_dir, prefix, size_of_chunk, incremental = True, allow_fraction = True):

    input_fn = set()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    with open(input_fofn) as f:
        for l in f:
            l = l.strip()
            input_fn.add(l)
            
    existing_fn = set()
    current_chunk_fn = glob.glob(os.path.join(out_dir, "%s_[0-9]*.fofn" % prefix))

    if incremental:
        for fn in current_chunk_fn:
            with open(fn) as f:
                for l in f:
                    l = l.strip()
                    existing_fn.add(l)
    else:
        for fn in current_chunk_fn:
            os.remove(fn)

    if len(existing_fn) != 0:
        chunk_fn_base = [ os.path.basename( fn ) for fn in current_chunk_fn  ]
        sn_re = re.compile( r"%s_([0-9]*).fofn" % prefix )
        try:
            serial_number = [  int(sn_re.match( fn ).group(1)) for fn in chunk_fn_base ]
        except:
            print  "Can't parse the serial number of the chunked files"
            sys.exit(-1)
        max_sn = max(serial_number)
        assert max_sn > 0
    else:
        max_sn = 0

    new_sn = max_sn + 1
    new_fn = input_fn - existing_fn
    new_fn = list(new_fn)
    new_fn.sort()
    for  i in range(0, len(new_fn), size_of_chunk):
        out_fn = os.path.join(out_dir, "%s_%05d.fofn" % (prefix, new_sn))
        if i + size_of_chunk <= len(new_fn):
            with open(out_fn, "w") as out_f:
                for j in range(i, i + size_of_chunk):
                    out_f.write(new_fn[j]+"\n")
            new_sn += 1

        elif allow_fraction:
            with open(out_fn, "w") as out_f:
                for j in range(i, len(new_fn) ):
                    out_f.write(new_fn[j]+"\n")
            new_sn += 1

def blasr_align(self):

    q_fofn = self.query_fofn
    target_fa = self.target_fa
    target_sa = self.target_sa
    output_dir = self.parameters["mapping_data_dir"]
    q_sn = self.parameters["q_sn"]
    t_sn = self.parameters["t_sn"]
    out_fn = os.path.join( output_dir, "q%05d_t%05d.m4" % (q_sn, t_sn))
    script_fn = os.path.join( output_dir, "q%05d_t%05d.sh" % (q_sn, t_sn))
    config = self.parameters["config"]
    blasr_opt = config["blasr_opt"]
    sge_option_dm = config["sge_option_dm"]
    install_prefix = config["install_prefix"]

    #blasr_cmd = """blasr %s %s -sa %s -noSplitSubreads -bestn 16 -nCandidates 32 -maxScore -1000 -minMatch 12 -maxLCPLength 15 -nproc 16 -m 4 -out %s""" % (fn(q_fofn), fn(target_fa), fn(target_sa), out_fn)

    blasr_cmd = """blasr {query} {target} -sa {target_sa} {blasr_opt} -noSplitSubreads -m 4 -out {out_fn}"""
    blasr_cmd = blasr_cmd.format( query=fn(q_fofn), target=fn(target_fa), target_sa=fn(target_sa), blasr_opt = blasr_opt, out_fn=out_fn )
                                                                                                                     

    with open(script_fn,"w") as script_file:
        script_file.write("source {install_prefix}/bin/activate\n".format(install_prefix = install_prefix))
        script_file.write(blasr_cmd+"\n")
        script_file.write("touch %s" % fn(self.job_done))

    job_name = self.URL.split("/")[-1]
    job_name += str(uuid.uuid1())[:8]
    job_data = {"job_name": job_name,
                "cwd": os.getcwd(),
                "sge_option": sge_option_dm,
                "script_fn": script_fn }
    run_script(job_data, job_type = config["job_type"])

    wait_for_file( fn(self.job_done), task=self, job_name=job_name )

def qrm_align(self):

    q_fofn = self.query_fofn
    t_fofn = self.target_fofn
    output_dir = self.parameters["mapping_data_dir"]
    q_sn = self.parameters["q_sn"]
    t_sn = self.parameters["t_sn"]
    out_fn = os.path.join( output_dir, "q%05d_t%05d.m4" % (q_sn, t_sn))
    script_fn = os.path.join( output_dir, "q%05d_t%05d.sh" % (q_sn, t_sn))
    config = self.parameters["config"]
    #blasr_opt = config["blasr_opt"]
    qrm_opt = config["qrm_opt"]
    sge_option_dm = config["sge_option_dm"]
    install_prefix = config["install_prefix"]

    #blasr_cmd = """blasr %s %s -sa %s -noSplitSubreads -bestn 16 -nCandidates 32 -maxScore -1000 -minMatch 12 -maxLCPLength 15 -nproc 16 -m 4 -out %s""" % (fn(q_fofn), fn(target_fa), fn(target_sa), out_fn)

    qrm_cmd = """falcon_qrm.py {qrm_opt} {target} {query} > {out_fn}"""
    qrm_cmd = qrm_cmd.format( query=fn(q_fofn), target=fn(t_fofn), qrm_opt=qrm_opt, out_fn=out_fn )
                                                                                                                     

    with open(script_fn,"w") as script_file:
        script_file.write("source {install_prefix}/bin/activate\n".format(install_prefix = install_prefix))
        script_file.write(qrm_cmd+"\n")
        script_file.write("touch %s" % fn(self.job_done))

    job_name = self.URL.split("/")[-1]
    job_name += str(uuid.uuid1())[:8]
    job_data = {"job_name": job_name,
                "cwd": os.getcwd(),
                "sge_option": sge_option_dm,
                "script_fn": script_fn }
    run_script(job_data, job_type = config["job_type"])

    wait_for_file( fn(self.job_done), task=self, job_name=job_name )

def query_filter(self):
    #print self.parameters
    #print [fn(f) for f in self.inputs.values()]
    output_dir = self.parameters["mapping_data_dir"]
    q_sn = self.parameters["q_sn"]
    script_fn = os.path.join( output_dir, "qf%05d.sh" % q_sn)
    qf_fofn = os.path.join( output_dir, "qf%05d_input.fofn" % (q_sn, ) )
    install_prefix = config["install_prefix"]
    sge_option_qf = config["sge_option_qf"]
    length_cutoff_pr = config["length_cutoff_pr"]
    bestn = config["bestn"]

    with open(script_fn,"w") as script_file:
        script_file.write("source {install_prefix}/bin/activate\n".format(install_prefix = install_prefix))
        script_file.write("""find %s -name "q[0-9]*_t[0-9]*.m4" > %s\n""" % (output_dir, qf_fofn))
        script_file.write("""query_m4_filtering.py %s 1 0 %d %d %s\n""" % (qf_fofn, bestn, length_cutoff_pr, fn(self.qf_out) ))
        script_file.write("""touch %s\n""" % fn(self.job_done) )

    job_name = self.URL.split("/")[-1]
    job_name += str(uuid.uuid1())[:8]
    job_data = {"job_name": job_name,
                "cwd": os.getcwd(),
                "sge_option": sge_option_qf,
                "script_fn": script_fn }
    run_script(job_data, job_type = config["job_type"])

    wait_for_file( fn(self.job_done), task=self, job_name=job_name )

def get_preads(self):

    q_fofn_fn = fn( self.query_fa_fofn )
    t_fofn_fn = fn( self.target_fa_fofn )
    qm4_fofn_fn = fn( self.qm4_fofn)
    pa_dir = self.parameters["pa_dir"] 
    pa_out = self.pa_out 
    pa_chunk = self.parameters["pa_chunk"]
    config = self.parameters["config"]
    min_cov = config["min_cov"]
    max_cov = config["max_cov"]
    trim_align = config["trim_align"]
    trim_plr = config["trim_plr"]
    nproc = config["q_nproc"]
    tmp_dir = config["tmpdir"]
    bestn = config["bestn"]
    preassembly_num_chunk = config["preassembly_num_chunk"]
    install_prefix = config["install_prefix"]
    sge_option_pa = config["sge_option_pa"]

    
    script_fn = os.path.join( pa_dir, "pa%05d.sh" % pa_chunk)
    with open(script_fn,"w") as script_file:
        script_file.write("source {install_prefix}/bin/activate\n".format(install_prefix = install_prefix))
        script_file.write("""get_rdata.py %s %s %s %d %d %d %d %d %d %d | falcon_wrap.py > %s \n""" % (q_fofn_fn, t_fofn_fn, qm4_fofn_fn,
                                                                                          bestn, pa_chunk, preassembly_num_chunk,
                                                                                          min_cov, max_cov, trim_align, trim_plr,
                                                                                          fn(pa_out) ) )
        script_file.write("""touch %s\n""" % fn(self.pa_job_done) )

    job_name = self.URL.split("/")[-1]
    job_name += str(uuid.uuid1())[:8]
    job_data = {"job_name": job_name,
                "cwd": os.getcwd(),
                "sge_option": sge_option_pa,
                "script_fn": script_fn }
    run_script(job_data, job_type = config["job_type"])
    wait_for_file( fn(self.pa_job_done), task=self, job_name=job_name )
    
def get_config(config_fn):

    import ConfigParser

    config = ConfigParser.RawConfigParser()

    config.read(config_fn)
    
    job_type = "SGE"
    if config.has_option('General', 'job_type'):
        job_type = config.get('General', 'job_type')
    
    concurrent_jobs = 8
    if config.has_option('General', 'concurrent_jobs'):
        concurrent_jobs = config.getint('General', 'concurrent_jobs')

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

    qrm_opt = """ --min_len 500 --n_core 24 --d_core 3"""
    if config.has_option('General', 'qrm_opt'):
        blasr_opt = config.get('General', 'qrm_opt')

    tmpdir = "/tmp"
    if config.has_option('General', 'tmpdir'):
        tmpdir = config.get('General', 'tmpdir')

    big_tmpdir = "/tmp"
    if config.has_option('General', 'big_tmpdir'):
        big_tmpdir = config.get('General', 'big_tmpdir')
    
    use_CA_spec = ""
    if config.has_option('General', 'use_CA_spec'):
        use_CA_spec = config.get('General', 'use_CA_spec')

    if not config.has_option('General', 'SEYMOUR_HOME'):
        print """ SEYMOUR_HOME not found in the configuration file, quiver can not be run """
        SEYMOUR_HOME = None
    else:
        SEYMOUR_HOME = config.get('General', 'SEYMOUR_HOME')

    if config.has_option('General', 'target'):
        target = config.get('General', 'target')
        if target not in ["mapping", "pre_assembly", "falcon_asm"]:
            print """ Target has to be "mapping", or "pre_assembly" in this verison. You have an unknown target %s in the configuration file.  """ % target
            sys.exit(1)
    else:
        print """ No target specified, assuming a "pre-assembly" only target """
        target = "pre_assembly"

    preassembly_num_chunk = 12
    if config.has_option('General', 'preassembly_num_chunk'):
        preassembly_num_chunk = config.getint('General', 'preassembly_num_chunk')

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

    q_chunk_size = 2
    if config.has_option('General', 'q_chunk_size'):
        q_chunk_size = config.getint('General', 'q_chunk_size')

    t_chunk_size = 10 
    if config.has_option('General', 't_chunk_size'):
        t_chunk_size = config.getint('General', 't_chunk_size')

    SYMOURE_HOME = config.get("General", "SEYMOUR_HOME")
    hgap_config = {"input_fofn_fn" : input_fofn_fn,
                   "job_type" : job_type,
                   "concurrent_jobs" : concurrent_jobs,
                   "length_cutoff" : length_cutoff,
                   "length_cutoff_pr" : length_cutoff_pr,
                   "sge_option_dm": config.get('General', 'sge_option_dm'),
                   "sge_option_qf": config.get('General', 'sge_option_qf'),
                   "sge_option_pa": config.get('General', 'sge_option_pa'),
                   "sge_option_fca": config.get('General', 'sge_option_fca'),
                   "sge_option_qv": config.get('General', 'sge_option_qv'),
                   "bestn" : bestn,
                   "blasr_opt" : blasr_opt,
                   "qrm_opt" : qrm_opt,
                   "RQ_threshold" : RQ_threshold,
                   "tmpdir" : tmpdir,
                   "SEYMOUR_HOME" : SEYMOUR_HOME,
                   "target" : target,
                   "preassembly_num_chunk": preassembly_num_chunk,
                   "big_tmpdir": big_tmpdir,
                   "use_CA_spec": use_CA_spec,
                   "min_cov": min_cov,
                   "max_cov": max_cov,
                   "trim_align": trim_align,
                   "trim_plr": trim_plr,
                   "q_nproc": q_nproc,
                   "q_chunk_size": q_chunk_size,
                   "t_chunk_size": t_chunk_size}

    hgap_config["install_prefix"] = sys.prefix
    
    return hgap_config


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print "you need to specify a configuration file"
        print "example: HGAP.py HGAP_run.cfg"
        sys.exit(1)
    
    fasta_dir = os.path.abspath("./0-fasta_files")
    dist_map_dir = os.path.abspath("./1-dist_map-falcon")
    pa_dir = os.path.abspath("./2-preads-falcon")
    falcon_asm_dir  = os.path.abspath("./3-asm-falcon")
    script_dir = os.path.abspath("./scripts")
    sge_log_dir = os.path.abspath("./sge_log")

    for d in (dist_map_dir, fasta_dir, pa_dir, script_dir, falcon_asm_dir,  sge_log_dir):
        try:
            os.makedirs(d)
        except:
            pass

    config = get_config(sys.argv[1])
    concurrent_jobs = config["concurrent_jobs"]
    PypeThreadWorkflow.setNumThreadAllowed(concurrent_jobs, concurrent_jobs)
    wf = PypeThreadWorkflow()


    #### Task to convert bas.h5 and bax.h5 to fasta files, it will generates two fofn files for the queries and targets
    input_h5_fofn = makePypeLocalFile(os.path.abspath( config["input_fofn_fn"] ))
    query_fa_fofn = makePypeLocalFile( os.path.join( fasta_dir, "queries.fofn" ) )
    target_fa_fofn = makePypeLocalFile( os.path.join( fasta_dir, "targets.fofn" ) )
    fasta_dump_done = makePypeLocalFile(os.path.abspath( os.path.join( fasta_dir, "fasta_dump_done") ) )
    parameters = {"fasta_dir": fasta_dir,
                  "min_length": config["length_cutoff"],
                  "min_read_score": config["RQ_threshold"]}

    @PypeTask(inputs = {"input_fofn": input_h5_fofn},
              outputs = {"fasta_dump_done": fasta_dump_done, 
                         "target_fa_fofn": target_fa_fofn,
                         "query_fa_fofn":  query_fa_fofn},
              parameters = parameters,
              TaskType = PypeThreadTaskBase)
    def h5fofn_to_fasta(self):
        os.system("h5fofn_to_fasta.py %s %s --min_length 500 --min_seed_length %d --min_read_score %f" %\
                   (fn(self.input_fofn), 
                    self.parameters["fasta_dir"], 
                    self.parameters["min_length"],
                    self.parameters["min_read_score"]))
        os.system("""find %s -name "*_t.fa" | sort > %s""" % (self.parameters["fasta_dir"], fn(self.target_fa_fofn)))
        os.system("""find %s -name "*_q.fa" | sort > %s""" % (self.parameters["fasta_dir"], fn(self.query_fa_fofn)))
        os.system("touch %s" % fn(self.fasta_dump_done))

    wf.addTasks([h5fofn_to_fasta])
    #we need to force the execution of the graph at this point to ensure generating correct downstream graph 
    wf.refreshTargets([fasta_dump_done]) 

     
    #### Task to split the fofn file into small chunks for parallel processing
    split_fofn_done = makePypeLocalFile(os.path.abspath( os.path.join( dist_map_dir, "split_fofn_done") ) )

    @PypeTask(inputs = {"target_fa_fofn": target_fa_fofn,
                        "query_fa_fofn":  query_fa_fofn},
              outputs = {"split_fofn_done": split_fofn_done},   
              parameters = {"config":config, "dist_map_dir": dist_map_dir},
              TaskType = PypeThreadTaskBase)
    def split_fofn_task(self):
        query_chunk_size = self.parameters["config"]["q_chunk_size"]
        target_chunk_size = self.parameters["config"]["t_chunk_size"]
        split_fofn( fn(self.query_fa_fofn), self.parameters["dist_map_dir"], "query", query_chunk_size, 
                    incremental = True, allow_fraction = True)
        split_fofn( fn(self.target_fa_fofn), self.parameters["dist_map_dir"], "target", target_chunk_size, 
                    incremental = True, allow_fraction = True)

        os.system("touch %s" % fn(self.split_fofn_done))

    wf.addTasks([split_fofn_task])
    #we need to force the execution of the graph at this point to ensure generating correct downstream graph 
    wf.refreshTargets([split_fofn_done])


    #### According to the split fofn, generate the targets sequence files and suffix array database
    gather_targets_done = makePypeLocalFile(os.path.abspath( os.path.join( dist_map_dir, "gather_target_done") ) )

    @PypeTask(inputs = {"split_fofn_done": split_fofn_done},
              outputs = {"gather_targets_done": gather_targets_done},   
              parameters = {"config":config, "dist_map_dir": dist_map_dir},
              TaskType = PypeThreadTaskBase)
    def gather_targets(self):
        os.system("touch %s" % fn(self.gather_targets_done))

                
    wf.addTasks([gather_targets])

    #we need to force the execution of the graph at this point to ensure generating correct downstream graph 
    wf.refreshTargets([gather_targets_done])


    #### Submit the MxN mapping jobs
    q_re = re.compile( r"query_([0-9]*).fofn"  )
    t_re = re.compile( r"target_([0-9]*).fofn"  )
    query_fofns  = glob.glob(os.path.join(dist_map_dir, "query_[0-9]*.fofn"))
    query_fofns.sort()
    target_fofns  = glob.glob(os.path.join(dist_map_dir, "target_[0-9]*.fofn"))
    target_fofns.sort()

    URL_to_object = {}
    all_qf_out = {}

    for q_fofn in query_fofns:
        q_dir, q_basename = os.path.split(q_fofn)
        q_sn = int(q_re.match( q_basename ).group(1))
        q_fofn = makePypeLocalFile(q_fofn)
        if q_fofn.URL in URL_to_object:
            q_fofn = URL_to_object[q_fofn.URL]
        else:
            URL_to_object[q_fofn.URL] = q_fofn 

        mapping_data_dir = os.path.join( dist_map_dir, "q%05d_md" % q_sn)
        mapping_data_dir = os.path.abspath( mapping_data_dir )
        try:
            os.makedirs( mapping_data_dir )
        except:
            pass

        query_group_done = {} 

        for t_fofn in target_fofns:

            t_dir, t_basename = os.path.split(t_fofn)
            t_sn = int(t_re.match( t_basename ).group(1))
            t_fofn = makePypeLocalFile(t_fofn)

            job_done = os.path.join( mapping_data_dir, "q%05d_t%05d_done" % (q_sn, t_sn) ) 
            job_done = makePypeLocalFile( job_done)
            query_group_done["q%05d_t%05d_done" % (q_sn, t_sn)] = job_done

            if t_fofn.URL in URL_to_object:
                t_fofn = URL_to_object[t_fofn.URL]
            else:
                URL_to_object[t_fofn.URL] = t_fofn 

            inputs = { "query_fofn" : q_fofn, "target_fofn": t_fofn}
            outputs = { "job_done": job_done }
            parameters = { "mapping_data_dir": mapping_data_dir, "q_sn": q_sn, "t_sn": t_sn, "config":config }
        
            make_mapping_task = PypeTask(inputs = inputs, 
                                 outputs = outputs, 
                                 parameters = parameters, 
                                 TaskType = PypeThreadTaskBase,
                                 URL = "task://localhost/mapping_task_q%05d_t%05d" % (q_sn, t_sn) )
            mapping_task = make_mapping_task ( qrm_align )
            wf.addTask(mapping_task)

        qf_out = os.path.join( mapping_data_dir, "qf%05d.m4" % q_sn ) 
        qf_out = makePypeLocalFile( qf_out )

        all_qf_out["qf_out_%s" % q_sn] = qf_out

        job_done = os.path.join( mapping_data_dir, "qf%05d_done" % q_sn ) 
        job_done = makePypeLocalFile(job_done)

        all_qf_out["qf_done_%s" % q_sn] = job_done
        parameters = { "mapping_data_dir": mapping_data_dir, "q_sn": q_sn }
        make_qf_task = PypeTask(inputs = query_group_done,
                                outputs = {"qf_out": qf_out, "job_done": job_done},
                                TaskType = PypeThreadTaskBase,
                                URL = "task://localhost/qf_task_q%05d" % q_sn,
                                parameters = parameters)
        qf_task = make_qf_task(query_filter)
        wf.addTask(qf_task)

    #### END: Submit the MxN mapping jobs

    qm4_fofn = os.path.join( pa_dir, "m4_files.fofn") 
    qm4_fofn = makePypeLocalFile( qm4_fofn )
    @PypeTask(inputs = all_qf_out, 
              outputs = {"qm4_fofn": qm4_fofn},   
              TaskType = PypeThreadTaskBase)
    def gather_qm4(self):
        all_qf = [fn(o) for o in self.inputs.values()]
        all_qf.sort()
        with open( fn( self.qm4_fofn ),"w" ) as f:
            for m4f in all_qf:
                if m4f.endswith("m4"):
                    print >> f, m4f
    wf.addTask(gather_qm4)

    
    ### Submit pre-assembly jobs

    inputs = {"target_fa_fofn": target_fa_fofn,
              "query_fa_fofn":  query_fa_fofn,
              "qm4_fofn":qm4_fofn}
    all_pa_out = {}
    for pa_chunk in range(config["preassembly_num_chunk"]): #...
        pa_job_done = os.path.join( pa_dir, "pa_%05d_done" % pa_chunk ) 
        pa_job_done = makePypeLocalFile( pa_job_done )
        pa_out = os.path.join( pa_dir, "pread_%05d.fa" % pa_chunk )
        pa_out = makePypeLocalFile( pa_out )
        make_pread_task = PypeTask(inputs = inputs,
                                   outputs = {"pa_job_done": pa_job_done,
                                              "pa_out": pa_out},   
                                   parameters = {"config":config, "pa_chunk":pa_chunk, "pa_dir": pa_dir},
                                   URL = "task://localhost/pa_task_%05d" % pa_chunk,
                                   TaskType = PypeThreadTaskBase)
        all_pa_out["pa_out_%05d" % pa_chunk] = pa_out
        all_pa_out["pa_done_%05d" % pa_chunk] = pa_job_done
        pread_task = make_pread_task( get_preads )
        wf.addTask( pread_task )
        
    pr_fofn = os.path.join( pa_dir, "pr_fasta.fofn") 
    pr_fofn = makePypeLocalFile( pr_fofn )
    @PypeTask(inputs = all_pa_out, 
              outputs = {"pr_fofn": pr_fofn},   
              TaskType = PypeThreadTaskBase)
    def gather_pr(self):
        all_pr = [fn(o) for o in self.inputs.values()]
        all_pr.sort()
        with open( fn( self.pr_fofn ),"w" ) as f:
            for fqn in all_pr:
                print >> f, fqn
    wf.addTask(gather_pr)


    fca_done = os.path.join( falcon_asm_dir, "fca_done"  ) 
    fca_done = makePypeLocalFile(fca_done)
    @PypeTask(inputs = {"pr_fofn": pr_fofn}, 
              outputs = {"fca_done": fca_done},   
              parameters = {"config": config, "fca_dir": falcon_asm_dir},
              TaskType = PypeThreadTaskBase)

    def run_FCA(self):
        from pbcore.io import FastaReader
        fca_dir = self.parameters["fca_dir"]
        config = self.parameters["config"]
        install_prefix = config["install_prefix"]
        sge_option_fca = config["sge_option_fca"]

        fasta_fn = os.path.abspath(os.path.join(fca_dir, "preads.fa"))
        with open(fasta_fn, "w") as fa_f:
            with open(fn(self.pr_fofn)) as fa_fofn:
                for fan in fa_fofn:
                    fan = fan.strip()
                    fan_h = FastaReader(fan)
                    for r in fan_h:
                        if len(r.sequence) > 8000:
                            fa_f.write(">"+r.name+"\n")
                            fa_f.write( r.sequence+"\n")
    
        fca_cmd  = "cd %s\n" % fca_dir
        fca_cmd += "falcon_overlap.py --n_core 24 --d_core 24 preads.fa | grep -v none > preads.ovlp\n"
        fca_cmd += "falcon_asm.py preads.ovlp preads.fa\n"
        fca_cmd += "falcon_fixasm.py\n"
                                                                                                                         

        script_fn = os.path.join( fca_dir, "runfca.sh")
        with open(script_fn,"w") as script_file:
            script_file.write("source {install_prefix}/bin/activate\n".format(install_prefix = install_prefix))
            script_file.write(fca_cmd+"\n")
            script_file.write("touch %s" % fn(self.fca_done))

        job_name = self.URL.split("/")[-1]
        job_name += str(uuid.uuid1())[:8]
        job_data = {"job_name": job_name,
                    "cwd": os.getcwd(),
                    "sge_option": sge_option_fca,
                    "script_fn": script_fn }
        run_script(job_data, job_type = config["job_type"])
        wait_for_file( fn(self.fca_done), task=self, job_name=job_name )

    wf.addTask(run_FCA)

    with open("./workflow.dot","w") as f:
        print >>f, wf.graphvizShortNameDot

    if config["target"] == "mapping":
        wf.refreshTargets([qm4_fofn])
    elif config["target"] == "pre_assembly":
        wf.refreshTargets([pr_fofn])
    elif config["target"] == "falcon_asm":
        wf.refreshTargets() #all
    else:
        wf.refreshTargets() #all

    with open("./workflow.dot","w") as f:
        print >>f, wf.graphvizShortNameDot

