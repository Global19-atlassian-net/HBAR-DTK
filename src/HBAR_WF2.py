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

def wait_for_file(filename):
    while 1:
        time.sleep(10)
        if os.path.exists(filename):
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
    job_data = {"job_name": job_name,
                "cwd": os.getcwd(),
                "sge_option": sge_option_dm,
                "script_fn": script_fn }
    run_script(job_data, job_type = config["job_type"])

    wait_for_file( fn(self.job_done) )

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
    job_data = {"job_name": job_name,
                "cwd": os.getcwd(),
                "sge_option": sge_option_qf,
                "script_fn": script_fn }
    run_script(job_data, job_type = config["job_type"])

    wait_for_file( fn(self.job_done) )

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
        script_file.write("""get_preads.py %s %s %s %d %d %d %d %d %d %d %d %s %s\n""" % (q_fofn_fn, t_fofn_fn, qm4_fofn_fn,
                                                                                          bestn, pa_chunk, preassembly_num_chunk,
                                                                                          min_cov, max_cov, trim_align, trim_plr,
                                                                                          nproc, fn(pa_out), tmp_dir ))
        script_file.write("""touch %s\n""" % fn(self.pa_job_done) )

    job_name = self.URL.split("/")[-1]
    job_data = {"job_name": job_name,
                "cwd": os.getcwd(),
                "sge_option": sge_option_pa,
                "script_fn": script_fn }
    run_script(job_data, job_type = config["job_type"])
    wait_for_file( fn(self.pa_job_done) )
    
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
        if target not in ["mapping", "all", "draft_assembly", "pre_assembly"]:
            print """ Target has to be "mapping", "all", "draft_assembly" or "pre_assembly". You have an unknown target %s in the configuration file.  """ % target
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
    dist_map_dir = os.path.abspath("./1-dist_map")
    pa_dir = os.path.abspath("./2-preads")
    celera_asm_dir  = os.path.abspath("./3-CA")
    script_dir = os.path.abspath("./scripts")
    sge_log_dir = os.path.abspath("./sge_log")

    for d in (dist_map_dir, fasta_dir, pa_dir, script_dir, celera_asm_dir,  sge_log_dir):
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
                  "min_seed_length": config["length_cutoff_pr"], 
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

        targets_dir = os.path.join(self.parameters["dist_map_dir"],"targets/") 
        try:
            os.makedirs(targets_dir)
        except:
            pass

        for fofn in glob.glob(os.path.join(self.parameters["dist_map_dir"], "target_[0-9]*.fofn")):
            fofn_base = os.path.basename(fofn).split(".")[0]
            out_f_fa = os.path.join(targets_dir, fofn_base+".fa")
            out_f_sa = os.path.join(targets_dir, fofn_base+".sa")
            with open( out_f_fa, "w" ) as out_f:
                with open(fofn) as f:
                    for l in f:
                        l = l.strip()
                        with open(l) as fasta_f:
                            out_f.write(fasta_f.read())
            os.system("sawriter %s %s -blt 12" % ( out_f_sa, out_f_fa))

        os.system("touch %s" % fn(self.gather_targets_done))

                
    wf.addTasks([gather_targets])

    #we need to force the execution of the graph at this point to ensure generating correct downstream graph 
    wf.refreshTargets([gather_targets_done])


    #### Submit the MxN mapping jobs
    q_re = re.compile( r"query_([0-9]*).fofn"  )
    t_re = re.compile( r"target_([0-9]*).fa"  )
    query_fofns  = glob.glob(os.path.join(dist_map_dir, "query_[0-9]*.fofn"))
    query_fofns.sort()
    target_fa  = glob.glob(os.path.join(dist_map_dir, "targets", "target_[0-9]*.fa"))
    target_fa.sort()

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

        for t_fa in target_fa:
            t_dir, t_basename = os.path.split(t_fa)
            t_sa = os.path.join(t_dir, t_basename.split(".")[0] + ".sa")
            t_fa = makePypeLocalFile(t_fa)

            if t_fa.URL in URL_to_object:
                t_fa = URL_to_object[t_fa.URL]
            else:
                URL_to_object[t_fa.URL] = t_fa

            t_sa = makePypeLocalFile(t_sa)
        
            if t_sa.URL in URL_to_object:
                t_sa = URL_to_object[t_sa.URL]
            else:
                URL_to_object[t_sa.URL] = t_sa

            t_sn = int(t_re.match( t_basename ).group(1))


            job_done = os.path.join( mapping_data_dir, "q%05d_t%05d_done" % (q_sn, t_sn) ) 
            job_done = makePypeLocalFile( job_done)
            query_group_done["q%05d_t%05d_done" % (q_sn, t_sn)] = job_done

            inputs = { "query_fofn" : q_fofn, "target_fa": t_fa, "target_sa": t_sa, "gather_targets_done": gather_targets_done }
            outputs = { "job_done": job_done }
            parameters = { "mapping_data_dir": mapping_data_dir, "q_sn": q_sn, "t_sn": t_sn, "config":config }
        
            #for testing 
            #task_decorator = PypeTask(inputs = inputs, outputs = outputs, parameters = parameters, TaskType = PypeTaskBase )
            #task() 
            make_mapping_task = PypeTask(inputs = inputs, 
                                 outputs = outputs, 
                                 parameters = parameters, 
                                 TaskType = PypeThreadTaskBase,
                                 URL = "task://localhost/mapping_task_q%05d_t%05d" % (q_sn, t_sn) )
            mapping_task = make_mapping_task ( blasr_align )
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

    ca_done = os.path.join( celera_asm_dir, "ca_done"  ) 
    ca_done = makePypeLocalFile(ca_done)
    @PypeTask(inputs = {"pr_fofn": pr_fofn}, 
              outputs = {"ca_done": ca_done},   
              parameters = {"config": config, "ca_dir": celera_asm_dir},
              TaskType = PypeThreadTaskBase)

    def run_CA(self):
        from pbcore.io import FastaReader
        ca_dir = self.parameters["ca_dir"]
        config = self.parameters["config"]
        install_prefix = config["install_prefix"]
        sge_option_ca = config["sge_option_ca"]

        fastq_fn = os.path.abspath(os.path.join(ca_dir, "preads.fastq"))
        with open(fastq_fn, "w") as fq:
            with open(fn(self.pr_fofn)) as fa_fofn:
                for fan in fa_fofn:
                    fan = fan.strip()
                    fan_h = FastaReader(fan)
                    for r in fan_h:
                        fq.write("@"+r.name+"\n")
                        fq.write( r.sequence+"\n")
                        fq.write( "+" + "\n" )
                        fq.write( "".join( [chr(33+24)] * len(r.sequence) ) + "\n" )

        ca_cmd  = "cd %s\n" % ca_dir
        ca_cmd += "fastqToCA -technology sanger -type sanger -libraryname pr_long -reads %s/preads.fastq > preads.frg\n" % ca_dir
        ca_cmd += "runCA -d . -p asm -s {install_prefix}/etc/asm.spec  preads.frg\n".format( install_prefix = install_prefix )
                                                                                                                         

        script_fn = os.path.join( ca_dir, "runca.sh")
        with open(script_fn,"w") as script_file:
            script_file.write("source {install_prefix}/bin/activate\n".format(install_prefix = install_prefix))
            script_file.write(ca_cmd+"\n")
            script_file.write("touch %s" % fn(self.ca_done))

        job_name = self.URL.split("/")[-1]
        job_data = {"job_name": job_name,
                    "cwd": os.getcwd(),
                    "sge_option": sge_option_ca,
                    "script_fn": script_fn }
        run_script(job_data, job_type = config["job_type"])
        wait_for_file( fn(self.ca_done) )

    wf.addTask(run_CA)

    with open("./workflow.dot","w") as f:
        print >>f, wf.graphvizShortNameDot

    if config["target"] == "mapping":
        wf.refreshTargets([qm4_fofn])
    elif config["target"] == "pre_assembly":
        wf.refreshTargets([pr_fofn])
    else:
        wf.refreshTargets() #all

    with open("./workflow.dot","w") as f:
        print >>f, wf.graphvizShortNameDot

