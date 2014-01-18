Installing and running the HBAR-DTK scripts
===========================================

:Authors: 
    Jason Chin

:Version: 0.1.1 of 2013/04/30


Prerequisites:

* A resonable size linux cluster with SGE cluster setup ( You will have to hack
  the code for job submissions if you don't use SGE. ) If one really wants, the
  code and its dependencies can be probably run on OS X, but why?
* python2.7
* gcc 4.4.3
* git ( http://git-scm.com/ for installing the dependencies from github )
* blasr ( https://github.com/PacificBiosciences/blasr )
* Quiver from SMRTAnalysis > v1.4 ( from http://pacbiodevnet.com ) or
  GenomicConsensus ( https://github.com/PacificBiosciences/GenomicConsensus )
  This is for the final consensus with Quiver. It is not required for
  pre-assembly. 


Install Python
--------------

Make sure you are using python2.7. We can create a clean virtualenv and
activate it::

    $ export HBAR_HOME=/some/path/to/your/HBAR_ENV
    $ virtualenv -p /usr/bin/python2.7 $HBAR_HOME
    $ cd $HBAR_HOME
    $ . bin/activate

Next you want to install ``pbcore`` ( http://www.numpy.org ) library and its
dependencies. First install ``numpy``::

    $ pip install numpy==1.6.2

We need to complile ``libhdf5`` (
http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-1.8.9/src/ ) and install it
in the virtualenv::

    $ wget http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-1.8.9/src/hdf5-1.8.9.tar.gz
    $ tar zxvf hdf5-1.8.9.tar.gz
    $ cd hdf5-1.8.9
    $ ./configure --prefix=$HBAR_HOME --enable-cxx
    $ make install
    $ cd ..

Install ``h5py`` ( http://h5py.googlecode.com/files/h5py-2.0.1.tar.gz )::

    $ wget http://h5py.googlecode.com/files/h5py-2.0.1.tar.gz
    $ tar zxvf h5py-2.0.1.tar.gz
    $ cd h5py-2.0.1
    $ python setup.py build --hdf5=$HBAR_HOME
    $ python setup.py install


If you alread have these libraries installed in your current python environment,
you can and maybe you should skip some of these steps.

Install HBAR Python Libraries
-----------------------------

Next, install the PacBio python libraries::
    

    $ pip install git+https://github.com/PacificBiosciences/pbcore.git#pbcore
    $ pip install git+https://github.com/PacificBiosciences/pbdagcon.git#pbdagcon
    $ pip install git+https://github.com/PacificBiosciences/pbh5tools.git#pbh5tools
    $ pip install git+https://github.com/cschin/pypeFLOW.git#pypeflow
    $ pip install git+https://github.com/PacificBiosciences/HBAR-DTK.git#hbar-dtk

If you do a ``pip freeze``, this is what you will see::

    $ pip freeze
    h5py==2.0.1
    html5lib==0.95
    isodate==0.4.9
    matplotlib==1.2.0
    numpy==1.6.2
    pbcore==0.6.0
    pbtools.hbar-dtk==0.1.0
    pbtools.pbdagcon==0.2.0
    pbtools.pbh5tools==0.75.0
    pyparsing==1.5.7
    pypeflow==0.1.0
    rdfextras==0.4
    rdflib==3.4.0
    wsgiref==0.1.2

Install Other HBAR Prerequisites
--------------------------------

We need BLASR for the pre-assembly mapping. BLASR is included in the SMRT(R)
Analysis installation and is also available on github. You need to copy a blasr
binary into your ``$HBAR_HOME/bin``::

    $ cp blasr $HBAR_HOME/bin

Last, we need a copy of Celera Assembler for the assembly itself::

    $ wget http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-7.0/wgs-7.0-PacBio-Linux-amd64.tar.bz2
    $ tar jxvf wgs-7.0-PacBio-Linux-amd64.tar.bz2 -C $HBAR_HOME/bin/
    $ ln -sf $HBAR_HOME/bin/wgs-7.0/Linux-amd64/bin/* $HBAR_HOME/bin/
 
For the final step of polishing the assembly using Quiver, we need a SMRT
Analysis installation, plus Quiver from github. Quiver is available via a link
from PacBio DevNet or directly on github. Please follow the installation
instructions there.

Running HBAR_WF.py
=================

Warning
---------

- While the general strategy of HBAR will work for larger genome in principle.
  Special consideration should be taken to do the distributed computing
  efficiently.

Set up the environment
-----------------------

Make sure you have clean UNIX shell environment. (Please be sure you do not
have ``PYTHON_PATH`` environment variable and other random non-standard paths
in your ``PATH`` environment variable.) If your shell environment is clean, do::

    $ export PATH_TO_HBAR_ENV=/the_full_path_to_your_installation
    $ source $PATH_TO_HBAR_ENV/bin/activate

You can "deactivate" the ``HBAR_ENV`` by::
 
    $ deactivate

Prepare data, set up the configuration and run
----------------------------------------------

Prepare a working directory and create a file ``input.fofn`` that points to the
base files (``bas.h5`` files) for assembly. Let call this directory
``my_assembly``.  You also need to make sure the paths in the ``input.fofn``
file are absolute and not relative paths.

Here is an example of the ``input.fofn`` files::

    /mnt/data/m120803_022519_42141_c100388772550000001523034210251234_s1_p0.bas.h5
    /mnt/data/m120803_041200_42141_c100388772550000001523034210251235_s1_p0.bas.h5
    /mnt/data/m120803_055858_42141_c100388772550000001523034210251236_s1_p0.bas.h5
    /mnt/data/m120803_074648_42141_c100388772550000001523034210251237_s1_p0.bas.h5

Copy the example configuration to the working directory::

    $ cd my_assembly
    $ cp $PATH_TO_HBAR_ENV/etc/HBAR.cfg .

Here is the content of ``HBAR.cfg``::

    [General]
    # list of files of the initial bas.h5 files
    input_fofn = input.fofn

    # The length cutoff used for seed reads used for initial mapping
    length_cutoff = 4500

    # The length cutoff used for seed reads usef for pre-assembly
    length_cutoff_pr = 4500

    # The read quality cutoff used for seed reads
    RQ_threshold = 0.75

    # SGE job option for distributed mapping 
    sge_option_dm = -pe smp 8 -q fas

    # SGE job option for m4 filtering
    sge_option_mf = -pe smp 4 -q fas

    # SGE job option for pre-assembly
    sge_option_pa = -pe smp 16 -q fas

    # SGE job option for CA 
    sge_option_ca = -pe smp 4 -q fas

    # SGE job option for Quiver
    sge_option_qv = -pe smp 16 -q fas

    # SGE job option for "qsub -sync y" to sync jobs in the different stages
    sge_option_ck = -pe smp 1 -q fas 

    # blasr for initial read-read mapping for each chunck (do not specific the "-out" option). 
    # One might need to tune the bestn parameter to match the number of distributed chunks to get more optimized results 
    blasr_opt = -nCandidates 50 -minMatch 12 -maxLCPLength 15 -bestn 4 -minPctIdentity 70.0 -maxScore -1000 -nproc 4 -noSplitSubreads

    #This is used for running quiver
    SEYMOUR_HOME = /mnt/secondary/Smrtpipe/builds/Assembly_Mainline_Nightly_Archive/build470-116466/

    #The number of best alignment hits used for pre-assembly
    #It should be about the same as the final PLR coverage, slight higher might be OK.
    bestn = 36

    # target choices are "pre_assembly", "draft_assembly", "all"
    # "pre_assembly" : generate pre_assembly for any long read assembler to use
    # "draft_assembly": automatic submit CA assembly job when pre-assembly is done
    # "all" : submit job for using Quiver to do final polish
    target = draft_assembly

    # number of chunks for distributed mapping
    preassembly_num_chunk = 8 

    # number of chunks for pre-assembly. 
    # One might want to use bigger chunk data sizes (smaller dist_map_num_chunk) to 
    # take the advantage of the suffix array index used by blasr
    dist_map_num_chunk = 4

    # "tmpdir" is for preassembly. A lot of small files are created and deleted during this process. 
    # It would be great to use ramdisk for this. Set tmpdir to a NFS mount will probably have very bad performance.
    tmpdir = /tmp

    # "big_tmpdir" is for quiver, better in a big disk
    big_tmpdir = /tmp
    
    # various trimming parameters
    min_cov = 8
    max_cov = 64
    trim_align = 50
    trim_plr = 50

    # number of processes used by by blasr during the preassembly process
    q_nproc = 16 

Please change the various ``sge_option_*`` to the proper SGE queue for the SGE
cluster to run the code.

You should estimate the overall coverage and length distribution for putting in
the correct options in the configuration file.  You will need to decide a
length cutoff for the seeding reads. The optimum cutoff length will depend on
the distribution of the sequencing read lengths, the genome size and the
overall yield. The general guideline is the coverage of the seeding sequences
should be above 20x of the genome and the overall coverage should be at least
3x of the coverage of the seeding sequences. Start the Hierarchical Genome
Assembly Process b the assembly process by::

    $ HBAR_WF.py HBAR.cfg  

If you want to kill the jobs, you should kill the python process using
``kill`` command and using ``qdel`` for the SGE jobs submitted by the python
process. 

The spec file used by the Celera Assembler is at ``$HBAR_HOME/etc/asm.spec``.
In the future, this will be configurable using the configuration file.

How to choose length cutoff
===========================

Here is some code snippet that might be useful for helping to get some
educational guess for the length cutoff.  First, loading some module::

    from pbcore.io import FastaIO
    import numpy as np
    from math import exp, log

Read the input reads and fill-in the ``seq_length`` list::

    f = FastaIO.FastaReader("all_norm.fa")
    seq_lengths = []
    for r in f:
        seq_lengths.append(len(r.sequence))
    seq_lengths = np.array(seq_lengths)

If you have `matplotlib` installed, you can check the histogram with::

    h=hist(seq_lengths,bins=50,range=(0,10000))

Set the genome size::

    genome_size = 47000000

Generate various coverage information and Lander-Waterman statistics for
different length cutoff::


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


Here is an example of the output::

    cutoff  total_base      totol/seed      cov     contig_count    contig_len/genome_size
    3000    338653714       1.36    36.03   0.00    78472724961.07
    3500    303058009       1.52    32.24   0.00    2319098472.00
    4000    267546052       1.72    28.46   0.00    68664337.09
    4500    231380061       1.99    24.61   0.00    1905600.54
    5000    197848047       2.33    21.05   0.00    69912.18
    5500    166824314       2.76    17.75   0.00    3362.59
    6000    139984665       3.29    14.89   0.00    251.54
    6500    116457831       3.96    12.39   0.04    26.81
    7000    95601124        4.82    10.17   0.26    3.82
    7500    78065025        5.90    8.30    1.29    0.78
    8000    63354342        7.28    6.74    4.68    0.21
    8500    50775028        9.08    5.40    13.47   0.07
    9000    40410887        11.41   4.30    30.49   0.03
    9500    31313271        14.72   3.33    58.92   0.02


Pick read length cutoffs that satisfy:
1. The ratio of the total number bases to the long read bases is larger than 3.
2. Estimated Lander-Waterman contig number less than 0.25. 
3. The estimated Lander-Waterman contig size is larger than 0.25x of the genome size.

::

    print "recommended cutoff (total/seed > 3, LW contig # <0.25, LW contig length > 0.25x genome)"
    print "cutoff\ttotal_base\ttotol/seed\tcov\tcontig_count\tcontig_len/genome_size" 

    for l in coverage_array[ (coverage_array[...,2]>3) & (coverage_array[...,4]<0.25) & (coverage_array[...,5]>0.25),...]:
        print "%d\t%d\t%0.2f\t%0.2f\t%0.2f\t%0.2f" % tuple(l)


The output::

    recommended cutoff (total/seed > 3, LW contig # <0.25, LW contig length > 0.25x genome)
    cutoff  total_base      totol/seed      cov     contig_count    contig_len/genome_size
    6000    139984665       3.29    14.89   0.00    251.54
    6500    116457831       3.96    12.39   0.04    26.81


In this example, length cutoffs from 6000 to 6500 satisfy the criteria.

