Installing and running the HBAR-DTK scripts
===========================================

:Authors: 
    Jason Chin

:Version: 0.1 of 2013/03/20


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
    $ python setup.py build --hdf5=$PATH_TO_HBAR_ENV
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

We need BLASR for the pre-assembly mapping. BLASR is included in the SMRT�
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

    genome_size = 22000000

Generate various coverage information and Lander-Waterman statistics for
different length cutoff::

    total = sum(seq_lengths)
    coverage_array = []
    for x in range(1000,6000,200):
        psum = sum(seq_lengths[seq_lengths>x])
        coverage = 0.5 * psum / genome_size # we loss 50% bases after the pre-assembly step
        contig_count = coverage * genome_size / x * exp( -coverage )
        contig_length = (exp(coverage) - 1) * x /coverage
        print x, psum, 1.0*total/psum, coverage,  contig_count,  contig_length/genome_size
        coverage_array.append( [x, psum, 1.0*total/psum, coverage,  contig_count, contig_length/genome_size, total] )
    coverage_array = np.array( coverage_array )

Here is an example of the output::

    1000 1795319853 1.59036696399 40.8027239318 1.70888952447e-12 585175335024.0
    1200 1684447291 1.69504703368 38.2828929773 1.6603399309e-11 60228630377.9
    1400 1581616454 1.80525270636 35.9458285 1.38314648717e-10 7229892200.71
    1600 1482267608 1.92624959797 33.6879001818 1.08469461447e-09 921918470.561
    1800 1389901045 2.05425946996 31.5886601136 7.37735364615e-09 135549961.133
    2000 1303401914 2.19058860765 29.6227707727 4.44644042439e-08 22489899.8875
    2200 1222781623 2.33501823244 27.7904914318 2.36940417001e-07 4220470.3303
    2400 1148529428 2.48597668844 26.1029415455 1.10290324698e-06 906697.847461
    2600 1079625953 2.64463574265 24.5369534773 4.58148716713e-06 218269.737205
    2800 1015165114 2.81256452239 23.0719344091 1.73115062617e-05 57765.048563
    3000 952304202 2.99821987344 21.6432773182 6.32511666946e-05 15809.9850463
    3200 893305655 3.19623789239 20.30240125 0.000212617632036 4703.2787869
    3400 834953276 3.41961336768 18.9762108182 0.000704513978576 1419.41824388
    3600 778391254 3.66810054626 17.6907103182 0.00224330115006 445.771616184
    3800 721976332 3.95472435515 16.408553 0.0071050206494 140.745533976
    4000 666984099 4.28078778532 15.1587295227 0.0217607023724 45.9543870361
    4200 614681938 4.64503218248 13.9700440455 0.06269864062 15.949295444
    4400 563419683 5.06765643826 12.8049927955 0.17587804557 5.68574235479
    4600 515441048 5.53936748941 11.7145692727 0.457950334934 2.18362505682
    4800 468934016 6.08874017789 10.6575912727 1.14896665263 0.870326807227
    5000 425971954 6.7028295107 9.68118077273 2.66009788684 0.375902539986
    5200 384694332 7.42204172636 8.743053 5.90232032443 0.169397860342
    5400 346010358 8.25182633405 7.86387177273 12.3148528561 0.0811715437399
    5600 309468288 9.22620344221 7.03337018182 24.3693682419 0.0409989309392
    5800 276458509 10.3278332591 6.28314793182 44.5077352936 0.0224260452905

Pick read length cutoffs that satisfy:
1. The ratio of the total number bases to the long read bases is larger than 3.
2. Estimated Lander-Waterman contig number less than 0.25. 
3. The estimated Lander-Waterman contig size is larger than 0.25x of the genome size.

::

    for l in coverage_array[ (coverage_array[...,2]>3) & (coverage_array[...,4]<0.25) & (coverage_array[...,5]>0.25),...]:
        print " ".join([str(c) for c in l])

The output::

    3200.0 893305655.0 3.19623789239 20.30240125 0.000212617632036 4703.2787869 2855217384.0
    3400.0 834953276.0 3.41961336768 18.9762108182 0.000704513978576 1419.41824388 2855217384.0
    3600.0 778391254.0 3.66810054626 17.6907103182 0.00224330115006 445.771616184 2855217384.0
    3800.0 721976332.0 3.95472435515 16.408553 0.0071050206494 140.745533976 2855217384.0
    4000.0 666984099.0 4.28078778532 15.1587295227 0.0217607023724 45.9543870361 2855217384.0
    4200.0 614681938.0 4.64503218248 13.9700440455 0.06269864062 15.949295444 2855217384.0
    4400.0 563419683.0 5.06765643826 12.8049927955 0.17587804557 5.68574235479 2855217384.0

In this example, length cutoffs from 3200 to 4400 satisfy the criteria.

