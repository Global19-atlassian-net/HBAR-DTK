What is HBAR-DTK?
==================

Update: Jun 13, 2014

The basic of this code base does not changed much. A new workflow script
(HBAR_WF3.py) that use the latest branch 0.1.3 Falcon
(https://github.com/PacificBiosciences/FALCON).  for doing intial overlapping,
consensus and assembly is provided. The dependence is updated and the 2nd minor
version number is bumped up to reflect the changes. Please understanding this
code is highly experimental by its nature. If something may not work for you
out of the box, I encourage you to see if you can find a fix to fit your need
and find a better way to do things. Or, you can submit some of your problem to
the Issue tracking page such that I can address them if possible.


Original Date: June 11, 2013


HBAR-DTK stands for "Hierarchical-Based AssembleR Development ToolKit". The
code is derived from the prototype process for developing the "Hierachical
Genome Assembly Process (HGAP)" that PacBio(R) provides to the scientists who
uses PacBio RS(R) for genome assembly.  For the most update-to-date instruction
on using the PacBio's official HGAP software, one should check the
https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/HGAP link to
see how to run HGAP from the official release.

HBAR-DTK provides a small number of python scripts for testing and developing
hierarchical-based assembly algorithms. While it might be useful for people who
like to test it out on using PacBio(R) data for assemblies at different scales,
it is not meant to be used by end-users who desire for a "push-button"
bioinformatics workflow for genome assemblies. 

Here is the list of the functions of the major scripts:

    HBAR_WF2.py:

        A pypeflow driven workflow script to process data and submmitting jobs
        to a SGE cluster for every steps of the hierachical genome assembly
        process. This workflow is designed for larger genome assembly project 
        where the data may not come at once. It uses different way from the 
        ``HBAR_WF.py`` for chunking the input sequencing data and gathering
        the alignment results. The old alignment data will be re-used if new 
        inputs are added. 

    HBAR_WF.py:

        A pypeflow driven workflow script to process data and submmitting jobs
        to a SGE cluster for every steps of the hierachical genome assembly
        process

    filterM4Query.py:

        A short script for filtering and sorting the alignment output from
        blasr aligner so best hits of each reads to seeds are identified from
        mutiple chuncks of the output files. This code is called by the
        HBAR_WF.py.  It is not meant to be called directly from command line.

    generate_preassemble_reads.py:
        
        A script taking the output from blasr aligner and the read sequence
        data to generate preassembled reads that can be feed into an
        assembler. This code is called by the
        HBAR_WF.py.  It is not meant to be called directly from command line.

    h5fofn_to_fasta.py:
        
        Used by ``HBAR_WF2.py`` to convert ``bas.h5``, ``bax.h5``, ``fasta``
        or ``fastq`` for its internal ``fasta`` files.

    get_pread.py:

        a modified version of generate_preassemble_reads.py used by ``HBAR-WF2.py``

    simple_asm.py:

        A preliminary implementation for the first step of "Consistent Long-read
        Evidence Essembly pRocess" (CLEAR)


Some convenient scripts:

    tig-sense.py:
        
        A hack to use pbdagcon to get consensus of all contigs from Celera(R) Assembler
        untig stage

    tig-sense_p.py:
        
        A multiprocessing version of the tig-sense.py. It can use more CPU cores for
        the consensus tasks 

    CA_best_edge_to_GML.py:
        
        A simple script to convert Celera(R) Assembler's "best.edges" to a GML
        which can be used to feed into Gephi to check the topology of the best
        overlapping graph.

        Usage: python CA_best_edge_to_GML.py asm.gkp_store asm.tigStore best.edge output.gml 

    circulization.py:

        The script provides an ad-hoc way to circularize a circular genome. The
        logic used here: (1) check whether the ends of a contig are overlapped.
        (2) If so, trimming both ends according the overlapped alignment. While
        this is useful for some case, IT IS NOT REALLY A RELIABLE WAY TO
        CIRCULARIZE A GENOME. The proper way is to look into the overlapping
        graph to ensure one actually gets a circular overlapping graph for a
        contig before trimming the ends. Overlapped ends can be due to repeats.
        Please understand the subtlety of the troubles that the repeats can
        cause during the process of assembling a genome.  One can use the
        ``CA_best_edge_to_GML.py`` script to generate a GML file to check a
        contig can be indeed unambiguously circularized visually.
    
        Uasge: python circulization.py initial_contigs.fastq 20000 /tmp circulaized_contigs.fastq

               where "20000" is the length to align at the ends of each contigs and
               "/tmp" is the location for the temporary output used to run ``blasr`` 


Most of these code was written on-fly when it was necessary. More handy scripts
may be added into the repository in the near furture if they are useful to make
the manual HGAP more scalable and easier. Some of the code might need
singificant and proper refactoring work.  While we use the code constantly and
get good results in our own computing environment, we have limited experience
and testing on more diversed computing environments.  However, everyone is
welcome to download the code, modify it for installing in your system if it
helps on automating some of these bioinformatics tasks. No official support
from PacBio for using these scripts will be provided at this moment. However,
bug reports or improvements suggestions are welcome.

A more detailed installation note and usage can be found in the file HBAR_README.rst


--------------------------------------------------------------------------------------
<script>
(function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
(i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
})(window,document,'script','//www.google-analytics.com/analytics.js','ga');
ga('create', 'UA-13166584-23', 'github.com');
ga('send', 'pageview');
</script>
