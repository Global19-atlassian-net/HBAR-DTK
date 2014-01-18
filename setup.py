from setuptools import setup, Extension, find_packages
import sys

setup(
    name = 'pbtools.hbar-dtk',
    version='0.1.5',
    author='pbiDevNet',
    author_email='pbiDevNet@pacificbiosciences.com',
    license='LICENCES',
    packages = find_packages('src'), 
    namespace_packages = ['pbtools'],
    scripts = ["src/HBAR_WF.py", 
               "src/HBAR_WF2.py", 
               "src/HBAR_WF3.py", 
               "src/generate_preassemble_reads.py",
               "src/get_preads.py",
               "src/filterM4Query.py",
               "src/query_m4_filtering.py",
               "src/h5fofn_to_fasta.py",
               "src/tig-sense_p.py",
               "src/tig-sense.py",
               "src/CA_best_edge_to_GML.py",
               "src/ContaminationTrimmer.py",
               "src/circulization.py"],
    package_dir = {'':'src'},
    data_files = [('{}/etc'.format(sys.prefix),['etc/asm.spec']),
                  ('{}/etc'.format(sys.prefix),['etc/HBAR.cfg']),
                  ('etc', ['HBAR_README.rst']),
                  ('etc', ['LICENCES'])],
    zip_safe = False,
    install_requires=[
           "pbcore >= 0.6.1",
           "pypeflow >= 0.1.1",
           "pbtools.pbdagcon >= 0.2.1",
           "pbtools.pbh5tools >= 0.75.0"
    ]
    )
