CNVera
========
 
CNVera - tool for contig copy number estimation

CNVera could be compiled with GCC 8.2.

Dependencies
-------------

CNVera has no external dependencies, although to test CNVera against Magnolia you should export paths to blastn, makeblastdb, art_illumina, quake, sga and sga/src/bin, KAligner and DistanceEst.

For now CNVera could be tested.

Building and running examples
-----------------------------

1. To build CNVera and supportive script SVSim, type:
    
        make

2. To test Magnolia vs CNVera performance on E.coli dataset. use:

        ./release_script.sh

Full manual will be added further


Authors
-------
- Dmitrii Meleshko (St. Petersburg University of the Russian Academy of Sciences, Centre for Algorithmic Biotechnologies (SPbSU))
- Vera Bushmanova (St. Petersburg State University)
- Son Pham (UCSD)



Contacts
--------
Please report any bugs or feedback to meleshko.dmitrii@gmail.com. 
