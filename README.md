clusterFAST
===========

ClusterFAST is a tool for finding translocation in next generation sequencing data devloped by Haley Abel, PhD and Eric Duncavage, MD at Washinton University in St. Louis, MO.

ClusterFAST is a pipeline for highly specific detection of translocations from high-coverage targeted capture sequence data.  It detects translocation breakpoints with single base accuracy and provides assembled contigs for PCR validation.

ClusterFAST is meant for tx detection from specific targeted regions and requires a picard-style interval file containing the targets (e.g., ALK_MLL.txt).  

WORKFLOW:
1) First pass through cf program identifies read pairs involving ALK or MLL.	(output: breakpoints.1.txt)

2) Divide in half split and unmapped reads whose partners map to ALK or MLL; remap using novoalign.

3) Run cf again on remapped reads to find breakpoints supported by both paired and split reads.  Options -m1 and -m2 set the minimum number of paired (default=2) and split (default=1) reads, respectively.  The defaults are chosen for maximum sensitivity (sometimes needed when tumor cellularlity is low); however, increasing these minimums (e.g., to m1=3, m2=2) will speed the pipeline.  Output: reads supported by both paired and split reads will be written to breakp
oints.2.txt.

4) Reads supporting the above breakpoints are prepared for assembly and then assembled.  For completeness, the current pipeline uses all 3 of pindel, velvet, and phrap for assembly.  All 3 work well in general; however,cphrap does not consider paired reads and so should be the least reliable.  (Any subset of methods may be used here, or other methods substituted.)

5) Assembled contigs are remapped to the human genome using blat.  The pipeline expects to be given the node and port of a blat server.

6) All assembled contigs that coincide with the breakpoints in breakpoints.2.txt are reported, along with their genomic coordinate in results.txt.
