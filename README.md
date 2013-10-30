clusterFAST
===========

ClusterFAST is a tool for finding translocation in next generation sequencing data devloped by Haley Abel, PhD and Eric Duncavage, MD at Washinton University in St. Louis, MO.  A manuscript detailing the use of ClusterFAST is currently under review.  For question 

ClusterFAST is a pipeline for highly specific detection of translocations from high-coverage targeted capture sequence data.  It detects translocation breakpoints with single base accuracy and provides assembled contigs for PCR validation.  ClustFAST is implemented in Java for improved interoperability and can be run from the command line via a perl script.  ClusterFAST is meant for translocation detection from specific targeted regions and requires a picard-style interval file containing the targets (e.g., ALK_MLL.txt).  

WORKFLOW:
1) First pass through ClusterFAST program identifies read pairs involving ALK or MLL.	(output: breakpoints.1.txt)

2) Divide in half split and unmapped reads whose partners map to ALK or MLL; remap using novoalign.

3) Run cf again on remapped reads to find breakpoints supported by both paired and split reads.  Options -m1 and -m2 set the minimum number of paired (default=2) and split (default=1) reads, respectively.  The defaults are chosen for maximum sensitivity (sometimes needed when tumor cellularlity is low); however, increasing these minimums (e.g., to m1=3, m2=2) will speed the pipeline.  Output: reads supported by both paired and split reads will be written to breakp
oints.2.txt.

4) Reads supporting the above breakpoints are prepared for assembly and then assembled.  For completeness, the current pipeline uses all 3 of pindel, velvet, and phrap for assembly.  All 3 work well in general; however,cphrap does not consider paired reads and so should be the least reliable.  (Any subset of methods may be used here, or other methods substituted.)

5) Assembled contigs are remapped to the human genome using blat.  The pipeline expects to be given the node and port of a blat server.

6) All assembled contigs that coincide with the breakpoints in breakpoints.2.txt are reported, along with their genomic coordinate in results.txt.

EXAMPLES:

Running ClusterFAST:

The current configuration perl script assumes that ClusterFAST is being run on a multi-node cluster.  A blat server must be running before invoking clusterFAST.  This can typically be done by:
./gfServer start [node] [port] [reference.2bit]

perl run_cf.pl -b file.bam -o outdir -f ref.fa -t translocation_region_file.txt -i ref.novoindex  -n node -p port -m1 2 -m2 1 -r 0

where:
-b is the input sequence file in sorted BAM format
-o is the output file directory
-f is the reference sequence in FASTA format (example hg19)
-t is the transloction partner_file (example ALK_MLL.txt)
-i is the binary novoalign indexed reference file
-n server node where blat server is running
-p port on which blast server is running
-m1 is number of discordant pairs required
-m2 is number of split reads required
-r remove temp file (T/F)
