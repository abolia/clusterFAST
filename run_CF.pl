#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;

my ($inbam, $outdir, $targets, $minct1, $minct2, $dist, $node, $port, $novoindex, $rmtmp, $FASTA);

#$targets="/home/comp/gtac/habel/translocations/ALK_MLL.txt";
#$novoindex="/srv/cgs/data/gtac/hg19.novoindex";
$minct1=2;
$minct2=1;
$dist=50000;
$rmtmp=1;


#my $FASTA="/srv/cgs/data/gtac/hg19/hg19.fa";
my $SAMTOOLS="/srv/cgs/local/samtools/latest/bin/samtools";
my $NOVOALIGN="/srv/cgs/local/novocraft/latest/novoalign";
my $PHRAP="/srv/cgs/local/phrap/latest/phrap";
my $PINDEL="/srv/cgs/local/pindel/latest/bin/pindel";
my $SAM2PINDEL="/home/comp/gtac/duncavage/sam2pindel";
my $VELVET="/home/comp/gtac/habel/bin/velvet/velvet_1.2.07";
my $FASTX="/srv/cgs/local/fastx/latest/bin";
my $BLAT="/srv/cgs/local/blat/latest";

#my $scriptspath="/home/comp/gtac/habel/bin/FindTrans/scripts0604";
#my $cffile="/home/comp/gtac/habel/bin/FindTrans/cf0520/cf.jar";
my $scriptspath="/home/comp/gtac/habel/tx_pipeline/scripts0604";
#my $cffile="/home/comp/gtac/habel/tx_pipeline/cf.jar";
#my $cffile="/home/comp/gtac/habel/bin/FindTrans/cf0806/cf0806.jar";
my $cffile="/home/comp/gtac/habel/bin/FindTrans/cf1010/cf1010.jar"; # FindTrans0311 ";
my $partnerfile="/home/comp/gtac/habel/tx_pipeline/known_partners.txt";
my $help='';


GetOptions(
    "b|inbam=s"=> \$inbam,
    "o|outdir=s"=> \$outdir,
    "t|targets=s"=> \$targets,
    "f|fasta=s"=> \$FASTA,
    "m1|minct1=i"=> \$minct1,
    "m2|minct2=i"=> \$minct2,
    "d|distance=i"=>\$dist,
    "n|node=s"=> \$node,
    "p|port=s"=> \$port,
    "i|novoindex=s"=> \$novoindex,
    "r|rmtmp=i"=> \$rmtmp,
    "h|help"=>\$help,
    ) or die "Usage:  run_cf_1010.pl -b inbam -o outdir -f fasta -t targets -d distance -n node -p port -i novoindex -r rmtmp -m1 min_pairs -m2 min_splits\n";

if($help) {
	print "Usage:  run_cf_1010.pl -b inbam -o outdir -f fasta -t targets  -d distance -n node -p port -i novoindex -r rmtmp -m1 min_pairs -m2 min_splits\n";
	print "targets=picard style interval file containing the regions to search for SV\n";
	print "fasta=reference fasta file\nnovoindex=reference index for novoalign\n";
	print "node=blat server node\nport=blat server port\n";
	print "distance=max distance from target to search (default=50000)\n";
	print "min_pairs=minimum number of read pairs supporting breakpoint (default=2)\nmin_splits=minimum number of split reads supporting breakpoint (default=1)\n";
	print "rmtmp=remove temp directory (default=1, true)\n";
	exit(1);
}

if(! ( defined $inbam && defined $targets && defined $FASTA && defined $outdir && defined $novoindex && defined $node && defined $port)) {
	print STDERR "Usage:  run_cf_0806.pl -b inbam -o outdir -f fasta -t targets  -d distance -n node -p port -i novoindex -r rmtmp -m1 min_pairs -m2 min_splits\n";
	exit(1);
}
elsif (! (-e $inbam && -e $targets && -e $FASTA && -e $novoindex)) {
	print STDERR "Critical file missing.\n";
	exit(1);
}


make_empty_output();

my $cmd="java -Xmx8g -Xms6g -classpath $cffile $inbam $targets $minct1 $dist $outdir 1";
print STDOUT "$cmd\n";
system($cmd);

chdir("${outdir}/temp");

$cmd="perl ${scriptspath}/splitunmapped40_0314.pl < final.1.sam filter1.fq filter2.fq FQ";
print STDOUT "$cmd\n";
system($cmd);

if(-z "filter1.fq") {
    print STDOUT "No breakpoints found.\n";
    chdir("$outdir");
    clean_up();
    exit(1);
}

$cmd="$NOVOALIGN -o SAM -i 230 140 -r all -e 999 -c2 -d $novoindex -F STDFQ -f filter1.fq filter2.fq > novoout.2.sam";
print STDOUT "$cmd\n";
system($cmd);

$cmd="java -Xmx6g  -Xms4g -classpath $cffile novoout.2.sam $targets $minct2 $dist $outdir 2";
print STDOUT "$cmd\n";
system($cmd);

print STDOUT "Filtering breakpoints...\n";
filter_breakpoints();
print STDOUT "Done\n";

print STDOUT "Formatting input for assembly...\n";
reformat_sams();
print STDOUT "Done\n";

print STDOUT "Assembling contigs...\n";
assemble_contigs();
print STDOUT "Done\n";

print STDOUT "Preparing output files...\n";
make_output();
print STDOUT "Done\n";

system("mv pindel_out_* contigs.txt $outdir");
chdir("$outdir");
clean_up();




#########################


sub make_empty_output {
    system("touch ${outdir}/breakpoints.1.txt");
    system("touch ${outdir}/breakpoints.2.txt");
    system("touch ${outdir}/breakpoints.final.txt");
    system("touch ${outdir}/breakpoints.1end.flt.txt");
    system("touch ${outdir}/breakpoints.2end.flt.txt");
    system("touch ${outdir}/contigs.txt");
}

sub filter_breakpoints {
    
    system("R --vanilla < ${scriptspath}/consolidate_hits0520.R --args ${outdir}/breakpoints.1.txt ${outdir}/breakpoints.2.txt ${outdir}/breakpoints.final.txt");
    system("perl ${scriptspath}/filter_by_partner.0719.pl $partnerfile  ${outdir}/breakpoints.final.txt 1 cf 99 > ${outdir}/breakpoints.1end.flt.txt");
    system("perl ${scriptspath}/filter_by_partner.0719.pl $partnerfile  ${outdir}/breakpoints.final.txt 2 cf 99 > ${outdir}/breakpoints.2end.flt.txt");
    if(-z "${outdir}/breakpoints.1end.flt.txt") {
        print "No breakpoints found.\n";
        exit(1);
    }
}

sub reformat_sams {
    system("samtools view -bS final.singleendmatenm.sam | samtools sort - final.1.sorted");
	system("samtools view -bS 2ends.final.2.sam | samtools sort - final.2.sorted");
	system("samtools view -bS final.singleendnm.sam | samtools sort - final.3.sorted");
	system("rm temp.bam");
	system("samtools merge temp.bam final.1.sorted.bam final.2.sorted.bam final.3.sorted.bam");
	system("samtools sort -n -m 10000000 temp.bam final");
	if(-z "final.bam") {
		print "No reads for assembly.\n";
		exit(1);
	}
	system("samtools view final.bam | sort -k 1,1 -k 2,2n -u > final.unique.sam");
    system("perl ${scriptspath}/split_for_velvet.pl < final.unique.sam pair1.fq pair2.fq single.fq");
	#system("cat pair1.fq pair2.fq single.fq | ${FASTX}/fastq_to_fasta -Q33 > phrap.fa");
}

sub assemble_contigs {
	system("$SAM2PINDEL  final.unique.sam  final.unique.sam.pindel 270 H 0");
	system("$PINDEL -c ALL -f $FASTA -p final.unique.sam.pindel -o pindel_out > /dev/null");
	system("${VELVET}/velveth vtest 31 -fastq -shortPaired -separate pair1.fq pair2.fq -short single.fq > /dev/null");
	system("${VELVET}/velvetg vtest -cov_cutoff 3 -ins_length 300 -exp_cov 100  -ins_length_sd 60 > /dev/null");
	#system("$PHRAP phrap.fa -minmatch 20 > /dev/null");
	system(" perl ${scriptspath}/pindel2fasta_0520.pl < pindel_out_BP > pindel.fa");
	system(" perl ${scriptspath}/pindel2fasta_0520.pl < pindel_out_LI >> pindel.fa");
}

sub make_output {
	system("${BLAT}/gfServer status $node $port ");
	#system("${BLAT}/gfClient $node $port / phrap.fa.contigs phrap.psl -minScore=2 -minIdentity=90");
	system("${BLAT}/gfClient $node $port / vtest/contigs.fa velvet.psl -minScore=2 -minIdentity=90");
	system("${BLAT}/gfClient $node $port / pindel.fa pindel.psl -minScore=2 -minIdentity=90");
	#system("perl ${scriptspath}/filter_contigs.pl phrap.psl  {$outdir}/breakpoints.2.txt phrap.fa.contigs phrap > phrap.txt");
	system("perl ${scriptspath}/filter_contigs.pl velvet.psl ${outdir}/breakpoints.2.txt vtest/contigs.fa velvet > velvet.txt");
	system("perl ${scriptspath}/filter_contigs.pl pindel.psl  ${outdir}/breakpoints.2.txt pindel.fa pindel > pindel.txt");
	#system("cat pindel.txt velvet.txt phrap.txt > results.txt");
	system("cat pindel.txt velvet.txt  > contigs.txt");
}
    
sub clean_up {
    
    if ($rmtmp) {
        system("rm -R temp");
    }
}
