#! /usr/bin/perl
use warnings;
use threads; use Thread::Queue; #use threads::shared;  # not being used but allows variable to be shared by all threads e.g.  $test :shared = "Text shared by all threads";
use strict;
use File::Copy; 
use Getopt::Std;

# check that the -s and -p flags were set
our($opt_r, $opt_p, $opt_e,$opt_s,$opt_m); 
getopts('r:p:es:m:');
unless($opt_r & $opt_p){die "\n\nUSAGE: perl RNAseqPipe.pl -r readfiles.txt -p params.txt -s samplesheet.txt -m merged.gtf\n\n"}

my $genDir; 
my $genome; 
my $gtf;
my $samplesheet=$opt_s;
my $outDir="./";
my $process_q = Thread::Queue -> new();
my %program; 
my %params; 
my %overwrite;
my @bamlist;

#read in parameters
open PARAM, "<$opt_p" or die $!;
while(<PARAM>){
	chomp;
	if(/^#/){next}
	if(/RefGenomeDirectory\s+(\S+)/){$genDir=$1}
	if(/RefGenomeFasta\s+(\S+)/){$genome=$1}
	if(/RefGenomeGTF\s+(\S+)/){$gtf=$1}
	if(/OutputDirectory\s+(\S+)/){$outDir=$1}
	
	if(/(\S+)\t([^\t]+)\t(\S+)/){
		$program{$1}=$2;
		$overwrite{$1}=$3;
		$params{$1}=$';
	}
}
close PARAM;

#### Edit GTF to include gene_id annotation for exons
open GTF, "<$gtf" or die $!;
my $gtf2=$gtf;
unless($gtf2 =~ /edited\.g[tf]f$/){  #only do this if "edited" is not in the name of the gff file already
	$gtf2 =~ s/\.g[tf]f/\.edited\.gff/;
	open GTF2, ">$gtf2";

	print "\nAdding gene_id to exon annotations in $gtf\n";
	while(<GTF>){
	chomp;
		if (/^\S+\t\S+\texon/){
			m/Parent=([^\.]+)/;
			print GTF2 "$_;gene_id=$1\n";
		}else{print GTF2 "$_\n";}
	}
	close GTF;close GTF2;
}
print "\nGenome Directory:	$genDir\n";
print "Genome FASTA:	$genome\n";
print "Annotations file:	$gtf2\n";
print "cufflinks output will be saved to:	$outDir\n\n";

unless($opt_e){system "mkdir -p $outDir"}	

# format genome if not already done
unless((-e "$genDir/SAindex") && ($overwrite{genomeGen} eq "F")){
	my $generateGenome = "STAR --runMode genomeGenerate --genomeFastaFiles $genome --genomeDir $genDir --sjdbGTFfile $gtf $params{genomeGen}";  #--sjdbOverhang should be Readlength -1 than
	print "\nEXECUTING: $generateGenome\n";
	unless($opt_e){system $generateGenome}
}

open READ, "<$opt_r" or die $!; #read in readfiles
while(<READ>){
	chomp;
	if(/^sample\tname/){next}elsif(/^#/){next}
	elsif(/(\S+)\t(\S+)\t(\S+)\t(\S+\.fastq(\.gz)?)/){
		my $sample=$1;
		my $name=$2;
		my $group=$3;
		my $readfile=$4;
		my $stem= "$1_$3_$2";
		unless($opt_e){system "mkdir -p $stem"}
		
		push @bamlist, "./$stem/$stem.Aligned.sortedByCoord.out.bam";  # needed later for cuffmerge
				
####### Trimmomatic 
	unless((-e "./$stem/$stem.trimmed.fq.gz") && ($overwrite{trim} eq "F")){		
		my $trim= "$program{trim} $readfile ./$stem/$stem.trimmed.fq.gz $params{trim}";
		print "\nEXECUTING: $trim\n";
		unless($opt_e){system $trim}		
	}

####### STAR aligner 
	unless((-e "./$stem/$stem.Aligned.sortedByCoord.out.bam") && ($overwrite{align} eq "F")){
		
		my $align= "$program{align} $params{align} --genomeDir $genDir --outFileNamePrefix ./$stem/$stem. --readFilesIn ./$stem/$stem.trimmed.fq.gz";
		#my $align= "hisat2 -t -p 30 -x $genDir/Pf3D7 --dta-cufflinks -U $_ | samtools sort -\@30 -O BAM --reference $genome -o ./$`/$`_hisat2.bam"; #hisat2 aligner  3 min
		#my $align= "tophat -p 30 -G $gtf -o ./$` $genome $_"; #tophat aligner, 63 min
		
		print "\nEXECUTING: $align\n";
		unless($opt_e){system $align}

		my $index="samtools index ./$stem/$stem.Aligned.sortedByCoord.out.bam";
		print "\nEXECUTING: $index\n";
		unless($opt_e){system $index}
		}
	}
	else{print "\nMALFORMATED LINE in $opt_r: $_\n"}
}
close READ;	


####### Cufflinks
sub worker{   # this gets executed to run cufflinks for each item in the $process_q, 
  while ( my $stem = $process_q -> dequeue() ){
	my $cufflinks="$program{cufflinks} $params{cufflinks} -G $gtf --output-dir ./$stem ./$stem/$stem.Aligned.sortedByCoord.out.bam";
	print "EXECUTING: cufflinks ";
	
	unless($opt_e){
		system $cufflinks;
		print "Cufflinks of $stem is done!\n";
		threads->exit();
	}
  }
}

foreach (@bamlist){
 	/\.\/([^\/]+)\/[^\.]+.Aligned.sortedByCoord.out.bam/;
 	my $stem=$1; 
 	unless((-e "./$stem/genes.fpkm_tracking") && ($overwrite{cufflinks} eq "F")){
 		$process_q -> enqueue ( $stem );	    # adds each $stem to the queue
 	}
} 
print "\n";

$process_q -> end();
foreach(@bamlist){threads -> create ( \&worker )}	# creates number of threads equal to length of @bamfile. each thread executes 1 item of the queue
foreach my $thr (threads -> list()){$thr -> join()} # waits for all threads to have finished

########### Cuffmerge
 unless(((-e "$outDir/merged.gtf") && ($overwrite{cuffmerge} eq "F")) || $opt_m){

 	if(-e "$outDir/assemblies.txt"){system "rm -f $outDir/assemblies.txt"} #removes existing assemblies.txt file
 	foreach (@bamlist){  # recreate assemblies.txt
 		/\.\/([^\/]+)\/[^\.]+.Aligned.sortedByCoord.out.bam/;
 		system "echo ./$1/transcripts.gtf >>$outDir/assemblies.txt\n";
 	}
 	
 	my $cuffmerge ="$program{cuffmerge} $params{cuffmerge} -g $gtf -s $genome -o $outDir $outDir/assemblies.txt ";
 	print "\nEXECUTING: $cuffmerge\n";
 	unless($opt_e){system $cuffmerge}
}
 
########### Cuffdiff
 unless((-e "$outDir/genes.fpkm_tracking") && ($overwrite{cuffdiff} eq "F")){
 	my $cuffdiff; my $merged;
 	if($opt_m){$merged=$opt_m}else{$merged="$outDir/merged.gtf"}	# without a sample sheet, just provide the merged.gtf file.
	if($opt_s){
		$cuffdiff ="$program{cuffdiff} --use-sample-sheet -o $outDir $params{cuffdiff} $merged $opt_s"} #if using a sample sheet to indicate groups
	else{
		$cuffdiff ="$program{cuffdiff} -o $outDir $params{cuffdiff} $merged @bamlist"
	}
	print "\nEXECUTING: $cuffdiff\n";
	unless($opt_e){system $cuffdiff}
}






