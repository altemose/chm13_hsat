#Assembly_HSat2and3_v3.pl
#Nicolas Altemose 2022
#purpose: given a fasta file with a reference sequence and two input files with HSat2 and HSat3 specific 24-mers,
#outputs all regions in the reference likely to be HSat2 or HSat3, along with their strand orientation,
#as a BED 9 file
#note: merges all adjacent regions within 5 kb (same strand, same type)
#HSat2_kmers.txt and HSat3_kmers.txt must be in the same directory in which this script is executed
#the kmers in these files were defined using HSat2/3 HuRef reads identified in Altemose et al. PLoS Comp Bio 2014

use warnings;
use strict;
my $tic = time; #store system time for runtime calculation
print "\n\n";

my $usage = "perl Assembly_HSat2and3_v3.pl /path/to/reference.fasta";

#read in reference file name from command line
my $reffile = "";
if(defined $ARGV[0]){
	$reffile = $ARGV[0];
	chomp($reffile);
}else{
	die "no reference fasta file specified!\n\nusage: $usage\n\n";
}

#create output file name
my $outfile = "";
if($reffile=~/.+\/(.+)\..+?$/){
	$outfile = "$1.HSat2and3_Regions.bed";
}
else{
	$outfile = "$reffile.HSat2and3_Regions.bed";
}

#read in kmers and store in hashes (dictionaries)
print "reading in kmers\n";
my %hsat2kmers; #store hsat2-specific kmers here
my %hsat3kmers; #store hsat3-specific kmers here
my %allkmersF; #store all kmers from hsat2/3 here
my %allkmersR; #store reverse complements of all kmers from hsat2/3 here

open(IN2, "HSat2_kmers.txt") or die "ERROR: Could not open HSat2_kmers.txt. Put file in this directory.\n";
while(my $line = <IN2>){
	chomp($line);
	$line=~/^(\S+)\t(\S+)\t(\S+)$/ or die "ERROR $line\n";
	if($3<10){ #only keep HSat2 kmers found fewer than 10 times in HSat3 reads
		$hsat2kmers{$1}=$2;
	}
	$allkmersF{$1}=0;
	my $rc = reverse($1);
	$rc =~ tr/ACGT/TGCA/;
	$allkmersR{$rc}=0;
}
close IN2;

open(IN3, "HSat3_kmers.txt") or die "ERROR: Could not open HSat3_kmers.txt. Put file in this directory.\n";
while(my $line = <IN3>){
	chomp($line);
	$line=~/^(\S+)\t(\S+)\t(\S+)$/ or die "ERROR $line\n";
	if($3<10){#only keep HSat3 kmers found fewer than 10 times in HSat2 reads
		$hsat3kmers{$1}=$2;
	}
	$allkmersF{$1}=0;
	my $rc = reverse($1);
	$rc =~ tr/ACGT/TGCA/;
	$allkmersR{$rc}=0;
}
close IN3;

#read in reference fasta one entry at a time and scan across sequence counting k-mer matches in 1-kb bins
my $totlen=0;
print "counting kmers across genome\n";
print "\tchr\ttotal_length\n";

open(OUTA,'>tempfileA.bed');

if(1==1){#use true conditional to restrict scope of contained variables
local $/ = ">"; #change the record separator (locally) so that one fasta entry is read at a time, instead of one line at a time
open(FASTA,$reffile) or die "ERROR: Could not open reference file $reffile\n";	
my $junk = <FASTA>; #remove first > from beginning of file
while(my $frecord = <FASTA>){
	chomp($frecord);
	$frecord=~/^(\S+).*?\n(.*)$/s or die "ERROR: problem in $reffile parsing this:\n$frecord\n";
			#note: $1 can be modified to change which parts of header are kept
	my $chr = $1;
	my $chrseq = $2;
	$chrseq=~s/[\n\s]//g; #remove newlines and spaces from chromosome sequence
	$chrseq=uc($chrseq); #make DNA sequence all uppercase

	my $chrlen = length($chrseq);
	$totlen+=$chrlen;
	print "\t$chr\t$chrlen\n";

	for(my $n=0;$n<=length($chrseq)-1000;$n+=1000){ #loop through sequence in 1 kb non-overlapping bins
		my $stop = $n+1000;
		my $seq = substr($chrseq,$n,1000);
		my $rc = reverse($seq);
		$rc =~ tr/ACGT/TGCA/;
		
		#count CATTC on both strands. if not more than 15 found on either strand, skip this bin
		#this serves as a pre-screening step (threshold validated on HuRef reads from Altemose et al. 2014)
		my $Fcattc = $seq =~ s/CATTC/CATTC/g;
		my $Rcattc = $rc =~ s/CATTC/CATTC/g;
		next unless($Fcattc>=15 || $Rcattc>=15);
		
		#define forward orientation as that with the greater number of CATTCs
		my $strand='+';
		if($Rcattc>$Fcattc){
			$seq = $rc;
			$strand='-';
		}
		
		#begin scanning through 24-mers in bin and count matches to HSat2/3-specific sets
		my $kmercount2=0;
		my $kmercount3=0;
		for(my $m=0;$m<=1000-24;$m++){
			my $mer = substr($seq,$m,24);
			if(exists $hsat2kmers{$mer}){
				$kmercount2++;
			}elsif(exists $hsat3kmers{$mer}){
				$kmercount3++;
			}
			if($kmercount3>100){ #any bin with >100 HSat3 kmer matches gets called as HSat3 (threshold validated on HuRef reads from Altemose et al. 2014)
				print OUTA "$chr\t$n\t$stop\tHSat3\t0\t$strand\n"; #output bin coordinates in BED6 format
				last; #once the threshold is reached, stop counting and move to the next bin to save time
			}elsif($kmercount2>150){ #otherwise any bin with >150 HSat2 kmer matches gets called as HSat2 (threshold validated on HuRef reads from Altemose et al. 2014)
				print OUTA "$chr\t$n\t$stop\tHSat2\t0\t$strand\n"; #output bin coordinates in BED6 format
				last; #once the threshold is reached, stop counting and move to the next bin to save time
			}
		}#closes for m
	}#closes for n
}#closes while
}#closes if TRUE
close OUTA;
%hsat2kmers = ();
%hsat3kmers = ();


#merge adjacent bins on same strand and of same type within 5 kb
open(IN,"tempfileA.bed");
open(OUTC,'>tempfileC.bed');

my $line1 = <IN>;
chomp($line1);
$line1=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+0\s+(\S+)/ or die $line1;
	
	my $prevchr = $1;
	my $prevstart = $2;
	my $prevstop = $3;
	my $prevtype = $4;
	my $prevstrand = $5;

while(my $line = <IN>){
	chomp($line);
	$line=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+0\s+(\S+)/ or die $line;
	unless($1 eq $prevchr && $2 - $prevstop <= 5000 && $4 eq $prevtype && $5 eq $prevstrand){
		print OUTC "$prevchr\t$prevstart\t$prevstop\t$prevtype\t0\t$prevstrand\n";
		$prevstart = $2;
	}
	$prevchr = $1;
	$prevstop = $3;
	$prevtype = $4;
	$prevstrand = $5;
}
close IN;
print OUTC "$prevchr\t$prevstart\t$prevstop\t$prevtype\t0\t$prevstrand\n";
close OUTC;


#refine boundaries of HSat2/3 regions by scanning for the furthest reaching 24-mer within 1-kb to either side of each region's boundaries
print "refining boundaries\n";
open(INC,'tempfileC.bed');
#record all region coordinates for each chromosome in a hash
my %memory;
while(my $line=<INC>){
	chomp($line);
	if($line=~/^(\S+)\t(.+)$/){
		$memory{$1}.=",$2";
	}
}
close INC;

#scan through reference fasta file again one chromosome at a time, as above
open(OUTD,'>tempfileD.bed');
if(1==1){#use true conditional to restrict scope of contained variables
local $/ = ">"; #change the record separator (locally), as above
open(FASTA,$reffile) or die "ERROR: Could not open $reffile\n";	
my $junk = <FASTA>; #remove first > from beginning of file
while(my $frecord = <FASTA>){
	chomp($frecord);
	$frecord=~/^(\S+).*?\n(.*)$/s or die "ERROR: problem in $reffile parsing this:\n$frecord\n";
			#note: $1 can be modified to change which parts of header are kept
	my $chr = $1;
	my $chrseq=$2;
	$chrseq=~s/[\n\s]//g;
	$chrseq=uc($chrseq);
	my $chrlen = length($chrseq);
	
	if(exists $memory{$chr}){
		my $info = $memory{$chr};
		$info=~s/^\,//;
		my @lines = split(',',$info);
		foreach my $line(@lines){ #loop through each reported region for this chromosome
			if($line=~/^(\S+)\t(\S+)\t(\S+)\t\S+\t(\S+)/){
				my $start = $1-1000;
				if($start<0){
					$start=0;
				}
				my $stop = $2+1000;
				if($stop>=$chrlen){
					$stop=$chrlen-1;
				}
				my $type = $3;
				my $strand = $4;
				my $sub1 = substr($chrseq,$start,2000); #get 2-kb window centered at start coordinate of region and identify the leftmost 24-mer coordinate on the same strand
				if($strand eq '+'){
					for(my $n=0;$n<2000-24;$n++){
						my $mer = substr($sub1,$n,24);
						if(exists $allkmersF{$mer}){
							$start = $start+$n;
							last;
						}
					}
					my $sub2 = reverse(substr($chrseq,$stop-2000,2000)); #get 2-kb window centered at stop coordinate of region and identify the rightmost 24-mer coordinate on the same strand
					$sub2 =~ tr/ACGT/TGCA/;
					for(my $n=0;$n<2000-24;$n++){
						my $mer = substr($sub2,$n,24);
						if(exists $allkmersR{$mer}){
							$stop = $stop-$n;
							last;
						}
					}
				}else{
					for(my $n=0;$n<2000-24;$n++){
						my $mer = substr($sub1,$n,24);
						if(exists $allkmersR{$mer}){
							$start = $start+$n;
							last;
						}
					}
					my $sub2 = reverse(substr($chrseq,$stop-2000,2000));
					$sub2 =~ tr/ACGT/TGCA/;
					for(my $n=0;$n<2000-24;$n++){
						my $mer = substr($sub2,$n,24);
						if(exists $allkmersF{$mer}){
							$stop = $stop-$n;
							last;
						}
					}
				}
				my $len = $stop-$start;
				#if($len>=1000){
					if($type eq "HSat2"){
						print OUTD "$chr\t$start\t$stop\t$type\t0\t$strand\t$start\t$stop\t51,51,102\n";
					}elsif($type eq "HSat3"){
						print OUTD "$chr\t$start\t$stop\t$type\t0\t$strand\t$start\t$stop\t120,161,187\n";
					}
				#}
			}#closes if
		}#closes foreach
	}#closes if exists 
}#closes while
}#closes if TRUE
close OUTD;


#resolve any overlaps between regions by setting the boundary at the midpoint of the overlap
print "resolving overlaps\n";
open(IN, 'tempfileD.bed');
open(OUT,'>'.$outfile);

my $line2 = <IN>;
chomp($line2);
$line2=~/^(\S+)\t(\S+)\t(\S+)\t(\S+)\t\S+\t(\S+)\t\S+\t\S+\t\S+$/ or die "ERROR:$line2\n";
$prevstop = $3;
$prevstart = $2;
$prevchr = $1;
$prevtype = $4;
$prevstrand = $5;

while(my $line = <IN>){
	chomp($line);
	$line=~/^(\S+)\t(\S+)\t(\S+)\t(\S+)\t\S+\t(\S+)\t\S+\t\S+\t\S+$/   or die "ERROR:$line\n";
	my $chr = $1;
	my $start = $2;
	my $stop = $3;
	my $type = $4;
	my $strand = $5;
	
	next if($chr eq $prevchr && $prevstart<=$start && $prevstop>=$stop); #skip entries fully contained in the previous entry (in effect, this removes small inversions smaller than the 5 kb merging threshold)
	my $printstop=0;
	if($chr eq $prevchr && $prevstop>$start){
		$printstop = $prevstop-int(($prevstop-$start)/2);
		if($prevtype eq "HSat2"){
			print OUT "$prevchr\t$prevstart\t$printstop\t$prevtype\t0\t$prevstrand\t$prevstart\t$printstop\t51,51,102\n";
		}else{
			print OUT "$prevchr\t$prevstart\t$printstop\t$prevtype\t0\t$prevstrand\t$prevstart\t$printstop\t120,161,187\n";
		}
		$prevstart = $printstop;
	}else{
		$printstop = $prevstop;
		if($prevtype eq "HSat2"){
			print OUT "$prevchr\t$prevstart\t$printstop\t$prevtype\t0\t$prevstrand\t$prevstart\t$printstop\t51,51,102\n";
		}else{
			print OUT "$prevchr\t$prevstart\t$printstop\t$prevtype\t0\t$prevstrand\t$prevstart\t$printstop\t120,161,187\n";
		}
		$prevstart=$start;
	}
	$prevstop = $stop;
	$prevchr = $chr;	
	$prevtype = $type;
	$prevstrand = $strand;
}
		#print last line
		if($prevtype eq "HSat2"){
			print OUT "$prevchr\t$prevstart\t$prevstop\t$prevtype\t0\t$prevstrand\t$prevstart\t$prevstop\t51,51,102\n";
		}else{
			print OUT "$prevchr\t$prevstart\t$prevstop\t$prevtype\t0\t$prevstrand\t$prevstart\t$prevstop\t120,161,187\n";
		}
		
close OUT;

system("rm -f tempfile*"); #remove intermediate files
print "printed output to $outfile\n";


#calculate and display runtime
my $toc = time;
my $elapsed = $toc-$tic;
printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($elapsed / 3600), int(($elapsed % 3600) / 60), int($elapsed % 60));
