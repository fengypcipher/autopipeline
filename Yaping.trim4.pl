#!/usr/bin/perl
#use Bio::DB::Bam;
use strict;
use warnings;
#use Gtk2::Unique

# output DNA must >=35, quality cutoff 15  modifiedy by Yaping Feng on Feb 17, 2014
# debug when @means is empty 
#remove the last base before triming

my @Readsname;
my $NumFile=0;



my @line;
my $l;
my $number;
my %qtable;
#for (my $Fi=0; $Fi<$NumFile; $Fi++) {

my $inFile = $ARGV[0];
open (mF, $inFile);
my $coreName;
#print "$inFile\n";
if ($inFile =~ m/\/(.+)\..+/){
	$coreName=$1;
}
print "$coreName\n";
open (Fout, ">Trim/$coreName.trimed.fastq");
open (Fout2, ">Trim/$coreName.removedReads.txt");
my $wSize = 5;
my $cutoffM = 15;
my $indx=0;
while(<mF>){ 
  	my @means;
  	# if it is the line before the quality line
  	chomp(my $readID = $_);
  	#print "$readID\n";
  	chomp(my $dnaSeq = <mF>);
	$dnaSeq=substr($dnaSeq,0,length($dnaSeq)-1);
	my $tmpl=int length($dnaSeq)/5;
	$tmpl =  $tmpl*5;
  	$dnaSeq=substr($dnaSeq,0,$tmpl);
	my $tmp = <mF>;
  	chomp(my $l = <mF>);  # get the quality line
	$l = substr($l,0,$tmpl);
  	#if ($readID =~ m/M0i0482:9:000000000-A1DTU:1:1101:3556:10310/){
	#	print "$readID\n$dnaSeq\n$l\n";
	#}

    
	@line = split(//,$l); # divide in chars
    for(my $i = 0; $i <=$#line; $i+=$wSize){ # for each char
      my $total = 0;
 	  for (my $j=$i;$j<$i+$wSize;$j++){
 	  		$number = ord($line[$j])-33;
 	  		$total+=$number;
 	  		#if($indx%10000==0){
 	  			#print "$number..";
 	  		#}
 	  	}
 	 	my $mean = $total/$wSize; 
 	  	push(@means, $mean);
 	  
 	  	#if ($indx%10000==0){
 	 	# print Fout "$mean^^\t";
      	#}
     }
	#print join ("_",@means); 
    # print "\t$#means\n";
    
   my $indxm;
   my $indxn;
   if ($#means>0){
	$indxm = ($#line+1)/$wSize; ###the first window >cutoffM from 3'-end  
	#print "$indxm\n";  
	while ($means[$indxm-1]<=$cutoffM &&$indxm>=0){ 
      	  $indxm--;
        }
    	if ($indxm>=0){ ###all windows are trimmed if $indxm<0
	#print "$indxm******\n";
      	$indxn = 0;
     	 while ($means[$indxn]<=$cutoffM && $indxn<($#line+1)/$wSize){
      	$indxn++;
	#print "$indxn\n";
      }
    }		
   else { 
	$indxn=0;
	$indxm=0;  
    }  
    if ($indxn<$indxm){
       	my $cutPos3 = $indxn*$wSize;
    	my $cutLen = $wSize*($indxm-$indxn);
	  #print Fout "\n*****$indxm\t$cutPos3\t$cutLen\n";
     	my $trimDna = substr($dnaSeq, $cutPos3, $cutLen);
     	my $trimQs = substr($l,$cutPos3, $cutLen);
     	#print Fout "$readID\t$cutPos3\t$cutLen\n";
     	 my @rid = split(/\s+/, $readID);
     	if (length($trimDna)>=35 && $trimDna =~ /^[ATCG]+$/){ 
		print Fout "$rid[0]\n";
     	
     		#print Fout "$dnaSeq\n";
     		print Fout "$trimDna\n";
     		#print Fout "$l";
    		print Fout "$tmp";
     		print Fout "$trimQs\n";
	}
	else {
		print Fout2 "$readID too short\n";
	}	
  	}
  	else {
  		print Fout2 "$readID\n";
  	}
    }
    else {
	print Fout2 "$readID\tremoved\n";	
     }	  
}

close(Fout2);
close(Fout);
close(mF);	


#print "\n";
#while (my ($hkey, $hvalue) = each(%qtable)){
#     print "$hkey\t$hvalue\n";
#}
#print "#######\n";


sub reverseIUPAC {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/NABCDGHMNRSTUVWXYabcdghmnrstuvwxy/NTVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}
sub plasmidArray {

	

}
