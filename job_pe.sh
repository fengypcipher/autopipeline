# This is version2 of the automatic pipeline for short read alignment updated on Nov 1st 2013
# This version is updated to to compare the QC filter fastq files
# it first separate the reads as paired(common) and single end(orphan reads (uniq) and merge the single end reads at the end of paired read.
##in this version, don't cat *common-out and *unique-out into the same *cmp.fastq. Only use *common-out for pair end mapping

#!/bin/sh
### script to run a serial job using one core on htc using queue

### Set the number of cores (cpus) and memory that will be used for this job
###PBS -l nodes=1:ppn=1,pmem=8gb
#PBS -l nodes=1:ppn=4,pmem=4gb
###PBS -l nodes=1:ppn=12,pmem=4gb

### Specify "wallclock time" required for this job, hhh:mm:ss
#PBS -l walltime=100:00:00

### Declares the standard error stream of the job will be merged with the standard output stream of the job
#PBS -j oe 

### Request email when job begins and ends
### PBS -m bea

### Specify email address to use for notification.
### PBS -M ypfeng@waksman.rutgers.edu

### cd: set directory for job execution
cd $mycwd
PRE=/ingens/solid/analysis/yaping/pipeline2/Script
echo working dir for pinepeline is $mycwd
filename1=`echo $fname1|cut -d '.' -f 1`
filename2=`echo $fname2|cut -d '.' -f 1`
sample=`echo $filename1|cut -d '_' -f 1`
if [ ! -d "logs" ]; then        
    mkdir "logs"
fi	
 echo ------------------------------------------------------------------------------------ >> logs/mappe_$filename1.log;     
 echo Started alignment for $fname1 and $fname2  at  `date` >> logs/mappe_$filename1.log;
 
 echo Comparing fastq files if QC filtering is done at  `date` >> logs/mappe_$filename1.log;
	 
		if [ "$adapter" == "yes"  ]; then
				if [ "$tr" == "yes"  ]; then
					### both removed adapters and trimmed reads for mapping
				#	bwa mem  -t 4 $para genome/$reffname  Trim/$filename1.CutAdaptor.trimed.fastq Trim/$filename2.CutAdaptor.trimed.fastq > bwaout/$filename1.pe.sam
					perl $PRE/cmpfastq.pl Trim/$filename1.CutAdaptor.trimed.fastq Trim/$filename2.CutAdaptor.trimed.fastq 
					#cat Trim/$filename1.CutAdaptor.trimed.fastq-common.out Trim/$filename1.CutAdaptor.trimed.fastq-unique.out > Trim/$filename1.CutAdaptor.trimed.cmp.fastq
					cat Trim/$filename1.CutAdaptor.trimed.fastq-common.out > Trim/$filename1.CutAdaptor.trimed.cmp.fastq
					#cat Trim/$filename2.CutAdaptor.trimed.fastq-common.out Trim/$filename2.CutAdaptor.trimed.fastq-unique.out > Trim/$filename2.CutAdaptor.trimed.cmp.fastq
					cat Trim/$filename2.CutAdaptor.trimed.fastq-common.out > Trim/$filename2.CutAdaptor.trimed.cmp.fastq
					echo Comparing fastq for adapter removed and trimmed $filename1  and $filename2 at is done `date` >> logs/mappe_$filename1.log;
				else
					###  only adapter removed reads for mapping
				#	bwa mem  -t 4 $para genome/$reffname  CutAdaptor/$filename1.CutAdaptor.fastq CutAdaptor/$filename2.Cutadaptor.fastq > bwaout/$filename1.pe.sam
					perl $PRE/cmpfastq.pl CutAdaptor/$filename1.CutAdaptor.fastq CutAdaptor/$filename2.Cutadaptor.fastq 
					#cat CutAdaptor/$filename1.CutAdaptor.fastq-common.out CutAdaptor/$filename1.CutAdaptor.fastq-unique.out > CutAdaptor/$filename1.CutAdaptor.cmp.fastq
					mv CutAdaptor/$filename1.CutAdaptor.fastq-common.out > CutAdaptor/$filename1.CutAdaptor.cmp.fastq
					#cat CutAdaptor/$filename2.CutAdaptor.fastq-common.out CutAdaptor/$filename2.CutAdaptor.fastq-unique.out > CutAdaptor/$filename2.CutAdaptor.cmp.fastq
					mv CutAdaptor/$filename2.CutAdaptor.fastq-common.out > CutAdaptor/$filename2.CutAdaptor.cmp.fastq
					echo Comparing fastq for just adapter removed $filename1 and $filename2  is done at  `date` >> logs/mappe_$filename1.log;
				fi
		else
				### adapters not removed
				if [ "$tr" == "yes"  ]; then
					### only trimmed reads for mapping
				#	bwa mem  -t 4 $para genome/$reffname  Trim/$filename1.trimed.fastq Trim/$filename2.trimed.fastq > bwaout/$filename1.pe.sam
					perl $PRE/cmpfastq.pl  Trim/$filename1.trimed.fastq Trim/$filename2.trimed.fastq 
					#cat Trim/$filename1.trimed.fastq-common.out Trim/$filename1.trimed.fastq-unique.out > Trim/$filename1.trimed.cmp.fastq			
					cat Trim/$filename1.trimed.fastq-common.out > Trim/$filename1.trimed.cmp.fastq
					#cat Trim/$filename2.trimed.fastq-common.out Trim/$filename2.trimed.fastq-unique.out > Trim/$filename2.trimed.cmp.fastq			
					cat Trim/$filename2.trimed.fastq-common.out > Trim/$filename2.trimed.cmp.fastq
					echo Comparing for just trimmed $filename1  and $filename2 is done at  `date`
				else
					### only raw reads used for mapping	
				#	bwa mem  -t 4 $para genome/$reffname  rawdata/$fname1 rawdata/$fname2 > bwaout/$filename1.pe.sam
					echo As no QC processing done on $filename1  and $filename2  comparing fastq not required at  `date` >> logs/mappe_$filename1.log;
				fi 	
		fi
			
 
 fnbwa()
 {
 		if [ ! -d "bwaout" ]; then
   			echo bwaout directory does not exist, creating directory to store alignment output..... >> logs/mappe_$filename1.log;
   			mkdir "bwaout"
		fi			

		if [ -e "bwaout/$filename1.pe.sam" ]; then
			echo Alignment for this aligner is  already done, moving to the next step >> logs/mappe_$filename1.log;
		else
			### doing alignment
			if [ "$adapter" == "yes"  ]; then
				if [ "$tr" == "yes"  ]; then
					### both removed adapters and trimmed reads for mapping
				#	bwa mem  -t 4 $para genome/$reffname  Trim/$filename1.CutAdaptor.trimed.fastq Trim/$filename2.CutAdaptor.trimed.fastq > bwaout/$filename1.pe.sam
					bwa mem  $para genome/$reffname  Trim/$filename1.CutAdaptor.trimed.cmp.fastq Trim/$filename2.CutAdaptor.trimed.cmp.fastq > bwaout/$filename1.pe.sam
					echo Finished Read Mapping for adapter removed and trimmed $filename1  and $filename2 using BWA-mem at  `date` >> logs/mappe_$filename1.log;
				else
					###  only adapter removed reads for mapping
				#	bwa mem  -t 4 $para genome/$reffname  CutAdaptor/$filename1.CutAdaptor.fastq CutAdaptor/$filename2.Cutadaptor.fastq > bwaout/$filename1.pe.sam
					bwa mem  $para genome/$reffname  CutAdaptor/$filename1.CutAdaptor.cmp.fastq CutAdaptor/$filename2.CutAdaptor.cmp.fastq > bwaout/$filename1.pe.sam	
					echo Finished Read Mapping for just adapter removed $filename1 and $filename2 using BWA-mem at  `date` >> logs/mappe_$filename1.log;
				fi
			else
				### adapters not removed
				if [ "$tr" == "yes"  ]; then
					### only trimmed reads for mapping
				#	bwa mem  -t 4 $para genome/$reffname  Trim/$filename1.trimed.fastq Trim/$filename2.trimed.fastq > bwaout/$filename1.pe.sam
					bwa mem  $para genome/$reffname  Trim/$filename1.trimed.cmp.fastq Trim/$filename2.trimed.cmp.fastq > bwaout/$filename1.pe.sam
					echo Finished Read Mapping for just trimmed $filename1  and $filename2 using BWA-mem at  `date`
				else
					### only raw reads used for mapping	
				#	bwa mem  -t 4 $para genome/$reffname  rawdata/$fname1 rawdata/$fname2 > bwaout/$filename1.pe.sam
					bwa mem   $para genome/$reffname  rawdata/$fname1 rawdata/$fname2 > bwaout/$filename1.pe.sam
					echo Finished Read Mapping for raw $filename1  and $filename2 using BWA-mem at  `date` >> logs/mappe_$filename1.log;
				fi 	
			fi
			
	   fi
		echo Calculating read mapping statistics for using BWA-mem at  `date` >> logs/mappe_$filename1.log;
		if [ -e "bwaout/$filename1.pe.sorted.bam" ]; then
			echo Sorting and indexing of BAM file already performed. >> logs/mappe_$filename1.log;
		else
			samtools view -Suh bwaout/$filename1.pe.sam | samtools sort - bwaout/$filename1.pe.sorted
			samtools index bwaout/$filename1.pe.sorted.bam
		fi
		### no of read in the raw file  for pair 1
		 
		rawreads1="$(cat rawdata/$fname1 | echo $((`wc -l`/4)))" 
		echo reads in pair 1 is $rawreads1
		 
		### no of read in the raw file  for pair 2
		 
		rawreads2="$(cat rawdata/$fname2 | echo $((`wc -l`/4)))" 
		echo reads in pair 2 are $rawreads2
		  
		### no of trimmed reads in pair 1 and 	### no of trimmed reads in pair 2
		
		if [ "$adapter" == "yes"  ]; then
				if [ "$tr" == "yes"  ]; then
					### adapter removed and trimmed read statistics
					reads1="$(cat Trim/$filename1.CutAdaptor.trimed.cmp.fastq | echo $((`wc -l`/4)))" 
					reads2="$(cat Trim/$filename2.CutAdaptor.trimed.cmp.fastq | echo $((`wc -l`/4)))" 
					
				else
					### just adapter removed not trimmed
					reads1="$(cat CutAdaptor/$filename1.CutAdaptor.cmp.fastq | echo $((`wc -l`/4)))" 
					reads2="$(cat CutAdaptor/$filename2.CutAdaptor.cmp.fastq | echo $((`wc -l`/4)))" 
				fi
		else
			### adapter not removed
			if [ "$tr" == "yes"  ]; then
					### adapter not removed but just  trimmed read statistics
					reads1="$(cat Trim/$filename1.trimed.cmp.fastq | echo $((`wc -l`/4)))" 
					reads2="$(cat Trim/$filename2.trimed.cmp.fastq | echo $((`wc -l`/4)))" 
				else
					### just raw reads
					reads1=$rawreads1
					reads2=$rawreads2
				fi
		fi
		
		echo Trimmed reads in pair 1 are $reads1
		
		echo Trimmed reads in pair 2 are $reads2
		 
	   	### no of hits
 	   
 	    hits="$(samtools view -c -F 4 bwaout/$filename1.pe.sorted.bam)"
	   	echo alignment are $hits
	     
	    ### uniqly mapped read 
		uniqm="$(samtools view -F 4 bwaout/$filename1.pe.sorted.bam | cut -f 1| sort | uniq -c | wc -l)" 
		echo UNiq Mapped are $uniqm
		
		 ### uniqly hit read 
		uniqhit="$(samtools view -F 4 bwaout/$filename1.pe.sorted.bam | cut -f 1| sort | uniq -u | wc -l)" 
		echo Uniqly hit reads  $uniqhit
		
		### count no of unmapped reads	  
		unmapped="$(samtools view -c -f 4 bwaout/$filename1.pe.sorted.bam)" 
		echo unmapped are $unmapped
	     
	   	### count perfect match  hits
	    PMhits="$(samtools view bwaout/$filename1.pe.sorted.bam | grep -o "NM:i:0" | wc -l )"
	   	echo PM hits is $PMhits
		 
		### count perfect match uniq hits
	    PM="$(samtools view bwaout/$filename1.pe.sorted.bam | grep -w "NM:i:0" | cut -f1 | sort | uniq -u | wc -l )"
	   	echo PM uniq hit is $PM
	   	  
	   	### count 1 mismatcth hits
		oneMMhits="$(samtools view bwaout/$filename1.pe.sorted.bam | grep -o "NM:i:1" | wc -l )"  
		echo 1MM hits $oneMMhits
		
		### count 1 mismatcth uniq hits
		oneMM="$(samtools view bwaout/$filename1.pe.sorted.bam | grep -w "NM:i:1" | cut -f1 | sort | uniq -u | wc -l )"  
		echo 1MM uniq  hits $oneMM
		
		### count 2 mismatcth hits
		twoMMhits="$(samtools view bwaout/$filename1.pe.sorted.bam | grep -o "NM:i:2" | wc -l)"  
		echo 2MM hit $twoMMhits
		
		### count 2 mismatcth hits
		twoMM="$(samtools view bwaout/$filename1.pe.sorted.bam | grep -w "NM:i:2" | cut -f1 | sort | uniq -u | wc -l)"  
		echo 2MM hit $twoMM

		# echo -e "$align\t$filename\t$rawreads\t$reads\t$hits\t$PM\t$oneMM\t$twoMM\t$uniqm\t$mapped\t$unmapped" >> $align.mappingstat.txt; 
		echo -e "'BWA'\t$sample\t$rawreads1\t$rawreads2\t$reads1\t$reads2\t$hits\t$uniqm\t$uniqhit\t$unmapped\t$PMhits\t$PM\t$oneMMhits\t$oneMM\t$twoMMhits\t$twoMM"  >> $align.pe.mappingstat.txt;  
		echo Read mapping statistics for $filename using BWA-mem  is done at  `date` and stored in file "$align.mappingstat.txt" >> logs/mappe_$filename1.log;
		  
 }
 
 fnbowtie()
 {
 		if [ ! -d "bt2out" ]; then
   			echo btout directory does not exist, creating directory to store alignment output..... >> logs/mappe_$filename1.log;
   			mkdir "bt2out"
		fi	
		
		if [ -e "bt2out/$filename1.pe.sam" ]; then
			echo read mapping using this $align for $filename1 and $filename2 sample is  already done, moving to the next step >> logs/mappe_$filename1.log;
		else	
		### doing alignment
			if [ "$adapter" == "yes"  ]; then
				if [ "$tr" == "yes"  ]; then
					### both removed adapters and then trimmed reads used for mapping
		 			bowtie2 $para -x genome/RefGenome -1 Trim/$filename1.CutAdaptor.trimed.cmp.fastq -2 Trim/$filename2.CutAdaptor.trimed.cmp.fastq  -S bt2out/$filename1.pe.sam
		 			echo Finished Read Mapping using Bowtie2  for pair end adapter removed and trimmed  $filename1 and $filename2 sample at `date` >> logs/mappe_$filename1.log;
		 		else
		 			### only adapter removed reads used for mapping
		 			bowtie2 $para -x genome/RefGenome -1 CutAdaptor/$filename1.CutAdaptor.cmp.fastq -2 CutAdaptor/$filename2.CutAdaptor.cmp.fastq  -S bt2out/$filename1.pe.sam
		 			echo Finished Read Mapping using Bowtie2  for pair end adapter removed  $filename1 and $filename2 sample at `date` >> logs/mappe_$filename1.log;
		 		fi	
		 	else
		 		if [ "$tr" == "yes"  ]; then
					### only  trimmed reads used for mapping
		 			bowtie2  $para -x genome/RefGenome -1 Trim/$filename1.trimed.cmp.fastq -2 Trim/$filename2.trimed.cmp.fastq  -S bt2out/$filename1.pe.sam
		 			echo Finished Read Mapping using Bowtie2  for pair end  trimmed  $filename1 and $filename2 sample at `date` >> logs/mappe_$filename1.log;
		 		else
		 			### only raw  reads used for mapping
		 			bowtie2  $para -x genome/RefGenome -1 rawdata/$fname1 -2 rawdata/$fname2  -S bt2out/$filename1.pe.sam
		 			echo Finished Read Mapping using Bowtie2  for pair end raw $filename1 and $filename2 sample at `date` >> logs/mappe_$filename1.log;
		 		fi	
		 	fi	
		fi	
		  
		echo Calculating read mapping statistics  using BOWTIE at  `date` >> logs/mappe_$filename1.log;
		  
	     if [ -e "bt2out/$filename1.pe.sorted.bam" ]; then
	     	echo  Sorting and indexing of BAM file already performed. >> logs/mappe_$filename1.log;
		 else
		 	echo  Sorting and indexing BAM file. >> logs/mappe_$filename1.log;
		 	samtools view -Suh bt2out/$filename1.pe.sam | samtools sort - bt2out/$filename1.pe.sorted
		 	samtools index bt2out/$filename1.pe.sorted.bam	
		 fi
		 ### no of read in the raw file for pair 1
		 rawreads1="$(cat rawdata/$fname1 | echo $((`wc -l`/4)))" 
		 echo reads in pair 1 is $rawreads1
		 
		 ### no of read in the raw file  for pair 2
		 
		 rawreads2="$(cat rawdata/$fname2 | echo $((`wc -l`/4)))" 
		 echo reads in pair 2 are $rawreads2
		  
		  ### no of trimmed reads in pair 1 and 	### no of trimmed reads in pair 2
		
		if [ "$adapter" == "yes"  ]; then
				if [ "$tr" == "yes"  ]; then
					### adapter removed and trimmed read statistics
					reads1="$(cat Trim/$filename1.CutAdaptor.trimed.cmp.fastq | echo $((`wc -l`/4)))" 
					reads2="$(cat Trim/$filename2.CutAdaptor.trimed.cmp.fastq | echo $((`wc -l`/4)))" 
					
				else
					### just adapter removed not trimmed
					reads1="$(cat CutAdaptor/$filename1.CutAdaptor.cmp.fastq | echo $((`wc -l`/4)))" 
					reads2="$(cat CutAdaptor/$filename2.CutAdaptor.cmp.fastq | echo $((`wc -l`/4)))" 
				fi
		else
			### adapter not removed
			if [ "$tr" == "yes"  ]; then
					### adapter not removed but just  trimmed read statistics
					reads1="$(cat Trim/$filename1.trimed.cmp.fastq | echo $((`wc -l`/4)))" 
					reads2="$(cat Trim/$filename2.trimed.cmp.fastq | echo $((`wc -l`/4)))" 
				else
					### just raw reads
					reads1=$rawreads1
					reads2=$rawreads2
				fi
		fi
		  
		 echo Trimmed reads in pair 1 are $reads1
		 echo Trimmed reads in pair 2 are $reads2
		
		 ### no of hits
		 hits="$(samtools view -c -F 4 bt2out/$filename1.pe.sorted.bam)"
		 echo alignment is $hits
		 
		 ### uniqly  mapped reads
		 uniqm="$(samtools view bt2out/$filename1.pe.sorted.bam |cut -f1 | sort  | uniq -c | wc -l)" 
		 echo Uniq Mapped are $uniqm
		
		 ### uniq  hit reads
		 uniqhit="$(samtools view bt2out/$filename1.pe.sorted.bam |cut -f1 | sort  | uniq -u | wc -l)" 
		 echo Uniqly hit are $uniqhit
		 
		 ### count no of unmapped reads	
		 unmapped="$(samtools view -c -f 4 bt2out/$filename1.pe.sorted.bam)" 
		 echo unmapped are $unmapped
		 
		 ### count perfect matcth hits
		 PMhits="$(samtools view bt2out/$filename1.pe.sorted.bam | grep -o "NM:i:0" | wc -l )"
		 echo PM hits  $PMhits
		  
		 ### count perfect matcth uniq hits
		 PM="$(samtools view bt2out/$filename1.pe.sorted.bam | grep -w "NM:i:0" | cut -f1 | sort | uniq -u | wc -l )"
		 echo PM uniq hit $PM 
		
		 ### count 1 mismatcth hits
		 oneMMhits="$(samtools view bt2out/$filename1.pe.sorted.bam | grep -o "NM:i:1" | wc -l )"
		 echo 1MM hits $oneMMhits 
		 
		 ### count 1 mismatcth uniq hits
		 oneMM="$(samtools view bt2out/$filename1.pe.sorted.bam | grep -w "NM:i:1" | cut -f1 | sort | uniq -u | wc -l )"
		 echo 1MM hits $oneMM
		 
		 ### count 2 mismatcth hits
		 twoMMhits="$(samtools view bt2out/$filename1.pe.sorted.bam | grep -o "NM:i:2" | wc -l)"
		 echo 2MM hits $twoMMhits
		 
		  ### count 2 mismatcth uniq hits
		 twoMMhits="$(samtools view bt2out/$filename1.pe.sorted.bam | grep -w "NM:i:2" | cut -f1 | sort | uniq -u | wc -l)"
		 echo 2MM hits $twoMM
		  
		  
		  echo -e "'Bowtie2'\t$sample\t$rawreads1\t$rawreads2\t$reads1\t$reads2\t$hits\t$uniqm\t$uniqhit\t$unmapped\t$PMhits\t$PM\t$oneMMhits\t$oneMM\t$twoMMhits\t$twoMM"  >> $align.pe.mappingstat.txt; 
		  echo Read mapping statistics for $filename using BOWTIE2  is done at  `date` and stored in file "$align.mappingstat.txt;">> logs/mappe_$filename1.log;
		  
 }
  
 fntophat()
 {
 	  	#	echo Started Read Mapping for $filename using ***Tophat*** at  `date` >> logs/mappe_$filename1.log;

		if [ -e "$filename1.thout/accepted_hits.bam" ]; then
			echo read mapping using this $align for $filename1 and $filename2 sample is  already done, moving to the next step >> logs/mappe_$filename1.log;
		else			
	    #	for anno in  genome/*.gtf; do
		#		echo annotation file : $anno
				### doing alignment
				if [ "$adapter" == "yes"  ]; then
					if [ "$tr" == "yes"  ]; then
						### both adaptor removed and trimmed reads for mapping
		 				tophat $para -o $filename1.pe.thout genome/RefGenome Trim/$filename1.CutAdaptor.trimed.cmp.fastq Trim/$filename2.CutAdaptor.trimed.cmp.fastq 	
 						echo Finished Read Mapping of adaptor removed and trimmed sample using Tophat with parameter $para at  `date` >> logs/mappe_$filename1.log;
 					else
 						### just adapter removed reads
 						tophat $para -o $filename1.pe.thout genome/RefGenome CutAdaptor/$filename1.CutAdaptor.cmp.fastq CutAdaptor/$filename2.CutAdaptor.cmp.fastq 	
 						echo Finished Read Mapping of adaptor removed sample using Tophat with parameter $para at  `date` >> logs/mappe_$filename1.log;
 					fi		
 				else
 					if [ "$tr" == "yes"  ]; then
		 				### used just trimmed reads for mapping
		 				tophat $para -o $filename1.pe.thout genome/RefGenome Trim/$filename1.trimed.cmp.fastq Trim/$filename2.trimed.cmp.fastq 	
		 				echo Finished Read Mapping of trimmed sample using Tophat with with parameter $para at  `date` >> logs/mappe_$filename1.log;
 					else
 						### using raw reads for 
 						tophat $para -o $filename1.pe.thout genome/RefGenome rawdata/$fname1 rawdata/$fname2 	
 						echo Finished Read Mapping of raw sample using Tophat with parameter $para at  `date` >> logs/mappe_$filename1.log;
 					fi	
 			  fi
		#	done
			echo Finished Read Mapping using Tophat with parameter $para at  `date` >> logs/mappe_$filename1.log;
		fi
		
		echo Calculating read mapping statistics using Tophat at  `date` >> logs/mappe_$filename1.log;
		
		### no of read in the raw file  for pair 1
		 
		rawreads1="$(cat rawdata/$fname1 | echo $((`wc -l`/4)))" 
		echo reads in pair 1 is $rawreads1
		 
		### no of read in the raw file  for pair 2
		 
		rawreads2="$(cat rawdata/$fname2 | echo $((`wc -l`/4)))" 
		echo reads in pair 2 are $rawreads2
		  
		### no of trimmed reads in pair 1 and 	### no of trimmed reads in pair 2
		
		if [ "$adapter" == "yes"  ]; then
				if [ "$tr" == "yes"  ]; then
					### adapter removed and trimmed read statistics
					reads1="$(cat Trim/$filename1.CutAdaptor.trimed.cmp.fastq | echo $((`wc -l`/4)))" 
					reads2="$(cat Trim/$filename2.CutAdaptor.trimed.cmp.fastq | echo $((`wc -l`/4)))" 
					
				else
					### just adapter removed not trimmed
					reads1="$(cat CutAdaptor/$filename1.CutAdaptor.cmp.fastq | echo $((`wc -l`/4)))" 
					reads2="$(cat CutAdaptor/$filename2.CutAdaptor.cmp.fastq | echo $((`wc -l`/4)))" 
				fi
		else
			### adapter not removed
			if [ "$tr" == "yes"  ]; then
					### adapter not removed but just  trimmed read statistics
					reads1="$(cat Trim/$filename1.trimed.cmp.fastq | echo $((`wc -l`/4)))" 
					reads2="$(cat Trim/$filename2.trimed.cmp.fastq | echo $((`wc -l`/4)))" 
			else
					### just raw reads
					reads1=$rawreads1
					reads2=$rawreads2
			fi
		fi  
		
		echo Trimmed reads in pair 1 are $reads1
		echo Trimmed reads in pair 2 are $reads2
		  
		### no of hits
		hits="$(samtools view -c $filename1.pe.thout/accepted_hits.bam)"
		echo alignment is $hits
		 
		### uniq mapped read
		uniqm="$(samtools view $filename1.pe.thout/accepted_hits.bam | cut -f1 | sort | uniq -c | wc -l )" 
		echo Uniq Mapped are $uniqm
		
		### uniqly hit read
		uniqhit="$(samtools view $filename1.pe.thout/accepted_hits.bam | cut -f1 | sort | uniq -u | wc -l )" 
		echo Uniqly hit read $uniqhit
		
		### count no of unmapped reads
		#samtools view -c -f 4 	
		unmapped="$(samtools view -c -f 4 $filename1.pe.thout/unmapped.bam)" 
		echo unmapped are $unmapped
		
		### count perfect matcth hits
		PMhits="$(samtools view $filename1.pe.thout/accepted_hits.bam | grep -o "NM:i:0" | wc -l )"
		echo PM hits $PMhits
		  
		### count perfect matcth uniq hit
		PM="$(samtools view $filename1.pe.thout/accepted_hits.bam | grep -w "NM:i:0" | cut -f1 | sort | uniq -u | wc -l )"
		echo PM  uniq hit $PM  
		  
		### count 1 mismatcth hits
		oneMMhits="$(samtools view $filename1.pe.thout/accepted_hits.bam | grep -o "NM:i:1" | wc -l )"
		echo   1MM hits $oneMMhits
		
		### count 1 mismatcth uniq hit
		oneMM="$(samtools view $filename1.pe.thout/accepted_hits.bam | grep -w "NM:i:1" | cut -f1 | sort | uniq -u | wc -l )"
		echo   1MM uniq hits $oneMM
		
		### count 2 mismatcth hit 
		twoMMhits="$(samtools view $filename1.pe.thout/accepted_hits.bam | grep -o "NM:i:2" | wc -l)"
	 	echo 2MM hits $twoMMhits
	 	
	 	### count 2 mismatcth uniq hit 
		twoMM="$(samtools view $filename1.pe.thout/accepted_hits.bam | grep -w "NM:i:2" |  cut -f1 | sort | uniq -u | wc -l)"
	 	echo 2MM uniq hits $twoMM
	 	
		# echo -e "$align\t$filename\t$rawreads\t$reads\t$hits\t$PM\t$oneMM\t$twoMM\t$uniqm\t$mapped\t$unmapped" >> $align.mappingstat.txt; 
		  echo -e "'Tophat'\t$sample\t$rawreads1\t$rawreads2\t$reads1\t$reads2\t$hits\t$uniqm\t$uniqhit\t$unmapped\t$PMhits\t$PM\t$oneMMhits\t$oneMM\t$twoMMhits\t$twoMM"  >> $align.pe.mappingstat.txt; 
		  echo Read mapping statistics for $filename using Tophat is done at  `date` and stored in file "tophat.mappingstat.txt" >> logs/mappe_$filename1.log;
 } 
  
  case "$align" in 
	'BWA'|'bwa')		
		  
		fnbwa			 
		  
		;;
	'bowtie'|'BOWTIE')
	
		fnbowtie	  
			 
		;;
	'tophat'|'TOPHAT')
		
		fntophat
		
		;;
	'all'|'ALL')
		echo Started Read Mapping for $filename1 and $filename2 using BWA in All aligner option at  `date` >> logs/mappe_$filename1.log;
		fnbwa
		echo Finished Read Mapping for $filename1 and $filename2 using BWA in All aligner option at  `date` >> logs/mappe_$filename1.log;
		
		echo Started Read Mapping for $filename1 and $filename2 using Bowtie2 in All Aligner option at  `date` >> logs/mappe_$filename1.log;
		fnbowtie
		echo Finished Read Mapping for $filename1  and $filename2 using Bowtie2 in All aligner option at  `date` >> logs/mappe_$filename1.log;
	
		echo Started Read Mapping for $filename1 and $filename2 using Tophat in All aligner option at  `date` >> logs/mappe_$filename1.log;
		fntophat
		echo Finished Read Mapping for $filename1 and $filename2 using Tophat in All aligner option at  `date` >> logs/mappe_$filename1.log;
		;;
esac


# echo Finished Read Mapping for $filename using $align aligner at  `date` >> logs/mappe_$filename1.log;
echo Finish alignment for $fname1 and $fname2 with $align at  `date` >> logs/mappe_$filename1.log;

echo ------------------------------------------------------------------------------------ >> logs/mappe_$filename1.log;    
