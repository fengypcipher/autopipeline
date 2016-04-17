#!/bin/sh
### script to run a serial job using one core on htc using queue

### Set the number of cores (cpus) and memory that will be used for this job
### PBS -l nodes=1:ppn=1,pmem=8gb

#PBS -l nodes=1:ppn=4,pmem=4gb

### Specify "wallclock time" required for this job, hhh:mm:ss
#PBS -l walltime=40:00:00

### Declares the standard error stream of the job will be merged with the standard output stream of the job
#PBS -j oe 

### Request email when job begins and ends
#PBS -m bea

### Specify email address to use for notification.
### PBS -M bharatijadhav@waksman.rutgers.edu
##PBS -M ypfeng@waksman.rutgers.edu
### cd: set directory for job execution
cd $mycwd
echo working dir for pinepeline is $mycwd
filename=`echo $fname|cut -d '.' -f 1`
sample=`echo $filename|cut -d '_' -f 1`
if [ ! -d "logs" ]; then        
    mkdir "logs"
fi	
 echo ------------------------------------------------------------------------------------ >> logs/mapse_$filename.log;    
 
 echo Started alignment for $filename with $align aligner at  `date` >> logs/mapse_$filename.log;
 
 fnbwa()
 {
 		if [ ! -d "bwaout" ]; then
   			echo bwaout directory does not exist, creating directory to store alignment output..... >> logs/mapse_$filename.log;
   			mkdir "bwaout"
		fi			

		if [ -e "bwaout/$filename.se.sam" ]; then
			echo Alignment for $align is  already done, moving to the next step >> logs/mapse_$filename.log;
		else
			### doing alignment
			if [ "$adapter" == "yes"  ]; then
				if [ "$tr" == "yes"  ]; then
					### both removed adapters and trimmed reads for mapping
					bwa bwasw  $para genome/$reffname Trim/$filename.CutAdaptor.trimed.fastq > bwaout/$filename.se.sam
					echo Finished Read Mapping for adapter removed and trimmed $filename using BWA-bwasw at  `date` >> logs/mapse_$filename.log;
				else
					###  only adapter removed reads for mapping
					bwa bwasw $para genome/$reffname  CutAdaptor/$filename.CutAdaptor.fastq > bwaout/$filename.se.sam	
					echo Finished Read Mapping for just adapter removed $filename using BWA-bwasw at  `date` >> logs/mapse_$filename.log;
				fi
			else
				### adapters not removed
				if [ "$tr" == "yes"  ]; then
					### only trimmed reads for mapping
					bwa bwasw  $para genome/$reffname  Trim/$filename.trimed.fastq > bwaout/$filename.se.sam
					echo Finished Read Mapping for just trimmed $filename using BWA-bwasw at  `date`
				else
					### only raw reads used for mapping	
					bwa bwasw $para genome/$reffname  rawdata/$fname > bwaout/$filename.se.sam
					echo Finished Read Mapping for raw $filename using BWA-bwasw at  `date` >> logs/mapse_$filename.log;
				fi 	
			fi
			
	   fi
	 echo Calculating read mapping statistics for using BWA-mem at  `date` >> logs/mapse_$filename.log;
		if [ -e "bwaout/$filename.se.sorted.bam" ]; then
			echo Sorting and indexing of BAM file already performed. >> logs/mapse_$filename.log;
		else
			samtools view -Suh bwaout/$filename.se.sam | samtools sort - bwaout/$filename.se.sorted
			samtools index bwaout/$filename.se.sorted.bam
		fi
		### no of read in the raw file  for single end file
		 
		rawreads="$(cat rawdata/$fname | echo $((`wc -l`/4)))" 
		echo reads in  $rawreads		 
		
		### no of trimmed reads in pair 1 and 	### no of trimmed reads in pair 2
		
		if [ "$adapter" == "yes"  ]; then
				if [ "$tr" == "yes"  ]; then
					### adapter removed and trimmed read statistics
					reads="$(cat Trim/$filename.CutAdaptor.trimed.fastq | echo $((`wc -l`/4)))" 
					
				else
					### just adapter removed not trimmed
					reads="$(cat CutAdaptor/$filename.CutAdaptor.fastq | echo $((`wc -l`/4)))" 
				fi
		else
			### adapter not removed
			if [ "$tr" == "yes"  ]; then
					### adapter not removed but just  trimmed read statistics
					reads="$(cat Trim/$filename.trimed.fastq | echo $((`wc -l`/4)))" 
			else
					### just raw reads
					reads=$rawreads
			fi
		fi
		
		echo Trimmed reads in single end reads are $reads
	
	   	### no of hits
 	   
 	    hits="$(samtools view -c -F 4 bwaout/$filename.se.sorted.bam)"
	   	echo alignment are $hits
	     
	    ### uniqly mapped read 
		uniqm="$(samtools view -F 4 bwaout/$filename.se.sorted.bam | cut -f 1| sort | uniq -c | wc -l)" 
		echo UNiq Mapped are $uniqm
		
		 ### uniqly hit read 
		uniqhit="$(samtools view -F 4 bwaout/$filename.se.sorted.bam | cut -f 1| sort | uniq -u | wc -l)" 
		echo Uniqly hit reads  $uniqhit
		
		### count no of unmapped reads	  
		unmapped="$(samtools view -c -f 4 bwaout/$filename.se.sorted.bam)" 
		echo unmapped are $unmapped
	     
	   	### count perfect match  hits
	    PMhits="$(samtools view bwaout/$filename.se.sorted.bam | grep -o "NM:i:0" | wc -l )"
	   	echo PM hits is $PMhits
		 
		### count perfect match uniq hits
	    PM="$(samtools view bwaout/$filename.se.sorted.bam | grep -w "NM:i:0" | cut -f1 | sort | uniq -u | wc -l )"
	   	echo PM uniq hit is $PM
	   	  
	   	### count 1 mismatcth hits
		oneMMhits="$(samtools view bwaout/$filename.se.sorted.bam | grep -o "NM:i:1" | wc -l )"  
		echo 1MM hits $oneMMhits
		
		### count 1 mismatcth uniq hits
		oneMM="$(samtools view bwaout/$filename.se.sorted.bam | grep -w "NM:i:1" | cut -f1 | sort | uniq -u | wc -l )"  
		echo 1MM uniq  hits $oneMM
		
		### count 2 mismatcth hits
		twoMMhits="$(samtools view bwaout/$filename.se.sorted.bam | grep -o "NM:i:2" | wc -l)"  
		echo 2MM hit $twoMMhits
		
		### count 2 mismatcth hits
		twoMM="$(samtools view bwaout/$filename.se.sorted.bam | grep -w "NM:i:2" | cut -f1 | sort | uniq -u | wc -l)"  
		echo 2MM hit $twoMM

		echo -e "'BWA'\t$sample\t$rawreads\t$reads\t$hits\t$uniqm\t$uniqhit\t$unmapped\t$PMhits\t$PM\t$oneMMhits\t$oneMM\t$twoMMhits\t$twoMM"  >> $align.se.mappingstat.txt;  
		echo Read mapping statistics for $filename using BWA-mem  is done at  `date` and stored in file "$align.se.mappingstat.txt" >> logs/mapse_$filename.log;
		  
 }
 
  ####################################
  	# calling bowtie function		
  ####################################
  fnbowtie()
 {
 		if [ ! -d "bt2out" ]; then
   			echo btout directory does not exist, creating directory to store alignment output..... >> logs/mapse_$filename.log;
   			mkdir "bt2out"
		fi	
		
		if [ -e "bt2out/$filename.se.sam" ]; then
			echo read mapping using this $align for $filename sample is  already done, moving to the next step >> logs/mapse_$filename.log;
		else	
		### doing alignment
			if [ "$adapter" == "yes"  ]; then
				if [ "$tr" == "yes"  ]; then
					### both removed adapters and then trimmed reads used for mapping
		 			bowtie2 $para -x genome/RefGenome Trim/$filename.CutAdaptor.trimed.fastq  -S bt2out/$filename.se.sam
		 			echo Finished Read Mapping using Bowtie2  for adapter removed and trimmed  $filename sample at `date` >> logs/mapse_$filename.log;
		 		else
		 			### only adapter removed reads used for mapping
		 			bowtie2 $para -x genome/RefGenome CutAdaptor/$filename.CutAdaptor.fastq  -S bt2out/$filename.se.sam
		 			echo Finished Read Mapping using Bowtie2  for adapter removed  $filename sample at `date` >> logs/mapse_$filename.log;
		 		fi	
		 	else
		 		if [ "$tr" == "yes"  ]; then
					### only  trimmed reads used for mapping
		 			bowtie2  $para -x genome/RefGenome Trim/$filename.trimed.fastq  -S bt2out/$filename.se.sam
		 			echo Finished Read Mapping using Bowtie2  for  trimmed  $filename sample at `date` >> logs/mapse_$filename.log;
		 		else
		 			### only raw  reads used for mapping
		 			bowtie2  $para -x genome/RefGenome rawdata/$fname  -S bt2out/$filename.se.sam
		 			echo Finished Read Mapping using Bowtie2  for raw $filename sample at `date` >> logs/mapse_$filename.log;
		 		fi	
		 	fi	
		fi	
		  
		echo Calculating read mapping statistics  using BOWTIE at  `date` >> logs/mapse_$filename.log;
		  
	     if [ -e "bt2out/$filename.se.sorted.bam" ]; then
	     	echo  Sorting and indexing of BAM file already performed. >> logs/mapse_$filename.log;
		 else
		 	echo  Sorting and indexing BAM file. >> logs/mapse_$filename.log;
		 	samtools view -Suh bt2out/$filename.se.sam | samtools sort - bt2out/$filename.se.sorted
		 	samtools index bt2out/$filename.se.sorted.bam	
		 fi
		 ### no of raw reads
		 rawreads="$(cat rawdata/$fname | echo $((`wc -l`/4)))" 
		 echo reads in pair 1 is $rawreads
		  
		  ### no of trimmed reads in pair 1 and 	### no of trimmed reads in pair 2
		
		if [ "$adapter" == "yes"  ]; then
				if [ "$tr" == "yes"  ]; then
					### adapter removed and trimmed read statistics
					reads="$(cat Trim/$filename.CutAdaptor.trimed.fastq | echo $((`wc -l`/4)))" 
					
				else
					### just adapter removed not trimmed
					reads="$(cat CutAdaptor/$filename.CutAdaptor.fastq | echo $((`wc -l`/4)))" 
				fi
		else
			### adapter not removed
			if [ "$tr" == "yes"  ]; then
					### adapter not removed but just  trimmed read statistics
					reads="$(cat Trim/$filename.trimed.fastq | echo $((`wc -l`/4)))" 
				else
					### just raw reads
					reads=$rawreads
					
				fi
		fi
		  
		 echo Trimmed reads  are $reads
		
		
		 ### no of hits
		 hits="$(samtools view -c -F 4 bt2out/$filename.se.sorted.bam)"
		 echo alignment is $hits
		 
		 ### uniqly  mapped reads
		 uniqm="$(samtools view -F 4 bt2out/$filename.se.sorted.bam |cut -f1 | sort  | uniq -c | wc -l)" 
		 echo Uniq Mapped are $uniqm
		
		 ### uniq  hit reads
		 uniqhit="$(samtools view -F 4 bt2out/$filename.se.sorted.bam |cut -f1 | sort  | uniq -u | wc -l)" 
		 echo Uniqly hit are $uniqhit
		 
		 ### count no of unmapped reads	
		 unmapped="$(samtools view -c -f 4 bt2out/$filename.se.sorted.bam)" 
		 echo unmapped are $unmapped
		 
		 ### count perfect matcth hits
		 PMhits="$(samtools view bt2out/$filename.se.sorted.bam | grep -o "NM:i:0" | wc -l )"
		 echo PM hits  $PMhits
		  
		 ### count perfect matcth uniq hits
		 PM="$(samtools view bt2out/$filename.se.sorted.bam | grep -w "NM:i:0" | cut -f1 | sort | uniq -u | wc -l )"
		 echo PM uniq hit $PM 
		
		 ### count 1 mismatcth hits
		 oneMMhits="$(samtools view bt2out/$filename.se.sorted.bam | grep -o "NM:i:1" | wc -l )"
		 echo 1MM hits $oneMMhits 
		 
		 ### count 1 mismatcth uniq hits
		 oneMM="$(samtools view bt2out/$filename.se.sorted.bam | grep -w "NM:i:1" | cut -f1 | sort | uniq -u | wc -l )"
		 echo 1MM hits $oneMM
		 
		 ### count 2 mismatcth hits
		 twoMMhits="$(samtools view bt2out/$filename.se.sorted.bam | grep -o "NM:i:2" | wc -l)"
		 echo 2MM hits $twoMMhits
		 
		  ### count 2 mismatcth uniq hits
		 twoMMhits="$(samtools view bt2out/$filename.se.sorted.bam | grep -w "NM:i:2" | cut -f1 | sort | uniq -u | wc -l)"
		 echo 2MM hits $twoMM
		  
		  
		  echo -e "'Bowtie2'\t$sample\t$rawreads\t$reads1\t$hits\t$uniqm\t$uniqhit\t$unmapped\t$PMhits\t$PM\t$oneMMhits\t$oneMM\t$twoMMhits\t$twoMM"  >> $align.se.mappingstat.txt; 
		  echo Read mapping statistics for $filename using BOWTIE2  is done at  `date` and stored in file "$align.se.mappingstat.txt;">> logs/mapse_$filename.log;
		  
 }

 
  
   ####################################
  	# calling tophat function		
   ####################################
  
  
  fntophat()
 {

		if [ -e "$filename.thout/accepted_hits.bam" ]; then
			echo read mapping using this $align for $filename1 and $filename2 sample is  already done, moving to the next step >> logs/mapse_$filename.log;
		else			
				### doing alignment
				if [ "$adapter" == "yes"  ]; then
					if [ "$tr" == "yes"  ]; then
						### both adaptor removed and trimmed reads for mapping
		 				tophat $para -o $filename.thout genome/RefGenome Trim/$filename.CutAdaptor.trimed.fastq 	
 						echo Finished Read Mapping of adaptor removed and trimmed sample using Tophat with parameter $para at  `date` >> logs/mapse_$filename.log;
 					else
 						### just adapter removed reads
 						tophat $para -o $filename.thout genome/RefGenome CutAdaptor/$filename.CutAdaptor.fastq C
 						echo Finished Read Mapping of adaptor removed sample using Tophat with parameter $para at  `date` >> logs/mapse_$filename.log;
 					fi		
 				else
 					if [ "$tr" == "yes"  ]; then
		 				### used just trimmed reads for mapping
		 				tophat $para -o $filename.thout genome/RefGenome Trim/$filename.trimed.fastq 
		 				echo Finished Read Mapping of trimmed sample using Tophat with with parameter $para at  `date` >> logs/mapse_$filename.log;
 					else
 						### using raw reads for 
 						tophat $para -o $filename1.pe.thout genome/RefGenome rawdata/$fname	
 						echo Finished Read Mapping of raw sample using Tophat with parameter $para at  `date` >> logs/mapse_$filename.log;
 					fi	
 			  fi
		#	done
			echo Finished Read Mapping using Tophat with parameter $para at  `date` >> logs/mapse_$filename.log;
		fi
		
		echo Calculating read mapping statistics using Tophat at  `date` >> logs/mapse_$filename.log;
		
		### no of read in the raw file  for pair 1
		 
		rawreads="$(cat rawdata/$fname | echo $((`wc -l`/4)))" 
		echo reads in pair 1 is $rawreads
		 
		### no of trimmed reads in pair 1 and 	### no of trimmed reads in pair 2
		
		if [ "$adapter" == "yes"  ]; then
				if [ "$tr" == "yes"  ]; then
					### adapter removed and trimmed read statistics
					reads="$(cat Trim/$filename.CutAdaptor.trimed.fastq | echo $((`wc -l`/4)))" 
					
				else
					### just adapter removed not trimmed
					reads="$(cat CutAdaptor/$filename.CutAdaptor.fastq | echo $((`wc -l`/4)))" 
				fi
		else
			### adapter not removed
			if [ "$tr" == "yes"  ]; then
					### adapter not removed but just  trimmed read statistics
					reads="$(cat Trim/$filename.trimed.fastq | echo $((`wc -l`/4)))" 
			else
					### just raw reads
					reads=$rawreads
			fi
		fi  
		
		echo Trimmed reads are $reads
		  
		### no of hits
		hits="$(samtools view -c $filename.thout/accepted_hits.bam)"
		echo alignment is $hits
		 
		### uniq mapped read
		uniqm="$(samtools view $filename.thout/accepted_hits.bam | cut -f1 | sort | uniq -u | wc -l )" 
		echo Uniq Mapped are $uniqm
		
		### uniqly hit read
		uniqhit="$(samtools view $filename.thout/accepted_hits.bam | cut -f1 | sort | uniq -u | wc -l )" 
		echo Uniqly hit read $uniqhit
		
		### count no of unmapped reads	
		unmapped="$(samtools view -c $filename.thout/unmapped.bam)" 
		echo unmapped are $unmapped
		
		### count perfect matcth hits
		PMhits="$(samtools view $filename.thout/accepted_hits.bam | grep -o "NM:i:0" | wc -l )"
		echo PM hits $PMhits
		  
		### count perfect matcth uniq hit
		PM="$(samtools view $filename.thout/accepted_hits.bam | grep -w "NM:i:0" | cut -f1 | sort | uniq -u | wc -l )"
		echo PM  uniq hit $PM  
		  
		### count 1 mismatcth hits
		oneMMhits="$(samtools view $filename.thout/accepted_hits.bam | grep -o "NM:i:1" | wc -l )"
		echo   1MM hits $oneMMhits
		
		### count 1 mismatcth uniq hit
		oneMM="$(samtools view $filename.thout/accepted_hits.bam | grep -w "NM:i:1" | cut -f1 | sort | uniq -u | wc -l )"
		echo   1MM uniq hits $oneMM
		
		### count 2 mismatcth hit 
		twoMMhits="$(samtools view $filename.thout/accepted_hits.bam | grep -o "NM:i:2" | wc -l)"
	 	echo 2MM hits $twoMMhits
	 	
	 	### count 2 mismatcth uniq hit 
		twoMM="$(samtools view $filename.thout/accepted_hits.bam | grep -w "NM:i:2" |  cut -f1 | sort | uniq -u | wc -l)"
	 	echo 2MM uniq hits $twoMM
	 	
		  echo -e "'Tophat'\t$sample\t$rawreads\t$reads\t$hits\t$uniqm\t$uniqhit\t$unmapped\t$PMhits\t$PM\t$oneMMhits\t$oneMM\t$twoMMhits\t$twoMM"  >> $align.se.mappingstat.txt; 
		  echo Read mapping statistics for $filename using Tophat is done at  `date` and stored in file "tophat.se.mappingstat.txt" >> logs/mapse_$filename.log;
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
		echo Started Read Mapping for $filename using BWA in All aligner option at  `date` >> logs/mapse_$filename.log;
		fnbwa
		echo Finished Read Mapping for $filename using BWA in All aligner option at  `date` >> logs/mapse_$filename.log;
		
		echo Started Read Mapping for $filename using Bowtie2 in All Aligner option at  `date` >> logs/mapse_$filename.log;
		fnbowtie
		echo Finished Read Mapping for $filenameusing Bowtie2 in All aligner option at  `date` >> logs/mapse_$filename.log;
	
		echo Started Read Mapping for $filename using Tophat in All aligner option at  `date` >> logs/mapse_$filename.log;
		fntophat
		echo Finished Read Mapping for $filename using Tophat in All aligner option at  `date` >> logs/mapse_$filename.log;
		
		;;
  esac

# echo Finished Read Mapping for $filename using $align aligner at  `date` >> logs/mapse_$filename.log;
echo Finished alignment for $filename with $align at  `date` >> logs/mapse_$filename.log;

echo ------------------------------------------------------------------------------------ >> logs/mapse_$filename.log;    
