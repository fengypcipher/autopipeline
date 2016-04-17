# This is version2 of the automatic pipeline for short read alignment updated on Nov 1st 2013
# This version is updated to to compare the QC filter fastq files 
# To use cmpfastq.pl need to modify Yaping.trim3.pl and updated version of trim program Yaping.trim4.pl is used instead
#06-24-2014 add RNAseq primer ATGTCAGAATCTCGTATGCCGTCTTCTGCTTG as adaptor
#06-24-2014 add Illumina Single End PCR Primer GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAA  GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG
#08-06-2014 This cutadaptor include all adaptors found in Qinghua small RNA sequencing as follows 
#!/bin/sh
### script to run a serial job using one core on htc using queue

### Set the number of cores (cpus) and memory that will be used for this job
#PBS -l nodes=1:ppn=1,pmem=8gb

### Specify "wallclock time" required for this job, hhh:mm:ss
#PBS -l walltime=40:00:00

### Declares the standard error stream of the job will be merged with the standard output stream of the job
#PBS -j oe 

### Request email when job begins and ends
### PBS -m bea

### Specify email address to use for notification.
### PBS -M bharatijadhav@waksman.rutgers.edu

### cd: set directory for job execution
 cd $mycwd
 # echo working dir for pinepeline is $mycwd
 PRE=/ingens/solid/analysis/yaping/pipeline2/Script
 filename=`echo $fname|cut -d '.' -f 1`

 if [ ! -d "logs" ]; then        
    mkdir "logs"
 fi	
 
 if [ ! -d "fastqcPDF" ]; then        
    mkdir "fastqcPDF"
 fi	
 echo ------------------------------------------------------------------------------------ >> logs/prepro_$filename.log;    
 echo Started preprocessing for $fname at  `date` >> logs/prepro_$filename.log;
 echo Started FASTQC for $fname at `date` >> logs/prepro_$filename.log
 rawfile=rawdata/$filename'_fastqc'
# if [ ! -d "rawdata/$filename_fastqc" ]; then
	if [ ! -d "$rawfile" ]; then
	 fastqc rawdata/$fname
	 raw=rawdata/$filename'_fastqc'/fastqc_report.html
	 prince $raw -o fastqcPDF/$filename.rawqc.pdf
 else
 	echo fastqc is already done for this sample, moving to the next step >> logs/prepro_$filename.log;
 	
 fi
 echo Finisheded FASTQC for $filename at `date` >> logs/prepro_$filename.log;
 if [ "$adapter" == "yes" ]; then
 
 	if [ ! -d "CutAdaptor" ]; then
                echo CutAdaptor directory does not exist, creating directory..... >> logs/prepro_$filename.log;
                mkdir "CutAdaptor"
 	fi		
   		 
 	if [ ! -e "CutAdaptor/$filename.CutAdaptor.fastq" ] ; then
	 	echo Started cut adapt for $fname at `date` >> logs/prepro_$filename.log;
	 	cutadapt -O 5 -e 0.05  -b GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG  -b ACACTCTTTCCCTACACGACGCTCTTCCGATCT  -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT  -b CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT  -b ACACTCTTTCCCTACACGACGCTCTTCCGATCT  -b GATCGGAAGAGCACACGTCT  -b ACACTCTTTCCCTACACGACGCTCTTCCGATCT  -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT  -b GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT  -b ACACTCTTTCCCTACACGACGCTCTTCCGATCT  -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -b GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT  -b CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC  -b CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC  -b CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC  -b CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC  -b CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC  -b CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC  -b CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC  -b CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC  -b CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC  -b CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC  -b CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC  -b CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC  -b GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG  -b ACACTCTTTCCCTACACGACGCTCTTCCGATCT  -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT  -b CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT  -b ACACTCTTTCCCTACACGACGCTCTTCCGATCT  -b CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT  -b GATCGTCGGACTGTAGAACTCTGAAC  -b ACAGGTTCAGAGTTCTACAGTCCGAC  -b CAAGCAGAAGACGGCATACGA  -b TCGTATGCCGTCTTCTGCTTG  -b CAAGCAGAAGACGGCATACGA  -b AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA  -b CGACAGGTTCAGAGTTCTACAGTCCGACGATC  -b TCGGACTGTAGAACTCTGAAC  -b ACAGGTTCAGAGTTCTACAGTCCGACATG  -b CAAGCAGAAGACGGCATACGA  -b TCGTATGCCGTCTTCTGCTTG  -b CAAGCAGAAGACGGCATACGA  -b AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA  -b CCGACAGGTTCAGAGTTCTACAGTCCGACATG  -b CAAGCAGAAGACGGCATACGA  -b GTTCAGAGTTCTACAGTCCGACGATC  -b TCGTATGCCGTCTTCTGCTTGT  -b CAAGCAGAAGACGGCATACGA  -b AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA  -b CGACAGGTTCAGAGTTCTACAGTCCGACGATC  -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTCTTTCCCTACACGACGCTCTTCCGATCT  -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC  -b CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT  -b CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT  -b CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT  -b CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT  -b CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT  -b CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT  -b CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT  -b CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT  -b CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT  -b AAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT  -b CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT  -b CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -b AGATCGGAAGAGC -b ATGTCAGAATCTCGTATGCCGTCTTCTGCTTG -b GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG -b GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAA -b CAAGCAGAAGACGGCATACGAGATTCTGACAT -b ATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTGA -b AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCG -b AGATCGGAAGAGCACACGTCTGAACTCCAGT -b GGAAGAGCACACGTCTGAACTCCAGTCACGC -b GATCGTCGGACTGTAGAACTCTGAACGTGT -b CGGACTGTAGAACTCTGAACGTGTAGATCT -b ATATATATATATATATATATATATATATATATATATATATATATATATAT -b CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -b AGATGTGTATAAGAGACAG -b GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCC -b GGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATC -b "A{150}" -b "G{150}" -b "T{150}" -b "C{150}"  rawdata/$fname -o CutAdaptor/$filename.CutAdaptor.fastq >CutAdaptor/$filename.CutAdaptor.logs 
###This cutadaptor include all adaptors found in Qinghua small RNA sequencing as follows	 
##Illumina Multiplexing PCR Primer 2.01 AGATCGGAAGAGCACACGTCTGAACTCCAGT
##TruSeq Adapter, Index 5 (96% over 31bp) GGAAGAGCACACGTCTGAACTCCAGTCACGC       
##Illumina DpnII expression Sequencing Primer (96% over 30bp) GATCGTCGGACTGTAGAACTCTGAACGTGT
##Illumina DpnII expression Adapter 1 (96% over 26bp)   AATCGTCGGACTGTAGAACTCTGAACGTGT
##Illumina RNA PCR Primer (100% over 30bp) CGGACTGTAGAACTCTGAACGTGTAGATCT
###These adaptors are Nextera mate pairs:
##Circularized Duplicate Junction Adapter CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG 
##Circularized Single Junction Adapter CTGTCTCTTATACACATCT 
##Circularized Single Junction Adapter Reverse Complement AGATGTGTATAAGAGACAG 
##Read 1 External Adapter GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
##Read 2 External Adapter GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
##Illumina Single End PCR Primer 1 (100% over 50bp)  GGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATC
##Illumina Single End PCR Primer 1 (100% over 50bp)  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCC

		echo Finished cutadapt for $filename at `date` >> logs/prepro_$filename.log;	
 	else
 		echo adaptors are already removed for this sample, moving to the next step >> logs/prepro_$filename.log;
 	fi


	 if [ ! -d "CutAdaptor/$filename.CutAdaptor_fastqc" ]; then
		echo Started FASTQC for $filename.CutAdaptor.fastq at `date` >> logs/prepro_$filename.log;
		fastqc CutAdaptor/$filename.CutAdaptor.fastq
		prince CutAdaptor/$filename.CutAdaptor_fastqc/fastqc_report.html -o fastqcPDF/$filename.cutqc.pdf	
    	echo Finished FASTQC for $filename.CutAdaptor.fastq at `date` >> logs/prepro_$filename.log;
 	else
 		echo fastqc is already done for $filename.CutAdaptor.fastq sample, moving to the next step >> logs/prepro_$filename.log;
	fi
 else
	echo User chose not to remove adaptors this sample, moving to the next step >> logs/prepro_$filename.log;
 fi	 
 
 if [ "$tr" == "yes" ]; then
 
 	if [ ! -d "Trim" ]; then
   			echo Trim directory does not exist, creating directory..... >> logs/prepro_$filename.log;
   			mkdir "Trim"  		
 	fi
 	### trim reads for adapter removed reads
	if [ "$adapter" == "yes" ]; then
 		if [ ! -e "Trim/$filename.CutAdaptor.trimed.fastq" ] ; then
 			echo Started trimming for $filename.CutAdaptor.fastq at `date` >> logs/prepro_$filename.log;
			perl $PRE/Yaping.trim4.pl CutAdaptor/$filename.CutAdaptor.fastq 	
			echo Finished trimming for $filename.CutAdaptor.fastq at `date` >> logs/prepro_$filename.log;
 		else
 			echo The $filename sample is already trimmed, moving to the next step >> logs/prepro_$filename.log;
 		fi
 		if [ ! -d "Trim/$filename.CutAdaptor.trimed_fastqc" ]; then	
 			echo Started FASTQC for $filename.CutAdaptor.trimed.fastq at `date` >> logs/prepro_$filename.log;
			fastqc Trim/$filename.CutAdaptor.trimed.fastq	
 			prince Trim/$filename.CutAdaptor.trimed_fastqc/fastqc_report.html -o fastqcPDF/$filename.cuttrimqc.pdf   
    		echo Finished FASTQC for $filename.CutAdaptor.trimed.fastq at `date` >> logs/prepro_$filename.log;
 		else
 			echo fastqc is already done for adapter removed and trimmed $filename moving to the next step >> logs/prepro_$filename.log;
 		fi
 	### trimming without cutting adapters	
 	else
 		if [ ! -e "Trim/$filename.trimed.fastq" ] ; then
 			echo Started trimming for raw $filename at `date` >> logs/prepro_$filename.log;
			perl $PRE/Yaping.trim4.pl rawdata/$fname	
			echo Finished trimming for $fname at `date` >> logs/prepro_$filename.log;
 		else
 			echo The $filename sample is already trimmed, moving to the next step >> logs/prepro_$filename.log;
 		fi
 		if [ ! -d "Trim/$filename.trimed_fastqc" ]; then	
 			echo Started FASTQC for $filename.trimed.fastq at `date` >> logs/prepro_$filename.log;
			fastqc Trim/$filename.trimed.fastq	
 			prince Trim/$filename.trimed_fastqc/fastqc_report.html -o fastqcPDF/$filename.trimqc.pdf   
    		echo Finished FASTQC for $filename.trimed.fastq at `date` >> logs/prepro_$filename.log;
 		else
 			echo fastqc is already done for trimmed $filename  moving to the next step >> logs/prepro_$filename.log;
 		fi
 	fi		 
 else
 	### don't trim
	echo User chose not to trim this sample, moving to the next step >> logs/prepro_$filename.log;
 fi	
 
 echo Finished pre-processing of $filename at `date` >> logs/prepro_$filename.log;
 echo ------------------------------------------------------------------------------------ >> logs/prepro_$filename.log; 

 
