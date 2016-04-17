#!/bin/sh
#PBS -l nodes=1:ppn=1,pmem=8gb
#PBS -l walltime=40:00:00
#PBS -j oe 
#PBS -M ypfeng@waksman.rutgers.edu

cd $mycwd
# echo my current working dir is $mycwd
# ref=`echo $reffname|cut -d '.' -f 1`

if [ ! -d "logs" ]; then        
    mkdir "logs"
fi	
echo Started Building indexes for $align aligner using reference genome $reffname at  `date` >> logs/buildindex.log;

fnbwaidx()
{
		echo checking index for $align using bwa index for $reffname at  `date` >> logs/buildindex.log;
		if [ -e "genome/$reffname.bwt" ] && [ -e "genome/$reffname.ann" ] && [ -e "genome/$reffname.pac" ] && [ -e "genome/$reffname.sa" ] && [ -e "genome/$reffname.amb" ] ; then
			echo BWA Indexes already exist, no need to create index again >> logs/buildindex.log;
		else
			echo Started Building index for $align using bwa index for $reffname at  `date` >> logs/buildindex.log;
		#	bwa index $homedir/genome/$reffname
			bwa index genome/$reffname
			echo Finished Building index using bwa index for $reffname at  `date` >> logs/buildindex.log;
		fi
}

fnbowtieidx()
{
		echo checking index for $align using bowtie index for $reffname at  `date` >> logs/buildindex.log;
	#	if [ -e "$ref.1.bt2" ]  && [ -e "$ref.2.bt2" ] && [ -e "$ref.3.bt2" ] && [ -e "$ref.4.bt2" ]  && [ -e "$ref.rev.1.bt2" ] && [ -e "$ref.rev.2.bt2" ]; then
		if [ -e "genome/RefGenome.1.bt2" ]  && [ -e "genome/RefGenome.2.bt2" ] && [ -e "genome/RefGenome.3.bt2" ] && [ -e "genome/RefGenome.4.bt2" ]  && [ -e "genome/RefGenome.rev.1.bt2" ] && [ -e "genome/RefGenome.rev.2.bt2" ]; then
			echo Indexes already exist, no need to create index again >> logs/buildindex.log;
		else
			echo Started Building index for $align using bowtie2-build for $reffname at  `date` >> logs/buildindex.log;
			cd $mycwd/genome
		    bowtie2-build $reffname RefGenome
		    cd $mycwd
			# bowtie2-build $reffname $ref
			echo Finished Building index using bowtie index at  `date` >> logs/buildindex.log;
		fi
}

case "$align" in 	
	'BWA'|'bwa')				
			fnbwaidx
			;;
	'bowtie'|'BOWTIE' |'tophat' |'TOPHAT')
			fnbowtieidx
			;;
	'all'|'ALL')
			fnbwaidx
			fnbowtieidx
			echo Finished Building index using for all aligners  at  `date` >> logs/buildindex.log;
		;;
esac
echo ---------------------------------------------------------------------------------------------------------------- >> logs/buildindex.log;
