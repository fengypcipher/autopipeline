import re,os,sys,logging,time,datetime,subprocess;
import commands; 
from optparse import OptionParser
from argparse import Action
ALIGNER="BWA"
ADAPTER="yes"
TRIM="yes"
GENOME="./genome/RefGenome.fa"
def main():
    global ALIGNER,ADAPTER,TRIM,GENOME,PARA_MAP
    usage = "usage: autopipeline.pl [arguments]\n\
     \nOptional arguments: --help print help message --man print complete\
    documentation --aligner Aligner (default : BWA). User can specify BWA,\
    BOWTIE, TOPHAT or ALL to do alignment using all alingers. --adapter\
    Remove adapters from raw reads. --trim Trim low quality reads. --genome\
    Input reference genome with full path.\n\
    sampleSheet.txt has to be in working directory."  
    parser = OptionParser(usage) 
    parser.add_option("--aligner",dest="aligner",default="BWA",help="Aligner to\
align sequencing reads to refrence sequences. \n\
The default aligner is BWA but user can specify BOWTIE and internally it \
will use BOWTIE 2, TOPHAT or choose ALL. \nIf no parameter are specified \
after prompt, the pipeline uses specified aligner with all default\
settings.")
    parser.add_option("--man",dest="man",action="store_true",default=False,help="Print complete manual of the program")
    parser.add_option("--adapter",dest="adapter",default="yes",help="Remove \
adapters from raw reads and perform quality control check on the resulting reads. \
\nThe default option is 'yes' but user can specify 'no' to skip this step. \
\nIf this parameter is set as yes, the pipeline will remove adaptors from raw \
reads and generate report after performing quality check.")
    parser.add_option("--trim",dest="trim",default="yes",help="Trim the low quality reads. \
\nThe default option is 'yes' but user can specify 'no' to skip this step. \
\nIf this parameter is set as yes, the pipeline will remove the low quality \
reads and generate report after performing quality check.")
    parser.add_option("--genome",dest="genome",default="./genome/RefGenome.fa",help="Reference genomepathway and file name. \
\nThe default option is './genome/RefGenome.fa'. This is a mandatory option. \
\nYou have to give pipepline the correct reference genome path and file name.")
    parser.add_option("--para_map",dest="para_map",default="",help="parameters for mapping. \
\nThe default options is alinger defualt. This is an optional option.")
    (options, args) = parser.parse_args()
    ALIGNER=options.aligner
    ADAPTER=options.adapter
    TRIM=options.trim
    GENOME=options.genome
    PARA_MAP=options.para_map
    if options.man:
        print("SYNOPSIS:autopipeline.pl [arguments]\n\
    \nOptional arguments: --help print help message --man print complete\
    documentation --aligner Aligner (default : BWA). User can specify BWA,\
    BOWTIE, TOPHAT or ALL to do alignment using all alingers. --adapter\
    Remove adapters from raw reads. --trim Trim low quality reads. --genome\
    Input reference genome with full path. sampleSheet.txt has to be in working\
    directory.\n")
        print(parser.format_option_help())
        parser.exit()
if __name__ == '__main__':
    main()
    mydir = os.getcwd()
    print("testing:"+mydir)
    print("testing:"+ALIGNER)
    print("testing:"+GENOME)


##################################################
####index genome#####
if ("$GENOME"):
    #print refFile, aligner, mydir 

    ### submit a job to build indices
    jobqsub = 'qsub -v reffname='+GENOME+',align='+ALIGNER+',mycwd='+mydir+' -N job_buildIndex_'+GENOME+' /ingens/home/ypfeng/clear/Script/optionalindex.sh'
    proc = subprocess.Popen(jobqsub, shell=True, stdout=subprocess.PIPE)
    data = proc.stdout.readline() 
    #jobid_index = data
    a = data.split('.', 1)
    jobid_index = a[0]
    print jobid_index
    print 'jobid_indexr Submitted - Build Index'
else:
    print "Referece genome: reference genome does not exists, stopping processing pipeline!"
    exit (0)


################################################
####read samplesheet#####
# Function LoadSampleSheet
# Input: "file_name" is the name of the sample sheet file
# Output: a two dimensional list will be returned
# Call example: data = LoadSampleSheet("sampleSheet.txt");
# Specification: 1: sample sheet may contain comment lines, labeld with leading '#'
#                                2: each line must have 6 fieds
#                                3: the file should be tab-delimited
def LoadSampleSheet(file_name):
        try:
                fh = open(file_name, 'r')
        except:
                print "LoadSampleSheet: failed to open file %s." % file_name
        sample_info = []
        for line in fh:
                if line[0] == '#':
                        continue
                info = line.split()
                if(len(info) != 6):
                        print "LoadSampleSheet: error in table size, file may be corrupted."
                sample_info.append(info)
        fh.close()
        return sample_info


#################################################
####call mapping#####
 
# under the same directory
####fastq preprocessing
slist = LoadSampleSheet("sampleSheet.txt");
for i in range(0, len(slist)):
        if slist[i][0]== 'p':
		#print 'depend=afterok:'+jobid_index+''
		jobqsubpp1 = 'qsub -W depend=afterok:'+jobid_index+' -v fname='+slist[i][1]+',mycwd='+mydir+',adapter='+ADAPTER+',tr='+TRIM+' -N prepro_'+slist[i][1]+' /ingens/solid/analysis/yaping/pipeline2/Script/job_prepro.sh'
		proc = subprocess.Popen(jobqsubpp1, shell=True, stdout=subprocess.PIPE)
   		datapp1 = proc.stdout.readline()
   		a = datapp1.split('.', 1)
   		jobid_pp1 = a[0] 
		print ''+jobid_pp1+' submitted!'
		jobqsubpp2 = 'qsub -W depend=afterok:'+jobid_index+' -v fname='+slist[i][2]+',mycwd='+mydir+',adapter='+ADAPTER+',tr='+TRIM+' -N prepro_'+slist[i][2]+' /ingens/solid/analysis/yaping/pipeline2/Script/job_prepro.sh'
                proc = subprocess.Popen(jobqsubpp2, shell=True, stdout=subprocess.PIPE)
                datapp2 = proc.stdout.readline()
                a = datapp2.split('.', 1)
		jobid_pp2 = a[0]
		print ''+jobid_pp2+' submitted!'
		####mapping		  
		jobsubpe = 'qsub -W depend=afterok:'+jobid_pp1+':'+jobid_pp2+' -v fname1='+slist[i][1]+',fname2='+slist[i][2]+',align='+ALIGNER+',reffname='+GENOME+',mycwd='+mydir+',para='+PARA_MAP+',adapter='+ADAPTER+',tr='+TRIM+' -N mappe_'+slist[i][1]+' /ingens/solid/analysis/yaping/pipeline2/Script/job_pe.sh'
		proc = subprocess.Popen(jobsubpe, shell=True, stdout=subprocess.PIPE)
                datape = proc.stdout.readline()
		a = datape.split('.', 1)
                jobid_pe = a[0]
		print ''+jobid_pe+' submitted!'
	else:
   		print ''+jobid_pe+' preprocessing stopped!'
    		exit (0)




