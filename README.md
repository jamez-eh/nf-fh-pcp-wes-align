# nf-fh-pcp-wes-align
Alignment of exome data

## How To Run:
Make sure that Docker/Singularity and Nextflow are installed on your machine. On the fred hutch server we can use ```module load nextflow```

It is always a good idea to make sure that you have the most up to date pipeline: ```nextflow pull jamez-eh/nf-fh-pcp-wes-align```
Can then clone the repo, alter the parameters as necessary and then run ```bash run.sh```

Otherwise can be run by just downloading the ```run_from_git.sh``` and ```nextflow.config``` and then running
```bash run_from_git.sh```

## input format
   
   ### input_csv:
   ```sampleID,kitID,type,patientID,R1,R2```
   
   patientID and kitID can be left blank for this pipeline, but must be kept to maintain consistency with other pipelines. (it is probably ok if the columns don't exist in the manifest, but I have not tested this properly)
   The type column is important to indicate which samples will have the host reference removal step performed. A value of Tumor, or Normal will avoid the host reference removal, but a value of PDX will indicate a sample that requires this step.
   
   example:
   sample_1,,Tumor,,sample_1.R1,sample_1.R2


## Details:
   
    
   Pipeline does automatic subtraction of host exome if type of sample is specified as 'PDX' in ```input_csv```. If type is specified as 'Tumor' or 'Normal' no subtraction is performed. This subtraction is done by combining the reference and the host reference provided as parameters before aligning the PDX samples. The host reads can be removed via a sed command.

   Alignment is performed with BWA-MEM. 
   Duplicates are marked/removed and base recalibration is performed with the supplied ```rear``` file as high quality snp set. 
Post processing of bams details can be found here https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
   Does not do indel realignment since it was removed from GATK best practices.
   
   An additional quality filter is applied for a minimum base quality. The default is 20, but can be changed with the parameter ```--base_quality```

## Example Files:

   ### run.sh
   
   Here is an example of a run script with all the required files passed as parameters. This run script assumes that you have cloned this git repository and are running the pipeline from the repo.
   Alternatively you can run directly from the repo as shown in run_from_git.sh.
   Run with ```bash run.sh``` 

   ### run_from_git.sh
   
   Here is an example of running the pipeline without first cloning the repo.

   ## nextflow.config

   Here is an example of the nextflow configuration file, you will need to alter the AWS parameters.
   You will need to specify the location of your nextflow.config file in one the run script above 
   ``` -c location/of/nextflow.config ```

   The config contains three profiles, each is a different configuration for running the pipeline. The profile is specified in the ```run.sh```
   hpc : runs on slurm cluster
   aws : runs on aws batch
   local: runs on a single (large) node

   ### running on aws batch
   If you place ALL files (including references) on s3 and replace the necessary paths in both the ```input_csv``` and the ```run.sh``` then you can use the aws profile.
   Currently the provided queue is cpu-spot-50, but this queue may not be optimal for the job. Read more about AWS batch queues here : https://sciwiki.fredhutch.org/scicomputing/compute_cloud/

   

## Required Parameters

   The pipeline has several required files that are required
   
   ### pdx_reference: 
   Fasta file used for subtraction for PDX labelled files (see above).
   Only required if files with ```type``` PDX are included in the run.

   ### reference:
   Fasta file for main alignment.

   ### reference_index:
   Index generated by samtools for ```reference```.

   ### reference_dict:
   Dictionary for ```reference```. 
   
   ### rear:
   High quality SNP vcf used for base recalibration. Details see: https://gatk.broadinstitute.org/hc/en-us/articles/360042477672-BaseRecalibrator
   
   ### rear_index:
   Index file for rear
   
 
```
    --pdx_reference $BASE_BUCKET/references/mm10/GRCm38.primary_assembly.genome.fa \
    --reference $BASE_BUCKET/references/hg38/Homo_sapiens_assembly38.fa \
    --reference_index $BASE_BUCKET/references/hg38/Homo_sapiens_assembly38.fa.fai \
    --reference_dict $BASE_BUCKET/references/hg38/Homo_sapiens_assembly38.dict \
    --rear $BASE_BUCKET/references/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --rear_index $BASE_BUCKET/references/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi \
```


### Additional parameters.
       #### output_folder:
       Results output location

       #### ref_name:
       hg38 or hg19, used to download funcotating data. Details: https://gatk.broadinstitute.org/hc/en-us/articles/360037224432-Funcotator

       ### base_quality:
       minimum base quality cutoff.


```
    --ref_name hg38 \
    --output_folder /fh/scratch/delete90/nelson_p/james/bams \
    --base_quality 20 \
```