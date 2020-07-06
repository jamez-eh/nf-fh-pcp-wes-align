# nf-fh-pcp-wes-align
Alignment of exome data

## input format
   
   ### input_csv:
   ```sampleID,kitID,type,patientID,R1,R2```
   
   patientID and kitID can be left blank for this pipeline, but must be kept to maintain consistency with other pipelines.
   
   example:
   sample_1,,Tumor,,sample_1.R1,sample_1.R2

## Details:
   
    
   Pipeline does automatic subtraction of host exome if type of sample is specified as 'PDX' in ```input_csv```. If type is specified as 'Tumor' or 'Normal' no subtraction is performed. This subtraction is done by combining the reference and the host reference provided as parameters before aligning the PDX samples. The host reads can be removed via a sed command.

   Alignment is performed with BWA-MEM. 
   Duplicates are marked/removed and base 
Post processing of bams details can be found here https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
   Does not do indel realignment since it was removed from GATK best practices.
   

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

   

   