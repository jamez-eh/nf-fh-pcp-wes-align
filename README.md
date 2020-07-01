# nf-fh-pcp-wes-align
Alignment of exome data

## input format
   
   ### input_csv:
   ```sampleID,kitID,type,patientID,R1,R2```
   
   patientID and kitID can be left blank for this pipeline, but must be kept to maintain consistency with other pipelines.
   
   example:
   sample_1,,Tumor,,sample_1.R1,sample_1.R2

## Details:
   
   Pipeline does automatic subtraction of xenograph genome specified if params.pdx = true.
   Does not do indel realignment since it was removed from GATK best practices.

## Example Files:

   ### run.sh

   Here is an example of a run script with all the required files passed as parameters

   ## nextflow.config

   Here is an example of the nextflow configuration file, you will need to alter the AWS parameters.
   You will need to specify the location of your nextflow.config file in the run script above 
   ``` -c location/of/nextflow.config ```
   