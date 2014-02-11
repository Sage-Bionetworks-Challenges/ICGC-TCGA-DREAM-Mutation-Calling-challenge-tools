ICGC-TCGA-DREAM-Mutation-Calling-challenge-tools
================================================

Tools for participants in the ICGC-TCGA DREAM Mutation Calling challenge

You will find here the tools and step-by-step instructions for uploading result files 
for participating in the The ICGC-TCGA DREAM Genomic Mutation Calling Challenge 
(referred herein as The Challenge). You must first 
[join the Challenge](https://www.synapse.org/#!Synapse:syn2298928/wiki/60977) and 
[be approved by the ICGC DACO](https://www.synapse.org/#!Synapse:syn2298928/wiki/60980) 
before you can submit.


# Submission steps

1. Obtain and install submission tool
2. Create a private working project
3. Run the validator on your VCF file
4. Ensure that there are no Germline calls present in the VCF
5. Upload and submit the entity to the evaluation


## Obtain and install submission tool
We provide a command line program to help validate and submit files to the contest.
First, you will need to install the program.
There are a few dependencies that are needed in order for the submission program to work. 
These package help to parse VCF files and contact the Synapse servers.

You can use the standard Python ‘pip’ installer to get these dependencies
```
pip install PyVCF
pip install synapseclient
```

or, to install in your home directory:
```
pip install --user PyVCF
```

Get a copy of the Mutation Calling contestant toolkit: 
https://github.com/Sage-Bionetworks/ICGC-TCGA-DREAM-Mutation-Calling-challenge-tools


```
git clone https://github.com/Sage-Bionetworks/ICGC-TCGA-DREAM-Mutation-Calling-challenge-tools.git
```


##Create a private working project
If you don’t already have a personal working project, go to the front page of Synapse and 
use the ‘Create Project’ widget. Please make sure that this folder is not visible to the 
public. There is more [info available](https://www.synapse.org/#!ProjectsHome:0) if you 
want to find out about Synapse projects.


##Run the validator on your VCF file


The ‘dream_submit’ program includes a ‘validate’ function, to check your VCF files for 
valid formatting.

```
dream_submit validate VarScan.HCC1143.mix1.n60t40.Somatic.vcf.gz 
```


####Submission format
VCF format, compressed with gzip. It should be able to be indexed using tabix 
( http://samtools.sourceforge.net/tabix.shtml )
For evaluation, the tumor source file needs to be identified. This is in the form of 
the UUID assigned to the BAM file from CGHub. 

The submission format is VCF 4.1, which you can find the spec for 
[online](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41)

For all submissions
* The file must end with .vcf or .vcf.gz
* All variants must be marked somatic as per methods defined in the spec (preferably using the SOMATIC tag in the INFO field).
* No whitespaces in fields, fields must be tab-delimited (this is also in the spec, but 
it's a common reason for the vcf parser to fail).


For SVs:

* If END is not specified in the INFO column, END is assumed to be the same as POS
(i.e. a single breakend)
* If a variant is marked 'IMPRECISE' it must have confidence intervals defined by 
'CIPOS' and 'CIEND' (the latter only if the 'END' field is used in the INFO column).


##Ensure that there are no Germline calls present in the VCF

Please do not include Germline mutations in your entry. This challenge is centered on 
the prediction of cancer level event. In the case of real patient genome data, germline
mutation calls represent a possible identifiable genomic characteristic, and should be 
kept private in compliance to the data usage agreements. 
The VCF validator function in the dream_submit program should identify germline calls. 

Note: There is a [Python Bug](http://bugs.python.org/issue17666) that breaks compressed 
VCF parsing in Python version 2.7.3 and 2.7.4. Please upgrade if the 'dream_submit' 
script fails during VCF validation.

##SynapseClient Login
You can either pass you name and password to the dream_submit program via the command 
line, or to can follow the instructions setting up an 
[authentication config file](https://www.synapse.org/#!Synapse:syn1768504/wiki/56068)

You can also visit [the settings page](https://www.synapse.org/#!Settings:0) to get an API 
key, then set the environmental variables:

```
export SYNAPSE_APIKEY=(really long string you copied from the setting page)
export SYNAPSE_EMAIL=(your email address)
```


##Upload and submit the entity to the evaluation

There are two leaderboards for the challenge: one for SNV calls and one for structural 
variant calls. Submit separate VCF files to these leaderboards using the ‘dream_submit’ 
program in the command line submission tool.

To see commands to submit a structural variant VCF:
```
dream_submit structural_variant -h
```
To see commands to submit a SNV VCF:
```
dream_submit snv -h
```
The arguments to the SNV submission command are shown below. Arguments for the structural 
variant submission are identical:

```
dream_submit snv -h
usage: dream_submit submit [-h] [--vcf VCF] [--entity ENTITY] --name NAME
                           --team-name TEAM_NAME [--project-id PROJECT_ID]
                           [--test]
                           sample_name

positional arguments:
  sample_name           Sample name (or a UUID of the tumor file) of the tumor data source.

optional arguments:
  -h, --help            show this help message and exit
  --vcf VCF             Path to the VCF file to be submitted, if not already uploaded to Synapse.
  --entity ENTITY       Synapse Entity ID of VCF file to be submitted, if already uploaded to Synapse.
  --name NAME           Name of the submission (REQUIRED).
  --team-name TEAM_NAME
                        Name Team (REQUIRED).
  --project-id PROJECT_ID
                        The SYN id of your personal private working directory.
```


Example usage:
```
dream_submit  --user (your user name) --password (your password) snv synthetic.challenge.set1.v2 --vcf VarScan.Somatic.vcf.gz --project-id (your private synapse project) --name "Varscan_test" --team-name "My Dream Team"
```

Alternative, if you have set up an [authentication config file](https://www.synapse.org/#!Synapse:syn1768504/wiki/56068):
```
dream_submit structural_variant synthetic.challenge.set1.v2 --vcf Breakdancer.Somatic.vcf.gz --project-id (your private synapse project) --name "breakdancer_test" --team-name "My Dream Team"
```

Or if your VCF is already stored on Synapse:
```
dream_submit structural_variant synthetic.challenge.set1.v2 --name "MyMutationPrediction" --entity syn123456 --team-name "My Dream Team" 
```



##Viewing the results
The results can be viewed at syn2177211. Note, this page is updated periodically, and your 
entry won't instantaneously appear.


#Additional helper functions

##List samples

You can use the dream_submit program to get a listing of the different source files that are 
available for processing for each sample with the command:

```
dream_submit list-sample
```

This list of source files will be mirrored at syn2280639


##List challenges

You can view a list of the active challenges with the command

```
dream_submit list-challenge
```

Which will produce a table like:

```
SampleName	snv Evaluation	sv Evaluation
synthetic.challenge.set1.v2	2177213	2294712
```

You don't need to use the Evaluation id numbers, but you can use the SampleName to help 
guide your submission to the correct contest.

