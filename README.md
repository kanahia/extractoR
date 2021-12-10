# extractoR
Script extracts sequences based on PCR primers from the transcriptome,
generate fasta files which are later on mapped to the genome,
converted to bed files and ready for visualization in genome browser.

args1 = path to the template

args2 = path to the transcriptome

args3 = path to the output dir

args4 = path to STAR genome (index)

./extractoR_pipe.sh args1 args2 args3 args4 

Mandatory fields that have to be provided in the template form:

1. Gene name

2. Primer sequences

3. Ensembl gene ID
