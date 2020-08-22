#!/bin/bash

path_template=$1
path_transcriptome=$2
path_out=$3
path_genome=$4

Rscript extractoR.R $path_template $path_transcriptome $path_out

cd $path_out/output_dir
mkdir star_probes
for files in ./*.fa; do 
	name=`basename $files`
	fasta=$name
	STAR --runThreadN 10 --genomeDir $path_genome --readFilesIn ./${fasta} --outFileNamePrefix ./star_probes/${name%.fa}.star --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts 
done

#convert bam to bed (intron+exons and only exons)
mkdir probes_bed

for x in $path_out/output_dir/star_probes/*.bam ; do
	echo "print current:$x";
	bedtools bamtobed -i "$x" > "${x%.bam}.bed"
	bedtools bamtobed -split -i "$x" > "${x%.bam}.split.bed";
done
echo "done"

mkdir $path_out/output_dir/probes_fasta
mv $path_out/output_dir/*fa probes_fasta
mv star_probes/*.bed $path_out/output_dir/probes_bed

