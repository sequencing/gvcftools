bin_dir=../../bin
data_dir=data

$bin_dir/twins --ref /illumina/scratch/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
--twin1 $data_dir/twin1.vcf.gz \
--twin2 $data_dir/twin2.vcf.gz \
--no-variable-metadata
