
#!/bin/bash 

### Change these settings as needed ###

# Read file information
read_dir="/g/steinmetz/project/IESY/sequencing/Data/PacBio/XEMBL.20161504.PACBIO_DATA"
read_file="XEMBL_20160408_RS42150_PL100065328A-1_H01-1_lbc4.fastq"
read_ID="JS599"

# Reference file information
ref_dir="/g/steinmetz/genome/Saccharomyces_cerevisiae/"
ref_file="S288c_SynIXR/SynIXR_S299C_combined.fsa"

# Results information
result_dir="/g/steinmetz/project/IESY/sequencing/Results/PacBio/"

# Append sample ID (well information etc.)
sample_ID_append="-F3-lowCoverage"

# Aligner arguments
bestn="2"
nproc="16" #number of spawned processes

# Samtools arguments
mem_alloc="20000000000" #in bytes

# Chromosome ID for reads of interest
chr_ID="gi|346228209|gb|JN020955.1|" #SynIXR in this case

# Assembler arguments
genome_size="155280" #in bytes
seq_lib_type="pacbio-raw" #library input type
use_grid="0" #attach to grid
non_default=true #true will invoke the non-default settings below (useful for low coverage assemblies)
corMhapSensitivity="high"
corMinCoverage="2"
errorRate="0.035"
minOverlapLength="100"
corMaxEvidenceErate="0.3"

### Run pipeline ###

# Read input full path
read_in="$read_dir/Demux/$read_ID/$read_file"

# Reference input full path
ref_in=$ref_dir$ref_file

# Sample ID
sample_ID=$read_ID$sample_ID_append

# Aligned reads output
align_append="_PacBio_all_reads_aligned"
align_out=$sample_ID$align_append

# Extracted reads output
extract_append="_PacBio_SynIXR_only_sorted"

# Switch between default and non-default settings for the assembly
if [ $non_default == true ]; 
    then
        low_coverage_settings="corMhapSensitivity=$corMhapSensitivity corMinCoverage=$corMinCoverage \
        errorRate=$errorRate minOverlapLength=$minOverlapLength corMaxEvidenceErate=$corMaxEvidenceErate"
        echo "Using non-default settings for the assembly"
    else
        low_coverage_settings=""
        echo "Using default settings for the assembly"
fi

# Create results directories
cd $result_dir
mkdir ./$sample_ID
mkdir ./$sample_ID/Ref_align
mkdir ./$sample_ID/Ref_align/All_reads
mkdir ./$sample_ID/Canu_assembly
mkdir ./$sample_ID/Canu_assembly/SynIXR_only

# Align reads using blasr
cd $result_dir$sample_ID/Ref_align/All_reads
~/.linuxbrew/bin/blasr $read_in $ref_in -bestn $bestn -sam -nproc $nproc -out $align_out.sam & wait

echo "Completed collecting reads"

# Convert to bam, sort and index
cd $result_dir$sample_ID/Ref_align/All_reads
/g/software/bin/samtools view -b -S $align_out.sam > $align_out.bam
/g/software/bin/samtools sort -m $mem_alloc $align_out.bam -o ${align_out}_sorted.bam
/g/software/bin/samtools index ${align_out}_sorted.bam

# Extract reads mapping to SynIXR
/g/software/bin/samtools view -b ${align_out}_sorted.bam "$chr_ID" > $sample_ID$extract_append.bam

# Convert to fasta
/g/software/bin/samtools view $sample_ID$extract_append.bam | awk '{OFS="\t"; print ">"$1"\n"$10}' - > \
$sample_ID$extract_append.fasta

# Set JAVA environment and run assembler
cd /g/steinmetz/project/IESY/Tools/environ
. set_java-1.8 

cd /g/steinmetz/project/IESY/Tools/canu/Linux-amd64/bin/

./canu -p JS599-F3-gDNA_SynIXR_only -d "$result_dir$sample_ID/Canu_assembly/SynIXR_only" \
genomeSize=$genome_size -$seq_lib_type "$result_dir$sample_ID/Ref_align/All_reads/$sample_ID$extract_append.fasta" \
useGrid=$use_grid $low_coverage_settings

echo "Completed assembly"
