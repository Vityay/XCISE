# XCISE
Calling X-Chromosome Inactivation from Single-cell Expression data (XCISE)

## Prerequisites:
* Samtools (https://github.com/samtools/samtools)

* Bcftools (https://github.com/samtools/bcftools)

* STAR aligner (https://github.com/alexdobin/STAR)

## Example of analysis 10x Chromium scRNA-seq data using sample 133C of IPF cell atlas (GSE136831)

### Step 1. Download all data needed

#### Step 1.1. Get genome reference and gene build

`wget http://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`

`wget http://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz`

#### Step 1.2 Unpack them

`gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`

`gunzip Homo_sapiens.GRCh38.109.gtf.gz`

#### Step 1.3. Index genome for STAR aligner

`STAR --runMode genomeGenerate --genomeDir GRCh38_STAR --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.109.gtf --sjdbOverhang 97 --runThreadN 16`

#### Step 1.4. Get scRNA-seq data for sample 133C

`wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/070/SRR21407770/SRR21407770_1.fastq.gz`

`wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/070/SRR21407770/SRR21407770_2.fastq.gz`

### Step 2. Aligning

#### Step 2.1. Copy or link folder 'whitelists' from STAR installation to current directory

#### Step 2.2. Run STAR alignment 

`STAR --genomeDir GRCh38_STAR --readFilesIn SRR21407770_2.fastq.gz SRR21407770_1.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --soloType Droplet --soloCBwhitelist whitelists/737K-august-2016.txt --outFileNamePrefix 133C --outSAMmapqUnique 50 --outSAMstrandField intronMotif --outSAMattrRGline 'ID:133C SM:133C LB:133C PL:Illumina' --runThreadN 16`

#### Step 2.3. Index BAM file

`samtools index -@ 16 133CAligned.sortedByCoord.out.bam`

### Step 3. Calling variants

#### Step 3.1. Run BCFtools 

`bcftools mpileup -a AD -f Homo_sapiens.GRCh38.dna.primary_assembly.fa -q 20 -d 1000000 -r X 133CAligned.sortedByCoord.out.bam | bcftools call -mv -Oz -o 133C_variants.vcf.gz`

#### Step 3.2 Normalize VCF file

`bcftools norm -m-any --check-ref -w -f Homo_sapiens.GRCh38.dna.primary_assembly.fa 133C_variants.vcf.gz -Oz -o 133C_variants_normalized.vcf.gz`

#### Step 3.3. Index VCF file

`bcftools index 133C_variants_normalized.vcf.gz`

#### Step 3.4. Overlap variants that are not common in population

`bcftools isec -i'QUAL>=200' 133C_variants_normalized.vcf.gz GRCh38_Common_Xlinked_SNVs.vcf.gz -p 133C_common_SNVs`

#### Step 3.5. Reformat VCF file for use with STAR aligner in WASP mode

`perl vcf4wasp.pl 133C_common_SNVs/0002.vcf >133C_heterozygous_SNVs.vcf`

### Step 4. Quantify allele-specific expression

#### Step 4.1. Run STAR alignment 

`STAR --genomeDir GRCh38_STAR --readFilesIn SRR21407770_2.fastq.gz SRR21407770_1.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --soloType Droplet --soloCBwhitelist whitelists/737K-august-2016.txt --outFileNamePrefix 133C_WASP --varVCFfile 133C_heterozygous_SNVs.vcf --waspOutputMode SAMtag --outSAMattributes vA vG NH HI AS NM nM CB UB --outSAMmapqUnique 50 --outSAMstrandField intronMotif --outSAMattrRGline 'ID:133C SM:133C LB:133C PL:Illumina' --runThreadN 16`

#### Step 4.2 Index Alignment

`samtools index -@ 16 133C_WASPAligned.sortedByCoord.out.bam`

### Step 5. Run XCISE

`perl xcise.pl -o 133C -s 133C_heterozygous_SNVs.vcf -b 133C_WASPAligned.sortedByCoord.out.bam`


## Example of Smart-seq scRNA-seq data analysis using sample (E-MTAB-6385)

### Step 1. Download all data needed

#### Step 1.1. Get genome reference and gene build

`wget http://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz`

`wget http://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz`

#### Step 1.2 Unpack them

`gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz`

`gunzip Mus_musculus.GRCm39.110.gtf.gz`

####Step 1.3. Index genome for STAR aligner

`STAR --runMode genomeGenerate --genomeDir GRCm39_STAR --genomeFastaFiles Mus_musculus.GRCm39.dna.primary_assembly.fa --sjdbGTFfile Mus_musculus.GRCm39.110.gtf --sjdbOverhang 42 --runThreadN 16`

#### Step 1.4. Get scRNA-seq data for sample 133C

`for i in {2717189..2717786}; do
    wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR271/ERR${i}/*.fastq.gz
done` 

### Step 2. Aligning

#### Step 2.1. Prepare manifest file to associate each file with cellID
We created a tab-delinmited text file E-MTAB-6385_manifest.tsv with 3 columns with top 10 lines looking like this:

`Fib_BxC_p1_a10_R1.batch1.fastq.gz	-	p1_a10`

`Fib_BxC_p1_a10_R1.fastq.gz	-	p1_a10`

`Fib_BxC_p1_a11_R1.batch1.fastq.gz	-	p1_a11`

`Fib_BxC_p1_a11_R1.fastq.gz	-	p1_a11`

`Fib_BxC_p1_a12_R1.batch1.fastq.gz	-	p1_a12`

`Fib_BxC_p1_a12_R1.fastq.gz	-	p1_a12`

`Fib_BxC_p1_a1_R1.batch1.fastq.gz	-	p1_a1`

`Fib_BxC_p1_a1_R1.fastq.gz	-	p1_a1`

`Fib_BxC_p1_a2_R1.batch1.fastq.gz	-	p1_a2`

`Fib_BxC_p1_a2_R1.fastq.gz	-	p1_a2`

#### Step 2.2. Run STAR alignment

`STAR --genomeDir GRCm39_STAR --readFilesManifest E-MTAB-6385_manifest.tsv --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --soloType SmartSeq --bamRemoveDuplicatesType UniqueIdentical --soloUMIdedup NoDedup --outFileNamePrefix Fib_BxC --outSAMmapqUnique 50 --outSAMstrandField intronMotif --runThreadN 16`

#### Step 2.3. Index BAM file

`samtools index -@ 16 Fib_BxCAligned.sortedByCoord.out.bam`

### Step 3. Calling variants

#### Step 3.1. Run BCFtools

`bcftools mpileup -a AD -f Mus_musculus.GRCm39.dna.primary_assembly.fa -q 20 -d 1000000 -r X Fib_BxCAligned.sortedByCoord.out.bam | bcftools call -mv -Oz -o Fib_BxC_variants.vcf.gz`

#### Step 3.2 Normalize VCF file

`bcftools norm -m-any --check-ref -w -f Mus_musculus.GRCm39.dna.primary_assembly.fa Fib_BxC_variants.vcf.gz -Oz -o Fib_BxC_variants_normalized.vcf.gz`

#### Step 3.3. Index VCF file

`bcftools index Fib_BxC_variants_normalized.vcf.gz`

#### Step 3.4. Overlap variants that are not common in population *

`bcftools isec -i'QUAL>=200' Fib_BxC_variants_normalized.vcf.gz GRCm39_Common_Xlinked_SNVs.vcf.gz -p Fib_BxC_common_SNVs`

#### Step 3.5. Reformat VCF file for use with STAR aligner in WASP mode

`perl vcf4wasp.pl Fib_BxC_common_SNVs/0002.vcf >Fib_BxC_heterozygous_SNVs.vcf`

### Step 4. Quantify allele-specific expression

#### Step 4.1. Run STAR alignment

`STAR --genomeDir GRCm39_STAR --readFilesManifest E-MTAB-6385_manifest.tsv --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --soloType SmartSeq --bamRemoveDuplicatesType UniqueIdentical --soloUMIdedup NoDedup --outFileNamePrefix Fib_BxC_WASP --varVCFfile Fib_BxC_heterozygous_SNVs.vcf --waspOutputMode SAMtag --outSAMattributes vA vG RG --outSAMmapqUnique 50 --outSAMstrandField intronMotif --runThreadN 16`

#### Step 4.2 Index alignment

`samtools index -@ 16 Fib_BxC_WASPAligned.sortedByCoord.out.bam`

### Step 5. Run XCISE

`perl xcise.pl -o Fib_BxC -s Fib_BxC_heterozygous_SNVs.vcf -b Fib_BxC_WASPAligned.sortedByCoord.out.bam`

## Notes:

- When using multiple file pairs in 10x processing, join them with commas during STAR alignment steps ( 2.2 and 5.1 ), e.g.:

`STAR --genomeDir GRCh38_STAR --readFilesIn SRR21407781_2.fastq.gz,SRR21407782_2.fastq.gz SRR21407781_1.fastq.gz,SRR21407782_1.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --soloType Droplet --soloCBwhitelist whitelists/737K-august-2016.txt --outFileNamePrefix 098C --outSAMmapqUnique 50 --outSAMstrandField intronMotif --outSAMattrRGline "ID:098C SM:098C LB:098C PL:Illumina" --runThreadN 16`

- Here you first need to create a file with X-linked SNVs common in human population. We supply such a file with this code (), but you could also create it from the
latest version of gnomAD database (https://gnomad.broadinstitute.org) by downloading VCF file for X chromosome and running extraction of common SNVs:

`perl -e 'print "##fileformat=VCF\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\n"' >GRCh38_Common_Xlinked_SNVs.vcf`

`bcftools view -e 'AF_non_cancer<0.01' gnomad.genomes.v3.1.2.sites.chrX.vcf.bgz | grep 'allele_type=snv' | perl -pe 's/^chr// ' | cut -s -f1-5 | perl -pe 's/([GATC])\n$/$1\t300\n/' >>Common_Xlinked_SNVs.vcf`

`bgzip GRCh38_Common_Xlinked_SNVs.vcf`

`bcftools index GRCh38_Common_Xlinked_SNVs.vcf.gz`

- A similar task for common SNVs in Mouse inbred strains:
  
`perl -e 'print "##fileformat=VCF\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\n"' >GRCm39_Common_Xlinked_SNVs.vcf`

`bcftools view -r X mgp_REL2021_snps.rsID.vcf.gz | cut -s -f1-5 | perl -pe 's/([GATC])\n$/$1\t300\n/' >>GRCm39_Common_Xlinked_SNVs.vcf`

`bgzip Common_Xlinked_SNVs.vcf`

`bcftools index GRCm39_Common_Xlinked_SNVs.vcf.gz`

- When using multiple file pairs in 10x processing, join them with commas during STAR alignment steps ( 2.2 and 4.1 ), e.g.:
  
`STAR --genomeDir GRCh38_STAR --readFilesIn SRR21407781_2.fastq.gz,SRR21407782_2.fastq.gz SRR21407781_1.fastq.gz,SRR21407782_1.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --soloType Droplet --soloCBwhitelist whitelists/737K-august-2016.txt --outFileNamePrefix 098C --outSAMmapqUnique 50 --outSAMstrandField intronMotif --outSAMattrRGline "ID:098C SM:098C LB:098C PL:Illumina" --runThreadN 16`
