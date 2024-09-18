# snoRNA-chimeric
snoRNA-chimeric pipeline is designed to identify snoRNA-focused chimeric interactions.
## Description/methods:
* Demultiplexes reads using inline barcodes
* Trims adapters with cutadapt
* maps to repeat elements and genome with STAR to remove nonchimeric reads
* map filtered reads to snoRNA reference with bowtie to identify snoRNA chimeric reads
* filter the mapped chimeric reads including snoRNAs with given criteria (>18nt) with "find_putitavie_target"
* map the reads passed the cirteira to snoRNA again to remove snoRNA:snoRNA chimeirc reads
* map the non-snoRNA fragment from the reads to target RNA sequence(rRNA,tRNA,snRNA)
* map the unmapped reads to hg19 to identify other interactions
* identify the high-confident chimeric targets by filtering the mapped reads via given criteira (>16nt, gap between snoRNAs less than 4nt) with "identify_chimeric_reads"

## Installation:
### The pipeline has been tested using the following softwares and their versions:
* anytree=2.8.0
* bowtie=1.3.0
* clipper=2.0.0
* cutadapt=3.2
* python=2.7.16
pysam=0.16.0.1
numpy=1.19.5
pandas=1.2.1
* samtools=1.11
* star=STAR_2.4.0j
* ucsc-tools=377
* umi_tools=1.0.0
* 

## Workflow in commandline

### sno-chimeras
The following code shows how snoRNA chimeric reads were identified step by step:

```shell
$ cat snorna.chimeras.pipeline.sh

#!/usr/bin/env bash
set -e
source /storage/vannostrand/software/chimeras/venv/environment.sh
cd /storage/vannostrand/software/chimeras/tests/mouse
umi_tools extract \
  --random-seed 1 \
  --bc-pattern NNNNNNNNNN \
  --stdin /storage/vannostrand/software/chimeras/tests/snorna/S6_R1_001.fastq.gz \
  --log /dev/null | cutadapt \
  -j 32 \
  --match-read-wildcards \
  --times 3 \
  -e 0.1 \
  --quality-cutoff 6 \
  -m 28 \
  -a file:/storage/vannostrand/software/eclip/data/se.2.round.adapters.fasta \
  -O 1 \
  - \
  --report minimal \
  2> S6.cut.adapt.log | cutadapt \
  -j 32 \
  -u \
  -10 \
  -o S6.trim.fastq.gz \
  - \
  --report minimal \
  >> S6.cut.adapt.log \
  2> /dev/null
[ -d S6.repeat.map ] || mkdir S6.repeat.map
STAR \
  --runMode alignReads \
  --runThreadN 32 \
  --alignEndsType EndToEnd \
  --genomeDir /storage/vannostrand/genomes/repbase/species_specific/mus_musculus_repbase_v2 \
  --genomeLoad LoadAndRemove \
  --outBAMcompression 10 \
  --outFileNamePrefix S6.repeat.map/ \
  --outFilterMultimapNmax 100 \
  --outFilterMultimapScoreRange 1 \
  --outFilterScoreMin 10 \
  --outFilterType BySJout \
  --outReadsUnmapped Fastx \
  --outSAMattrRGline ID:foo \
  --outSAMattributes All \
  --outSAMmode Full \
  --outSAMtype BAM Unsorted \
  --outSAMunmapped None \
  --outStd Log \
  --readFilesCommand zcat \
  --readFilesIn S6.trim.fastq.gz
mv S6.repeat.map/Log.final.out S6.mask.repeat.map.log
mv S6.repeat.map/Unmapped.out.mate1 S6.repeat.unmap.fastq
rm -f S6.repeat.map
[ -d S6.genome.map ] || mkdir S6.genome.map
STAR \
  --runMode alignReads \
  --runThreadN 32 \
  --alignEndsType EndToEnd \
  --genomeDir /storage/vannostrand/genomes/mm10/yeolab/star_2_4_0i_gencode15_sjdb \
  --genomeLoad LoadAndRemove \
  --outBAMcompression 10 \
  --outFileNamePrefix S6.genome.map/ \
  --outFilterMultimapNmax 1 \
  --outFilterMultimapScoreRange 1 \
  --outFilterScoreMin 10 \
  --outFilterType BySJout \
  --outReadsUnmapped Fastx \
  --outSAMattrRGline ID:foo \
  --outSAMattributes All \
  --outSAMmode Full \
  --outSAMtype BAM Unsorted \
  --outSAMunmapped None \
  --outStd Log \
  --readFilesIn S6.repeat.unmap.fastq
mv S6.genome.map/Log.final.out S6.mask.genome.map.log
fastq_to_fasta S6.genome.map/Unmapped.out.mate1 S6.mask.fasta
rm -r S6.genome.map
bowtie2-build --quiet --offrate 2 /storage/vannostrand/software/chimeras/data/mm10/20220224/snoRNA.fa snoRNA.bowtie2.index 2> /dev/null
bowtie2 \
  -D 20 \
  -R 3 \
  -N 0 \
  -L 16 \
  --local \
  --norc \
  -i S,1,0.50 \
  --score-min L,16,0 \
  --ma 1 \
  --np 0 \
  --mp 2,2 \
  --rdg 5,1 \
  --rfg 5,1 \
  -p 32 \
  --no-unal \
  -x snoRNA.bowtie2.index \
  -f \
  -U S6.mask.fasta \
  -a \
  -S S6.mask.map.to.snoRNA.sam \
  2> S6.mask.map.to.snoRNA.log
split_source_map_sam S6.mask.map.to.snoRNA.sam S6.mask.map.to.snoRNA.AA.bam
find_putative_target S6.mask.map.to.snoRNA.AA.bam S6.mask.map.to.snoRNA.AA.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.AT.bam S6.mask.map.to.snoRNA.AT.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.AG.bam S6.mask.map.to.snoRNA.AG.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.AC.bam S6.mask.map.to.snoRNA.AC.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.AN.bam S6.mask.map.to.snoRNA.AN.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.TA.bam S6.mask.map.to.snoRNA.TA.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.TT.bam S6.mask.map.to.snoRNA.TT.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.TG.bam S6.mask.map.to.snoRNA.TG.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.TC.bam S6.mask.map.to.snoRNA.TC.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.TN.bam S6.mask.map.to.snoRNA.TN.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.GA.bam S6.mask.map.to.snoRNA.GA.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.GT.bam S6.mask.map.to.snoRNA.GT.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.GG.bam S6.mask.map.to.snoRNA.GG.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.NT.bam S6.mask.map.to.snoRNA.NT.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.GC.bam S6.mask.map.to.snoRNA.GC.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.CA.bam S6.mask.map.to.snoRNA.CA.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.NG.bam S6.mask.map.to.snoRNA.NG.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.GN.bam S6.mask.map.to.snoRNA.GN.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.CT.bam S6.mask.map.to.snoRNA.CT.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.CG.bam S6.mask.map.to.snoRNA.CG.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.CC.bam S6.mask.map.to.snoRNA.CC.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.CN.bam S6.mask.map.to.snoRNA.CN.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.NA.bam S6.mask.map.to.snoRNA.NA.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.NC.bam S6.mask.map.to.snoRNA.NC.target.fasta snoRNA False
find_putative_target S6.mask.map.to.snoRNA.NN.bam S6.mask.map.to.snoRNA.NN.target.fasta snoRNA False
cat S6.mask.map.to.snoRNA.[ATGCN][ATGCN].target.fasta > S6.mask.map.to.snoRNA.target.fasta
rm S6.mask.map.to.snoRNA.[ATGCN][TGCN].target.fasta || true
cat S6.mask.map.to.snoRNA.[ATGCN][ATGCN].target.csv > S6.mask.map.to.snoRNA.target.tmp.csv
grep -m 1 name S6.mask.map.to.snoRNA.target.tmp.csv > S6.mask.map.to.snoRNA.target.csv
grep -v name S6.mask.map.to.snoRNA.target.tmp.csv >> S6.mask.map.to.snoRNA.target.csv
rm S6.mask.map.to.snoRNA.target.tmp.csv
rm S6.mask.map.to.snoRNA.[ATGCN][ATGCN].target.csv || true
rm S6.mask.map.to.snoRNA.[ATGCN][ATGCN].bam || true
bowtie2 \
  -D 20 \
  -R 3 \
  -N 0 \
  -L 16 \
  --local \
  --norc \
  -i S,1,0.50 \
  --score-min L,12,0 \
  --ma 1 \
  --np 0 \
  --mp 2,2 \
  --rdg 5,1 \
  --rfg 5,1 \
  -p 32 \
  --no-unal \
  -x snoRNA.bowtie2.index \
  --un S6.mask.map.to.snoRNA.true.target.fasta \
  -f \
  -U S6.mask.map.to.snoRNA.target.fasta \
  -a \
  -S S6.mask.map.to.snoRNA.target.back.map.sam \
  2> S6.mask.map.to.snoRNA.target.back.map.log
bowtie-build --quiet --offrate 2 /storage/vannostrand/software/chimeras/data/mm10/20220224/rRNA.fa rRNA.bowtie.index 2> /dev/null
bowtie \
  -a \
  --best \
  --strata \
  -e 35 \
  -q \
  -l 8 \
  -n 1 \
  -p 32 \
  -x rRNA.bowtie.index \
  --no-unal \
  --norc \
  -f S6.mask.map.to.snoRNA.true.target.fasta \
  --sam S6.snoRNA.rRNA.sam \
  --un S6.snoRNA.rRNA.unmap.fasta \
  2> S6.snoRNA.rRNA.log
bowtie-build --quiet --offrate 2 /storage/vannostrand/software/chimeras/data/mm10/20220224/scaRNA.fa scaRNA.bowtie.index 2> /dev/null
bowtie \
  -a \
  --best \
  --strata \
  -e 35 \
  -q \
  -l 8 \
  -n 1 \
  -p 32 \
  -x scaRNA.bowtie.index \
  --no-unal \
  --norc \
  -f S6.mask.map.to.snoRNA.true.target.fasta \
  --sam S6.snoRNA.scaRNA.sam \
  --un S6.snoRNA.scaRNA.unmap.fasta \
  2> S6.snoRNA.scaRNA.log
bowtie-build --quiet --offrate 2 /storage/vannostrand/software/chimeras/data/mm10/20220224/snRNA.fa snRNA.bowtie.index 2> /dev/null
bowtie \
  -a \
  --best \
  --strata \
  -e 35 \
  -q \
  -l 8 \
  -n 1 \
  -p 32 \
  -x snRNA.bowtie.index \
  --no-unal \
  --norc \
  -f S6.mask.map.to.snoRNA.true.target.fasta \
  --sam S6.snoRNA.snRNA.sam \
  --un S6.snoRNA.snRNA.unmap.fasta \
  2> S6.snoRNA.snRNA.log
bowtie-build --quiet --offrate 2 /storage/vannostrand/software/chimeras/data/mm10/20220224/tRNA.fa tRNA.bowtie.index 2> /dev/null
bowtie \
  -a \
  --best \
  --strata \
  -e 35 \
  -q \
  -l 8 \
  -n 1 \
  -p 32 \
  -x tRNA.bowtie.index \
  --no-unal \
  --norc \
  -f S6.mask.map.to.snoRNA.true.target.fasta \
  --sam S6.snoRNA.tRNA.sam \
  --un S6.snoRNA.tRNA.unmap.fasta \
  2> S6.snoRNA.tRNA.log
compile_unmapped_putative_target S6 snoRNA rRNA,scaRNA,snRNA,tRNA S6.snoRNA.RNA.unmap.fasta
[ -d S6.snoRNA.mm10 ] || mkdir S6.snoRNA.mm10
STAR \
  --alignEndsType EndToEnd \
  --genomeDir /storage/vannostrand/genomes/mm10/yeolab/star_2_4_0i_gencode15_sjdb \
  --genomeLoad LoadAndRemove \
  --outBAMcompression 10 \
  --outFileNamePrefix S6.snoRNA.mm10/ \
  --outFilterMatchNminOverLread 0.66 \
  --outFilterMultimapNmax 1 \
  --outFilterMultimapScoreRange 1 \
  --outFilterScoreMin 10 \
  --outFilterScoreMinOverLread 0.66 \
  --outFilterType BySJout \
  --outReadsUnmapped Fastx \
  --outSAMattrRGline ID:foo \
  --outSAMattributes Standard \
  --outSAMmode Full \
  --outSAMtype SAM \
  --outSAMunmapped None \
  --outStd Log \
  --readFilesIn S6.snoRNA.RNA.unmap.fasta \
  --runMode alignReads \
  --runThreadN 32
mv S6.snoRNA.mm10/Aligned.out.sam S6.snoRNA.mm10.sam
mv S6.snoRNA.mm10/Log.final.out S6.snoRNA.mm10.log
rm -r S6.snoRNA.mm10
sam_to_bam S6.snoRNA.rRNA.sam S6.snoRNA.rRNA.bam True
sam_to_bam S6.snoRNA.scaRNA.sam S6.snoRNA.scaRNA.bam True
sam_to_bam S6.snoRNA.tRNA.sam S6.snoRNA.tRNA.bam True
sam_to_bam S6.snoRNA.snRNA.sam S6.snoRNA.snRNA.bam True
sam_to_bam S6.snoRNA.mm10.sam S6.snoRNA.mm10.bam True
identify_chimeric_read S6.snoRNA.rRNA.bam snoRNA.rRNA.scaRNA.snRNA.tRNA.chimeras.count.tsv S6.snoRNA.rRNA.chimeras.bam
identify_chimeric_read S6.snoRNA.scaRNA.bam snoRNA.rRNA.scaRNA.snRNA.tRNA.chimeras.count.tsv S6.snoRNA.scaRNA.chimeras.bam
identify_chimeric_read S6.snoRNA.snRNA.bam snoRNA.rRNA.scaRNA.snRNA.tRNA.chimeras.count.tsv S6.snoRNA.snRNA.chimeras.bam
identify_chimeric_read S6.snoRNA.tRNA.bam snoRNA.rRNA.scaRNA.snRNA.tRNA.chimeras.count.tsv S6.snoRNA.tRNA.chimeras.bam
identify_chimeric_read S6.snoRNA.mm10.bam snoRNA.rRNA.scaRNA.snRNA.tRNA.chimeras.count.tsv S6.snoRNA.mm10.chimeras.bam
bam_to_bw S6.snoRNA.rRNA.chimeras.bam
bam_to_bw S6.snoRNA.scaRNA.chimeras.bam
bam_to_bw S6.snoRNA.snRNA.chimeras.bam
bam_to_bw S6.snoRNA.tRNA.chimeras.bam
bam_to_bw S6.snoRNA.mm10.chimeras.bam
clipper \
  --species mm10 \
  --processors 32 \
  --bam S6.snoRNA.mm10.chimeras.bam \
  --outfile S6.snoRNA.mm10.chimeras.peak.bed
sno_chimeras_summary \
  --stag snoRNA \
  --tags rRNA,scaRNA,snRNA,tRNA,mm10 \
  --ids S6, \
  --output snoRNA.rRNA.scaRNA.snRNA.tRNA.mm10.chimeras.report.html \
  --count snoRNA.rRNA.scaRNA.snRNA.tRNA.chimeras.count.tsv \
  --gtag mm10
rm S6.*.trim.fastq.gz || true
rm S6.*.target.* || true
rm S6.*.fastq || true
rm S6.*.tsv || true
rm S6.*.duplicated.csv || true
rm S6.*.sam || true
rm S6.[ATGCN][ATGCN].bam || true
rm S6.candidate.bam* || true
rm -r S6*.repeat.map || true
rm -r S6*.genome.map || true
rm -r S6*.genome.length.tsv || true
rm -r S6*.bw.*.bg || true
rm -r S6*.snoRNA.mm10 || true
tail -n +1 S6*.log > map.metric.and.log
rm S6*.log || true
