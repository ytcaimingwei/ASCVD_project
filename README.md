# ASCVD_project
Gut bacteria derived odd chain fatty acid modulates cholesterol homeostasis and alleviates atherosclerosis


I. Megatnomic analysis

1.	Quality control of metagenomic reads

sickle se -t sanger -q 25 -f $entry -o $entry.t.fq

2.	Remove human reads

bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=/lustre/home/mwcai/data/remove_human_reads_from_MG_MT qtrim=rl trimq=10 untrim -Xmx23g in=$entry outu=$entry.clean.fq outm=$entry.human.fq

3.	Fastq to fasta format
seqtk seq -a  $entry > $entry.fa

4.	Taxonomic and abundance of the high-quality reads
metawrap kraken2 -t 32 -o kraken2_SRR11461968 SRR11461968.fastq.t.fq.clean.fq.fa

5.	Assessment of PA synthesis gene abundance in the two published cohorts
coverm contig --single $entry -r mutA.fa --min-read-percent-identity 95 --min-read-aligned-percent 50  -o $entry.coverm.mutA






II. Identification of PA synthesis genes in Unified Human Gastrointestinal Genome (UHGG) collection

blastp -db gene.fasta.faa -query NBT_total.faa -outfmt 6 -out NBT_total.faa_1E3_50.blastp -evalue 1E-3 -num_threads 64 -qcov_hsp_perc 50
