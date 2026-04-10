.PHONY: clean all

all: gene_enrichment_plot.png

tss_start_sites.bed: parse_annotation.py human_gene_annotation.tsv
	python parse_annotation.py

promoters.bed: tss_start_sites.bed hg38.chrom.sizes
	bedtools slop -i tss_start_sites.bed -g hg38.chrom.sizes -l 500 -r 0 -s > promoters.bed

promoters.fa: promoters.bed hg38.fa
	bedtools getfasta -fi hg38.fa -bed promoters.bed -s -nameOnly > promoters.fa

cre_hits.gff: promoters.fa
	dreg -sequence promoters.fa -pattern "TGACGTCA" -outfile cre_hits.gff -rformat gff

gene_list.txt: cre_hits.gff parse_hits.py
	python parse_hits.py

gene_enrichment_plot.png: gene_list.txt run_go.R
	Rscript run_go.R

human_gene_annotation.tsv: human_gene_annotation.tsv.gz
	gunzip -k human_gene_annotation.tsv.gz

hg38.fa: hg38.fa.gz
	gunzip -k hg38.fa.gz

hg38.chrom.sizes:
	wget -O hg38.chrom.sizes https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes


clean:
	rm -f tss_start_sites.bed promoters.bed promoters.fa cre_hits.gff gene_list.txt gene_enrichment_plot.png hg38.chrom.sizes hg38.fa hg38.fa.fai
