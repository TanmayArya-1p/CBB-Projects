run: gene_enrichment_plot.pdf

hg38.fa:
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
	gunzip hg38.fa.gz

promoters.bed: human_gene_annotation.tsv.gz parse_annotation.py
	python parse_annotation.py

promoters.fasta: hg38.fa promoters.bed
	bedtools getfasta -fi hg38.fa -bed promoters.bed -s -nameOnly > promoters.fasta

cre_hits.gff: promoters.fasta
	dreg -sequence promoters.fasta -pattern "TGACGTCA" -outfile cre_hits.gff -rformat gff

gene_list.txt: cre_hits.gff parse_hits.py
	python parse_hits.py
	sed -i 's/(+)//g; s/(-)//g' gene_list.txt

gene_enrichment_plot.pdf: gene_list.txt run_go.R
	Rscript run_go.R

clean:
	rm -f promoters.bed promoters.fasta cre_hits.gff gene_list.txt gene_enrichment_plot.pdf
