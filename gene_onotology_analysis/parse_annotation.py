import gzip

with gzip.open('human_gene_annotation.tsv.gz', 'rt') as fin, open('promoters.bed', 'w') as fout:
    next(fin)
    for line in fin:
        cols = line.strip().split('\t')
        if len(cols) < 8:
            continue
        chrom = "chr" + cols[4]
        if chrom == "chrMT":
            chrom = "chrM"

        strand = "+" if cols[5] == "1" else "-"
        tss = int(cols[7])
        gene_name = cols[6]

        if strand == "+":
            start = max(0, tss - 1000)
            end = tss
        else:
            start = tss
            end = tss + 1000

        fout.write(f"{chrom}\t{start}\t{end}\t{gene_name}\t.\t{strand}\n")
