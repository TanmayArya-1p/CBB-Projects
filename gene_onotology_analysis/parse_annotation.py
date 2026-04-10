valid_chroms = set({str(i) for i in range(1, 23)} | {'X', 'Y', 'MT'})

with open('human_gene_annotation.tsv', 'r') as fin, open('tss_start_sites.bed', 'w') as fout:
    for line in fin.readlines()[1:]:
        cols = line.strip().split()
        if len(cols) < 8:
            continue

        raw_chrom = cols[4]
        if raw_chrom not in valid_chroms:
            continue

        chrom = "chr" + raw_chrom
        if chrom == "chrMT":
            chrom = "chrM"

        strand = "+" if cols[5] == "1" else "-"
        tss = int(cols[7])
        gene_name = cols[6]

        fout.write(f"{chrom}\t{tss}\t{tss + 1}\t{gene_name}\t.\t{strand}\n")
