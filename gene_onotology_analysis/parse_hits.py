unq_genes = set()

with open('cre_hits.gff', 'r') as f:
    for line in f.readlines():
        if line.startswith('#') or not line.strip():
            continue
        gene = line.split('\t')[0]
        unq_genes.add(gene)

with open('gene_list.txt', 'w') as f:
    for gene in sorted(unq_genes):
        f.write(f"{gene}\n")
