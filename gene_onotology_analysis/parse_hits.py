unq_genes = set()

with open('cre_hits.gff', 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        gene = line.split('\t')[0]
        gene = gene.replace('(+)', '').replace('(-)', '')
        unq_genes.add(gene)

with open('gene_list.txt', 'w') as f:
	f.writelines([i + "\n" for i in sorted(unq_genes)])
