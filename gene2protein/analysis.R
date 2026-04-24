.libPaths(c("/usr/local/lib/R/site-library", "/usr/lib/R/library", .libPaths()))
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(IRanges)

edb <- EnsDb.Hsapiens.v86
sftpa1_transcripts <- transcripts(
  edb,
  filter      = GeneNameFilter("SFTPA1"),
  return.type = "DataFrame"
)

print(sftpa1_transcripts)
write.csv(
  as.data.frame(sftpa1_transcripts),
  file      = "transcripts.csv",
  row.names = FALSE
)

sftpa1_proteins <- proteins(
  edb,
  filter = GeneNameFilter("SFTPA1")
)

print(sftpa1_proteins)
write.csv(
  as.data.frame(sftpa1_proteins),
  file      = "proteins.csv",
  row.names = FALSE
)

prot_ids     <- sftpa1_proteins$protein_id
prot_lengths <- nchar(sftpa1_proteins$protein_sequence)

prot_ranges <- IRanges(
  start = rep(1L, length(prot_ids)),
  end   = prot_lengths,
  names = prot_ids
)

p2g_result <- proteinToGenome(prot_ranges, edb)



p2g_df_list <- lapply(names(p2g_result), function(pid) {
  gr <- p2g_result[[pid]]
  if (is.null(gr) || length(gr) == 0) {
    return(NULL)
  }
  df                  <- as.data.frame(gr)
  df$query_protein_id <- pid
  df
})

p2g_df_list <- Filter(Negate(is.null), p2g_df_list)

if (length(p2g_df_list) > 0) {
  p2g_combined        <- do.call(rbind, p2g_df_list)
  rownames(p2g_combined) <- NULL

  write.csv(
    p2g_combined,
    file      = "protein_to_genome.csv",
    row.names = FALSE
  )
  cat("total exon rows exported:", nrow(p2g_combined), "\n")
} else {
  cat("no mappings found")
}
