library(tidyverse)

# Wrapped callback function to collect unique values of a column of a chunk
collect_id <- function(index, output_file, max_size = 1000000) {
    # index is the index column of the table that we want to get
    # chunk & pos are fixed arguments of the callback function
    function(chunk, pos) {
        ids <- c(ids, chunk[[index]])
        ids <<- unique(ids)
        # if ids grow too large, write the partial result and reset the vector
        if (object.size(ids) > max_size) {
            write.table(as_tibble(ids), output_file, sep = "\t",
                        col.names = F, row.names = F, quote = F, append = T)
            ids <<- c()
        }
    }
}

# --- Collect IDs from AlphaMissense_hg38.tsv.gz
ids <- c()
output_file <- "AM_uniprot_id.tsv"
read_tsv_chunked("AlphaMissense_hg38.tsv.gz",
                 DataFrameCallback$new(collect_id(6, output_file)),
                 comment = "#", chunk_size = 1000, show_col_types = F)
# write the last chunk to output
write.table(as_tibble(ids), output_file, sep = "\t",
            col.names = F, row.names = F, quote = F, append = T)
# free up memory
rm(ids)

# --- Collect IDs from AlphaMissense_isoforms_hg38.tsv.gz
ids <- c()
output_file <- "AM_isoform_ens_id.tsv"
read_tsv_chunked("AlphaMissense_isoforms_hg38.tsv.gz",
                 DataFrameCallback$new(collect_id(6, output_file)),
                 comment = "#", chunk_size = 1000, show_col_types = F)
write.table(as_tibble(ids), output_file, sep = "\t",
            col.names = F, row.names = F, quote = F, append = T)
rm(ids)

# --- Collect IDs from AlphaMissense_isoforms_aa_substitutions.tsv.gz
ids <- c()
output_file <- "AM_isoform_aa_ens_id.tsv"
read_tsv_chunked("AlphaMissense_isoforms_aa_substitutions.tsv.gz",
                 DataFrameCallback$new(collect_id(1, output_file)),
                 comment = "#", chunk_size = 1000, show_col_types = F)
write.table(as_tibble(ids), output_file, sep = "\t",
            col.names = F, row.names = F, quote = F, append = T)
rm(ids)


# --- Querying the IDs
ensembl_canonical <- read_tsv("canonical_transcripts.tsv", col_names = T)
am_uniprot <- unique(read_tsv("AM_uniprot_id.tsv", col_names = F))
am_iso_ens <- unique(read_tsv("AM_isoform_ens_id.tsv", col_names = F))
### am_iso_aa is a superset of am_iso_ens, we don't need to use am_iso_ens.
am_iso_aa_ens <- unique(read_tsv("AM_isoform_aa_ens_id.tsv", col_names = F))
idmapping <- read_tsv("idmapping_1.tsv.gz")
### idmapping obtained from Uniprot tool, 49328 hits.

idmapping <- dplyr::rename(idmapping, "uniprot_ID" = "From", "ID" = "To")
am_mapped_to_ensembl_canonical <- ensembl_canonical %>% inner_join(idmapping, by = "ID")
am_iso_aa_ens <- dplyr::rename(am_iso_aa_ens, "ID" = "X1")
missing_from_am <- ensembl_canonical %>% filter(!ID %in% idmapping$ID)
### 1569 IDs in ensembl canonical but not in the main AM.
print(missing_from_am %>% filter(ID %in% am_iso_aa_ens$ID))
### 862 out of 1569 found in am_iso_aa
missing_from_am <- missing_from_am %>% mutate(ID2 = sub("\\.\\d+$", "", ID))
am_iso_aa_ens <- am_iso_aa_ens %>% mutate(ID2 = sub("\\.\\d+$", "", ID))
print(missing_from_am %>% filter(ID2 %in% am_iso_aa_ens$ID2))
### 1181 out of 1569 found in am_iso_aa if using stable IDs without version.
ids_to_add <- missing_from_am %>% filter(ID2 %in% am_iso_aa_ens$ID2)
ids_to_add <- ids_to_add %>% mutate(uniprot_ID = "") %>% select("ID", "uniprot_ID")

am_mapped_to_ensembl_canonical <- rows_insert(am_mapped_to_ensembl_canonical, ids_to_add)
### Combine all IDs

# Ensembl IDs that match to >1 Uniprot ID
multiple_match <- am_mapped_to_ensembl_canonical %>% 
    group_by(ID) %>%
    filter(n() > 1) %>%
    summarise(uniprot_IDs = list(unique(uniprot_ID))) %>%
    unnest(cols = uniprot_IDs)
print(multiple_match)
### 5 IDs flagged. How did OP decide which uniprot_ID to use for each of these?
obsolete_ids <- c("Q32Q52", "Q8IXS6", "Q5VZT2", "Q9UPP5", "Q6ZW33")
am_mapped_to_ensembl_canonical <- am_mapped_to_ensembl_canonical %>%
    filter(!uniprot_ID %in% obsolete_ids)

# Uniprot IDs that did not map using Uniprot tool
unmatched_uniprot <- am_uniprot %>%
    anti_join(idmapping, by = c("X1" = "uniprot_ID"))
print(unmatched_uniprot)
### 215 unmatched
### Not sure how OP identified the 25 perfectly fine proteins, skipping.
am_mapped_to_ensembl_canonical <- am_mapped_to_ensembl_canonical %>%
    filter(!uniprot_ID %in% unmatched_uniprot)

# Check for valid transcripts
library(Biostrings)
fasta <- readDNAStringSet("Homo_sapiens.GRCh38.cds.all.fa.gz")
ref_df <- tibble(ID = names(fasta), sequence = as.character(fasta))
ref_df <- ref_df %>% mutate(ID = str_extract(ID, "ENST\\d+\\.\\d+"))

is_invalid_seq <- function(seq) {
    return(any(
        # not divisible by 3
        nchar(seq) %% 3 != 0 |
            # contains N
            grepl("N", seq) |
            # does not start with start codon
            !(startsWith(seq, "ATG") | startsWith(seq, "CTG") | startsWith(seq, "TTG")) |
            # does not end with stop codon
            !(endsWith(seq, "TGA") | endsWith(seq, "TAG") | endsWith(seq, "TAA"))
    ))
}

invalid_am_ids <- c()
for (id in am_mapped_to_ensembl_canonical$ID) {
    seq <- ref_df %>% filter(ID == id) %>% pull(sequence)
    if (is_invalid_seq(seq)) {
        invalid_am_ids <- c(invalid_am_ids, id)
    }
}
print(invalid_am_ids)
### 85 transcripts are invalid.
invalid_am_ids <- c(invalid_am_ids, "ENST00000569103.2")
### 1 transcript is a non-functional IG according to OP. Added to invalids.

am_mapped_to_ensembl_canonical <- am_mapped_to_ensembl_canonical %>% filter(!(ID %in% invalid_am_ids))
### Final list of valid IDs.


# --- Writing output

# AlphaMissense_hg38.tsv.gz -> am_mapped_to_ensembl_canonical -> am_output.tsv
cat("ID\tchange\tscore\n", file = "am_output.tsv")

write_chunk <- function(chunk, pos) {
    chunk <- unique(chunk)
    colnames(chunk) <- c("CHROM", "POS", "REF", "ALT", "genome",
                         "uniprot_ID", "old_ID", "change", "score", "am_class")
    # Filter for uniprot IDs that are in the mapping
    chunk <- chunk %>% filter(uniprot_ID %in% am_mapped_to_ensembl_canonical$uniprot_ID)
    if (nrow(chunk) > 0) {
        # Map to the updated Ensembl IDs
        chunk <- left_join(chunk, am_mapped_to_ensembl_canonical, by = c("uniprot_ID" = "uniprot_ID"))
        subchunk <- select(chunk, "ID", "change", "score")
        subchunk <- unique(subchunk)
        write.table(subchunk, "am_output.tsv",
                    col.names = F, row.names = F, quote = F, append = T)
    }
}
read_tsv_chunked("AlphaMissense_hg38.tsv.gz",
                 DataFrameCallback$new(write_chunk),
                 comment = "#", chunk_size = 1000, show_col_types = F)


# AlphaMissense_isoforms_aa_substitutions.tsv.gz -> am_mapped_to_ensembl_canonical -> am_iso_aa_output.tsv
cat("ID\tchange\tscore\n", file = "am_iso_aa_output.tsv")

write_chunk <- function(chunk, pos) {
    chunk <- unique(chunk)
    colnames(chunk) <- c("ID", "change", "score", "am_class")
    chunk <- chunk %>% filter(ID %in% am_mapped_to_ensembl_canonical$ID)
    if (nrow(chunk) > 0) {
        subchunk <- select(chunk, "ID", "change", "score")
        subchunk <- unique(subchunk)
        write.table(subchunk, "am_iso_aa_output.tsv",
                    col.names = F, row.names = F, quote = F, append = T)
    }
}
read_tsv_chunked("AlphaMissense_isoforms_aa_substitutions.tsv.gz",
                 DataFrameCallback$new(write_chunk),
                 comment = "#", chunk_size = 1000, show_col_types = F)
