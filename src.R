library(tidyverse)

# Wrapped callback function to collect unique values of a column of a chunk
collect_id <- function(index, output_file, max_size = 1000000) {
    # index is the index column of the table that we want to get
    # chunk & pos are fixed arguments of the callback function
    function(chunk, pos) {
        # ids is a global vector hence manipulated using <<-
        ids <- c(ids, chunk[[index]])
        ids <<- unique(ids)
        # if ids grow too large, write the partial result and reset the vector
        # there is a small chance a repeated value will slip through if the
        # value was already saved in the file but is unique in the ids vector
        # but without this vector memory reset, RAM will be an issue again
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
                 comment = "#", chunk_size = 10000, show_col_types = F)
# write the last chunk to output
write.table(as_tibble(ids), output_file, sep = "\t",
            col.names = F, row.names = F, quote = F, append = T)
# free up memory
rm(ids)

# --- Collect IDs from AlphaMissense_aa_substitutions.tsv.gz
ids <- c()
output_file <- "AM_aa_uniprot_id.tsv"
read_tsv_chunked("AlphaMissense_aa_substitutions.tsv.gz",
                 DataFrameCallback$new(collect_id(1, output_file)),
                 comment = "#", chunk_size = 10000, show_col_types = F)
write.table(as_tibble(ids), output_file, sep = "\t",
            col.names = F, row.names = F, quote = F, append = T)
rm(ids)

# --- Collect IDs from AlphaMissense_isoforms_hg38.tsv.gz
ids <- c()
output_file <- "AM_isoform_ens_id.tsv"
read_tsv_chunked("AlphaMissense_isoforms_hg38.tsv.gz",
                 DataFrameCallback$new(collect_id(6, output_file)),
                 comment = "#", chunk_size = 10000, show_col_types = F)
write.table(as_tibble(ids), output_file, sep = "\t",
            col.names = F, row.names = F, quote = F, append = T)
rm(ids)

# --- Collect IDs from AlphaMissense_isoforms_aa_substitutions.tsv.gz
ids <- c()
output_file <- "AM_isoform_aa_ens_id.tsv"
read_tsv_chunked("AlphaMissense_isoforms_aa_substitutions.tsv.gz",
                 DataFrameCallback$new(collect_id(1, output_file)),
                 comment = "#", chunk_size = 10000, show_col_types = F)
write.table(as_tibble(ids), output_file, sep = "\t",
            col.names = F, row.names = F, quote = F, append = T)
rm(ids)


# --- Harmonizing the IDs
ensembl_canonical <- read_tsv("canonical_transcripts.tsv", col_names = T)
am_uniprot <- unique(read_tsv("AM_uniprot_id.tsv", col_names = F))
am_aa_uniprot <- unique(read_tsv("AM_aa_uniprot_id.tsv", col_names = F))
am_iso_ens <- unique(read_tsv("AM_isoform_ens_id.tsv", col_names = F))
### am_iso_aa is a superset of am_iso_ens, we don't need to use am_iso_ens.
am_iso_aa_ens <- unique(read_tsv("AM_isoform_aa_ens_id.tsv", col_names = F))
idmapping1 <- read_tsv("idmapping_1.tsv.gz")
idmapping2 <- read_tsv("idmapping_2.tsv.gz")
### idmapping1 obtained from Uniprot tool (AM_uniprot_id.tsv), 49328 hits.
### idmapping2 obtained from Uniprot tool, (AM_aa_uniprot_id.tsv), 50550 hits.

idmapping1 <- dplyr::rename(idmapping1, "uniprot_ID" = "From", "ID" = "To")
idmapping2 <- dplyr::rename(idmapping2, "uniprot_ID" = "From", "ID" = "To")
idmapping1_canonical <- ensembl_canonical %>% inner_join(idmapping1, by = "ID")
idmapping2_canonical <- ensembl_canonical %>% inner_join(idmapping2, by = "ID")
am_mapped_to_ensembl_canonical <- bind_rows(idmapping1_canonical, idmapping2_canonical) %>% distinct()
### Combine all IDs from AM_uniprot and AM_aa_uniprot, 18407 IDs.

am_iso_aa_ens <- dplyr::rename(am_iso_aa_ens, "ID" = "X1")
missing_from_am <- ensembl_canonical %>% filter(!ID %in% am_mapped_to_ensembl_canonical$ID)
print(missing_from_am)
### 1288 IDs in ensembl canonical but not in the main AM.
print(missing_from_am %>% filter(ID %in% am_iso_aa_ens$ID))
### 752 out of 1288 found in AM_iso_aa.
missing_from_am <- missing_from_am %>% mutate(ID2 = sub("\\.\\d+$", "", ID))
am_iso_aa_ens <- am_iso_aa_ens %>% mutate(ID2 = sub("\\.\\d+$", "", ID))
print(missing_from_am %>% filter(ID2 %in% am_iso_aa_ens$ID2))
### 1004 out of 1288 found in AM_iso_aa if using stable IDs without version.
ids_to_add <- missing_from_am %>% filter(ID2 %in% am_iso_aa_ens$ID2)
ids_to_add <- ids_to_add %>% mutate(uniprot_ID = "") %>% select("ID", "uniprot_ID")

am_mapped_to_ensembl_canonical <- rows_insert(am_mapped_to_ensembl_canonical, ids_to_add)
am_mapped_to_ensembl_canonical <- am_mapped_to_ensembl_canonical %>% mutate(ID2 = sub("\\.\\d+$", "", ID))
### Combine all IDs, 19411 IDs in total.

# --- Filtering the IDs
# Ensembl IDs that match to >1 Uniprot ID
multiple_match <- am_mapped_to_ensembl_canonical %>% 
    group_by(ID) %>%
    filter(n() > 1) %>%
    summarise(uniprot_IDs = list(unique(uniprot_ID))) %>%
    unnest(cols = uniprot_IDs)
print(multiple_match)
### 6 IDs flagged.
obsolete_ids <- c("A0A140G945", "Q32Q52", "Q8IXS6", "Q5VZT2", "Q9UPP5", "Q6ZW33")
### Obsolete ids identified from looking up Uniprot entries.
am_mapped_to_ensembl_canonical <- am_mapped_to_ensembl_canonical %>%
    filter(!uniprot_ID %in% obsolete_ids)

# Uniprot IDs that did not map using Uniprot tool
idmapping <- bind_rows(idmapping1, idmapping2) %>% unique()
unmatched_uniprot <- am_uniprot %>%
    anti_join(idmapping, by = c("X1" = "uniprot_ID"))
print(unmatched_uniprot)
### 215 unmatched. Not sure how OP identified 25 perfectly fine proteins.
### skipping.
am_mapped_to_ensembl_canonical <- am_mapped_to_ensembl_canonical %>%
    filter(!uniprot_ID %in% unmatched_uniprot)

# Check for valid transcripts via CDS
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
print(length(invalid_am_ids))
### 97 transcripts are invalid.
invalid_am_ids <- c(invalid_am_ids, "ENST00000569103.2")
### 1 transcript is a non-functional IG according to OP. Added to invalids.

am_mapped_to_ensembl_canonical <- am_mapped_to_ensembl_canonical %>% filter(!(ID %in% invalid_am_ids))
write.table(am_mapped_to_ensembl_canonical, "am_mapped_to_ensembl_canonical.txt", sep = "\t", quote = F, row.names = F)


# --- Writing output (intermediate) to prevent redundant processing later during length check

# AlphaMissense_hg38.tsv.gz + AlphaMissense_isoforms_hg38.tsv.gz
# -> filtered through `am_mapped_to_ensembl_canonical`
# -> am.tsv
cat("ID\tchange\tscore\n", file = "am.tsv")

write_chunk <- function(chunk, pos) {
    chunk <- unique(chunk)
    colnames(chunk) <- c("CHROM", "POS", "REF", "ALT", "genome",
                         "uniprot_ID", "old_ID", "change", "score", "am_class")
    chunk <- chunk %>% filter(uniprot_ID %in% am_mapped_to_ensembl_canonical$uniprot_ID)
    if (nrow(chunk) > 0) {
        chunk <- left_join(chunk, am_mapped_to_ensembl_canonical, by = "uniprot_ID")
        chunk <- select(chunk, "ID", "change", "score")
        chunk <- unique(chunk)
        write.table(chunk, "am.tsv", sep = "\t",
                    col.names = F, row.names = F, quote = F, append = T)
    }
}
read_tsv_chunked("AlphaMissense_hg38.tsv.gz",
                 DataFrameCallback$new(write_chunk), col_names = F,
                 comment = "#", chunk_size = 10000, show_col_types = F)

write_chunk <- function(chunk, pos) {
    chunk <- unique(chunk)
    colnames(chunk) <- c("CHROM", "POS", "REF", "ALT", "genome",
                         "old_ID", "change", "score", "am_class")
    chunk <- chunk %>% mutate(ID2 = sub("\\.\\d+$", "", old_ID))
    # match via IDs without version
    chunk <- chunk %>% filter(ID2 %in% am_mapped_to_ensembl_canonical$ID2)
    if (nrow(chunk) > 0) {
        chunk <- left_join(chunk, am_mapped_to_ensembl_canonical, by = "ID2")
        chunk <- select(chunk, "ID", "change", "score")
        chunk <- unique(chunk)
        write.table(chunk, "am.tsv", sep = "\t",
                    col.names = F, row.names = F, quote = F, append = T)
    }
}
read_tsv_chunked("AlphaMissense_isoforms_hg38.tsv.gz",
                 DataFrameCallback$new(write_chunk), col_names = F,
                 comment = "#", chunk_size = 10000, show_col_types = F)


# AlphaMissense_aa_substitutions.tsv.gz + AlphaMissense_isoforms_aa_substitutions.tsv.gz
# -> filtered through `am_mapped_to_ensembl_canonical`
# -> am_full.tsv
cat("ID\tchange\tscore\n", file = "am_full.tsv")

write_chunk <- function(chunk, pos) {
    chunk <- unique(chunk)
    colnames(chunk) <- c("uniprot_ID", "change", "score", "am_class")
    chunk <- chunk %>% filter(uniprot_ID %in% am_mapped_to_ensembl_canonical$uniprot_ID)
    if (nrow(chunk) > 0) {
        chunk <- left_join(chunk, am_mapped_to_ensembl_canonical, by = "uniprot_ID")
        chunk <- select(chunk, "ID", "change", "score")
        chunk <- unique(chunk)
        write.table(chunk, "am_full.tsv", sep = "\t",
                    col.names = F, row.names = F, quote = F, append = T)
    }
}
read_tsv_chunked("AlphaMissense_aa_substitutions.tsv.gz",
                 DataFrameCallback$new(write_chunk), col_names = F,
                 comment = "#", chunk_size = 10000, show_col_types = F)

write_chunk <- function(chunk, pos) {
    chunk <- unique(chunk)
    colnames(chunk) <- c("old_ID", "change", "score", "am_class")
    chunk <- chunk %>% mutate(ID2 = sub("\\.\\d+$", "", old_ID))
    chunk <- chunk %>% filter(ID2 %in% am_mapped_to_ensembl_canonical$ID2)
    if (nrow(chunk) > 0) {
        chunk <- left_join(chunk, am_mapped_to_ensembl_canonical, by = "ID2")
        chunk <- select(chunk, "ID", "change", "score")
        chunk <- unique(chunk)
        write.table(chunk, "am_full.tsv", sep = "\t",
                    col.names = F, row.names = F, quote = F, append = T)
    }
}
read_tsv_chunked("AlphaMissense_isoforms_aa_substitutions.tsv.gz",
                 DataFrameCallback$new(write_chunk), col_names = F,
                 comment = "#", chunk_size = 10000, show_col_types = F)


# --- Check the lengths
library(tidyverse)
am_map <- read_tsv("am_mapped_to_ensembl_canonical.txt")
am_map <- am_map %>% mutate("amino_acids" = -1, "bases" = -1)
### store the lengths

am_length_aa <- am_map %>% select("ID", "amino_acids")
collect_aa <- function(chunk, pos) {
    chunk <- unique(chunk)
    colnames(chunk) <- c("ID", "change", "score")
    chunk <- chunk %>% filter(ID %in% am_map$ID)
    if (nrow(chunk) > 0) {
        chunk$aa_no <- chunk$change %>% str_extract("\\d+")
        chunk <- chunk %>% group_by(ID) %>% summarize(amino_acids = max(as.numeric(aa_no)))
        # update am_length_aa only if the value in chunk is larger for the particular ID
        # because we are reading in chunks, we don't immediately have the largest length
        am_length_aa <<- am_length_aa %>%
                        left_join(chunk, by = "ID") %>%
                        mutate(amino_acids = pmax(amino_acids.y, amino_acids.x, na.rm = TRUE)) %>%
                        select(ID, amino_acids)
    }
    return(invisible(NULL))
}
read_tsv_chunked("am.tsv.gz",
                 DataFrameCallback$new(collect_aa), col_names = F,
                 comment = "#", chunk_size = 10000, show_col_types = F)
read_tsv_chunked("am_full.tsv.gz",
                 DataFrameCallback$new(collect_aa), col_names = F,
                 comment = "#", chunk_size = 10000, show_col_types = F)

am_map$amino_acids <- am_length_aa$amino_acids

library(Biostrings)
fasta <- readDNAStringSet("Homo_sapiens.GRCh38.cds.all.fa.gz")
ref_df <- tibble(ID = names(fasta), sequence = as.character(fasta))
ref_df <- ref_df %>% mutate(ID = str_extract(ID, "ENST\\d+\\.\\d+"))

am_length_bases <- am_map %>% select("ID", "bases")
collect_bases <- function(chunk, pos) {
    chunk <- unique(chunk)
    colnames(chunk) <- c("ID", "change", "score")
    chunk <- chunk %>% filter(ID %in% am_map$ID)
    if (nrow(chunk) > 0) {
        chunk <- chunk %>% filter(ID %in% am_map$ID)
        chunk <- unique(chunk %>% select(ID))
        cds <- c()
        for (id in chunk$ID) {
            seq <- ref_df %>% filter(ID == id) %>% pull(sequence)
            seq <- sub("(TGA|TAG|TAA)$", "", seq)
            cds <- c(cds, nchar(seq))
        }
        chunk$bases <- cds
        am_length_bases <<- am_length_bases %>%
                            left_join(chunk, by = "ID") %>%
                            mutate(bases = pmax(bases.y, bases.x, na.rm = TRUE)) %>%
                            select(ID, bases)
        
    }
    return(invisible(NULL))
}
read_tsv_chunked("am.tsv.gz",
                 DataFrameCallback$new(collect_bases), col_names = F,
                 comment = "#", chunk_size = 10000, show_col_types = F)
read_tsv_chunked("am_full.tsv.gz",
                 DataFrameCallback$new(collect_bases), col_names = F,
                 comment = "#", chunk_size = 10000, show_col_types = F)


am_map$bases <- am_length_bases$bases
am_map <- am_map %>% mutate(codons = bases / 3)

mismatch_length_am <- am_map %>% filter(amino_acids != codons)
### 774 out of 19308 have mismatched lengths
matching_length_am <- am_map %>% filter(amino_acids == codons)


# --- Write final output where lengths are matching
write_chunk_filtered <- function(output_file) {
    function(chunk, pos) {
        colnames(chunk) <- c("ID", "change", "score")
        chunk <- chunk %>% filter(ID %in% matching_length_am$ID)
        if (nrow(chunk) > 0) {
            write.table(chunk, output_file, sep = "\t",
                        col.names = F, row.names = F, quote = F, append = T)
        }
    }
}
cat("ID\tchange\tscore\n", file = "am_filtered.tsv")
read_tsv_chunked("am.tsv.gz",
                 DataFrameCallback$new(write_chunk_filtered("am_filtered.tsv")),
                 col_names = F, chunk_size = 10000, show_col_types = F)

cat("ID\tchange\tscore\n", file = "am_full_filtered.tsv")
read_tsv_chunked("am_full.tsv.gz",
                 DataFrameCallback$new(write_chunk_filtered("am_full_filtered.tsv")),
                 col_names = F, chunk_size = 10000, show_col_types = F)
