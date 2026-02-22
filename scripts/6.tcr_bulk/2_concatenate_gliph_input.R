# 1 open required packages, define functions ----
library(dplyr)
library(stringr)

outdir = "/path/to/gliph_input"

underscore_to_dash = function(text) {
  text = stringr::str_replace_all(text, "_", "-")
  return(text)
}

transform_trb = function(text) {
  text = gsub("\\([^()]*\\)", "", text) # remove parentheses
  text = gsub("\\*00", "", text) # remove *00
  text = gsub(",", "-or", text) # replace comma with -or
  text = gsub("V0", "V", text) # remove 0 from V
  text = gsub("J0", "J", text) # remove 0 from J
  text = gsub("-0", "-", text) # remove 0 from -
  text = gsub("_0", "-", text)
  text = gsub("or0", "or", text) # remove 0 from or
  text = gsub("\\*0", "-", text) # remove 0 from *
  text = gsub("-(\\d)-\\d$", "-\\1", text) # remove last digit after -
  text = gsub("TRBV(\\d+)", "TRBV\\1", text)
  text = gsub("TRBJ(\\d+)", "TRBJ\\1", text)
  text = gsub("/0", "-or", text)
  text = gsub("/", "-or", text)
  text = gsub("orTRBV", "or", text)
  text = gsub("orTRBJ", "or", text)
  text = gsub("-or12-4-or12-4", "-or12-4", text)
  text = gsub("TRBV20-1-or20-or9-2", "TRBV20-1-or9-2", text)
  text = gsub("TRBV25-1-or25-or9-2", "TRBV25-1-or9-2", text)
  text = gsub("TRBV25-1-or25-or9-2-or6-4", "TRBV25-1-or9-2-or6-4", text)
  text = gsub("6-2-or6-3-or6-2", "6-2-or6-3", text)
  text = gsub("(-or\\d+)(-or\\1)+", "\\1", text)
  return(text)
}

normalize_tr_genes = function(text) {
  text = stringr::str_replace_all(text, "TCR", "TR")
  text = transform_trb(text)
  return(text)
}

standardize_gliph_df = function(
  df,
  cdr3_col,
  trbv_col,
  trbv_ties_col,
  trbj_col,
  trbj_ties_col,
  count_col,
  subject_condition,
  cdr3a_value = NA
) {
  out = df
  out$TRBV = ifelse(out[[trbv_col]] != "", out[[trbv_col]], out[[trbv_ties_col]])
  out$TRBJ = ifelse(out[[trbj_col]] != "", out[[trbj_col]], out[[trbj_ties_col]])
  out$TRBV = normalize_tr_genes(out$TRBV)
  out$TRBJ = normalize_tr_genes(out$TRBJ)
  out$CDR3a = cdr3a_value
  out$subject.condition = subject_condition

  out = out %>%
    dplyr::transmute(
      CDR3b = .data[[cdr3_col]],
      TRBV = TRBV,
      TRBJ = TRBJ,
      CDR3a = CDR3a,
      subject.condition = subject.condition,
      count = .data[[count_col]]
    ) %>%
    dplyr::distinct()

  return(out)
}

# 2 open input files ----

## tigerdb known antigen (4146 CDR3s) ----
tiger = read.delim("/path/to/tigerdb/VDJdb_MinimalScoreConfidence3_VersionII.tsv")

## ped cancer ----
### intercept cancer child (50572 CDR3s) ----
intercept = read.table("/pat/to/intercept_tcrs.txt") # no normal, no lfs, no singles

## adult cancer ----
### mdavis cancer adult NSCLC (19596 CDR3s) ----
mdavis = read.delim("/path/to/mdavis/MarkDavis_TumorEnrichedTRB-CDR3_v4.tsv")

### hamm adult cancer (490421) ----
hamm_adult_cancer = read.csv("/path/to/reference_tcrs/hamm_data_cancer.csv") %>%
  dplyr::select(aminoAcid, `count..templates.reads.`, vGeneName, vGeneNameTies, jGeneName, jGeneNameTies, filename)

## ped healthy ----
### mitchell ped healthy (727915) ----
mitch = read.csv("/path/to/mitchell_data.csv") %>%
  dplyr::select(aminoAcid, `count..templates.reads.`, vGeneName, vGeneNameTies, jGeneName, jGeneNameTies, filename)

### emerson ped healthy (87500) ----
emerson_child = read.csv("/path/to/emerson_child_tcrs.csv") %>%
  dplyr::select(amino_acid, seq_reads, v_resolved, v_gene_ties, j_resolved, j_gene_ties, sample_name)

### wherry ped healthy ----
wherry_ped = read.csv("/path/to/wherry_ped_tcrs/Rearrangements_all.csv") %>%
  dplyr::select(sample_name, amino_acid, v_trim, j_trim, templates)

## adult healthy ----
### emerson adult healthy (483003) ----
emerson_adult = read.csv("/path/to/emerson_adult_tcrs.csv") %>%
  dplyr::select(amino_acid, seq_reads, v_resolved, v_gene_ties, j_resolved, j_gene_ties, sample_name)

### hamm adult healthy (1100568) ----
hamm_adult_healthy = read.csv("/path/to/hamm_data_ctr.csv") %>%
  dplyr::select(aminoAcid, `count..templates.reads.`, vGeneName, vGeneNameTies, jGeneName, jGeneNameTies, filename)

## anbl + profyle ----
anbl_PBMC_tumor = list.files(
  "/path/to/mixcr4/output/",
  pattern = "clones_TRB.tsv",
  full.names = TRUE
) %>%
  lapply(function(file) {
    data = read.table(file, sep = "\t", header = TRUE)
    data$filename = basename(file) # Adds a new column with the filename
    return(data)
  }) %>%
  bind_rows() %>%
  dplyr::filter(!grepl("_", aaSeqCDR3)) %>%
  dplyr::filter(!grepl("\\*", aaSeqCDR3))

anbl_PBMC_tumor$patient = case_when(
  grepl("-SM-", anbl_PBMC_tumor$filename) ~ substr(anbl_PBMC_tumor$filename, 9, 14),
  grepl("PRO", anbl_PBMC_tumor$filename) ~ substr(anbl_PBMC_tumor$filename, 8, 13),
  TRUE ~ substr(anbl_PBMC_tumor$filename, 10, 15)
)

anbl_PBMC_tumor$patient = sub("-$", "", anbl_PBMC_tumor$patient)

anbl_PBMC_tumor$timepoint = case_when(
  grepl("-SM-", anbl_PBMC_tumor$filename) ~ substr(anbl_PBMC_tumor$filename, 16, 17),
  grepl("PRO", anbl_PBMC_tumor$filename) ~ substr(anbl_PBMC_tumor$filename, 14, 16),
  TRUE ~ substr(anbl_PBMC_tumor$filename, 17, 18)
)

anbl_PBMC_tumor$timepoint = case_when(
  anbl_PBMC_tumor$timepoint == "T0" ~ "pTx",
  anbl_PBMC_tumor$timepoint == "T1" ~ "pI4",
  anbl_PBMC_tumor$timepoint == "T2" ~ "ePC",
  anbl_PBMC_tumor$timepoint == "T3" ~ "pI3",
  anbl_PBMC_tumor$timepoint == "T4" ~ "pSu",
  anbl_PBMC_tumor$timepoint == "T5" ~ "Pgn",
  anbl_PBMC_tumor$timepoint == "T6" ~ "Rel",
  anbl_PBMC_tumor$timepoint == "T7" ~ "EOI",
  anbl_PBMC_tumor$timepoint == "T8" ~ "pI2",
  anbl_PBMC_tumor$timepoint == "T9" ~ "M72",
  anbl_PBMC_tumor$timepoint == "DX" ~ "tDX",
  anbl_PBMC_tumor$timepoint == "PT" ~ "tPT",
  anbl_PBMC_tumor$timepoint == "T1-" ~ "PT1",
  anbl_PBMC_tumor$timepoint == "T2-" ~ "PT2",
  TRUE ~ anbl_PBMC_tumor$timepoint
)

# 3 reformat to GLIPH2 input format ----
## tigerdb ----
tiger = tiger %>% dplyr::distinct()
tiger$number = 1:nrow(tiger)
tiger$tigernumber = paste0("tiger", tiger$number)

tigerdb = data.frame(
  CDR3b = tiger$CDR3b,
  TRBV = tiger$TRBV,
  TRBJ = tiger$TRBJ,
  CDR3a = NA,
  HLA = tiger$MHC,
  subject.condition = paste0(
    tiger$tigernumber, ":", tiger$Epitope.species, "_", tiger$Mutation, "_", tiger$Epitope.gene.protein.name
  ),
  count = 1
)
tigerdb$subject.condition = gsub("_NA", "", tigerdb$subject.condition)

## anbl ----
anbl = data.frame(
  CDR3b = anbl_PBMC_tumor$aaSeqCDR3,
  TRBV = anbl_PBMC_tumor$allVHitsWithScore,
  TRBJ = anbl_PBMC_tumor$allJHitsWithScore,
  CDR3a = NA,
  subject.condition = case_when(
    grepl("CHP-P", anbl_PBMC_tumor$filename) ~
      paste0("anbl", anbl_PBMC_tumor$patient, ":", "NB_", anbl_PBMC_tumor$timepoint, "-Tumor"),
    grepl("PRO", anbl_PBMC_tumor$filename) ~
      paste0("pro", anbl_PBMC_tumor$patient, ":", "NB_", anbl_PBMC_tumor$timepoint, "-PBMC"),
    grepl("SM-P", anbl_PBMC_tumor$filename) ~
      paste0("anbl", anbl_PBMC_tumor$patient, ":", "NB_", anbl_PBMC_tumor$timepoint, "-cfDNA"),
    TRUE ~
      paste0("anbl", anbl_PBMC_tumor$patient, ":", "NB_", anbl_PBMC_tumor$timepoint, "-PBMC")
  ),
  count = anbl_PBMC_tumor$readCount
)

## mitchell ----
mitchell_subject = paste0("mitchel", substr(mitch$filename, 1, 6))
mitchell_subject_condition = paste0(mitchell_subject, ":MitchelChild")
mitchell = standardize_gliph_df(
  df = mitch,
  cdr3_col = "aminoAcid",
  trbv_col = "vGeneName",
  trbv_ties_col = "vGeneNameTies",
  trbj_col = "jGeneName",
  trbj_ties_col = "jGeneNameTies",
  count_col = "count..templates.reads.",
  subject_condition = mitchell_subject_condition
)

## emerson ----
emerson_c_subject_condition = paste0("emerson", emerson_child$sample_name, ":EmersonChild")
emerson_c = standardize_gliph_df(
  df = emerson_child,
  cdr3_col = "amino_acid",
  trbv_col = "v_resolved",
  trbv_ties_col = "v_gene_ties",
  trbj_col = "j_resolved",
  trbj_ties_col = "j_gene_ties",
  count_col = "seq_reads",
  subject_condition = emerson_c_subject_condition
)

emerson_a_subject_condition = paste0("emerson", emerson_adult$sample_name, ":EmersonAdult")
emerson_a = standardize_gliph_df(
  df = emerson_adult,
  cdr3_col = "amino_acid",
  trbv_col = "v_resolved",
  trbv_ties_col = "v_gene_ties",
  trbj_col = "j_resolved",
  trbj_ties_col = "j_gene_ties",
  count_col = "seq_reads",
  subject_condition = emerson_a_subject_condition
)

## hamm ----
hamm_cancer_subject = paste0("hamm", str_extract(hamm_adult_cancer$filename, "(?<=_)[^_]+(?=\\.tsv)"))
hamm_cancer_subject_condition = paste0(hamm_cancer_subject, ":HammCancer")
hamm_cancer = standardize_gliph_df(
  df = hamm_adult_cancer,
  cdr3_col = "aminoAcid",
  trbv_col = "vGeneName",
  trbv_ties_col = "vGeneNameTies",
  trbj_col = "jGeneName",
  trbj_ties_col = "jGeneNameTies",
  count_col = "count..templates.reads.",
  subject_condition = hamm_cancer_subject_condition
)

hamm_healthy_subject = paste0("hamm", str_extract(hamm_adult_healthy$filename, "(?<=_)[^_]+(?=\\.tsv)"))
hamm_healthy_subject_condition = paste0(hamm_healthy_subject, ":HammHealthy")
hamm_healthy = standardize_gliph_df(
  df = hamm_adult_healthy,
  cdr3_col = "aminoAcid",
  trbv_col = "vGeneName",
  trbv_ties_col = "vGeneNameTies",
  trbj_col = "jGeneName",
  trbj_ties_col = "jGeneNameTies",
  count_col = "count..templates.reads.",
  subject_condition = hamm_healthy_subject_condition
)

## mdavis ----
mdv = mdavis %>%
  dplyr::mutate(
    TRBV = sapply(TRBV, underscore_to_dash),
    TRBJ = sapply(TRBJ, underscore_to_dash)
  ) %>%
  dplyr::distinct()
mdv$number = 1:nrow(mdv)
mdv$subject.condition = paste0("mdavis", mdv$number, ":MDavis")
mdv = mdv %>% dplyr::select(CDR3b, TRBV, TRBJ, CDR3a, subject.condition, count)

## intercept ----
colnames(intercept) = c("CDR3b", "TRBV", "TRBJ", "CDR3a", "subject.condition", "count")
intercept$subject.condition = gsub("_TCR_", "_", intercept$subject.condition)
intercept$patient = gsub(".*_(.*?)_.*", "\\1", intercept$subject.condition)
intercept$patient = case_when(
  grepl("CHP_", intercept$subject.condition) ~ paste0("chp", intercept$patient),
  TRUE ~ intercept$patient
)
intercept$timepoint = gsub(".*_X(.*?)\\-.*", "\\1", intercept$subject.condition)
intercept$sampletype = gsub(".*-", "", intercept$subject.condition)
intercept$cancertype = gsub("^(.*?):.*", "\\1", intercept$subject.condition)
intercept$subject.condition = paste0(intercept$patient, ":", intercept$cancertype, "_X", intercept$timepoint, "-", intercept$sampletype)
intercept = intercept %>%
  dplyr::select(CDR3b, TRBV, TRBJ, CDR3a, subject.condition, count) %>%
  dplyr::distinct()

## wherry ----
wherry = wherry_ped
colnames(wherry) = c("subject.condition", "CDR3b", "TRBV", "TRBJ", "count")
wherry$CDR3a = NA
wherry$subject.condition = paste0("wherry", substr(wherry$subject.condition, 8, 13), ":WherryChild")
wherry = wherry %>%
  dplyr::select(CDR3b, TRBV, TRBJ, CDR3a, subject.condition, count) %>%
  dplyr::distinct() %>%
  dplyr::filter(TRBV != "unknown")
wherry$TRBV = normalize_tr_genes(wherry$TRBV)
wherry$TRBJ = normalize_tr_genes(wherry$TRBJ)

# 4 get specific number of TCRs ----
set.seed(555)

hamm_healthy_final = hamm_healthy[sample(nrow(hamm_healthy), 134516), ]
emerson_a_final = emerson_a[sample(nrow(emerson_a), 134516), ]

mitchell_final = mitchell[sample(nrow(mitchell), 38212), ]
emerson_c_final = emerson_c[sample(nrow(emerson_c), 38212), ]

hamm_cancer_final = hamm_cancer[sample(nrow(hamm_cancer), 117830), ]
mdv_final = mdv

intercept_final = intercept
intercept_final$patient = substr(intercept_final$subject.condition, 1, 6)
intercept_final = intercept_final %>% dplyr::select(-patient)

# 5 combine all tcr dataframes ----
df_list = list(
  tigerdb, anbl, hamm_healthy_final, emerson_a_final, mitchell_final,
  emerson_c_final, hamm_cancer_final, mdv_final, intercept_final, wherry
)

df_list = lapply(df_list, function(df) {
  dplyr::select(df, CDR3b, TRBV, TRBJ, CDR3a, subject.condition, count)
})

gliph_input = do.call(rbind, df_list)
gliph_input = dplyr::filter(gliph_input, !is.na(TRBV))
gliph_input = unique(gliph_input)
gliph_input = gliph_input %>% dplyr::filter(TRBV != "")

gliph_input = gliph_input %>%
  dplyr::mutate(
    TRBV = sapply(TRBV, transform_trb),
    TRBJ = sapply(TRBJ, transform_trb)
  )

gliph_input$TRBJ = ifelse(gliph_input$TRBJ == "TRBJNA", NA, gliph_input$TRBJ)
gliph_input$TRBJ = ifelse(gliph_input$TRBJ == "TRBJ", NA, gliph_input$TRBJ)
gliph_input$TRBJ = ifelse(gliph_input$TRBJ == "", NA, gliph_input$TRBJ)
gliph_input$TRBV = gsub("TRBVTRBVA", "TRBVA", gliph_input$TRBV)

table(gliph_input$TRBV)
table(gliph_input$TRBJ)

## indicate which intercept patients have NB ----
intercept_nb = c("chp384", "chp385", "chp418", "chp379", "chp346")
gliph_input$subject.condition = ifelse(
  substr(gliph_input$subject.condition, 1, 6) %in% intercept_nb,
  str_replace_all(gliph_input$subject.condition, "Solidtumors", "NBSolidtumors"),
  gliph_input$subject.condition
)

# 7 save final gliph input as RDS and .txt file ----
saveRDS(gliph_input, file.path(outdir, "gliph_input.rds"))
write.table(
  gliph_input,
  file = file.path(outdir, "gliph_input.rds"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
