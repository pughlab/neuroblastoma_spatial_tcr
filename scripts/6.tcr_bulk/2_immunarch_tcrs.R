# Bulk TCR analysis from MiXCR clones output using immunarch.
#
# Computes per-sample Hill diversity (Q = 1) for PBMC and tumor compartments
# across TRA / TRB / TRG chains, tests longitudinal change with a linear mixed
# model (Value ~ timepoint + (1 | patient)), Holm-adjusts pairwise contrasts,
# and tracks four T-cell clonotypes across samples with ggalluvial.
#
# Per-figure CSVs are written to ../../data/processed/ and the corresponding
# PDFs to ../../figures/.

suppressPackageStartupMessages({
  library(immunarch)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(ggsignif)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(tibble)
  library(tidyr)
  library(ggalluvial)
})

# ---- Paths ----------------------------------------------------------------
# `input_repdir` should point to the MiXCR `clones_TR{A,B,G}.tsv` outputs
# from Step 6.1; per-sample file names are expected to already use the
# `PatientN` identifiers. When the directory is absent, the figure-generation
# half of the script reads the CSVs already in `outdir_csv`.
input_repdir <- "/path/to/xenium/mixcr4/analyze_output/"

outdir_csv   <- "../../data/processed"
outdir_plots <- "../../figures"

dir.create(outdir_csv,   showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_plots, showWarnings = FALSE, recursive = TRUE)

run_compute_stage <- dir.exists(input_repdir)

# ---- Plot theme -----------------------------------------------------------
theme_pub <- function() {
  theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "none",
      axis.text  = element_text(size = 6),
      axis.title = element_text(size = 6),
      plot.title = element_text(size = 6),
      axis.line.x = element_line(linewidth = 0.5),
      axis.line.y = element_line(linewidth = 0.5),
      axis.ticks  = element_line(linewidth = 0.5)
    )
}

# ---- Sample metadata ------------------------------------------------------
# Parses patient ID, TCR chain, treatment timepoint, project, sample type
# (tumor / PBMC / cfDNA) and treatment arm from the MiXCR sample name.
add_metadata_columns <- function(meta_df) {
  meta_df %>%
    mutate(
      patient = str_extract(Sample, "(?:ANBL|PRO)-[A-Z]+-([A-Za-z0-9]+)") %>%
        str_replace(".*-([A-Za-z0-9]+)$", "\\1"),
      chain = case_when(
        grepl("TRB",  Sample) ~ "TRB",
        grepl("TRAD", Sample) ~ "TRAD",
        grepl("TRA",  Sample) ~ "TRA",
        grepl("TRG",  Sample) ~ "TRG",
        grepl("TRD",  Sample) ~ "TRD",
        TRUE ~ NA_character_
      ),
      timepoint = case_when(
        grepl("-T0-",  Sample) ~ "pTx",
        grepl("-pTx-", Sample) ~ "pTx",
        grepl("-T1-",  Sample) ~ "pI4",
        grepl("-T2-",  Sample) ~ "ePC",
        grepl("-T3-",  Sample) ~ "pI3",
        grepl("-T4-",  Sample) ~ "pSu",
        grepl("-T5-",  Sample) ~ "Pgn",
        grepl("-T6-",  Sample) ~ "Rel",
        grepl("-T7-",  Sample) ~ "EOI",
        grepl("-T8-",  Sample) ~ "pI2",
        grepl("-T9-",  Sample) ~ "M72",
        grepl("-DX-",  Sample) ~ "tDX",
        grepl("-PT-",  Sample) ~ "tPT",
        TRUE ~ "NA"
      ),
      project = case_when(
        grepl("^ANBL", Sample) ~ "ANBL",
        grepl("^PRO",  Sample) ~ "PRO",
        TRUE ~ "unknown"
      ),
      # Compartment is encoded by the DNA-source token (-R- tumor,
      # -P- cfDNA plasma, -M- PBMC mononuclear).
      sampletype = case_when(
        project == "PRO"           ~ "PBMC",
        grepl("-R-DNA", Sample)    ~ "tumor",
        grepl("-P-DNA", Sample)    ~ "cfDNA",
        grepl("-M-DNA", Sample)    ~ "PBMC",
        TRUE                        ~ "unknown"
      ),
      arm = case_when(
        patient %in% c("Patient1", "Patient2")             ~ "A",
        patient %in% c("Patient3", "Patient4", "Patient5") ~ "B",
        TRUE                                                ~ "unknown"
      )
    )
}

# ---------------------------------------------------------------------------
# Repertoire diversity (Extended Data Fig 1a-f)
# ---------------------------------------------------------------------------

# Hill diversity (Q = 1) within one (compartment, chain) panel. Fits
# Value ~ timepoint + (1 | patient) and writes the per-sample diversity
# values and Holm-adjusted pairwise contrasts to two CSVs.
compute_diversity_csv <- function(data_list, meta_df, sampletype_keep, chain_keep,
                                  time_levels, data_csv, stats_csv) {
  filtered_meta <- meta_df %>%
    dplyr::filter(timepoint %in% time_levels,
                  sampletype == sampletype_keep,
                  chain == chain_keep) %>%
    mutate(timepoint = factor(timepoint, levels = time_levels))
  filtered_data <- data_list[filtered_meta$Sample]
  if (length(filtered_data) == 0) {
    message(sprintf("No %s/%s samples -- skipping %s.",
                    sampletype_keep, chain_keep, data_csv))
    return(invisible(NULL))
  }

  div_hill <- repDiversity(filtered_data, "hill")
  plot_data_div <- div_hill %>%
    as.data.frame() %>%
    left_join(filtered_meta %>% dplyr::select(Sample, patient, timepoint),
              by = "Sample") %>%
    dplyr::filter(Q == 1)

  write.csv(plot_data_div, file.path(outdir_csv, data_csv), row.names = FALSE)

  if (length(unique(plot_data_div$patient)) >= 2 &&
      nrow(plot_data_div) >= length(time_levels) + 2) {
    model_div    <- lmer(Value ~ timepoint + (1 | patient), data = plot_data_div)
    emm_div      <- emmeans(model_div, ~ timepoint)
    pairs_result <- as.data.frame(pairs(emm_div, adjust = "holm"))
    pairs_result$comp_key <- vapply(
      strsplit(as.character(pairs_result$contrast), " - ", fixed = TRUE),
      function(x) paste(sort(x), collapse = " - "),
      character(1)
    )
    write.csv(pairs_result, file.path(outdir_csv, stats_csv), row.names = FALSE)
  } else {
    message(sprintf("Not enough samples to fit LMM for %s -- skipping %s.",
                    data_csv, stats_csv))
  }
  invisible(plot_data_div)
}

# Renders the diversity panel from the CSVs written by `compute_diversity_csv`.
plot_diversity_from_csv <- function(data_csv, stats_csv, time_levels, colors,
                                    title_text, y_label = "Hill Q=1 diversity",
                                    with_boxplot = TRUE,
                                    width_cm = 6, height_cm = 4.2,
                                    out_pdf) {
  data_path  <- file.path(outdir_csv, data_csv)
  stats_path <- file.path(outdir_csv, stats_csv)
  if (!file.exists(data_path)) {
    warning("Missing ", data_path, " -- skipping ", out_pdf)
    return(invisible(NULL))
  }
  plot_data_div <- read.csv(data_path, check.names = FALSE) %>%
    mutate(timepoint = factor(timepoint, levels = time_levels))

  p <- ggplot(plot_data_div, aes(x = timepoint, y = Value)) +
    geom_line(aes(group = patient), alpha = 0.3, color = "gray50",
              linewidth = 0.2)
  if (with_boxplot) {
    p <- p + geom_boxplot(aes(fill = timepoint),
                          alpha = 0.7, outlier.shape = NA,
                          width = 0.5, linewidth = 0.2)
  }
  p <- p +
    geom_jitter(aes(color = timepoint),
                width = 0.15, height = 0, size = 0.4, alpha = 0.9) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    coord_cartesian(ylim = c(0, NA)) +
    labs(x = "Timepoint", y = y_label, title = title_text) +
    theme_pub()

  if (file.exists(stats_path)) {
    stats_df    <- read.csv(stats_path, check.names = FALSE)
    comparisons <- combn(time_levels, 2, simplify = FALSE)
    comp_keys   <- vapply(comparisons,
                          function(x) paste(sort(x), collapse = " - "),
                          character(1))
    pvals       <- stats_df$p.value[match(comp_keys, stats_df$comp_key)]
    p_labels    <- paste0("adj p=", format.pval(pvals, digits = 2, eps = 1e-3))

    ymax_div <- max(plot_data_div$Value, na.rm = TRUE)
    y_pos <- if (length(time_levels) == 2) {
      ymax_div * 0.90
    } else {
      # Stagger the six pairwise comparison bars across two y-tiers.
      c(ymax_div * 0.75, ymax_div * 0.83, ymax_div * 0.92,
        ymax_div * 0.50, ymax_div * 0.58, ymax_div * 0.67)
    }
    p <- p + geom_signif(
      comparisons = comparisons,
      annotations = p_labels,
      y_position  = y_pos,
      tip_length  = 0,
      textsize    = 1.8,
      size        = 0.3
    )
  }

  ggsave(file.path(outdir_plots, out_pdf), p,
         width = width_cm, height = height_cm, units = "cm")
  invisible(p)
}

# ---------------------------------------------------------------------------
# Tumor DX -> PT arm contingency (Extended Data Fig 1d,e,f panel)
# ---------------------------------------------------------------------------

# For each chain, tabulates whether each patient's TRA / TRB / TRG diversity
# rises or falls from DX to PT, split by treatment arm, and computes a 2 x 2
# Fisher's exact test.
compute_arm_contingency_csv <- function(tumor_csv_map) {
  rows <- list()
  for (chain_keep in names(tumor_csv_map)) {
    info      <- tumor_csv_map[[chain_keep]]
    data_path <- file.path(outdir_csv, info$data)
    if (!file.exists(data_path)) next

    df <- read.csv(data_path, check.names = FALSE) %>%
      mutate(arm = case_when(
        patient %in% c("Patient1", "Patient2")             ~ "A",
        patient %in% c("Patient3", "Patient4", "Patient5") ~ "B",
        TRUE                                                ~ "unknown"
      )) %>%
      pivot_wider(id_cols = c(patient, arm),
                  names_from = timepoint, values_from = Value) %>%
      dplyr::filter(!is.na(tDX) & !is.na(tPT)) %>%
      mutate(direction = ifelse(tPT > tDX, "increase", "decrease"))

    counts    <- table(df$arm, df$direction)
    arm_a_inc <- if ("A" %in% rownames(counts) && "increase" %in% colnames(counts)) counts["A","increase"] else 0
    arm_a_dec <- if ("A" %in% rownames(counts) && "decrease" %in% colnames(counts)) counts["A","decrease"] else 0
    arm_b_inc <- if ("B" %in% rownames(counts) && "increase" %in% colnames(counts)) counts["B","increase"] else 0
    arm_b_dec <- if ("B" %in% rownames(counts) && "decrease" %in% colnames(counts)) counts["B","decrease"] else 0
    mat <- matrix(c(arm_a_inc, arm_a_dec, arm_b_inc, arm_b_dec),
                  nrow = 2, byrow = TRUE,
                  dimnames = list(c("A","B"), c("increase","decrease")))
    fisher_p <- tryCatch(fisher.test(mat)$p.value,
                         error = function(e) NA_real_)

    rows[[chain_keep]] <- data.frame(
      Chain            = info$chain_display,
      `Arm A increase` = arm_a_inc,
      `Arm A decrease` = arm_a_dec,
      `Arm B increase` = arm_b_inc,
      `Arm B decrease` = arm_b_dec,
      `Fisher's P`     = fisher_p,
      check.names      = FALSE
    )
  }
  out <- do.call(rbind, rows)
  write.csv(out, file.path(outdir_csv, "ext_fig1_arm_contingency.csv"),
            row.names = FALSE)
  invisible(out)
}

# ---------------------------------------------------------------------------
# Clonotype tracking (Fig 5b/c/f/g)
# ---------------------------------------------------------------------------

# Tracks one target CDR3 amino-acid sequence across the samples matching
# `sample_filter_regex`, then writes a wide proportion table and a long-form
# table of per-sample clone counts and proportions used for stratum labels.
compute_clonotype_tracking_csv <- function(data_list, sample_filter_regex,
                                           chain_keep, target_aa,
                                           out_csv, out_counts_csv) {
  chain_pat <- paste0("\\.clones_", chain_keep, "$")
  keys <- grep(chain_pat, names(data_list), value = TRUE)
  keys <- keys[grepl(sample_filter_regex, keys)]
  if (length(keys) == 0) {
    warning("No samples match ", sample_filter_regex, " on ", chain_keep)
    return(invisible(NULL))
  }
  sub_data <- data_list[keys]
  names(sub_data) <- sub(chain_pat, "", names(sub_data))

  tc <- trackClonotypes(sub_data, target_aa, .col = "aa")
  write.csv(tc, file.path(outdir_csv, out_csv), row.names = FALSE)

  counts_df <- bind_rows(lapply(names(sub_data), function(s) {
    df  <- sub_data[[s]]
    hit <- df[df$CDR3.aa == target_aa, , drop = FALSE]
    if (nrow(hit) == 0) {
      return(data.frame(Sample = s, Clones = 0L, Proportion = 0))
    }
    data.frame(
      Sample     = s,
      Clones     = as.integer(hit$Clones[1]),
      Proportion = as.numeric(hit$Proportion[1])
    )
  }))
  write.csv(counts_df, file.path(outdir_csv, out_counts_csv), row.names = FALSE)
  invisible(tc)
}

# Renders the alluvial tracking plot for one target clonotype. The y-axis is
# clipped to `y_max` because the tracked proportions are small relative to the
# overall repertoire, and each stratum is annotated with the per-sample clone
# count and percentage.
plot_clonotype_tracking_from_csv <- function(out_csv, out_counts_csv, target_aa,
                                             fill_color = "#d7b5d8",
                                             y_max      = NA_real_,
                                             width_cm, height_cm, out_pdf) {
  csv_path    <- file.path(outdir_csv, out_csv)
  counts_path <- file.path(outdir_csv, out_counts_csv)
  if (!file.exists(csv_path)) {
    warning("Missing ", csv_path, " -- skipping ", out_pdf)
    return(invisible(NULL))
  }
  
  tc          <- read.csv(csv_path, check.names = FALSE)
  sample_cols <- setdiff(colnames(tc), "CDR3.aa")
  
  long_tc <- tc %>%
    tidyr::pivot_longer(cols = all_of(sample_cols),
                        names_to = "Sample", values_to = "Proportion") %>%
    mutate(Sample    = factor(Sample, levels = sample_cols),
           Clonotype = factor(CDR3.aa, levels = target_aa)) %>%
    dplyr::select(Sample, Clonotype, Proportion)
  
  # Per-sample "<clones>\n(<percentage>%)" labels for each alluvial stratum.
  if (file.exists(counts_path)) {
    counts_df <- read.csv(counts_path, check.names = FALSE)
    lab_df <- counts_df %>%
      mutate(lab = ifelse(Proportion > 0,
                          sprintf("%d\n(%.2f%%)",
                                  as.integer(Clones),
                                  100 * as.numeric(Proportion)),
                          ""))
  } else {
    lab_df <- long_tc %>%
      mutate(lab = ifelse(Proportion > 0,
                          sprintf("%.2f%%", 100 * Proportion),
                          ""))
  }
  lab_df$Sample    <- factor(lab_df$Sample, levels = sample_cols)
  lab_df$Clonotype <- factor(target_aa, levels = target_aa)
  
  df_plot <- long_tc %>%
    left_join(lab_df %>% dplyr::select(Sample, Clonotype, lab),
              by = c("Sample", "Clonotype")) %>%
    mutate(lab = ifelse(is.na(lab), "", lab))
  
  df_plot <- df_plot %>%
    mutate(Proportion = ifelse(Proportion <= 0, 1e-10, Proportion))
  
  p_alluv <- ggplot(df_plot,
                    aes(x = Sample, y = Proportion,
                        alluvium = Clonotype, stratum = Clonotype,
                        fill = Clonotype)) +
    ggalluvial::geom_alluvium(alpha = 0.6) +
    ggalluvial::geom_stratum(alpha = 1, color = "white", linewidth = 0) +
    geom_text(stat = "stratum",
              aes(label = lab),
              color = "#000000", size = 1.8) +
    scale_fill_manual(values = setNames(fill_color, target_aa)) +
    labs(x = "Sample", y = "Proportion",
         title = paste("Clonotype tracking:", target_aa)) +
    theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text   = element_text(size = 6),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      axis.title  = element_text(size = 6),
      plot.title  = element_text(size = 6),
      legend.position = "none"
    )
  
  if (!is.na(y_max)) {
    p_alluv <- p_alluv + coord_cartesian(ylim = c(0, y_max))
  }
  
  ggsave(file.path(outdir_plots, out_pdf), p_alluv,
         width = width_cm, height = height_cm, units = "cm")
  invisible(p_alluv)
}

# ---------------------------------------------------------------------------
# Per-figure configuration
# ---------------------------------------------------------------------------
PBMC_TIMEPOINTS  <- c("pTx", "pI4", "pSu", "ePC")
TUMOR_TIMEPOINTS <- c("tDX", "tPT")
batlow_colors_4  <- c("#d7b5d8", "#df65b0", "#dd1c77", "#980043")
tumor_colors_2   <- c("#c994c7", "#dd1c77")

pbmc_csv_map <- list(
  TRA = list(data = "ext_fig1a_data.csv", stats = "ext_fig1a_stats.csv",
             pdf  = "ext_fig1a.pdf"),
  TRB = list(data = "ext_fig1b_data.csv", stats = "ext_fig1b_stats.csv",
             pdf  = "ext_fig1b.pdf"),
  TRG = list(data = "ext_fig1c_data.csv", stats = "ext_fig1c_stats.csv",
             pdf  = "ext_fig1c.pdf")
)

tumor_csv_map <- list(
  TRA = list(data = "ext_fig1d_data.csv", stats = "ext_fig1d_stats.csv",
             pdf  = "ext_fig1d.pdf", chain_display = "\u03b1"),  # α
  TRB = list(data = "ext_fig1e_data.csv", stats = "ext_fig1e_stats.csv",
             pdf  = "ext_fig1e.pdf", chain_display = "\u03b2"),  # β
  TRG = list(data = "ext_fig1f_data.csv", stats = "ext_fig1f_stats.csv",
             pdf  = "ext_fig1f.pdf", chain_display = "\u03b3")   # γ
)

# `patients` lists the patient IDs whose samples enter each clonotype-tracking
# figure; samples are selected by matching "PatientN-" against the file name.
clonotype_map <- list(
  fig5b = list(
    target     = "CATWDRRKKLF",
    chain      = "TRG",
    patients   = c("Patient2", "Patient17", "Patient18",
                   "Patient23", "Patient26", "Patient27"),
    csv        = "fig5b.csv",
    counts_csv = "fig5b_counts.csv",
    pdf        = "fig5b.pdf",
    ymax       = 0.025,
    width      = 11, height = 10
  ),
  fig5c = list(
    target     = "CALWEVQELGKKIKVF",
    chain      = "TRG",
    patients   = c("Patient1", "Patient3", "Patient4", "Patient5"),
    csv        = "fig5c.csv",
    counts_csv = "fig5c_counts.csv",
    pdf        = "fig5c.pdf",
    ymax       = 0.1,
    width      = 9,  height = 9
  ),
  fig5f = list(
    target     = "CAAKQAGYSTLTF",
    chain      = "TRA",
    patients   = c("Patient2"),
    csv        = "fig5f.csv",
    counts_csv = "fig5f_counts.csv",
    pdf        = "fig5f.pdf",
    ymax       = 0.18,
    width      = 10, height = 10
  ),
  fig5g = list(
    target     = "CASSVIAETYEQYF",
    chain      = "TRB",
    patients   = c("Patient2"),
    csv        = "fig5g.csv",
    counts_csv = "fig5g_counts.csv",
    pdf        = "fig5g.pdf",
    ymax       = 0.3,
    width      = 10, height = 10
  )
)

build_patient_regex <- function(patient_ids) {
  paste0(patient_ids, "-", collapse = "|")
}

# ---------------------------------------------------------------------------
# Compute CSVs from the raw MiXCR clones output
# ---------------------------------------------------------------------------
if (run_compute_stage) {
  message("repLoad(", input_repdir, ")")
  anbl_repdata <- repLoad(input_repdir)
  samples <- names(anbl_repdata$data)
  data    <- anbl_repdata$data[samples]
  meta    <- anbl_repdata$meta[anbl_repdata$meta$Sample %in% samples, ]
  meta    <- add_metadata_columns(meta)

  # PBMC diversity (Ext Fig 1a/b/c)
  for (chain_keep in names(pbmc_csv_map)) {
    m <- pbmc_csv_map[[chain_keep]]
    compute_diversity_csv(data, meta, "PBMC", chain_keep,
                          PBMC_TIMEPOINTS, m$data, m$stats)
  }

  # Tumor diversity (Ext Fig 1d/e/f)
  for (chain_keep in names(tumor_csv_map)) {
    m <- tumor_csv_map[[chain_keep]]
    compute_diversity_csv(data, meta, "tumor", chain_keep,
                          TUMOR_TIMEPOINTS, m$data, m$stats)
  }

  # Arm contingency (Ext Fig 1 panel)
  compute_arm_contingency_csv(tumor_csv_map)

  # Clonotype tracking (Fig 5b/c/f/g)
  for (k in names(clonotype_map)) {
    info <- clonotype_map[[k]]
    compute_clonotype_tracking_csv(data,
                                   build_patient_regex(info$patients),
                                   info$chain,
                                   info$target,
                                   info$csv,
                                   info$counts_csv)
  }
} else {
  message("input_repdir not found, skipping recomputation: ", input_repdir)
}

# ---------------------------------------------------------------------------
# Render figures from CSVs
# ---------------------------------------------------------------------------
# PBMC diversity (Ext Fig 1a/b/c)
for (chain_keep in names(pbmc_csv_map)) {
  m <- pbmc_csv_map[[chain_keep]]
  plot_diversity_from_csv(
    data_csv     = m$data,
    stats_csv    = m$stats,
    time_levels  = PBMC_TIMEPOINTS,
    colors       = batlow_colors_4,
    title_text   = paste0("PBMC Hill Q=1 diversity (", chain_keep, ")"),
    with_boxplot = TRUE,
    width_cm     = 6,
    height_cm    = 4.2,
    out_pdf      = m$pdf
  )
}

# Tumor diversity (Ext Fig 1d/e/f)
for (chain_keep in names(tumor_csv_map)) {
  m <- tumor_csv_map[[chain_keep]]
  plot_diversity_from_csv(
    data_csv     = m$data,
    stats_csv    = m$stats,
    time_levels  = TUMOR_TIMEPOINTS,
    colors       = tumor_colors_2,
    title_text   = paste0("Tumor Hill Q=1 diversity (", chain_keep, ")"),
    with_boxplot = FALSE,
    width_cm     = 4,
    height_cm    = 4.2,
    out_pdf      = m$pdf
  )
}

# Clonotype tracking (Fig 5b/c/f/g)
for (k in names(clonotype_map)) {
  info <- clonotype_map[[k]]
  plot_clonotype_tracking_from_csv(
    out_csv        = info$csv,
    out_counts_csv = info$counts_csv,
    target_aa      = info$target,
    y_max          = info$ymax,
    width_cm       = info$width,
    height_cm      = info$height,
    out_pdf        = info$pdf
  )
}
