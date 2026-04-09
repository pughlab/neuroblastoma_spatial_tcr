# IMMUNARCH TCR ANALYSIS ----

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

# PATHS AND RUN CONFIGURATION ----
project_dir = "/path/to/xenium/"
input_repdir = "/path/to/xenium/mixcr4/analyze_output/"
outdir_csv = file.path(project_dir, "outdir")
outdir_plots = file.path(project_dir, "Rplots")
outdir_logs = file.path(project_dir, "logs")

dir.create(outdir_csv, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_logs, showWarnings = FALSE, recursive = TRUE)

# HELPERS ----
theme_pub = function() {
  theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 6),
      plot.title = element_text(size = 6),
      axis.line.x = element_line(linewidth = 0.5),
      axis.line.y = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5)
    )
}

add_metadata_columns = function(meta_df) {
  meta_df %>%
    mutate(
      patient = str_extract(Sample, "(?:ANBL|PRO)-[A-Z]+-([A-Za-z0-9]+)") %>%
        str_replace(".*-([A-Za-z0-9]+)$", "\\1"),

      chain = case_when(
        grepl("TRB", Sample) ~ "TRB",
        grepl("TRAD", Sample) ~ "TRAD",
        grepl("TRA", Sample) ~ "TRA",
        grepl("TRG", Sample) ~ "TRG",
        grepl("TRD", Sample) ~ "TRD"
      ),

      timepoint = case_when(
        grepl("-T0-", Sample) ~ "pTx",
        grepl("-pTx-", Sample) ~ "pTx",
        grepl("-T1-", Sample) ~ "pI4",
        grepl("-T2-", Sample) ~ "ePC",
        grepl("-T3-", Sample) ~ "pI3",
        grepl("-T4-", Sample) ~ "pSu",
        grepl("-T5-", Sample) ~ "Pgn",
        grepl("-T6-", Sample) ~ "Rel",
        grepl("-T7-", Sample) ~ "EOI",
        grepl("-T8-", Sample) ~ "pI2",
        grepl("-T9-", Sample) ~ "M72",
        grepl("-DX-", Sample) ~ "tDX",
        grepl("-PT-", Sample) ~ "tPT",
        grepl("-Ch1-", Sample) ~ "Ch1",
        grepl("-Ch2-", Sample) ~ "Ch2",
        grepl("-Ch3-", Sample) ~ "Ch3",
        grepl("-SC1-", Sample) ~ "SC1",
        grepl("-SC2-", Sample) ~ "SC2",
        grepl("-SC3-", Sample) ~ "SC3",
        grepl("-CT1-", Sample) ~ "CT1",
        grepl("-CT2-", Sample) ~ "CT2",
        grepl("-CT3-", Sample) ~ "CT3",
        grepl("-PT1-", Sample) ~ "PT1",
        grepl("-PT2-", Sample) ~ "PT2",
        TRUE ~ "NA"
      ),

      project = case_when(
        grepl("^ANBL", Sample) ~ "ANBL",
        grepl("^PRO", Sample) ~ "PRO",
        TRUE ~ "unknown"
      ),

      sampletype = case_when(
        project == "PRO" ~ "PBMC",
        grepl("^ANBL-CHP-[A-Z]{6}", Sample) ~ "tumor",
        grepl("^ANBL-SM-[A-Z]{6}", Sample) ~ "cfDNA",
        grepl("^ANBL-CHP-[0-9]{6}", Sample) ~ "PBMC",
        grepl("^ANBL-SM-[0-9]{6}", Sample) ~ "PBMC",
        TRUE ~ "unknown"
      ),

      arm = case_when(
        grepl("PBBCHF|PAZNRG|886457|907322", Sample) ~ "A",
        grepl("PBBKFP|PAZWZN|PBADJC|890689|911891|894374", Sample) ~ "B",
        TRUE ~ "unknown"
      )
    )
}

get_filtered_inputs = function(meta_df, data_list, sampletype_keep, chain_keep, time_levels) {
  filtered_meta = meta_df %>%
    filter(timepoint %in% time_levels) %>%
    filter(sampletype == sampletype_keep) %>%
    filter(chain == chain_keep) %>%
    mutate(timepoint = factor(timepoint, levels = time_levels))

  filtered_data = data_list[filtered_meta$Sample]

  list(meta = filtered_meta, data = filtered_data)
}

make_emm_labels = function(model_obj, comparisons_list, adjust_method = "holm") {
  emm = emmeans(model_obj, ~ timepoint)
  pairs_df = as.data.frame(pairs(emm, adjust = adjust_method))

  pairs_df$comp_key = vapply(
    strsplit(as.character(pairs_df$contrast), " - ", fixed = TRUE),
    function(x) paste(sort(x), collapse = " - "),
    character(1)
  )

  comparison_keys = vapply(
    comparisons_list,
    function(x) paste(sort(x), collapse = " - "),
    character(1)
  )

  pvals_ordered = pairs_df$p.value[match(comparison_keys, pairs_df$comp_key)]
  p_labels = paste0("adj p=", format.pval(pvals_ordered, digits = 2, eps = 1e-3))

  list(pairs_df = pairs_df, p_labels = p_labels)
}

fit_lmm = function(df, y_var) {
  formula_obj = as.formula(paste0(y_var, " ~ timepoint + (1 | patient)"))
  lmer(formula_obj, data = df)
}

make_metric_plot = function(plot_data, y_var, title_text, y_label, colors, y_lim = c(0, NA), with_boxplot = TRUE) {
  p = ggplot(plot_data, aes(x = timepoint, y = .data[[y_var]])) +
    geom_line(aes(group = patient), alpha = 0.3, color = "gray50", linewidth = 0.2)

  if (with_boxplot) {
    p = p + geom_boxplot(aes(fill = timepoint), alpha = 0.7, outlier.shape = NA, width = 0.5, linewidth = 0.2)
  }

  p = p +
    geom_jitter(aes(color = timepoint), width = 0.15, height = 0, size = 0.4, alpha = 0.9) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    coord_cartesian(ylim = y_lim) +
    labs(x = "Timepoint", y = y_label, title = title_text) +
    theme_pub()

  p
}

classify_clones = function(data_list, meta_df, thresholds = c(Small = 0.0001, Medium = 0.001, Large = 0.01, Hyperexpanded = 1)) {
  clone_classifications = list()

  for (sample_name in names(data_list)) {
    sample_data = data_list[[sample_name]]
    total_reads = sum(sample_data$Clones)

    sample_data = sample_data %>%
      mutate(
        Proportion = Clones / total_reads,
        Sample = sample_name,
        Clone_Type = case_when(
          Proportion >= thresholds["Large"] ~ "Hyperexpanded",
          Proportion >= thresholds["Medium"] ~ "Large",
          Proportion >= thresholds["Small"] ~ "Medium",
          TRUE ~ "Small"
        )
      )

    clone_classifications[[sample_name]] = sample_data
  }

  all_clones = bind_rows(clone_classifications) %>%
    left_join(meta_df %>% select(Sample, timepoint), by = "Sample") %>%
    mutate(Clone_Type = factor(Clone_Type, levels = c("Small", "Medium", "Large", "Hyperexpanded")))

  all_clones
}

summarize_clone_classes = function(clone_class_df) {
  clone_summary = clone_class_df %>%
    group_by(Sample, timepoint, Clone_Type) %>%
    summarise(
      n_clones = n(),
      total_reads = sum(Clones),
      mean_proportion = mean(Proportion),
      .groups = "drop"
    ) %>%
    arrange(timepoint, Sample, Clone_Type)

  timepoint_summary = clone_class_df %>%
    group_by(timepoint, Clone_Type) %>%
    summarise(
      n_clones = n(),
      n_samples = n_distinct(Sample),
      mean_proportion = mean(Proportion),
      sd_proportion = sd(Proportion),
      .groups = "drop"
    ) %>%
    arrange(timepoint, Clone_Type)

  list(clone_summary = clone_summary, timepoint_summary = timepoint_summary)
}

make_clonality_plot = function(filtered_data, filtered_meta, time_levels, colors) {
  imm_hom = repClonality(
    filtered_data,
    .method = "homeo",
    .clone.types = c(Small = 0.0001, Medium = 0.001, Large = 0.01, Hyperexpanded = 1)
  )

  plot_data_hom = imm_hom %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    pivot_longer(cols = -Sample, names_to = "Clone_Type", values_to = "Proportion") %>%
    left_join(filtered_meta %>% select(Sample, timepoint), by = "Sample") %>%
    mutate(
      timepoint = factor(timepoint, levels = time_levels),
      Clone_Type = case_when(
        grepl("Small", Clone_Type) ~ "Small",
        grepl("Medium", Clone_Type) ~ "Medium",
        grepl("Large", Clone_Type) ~ "Large",
        grepl("Hyperexpanded", Clone_Type) ~ "Hyperexpanded",
        TRUE ~ Clone_Type
      ),
      Clone_Type = factor(Clone_Type, levels = c("Small", "Medium", "Large", "Hyperexpanded"))
    )

  plot_data_hom_summary = plot_data_hom %>%
    group_by(timepoint, Clone_Type) %>%
    summarise(
      Mean_Proportion = mean(Proportion),
      SE = sd(Proportion) / sqrt(n()),
      .groups = "drop"
    ) %>%
    group_by(timepoint) %>%
    arrange(desc(Clone_Type)) %>%
    mutate(
      y_start = cumsum(lag(Mean_Proportion, default = 0)),
      y_end = cumsum(Mean_Proportion),
      y_pos = (y_start + y_end) / 2,
      label = ifelse(Mean_Proportion >= 0.1, paste0(round(Mean_Proportion, 3)), "")
    ) %>%
    ungroup()

  p_hom = ggplot(plot_data_hom_summary, aes(x = timepoint, y = Mean_Proportion, fill = Clone_Type)) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
    geom_text(aes(y = y_pos, label = label), position = position_identity(), size = 1.8, color = "white") +
    scale_fill_manual(values = rev(colors)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Timepoint", y = "Mean Proportion", title = "Clonal Homeostasis", fill = "Clone Type") +
    theme_pub()

  list(plot = p_hom, summary_df = plot_data_hom_summary)
}

run_compartment_analysis = function(
  compartment_label,
  sampletype_keep,
  time_levels,
  data_list,
  meta_df,
  colors,
  volume_ylim,
  diversity_ylim,
  with_boxplot_volume,
  with_boxplot_diversity,
  plot_width_cm,
  plot_height_cm,
  file_tag
) {
  inputs = get_filtered_inputs(
    meta_df = meta_df,
    data_list = data_list,
    sampletype_keep = sampletype_keep,
    chain_keep = "TRB",
    time_levels = time_levels
  )

  filtered_meta = inputs$meta
  filtered_data = inputs$data

  # VOLUME ----
  exp_vol = repExplore(filtered_data, .method = "volume")
  plot_data_vol = exp_vol %>%
    as.data.frame() %>%
    left_join(filtered_meta %>% select(Sample, patient, timepoint), by = "Sample")

  write.csv(plot_data_vol, file.path(outdir_csv, paste0(file_tag, "_volume_plot_data.csv")), row.names = FALSE)

  p_vol = make_metric_plot(
    plot_data = plot_data_vol,
    y_var = "Volume",
    title_text = paste0(compartment_label, " Repertoire Volume"),
    y_label = "Volume",
    colors = colors,
    y_lim = volume_ylim,
    with_boxplot = with_boxplot_volume
  )

  model_vol = fit_lmm(plot_data_vol, "Volume")
  capture.output(summary(model_vol), file = file.path(outdir_logs, paste0(file_tag, "_volume_model_summary.txt")))
  capture.output(anova(model_vol), file = file.path(outdir_logs, paste0(file_tag, "_volume_model_anova.txt")))

  comparisons = combn(time_levels, 2, simplify = FALSE)
  vol_stats = make_emm_labels(model_vol, comparisons, adjust_method = "holm")

  write.csv(vol_stats$pairs_df, file.path(outdir_csv, paste0(file_tag, "_volume_statistical_results.csv")), row.names = FALSE)

  if (length(time_levels) == 2) {
    y_positions = c(volume_ylim[2] * 0.9)
  } else {
    y_positions = seq(volume_ylim[2] * 0.5, volume_ylim[2] * 0.92, length.out = length(comparisons))
  }

  p_vol_stats = p_vol +
    geom_signif(
      comparisons = comparisons,
      annotations = vol_stats$p_labels,
      y_position = y_positions,
      tip_length = 0,
      textsize = 1.8,
      size = 0.3
    )

  ggsave(
    filename = file.path(outdir_plots, paste0(file_tag, "_volume.pdf")),
    plot = p_vol_stats,
    width = plot_width_cm,
    height = plot_height_cm,
    units = "cm"
  )

  # CLONALITY ----
  clonality = make_clonality_plot(filtered_data, filtered_meta, time_levels, colors)
  ggsave(
    filename = file.path(outdir_plots, paste0(file_tag, "_clonal_homeostasis.pdf")),
    plot = clonality$plot,
    width = ifelse(length(time_levels) > 2, 3.6, 2.5),
    height = 4.2,
    units = "cm"
  )

  clone_classification_results = classify_clones(filtered_data, filtered_meta)
  clone_summ = summarize_clone_classes(clone_classification_results)

  write.csv(
    clone_summ$timepoint_summary,
    file.path(outdir_csv, paste0(file_tag, "_clone_classification_results.csv")),
    row.names = FALSE
  )

  # DIVERSITY ----
  div_hill = repDiversity(filtered_data, "hill")
  plot_data_div = div_hill %>%
    as.data.frame() %>%
    left_join(filtered_meta %>% select(Sample, patient, timepoint), by = "Sample") %>%
    filter(Q == 1)

  write.csv(plot_data_div, file.path(outdir_csv, paste0(file_tag, "_div_plot_data.csv")), row.names = FALSE)

  p_div = make_metric_plot(
    plot_data = plot_data_div,
    y_var = "Value",
    title_text = paste0(compartment_label, " Shannon Diversity"),
    y_label = "Shannon diversity",
    colors = colors,
    y_lim = diversity_ylim,
    with_boxplot = with_boxplot_diversity
  )

  model_div = fit_lmm(plot_data_div, "Value")
  capture.output(summary(model_div), file = file.path(outdir_logs, paste0(file_tag, "_diversity_model_summary.txt")))
  capture.output(anova(model_div), file = file.path(outdir_logs, paste0(file_tag, "_diversity_model_anova.txt")))

  div_stats = make_emm_labels(model_div, comparisons, adjust_method = "holm")
  write.csv(div_stats$pairs_df, file.path(outdir_csv, paste0(file_tag, "_div_statistical_results.csv")), row.names = FALSE)

  ymax_div = max(plot_data_div$Value, na.rm = TRUE)
  if (length(time_levels) == 2) {
    y_positions_div = c(ymax_div * 0.9)
  } else {
    y_positions_div = seq(ymax_div * 0.5, ymax_div * 0.92, length.out = length(comparisons))
  }

  p_div_stats = p_div +
    geom_signif(
      comparisons = comparisons,
      annotations = div_stats$p_labels,
      y_position = y_positions_div,
      tip_length = 0,
      textsize = 1.8,
      size = 0.3
    )

  ggsave(
    filename = file.path(outdir_plots, paste0(file_tag, "_diversity.pdf")),
    plot = p_div_stats,
    width = plot_width_cm,
    height = plot_height_cm,
    units = "cm"
  )

  invisible(
    list(
      filtered_meta = filtered_meta,
      filtered_data = filtered_data,
      volume_data = plot_data_vol,
      diversity_data = plot_data_div
    )
  )
}

run_clonotype_tracking = function(data_list, target_aa, out_csv_name, out_plot_name) {
  tcr_keys = grep("\\.clones_TRG$", names(data_list), value = TRUE)
  keep_keys = tcr_keys[grepl("CHP.*(PAZWZN|PBBKFP|PAZNRG|PBADJC)", tcr_keys)]

  imm_sub = list(
    data = data_list[keep_keys],
    meta = NULL
  )

  names(imm_sub$data) = sub("\\.clones_TRG$", "", names(imm_sub$data))

  tc = trackClonotypes(
    imm_sub$data,
    target_aa,
    .col = "aa"
  )

  write.csv(tc, file.path(outdir_csv, out_csv_name), row.names = FALSE)

  p = vis(tc)

  lab_df = bind_rows(lapply(names(imm_sub$data), function(s) {
    df = imm_sub$data[[s]]
    hit = df[df$CDR3.aa == target_aa, , drop = FALSE]
    if (nrow(hit) == 0) {
      return(data.frame(Sample = s, lab = "", stringsAsFactors = FALSE))
    }
    data.frame(
      Sample = s,
      lab = sprintf("%d\n(%.2f)", as.integer(hit$Clones[1]), as.numeric(hit$Proportion[1])),
      stringsAsFactors = FALSE
    )
  }))

  p_alluv = p
  p_alluv$data = as.data.frame(p_alluv$data, check.names = FALSE) %>%
    left_join(lab_df, by = "Sample")

  # ALLUVIAL CUSTOMIZATION ----
  p_alluv = p_alluv +
    scale_fill_manual(values = setNames("#d7b5d8", target_aa)) +
    labs(x = "Timepoint", y = "Mean Proportion", title = "Clonotype tracking") +
    theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 6),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      axis.title = element_text(size = 6),
      plot.title = element_text(size = 6),
      legend.position = "none"
    ) +
    geom_text(
      stat = "stratum",
      aes(label = lab),
      color = "#000000",
      size = 1.8
    )

  p_alluv$layers[[1]]$aes_params$alpha = 0.6
  p_alluv$layers[[2]]$aes_params$alpha = 1
  p_alluv$layers[[2]]$aes_params$color = "white"
  p_alluv$layers[[2]]$aes_params$linewidth = 0

  ggsave(
    filename = file.path(outdir_plots, out_plot_name),
    plot = p_alluv,
    width = 11,
    height = 10,
    units = "cm"
  )

  invisible(tc)
}

# MAIN WORKFLOW ----
anbl_repdata = repLoad(input_repdir)

samples = names(anbl_repdata$data)
data = anbl_repdata$data[samples]
meta = anbl_repdata$meta[anbl_repdata$meta$Sample %in% samples, ]
meta = add_metadata_columns(meta)

batlow_colors_4 = c("#d7b5d8", "#df65b0", "#dd1c77", "#980043")
tumor_colors_2 = c("#c994c7", "#dd1c77")

# PBMC ANALYSES ----
pbmc_results = run_compartment_analysis(
  compartment_label = "PBMC",
  sampletype_keep = "PBMC",
  time_levels = c("pTx", "pI4", "pSu", "ePC"),
  data_list = data,
  meta_df = meta,
  colors = batlow_colors_4,
  volume_ylim = c(0, 1800),
  diversity_ylim = c(0, NA),
  with_boxplot_volume = TRUE,
  with_boxplot_diversity = TRUE,
  plot_width_cm = 6,
  plot_height_cm = 4.2,
  file_tag = "fig5_pbmc"
)

# TUMOR ANALYSES ----
tumor_results = run_compartment_analysis(
  compartment_label = "Tumor",
  sampletype_keep = "tumor",
  time_levels = c("tDX", "tPT"),
  data_list = data,
  meta_df = meta,
  colors = tumor_colors_2,
  volume_ylim = c(0, 1100),
  diversity_ylim = c(0, NA),
  with_boxplot_volume = FALSE,
  with_boxplot_diversity = FALSE,
  plot_width_cm = 4,
  plot_height_cm = 4.2,
  file_tag = "fig5_tumor"
)

# CLONOTYPE TRACKING ----
run_clonotype_tracking(
  data_list = data,
  target_aa = "CALWEVQELGKKIKVF",
  out_csv_name = "fig5j_clonotype_tracking_data_clone17.csv",
  out_plot_name = "clonetracking_17.pdf"
)