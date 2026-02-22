# LOAD PACKAGES AND PATHS ----
required_packages = c(
  "igraph",
  "dplyr",
  "stringr",
  "readr",
  "tibble",
  "ggplot2",
  "ggseqlogo",
  "janitor"
)

missing_packages = required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggseqlogo)
  library(janitor)
})

gliph_path = "/Users/yyj/Doc/mount/h4h/projects/pediatric/projects/dod/tcr_analysis/20251015/1_gliph"
gliph_fname = "gliph_output/anbl_vdj_cluster.csv"
output_path = "community_detection_output"
gliph_input_fname = "gliph_input/anbl_gliphinput_20251015.txt"

run_diagnostics = FALSE
run_example_plots = FALSE
example_community_id = 15231

group_col = c(
  "Viral" = "#001959",
  "Viral and bacterial" = "#0C325E",
  "Bacterial" = "#124560",
  "Cross-reactive" = "#185461",
  "Human" = "#26635F",
  "Adult healthy" = "#3B6C55",
  "Adult and pediatric healthy" = "#577646",
  "Pediatric healthy" = "#727E38",
  "Adult" = "#91852C",
  "Adult and pediatric" = "#B28D2E",
  "Pediatric" = "#D29243",
  "Adult cancer" = "#EC9961",
  "Adult and pediatric cancer" = "#FAA588",
  "Pediatric cancer" = "#FDB2AE",
  "Neuroblastoma" = "#FCBFD2",
  "ANBL" = "#F9CCF9"
)

# HELPERS ----
normalize_tokens = function(x) {
  str_replace_all(
    x,
    c(
      "T-cellmalignancies" = "T.cellmalignancies",
      "S-pneumoniae" = "S.pneumoniae",
      "_V" = "V"
    )
  )
}

safe_token_extract = function(x, split_pattern, index) {
  vapply(
    strsplit(x, split = split_pattern, perl = TRUE),
    function(parts) if (length(parts) >= index) parts[[index]] else NA_character_,
    character(1)
  )
}

collapse_unique_or_tbd = function(x) {
  x = x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    "TBD"
  } else {
    paste(sort(unique(x)), collapse = ",")
  }
}

build_edges_from_signatures = function(df) {
  specificity_signatures = unique(df$type)

  edge_tables = lapply(specificity_signatures, function(sig) {
    cdr3s = unique(df$Node[df$type == sig])
    if (length(cdr3s) <= 1) {
      return(NULL)
    }
    out = as.data.frame(t(combn(cdr3s, 2, simplify = TRUE)), stringsAsFactors = FALSE)
    colnames(out) = c("from", "to")
    out
  })

  edge_tables = edge_tables[!vapply(edge_tables, is.null, logical(1))]
  if (length(edge_tables) == 0) {
    return(tibble(from = character(), to = character()))
  }

  bind_rows(edge_tables) %>% distinct()
}

build_vertex_table = function(df) {
  df %>%
    group_by(Node) %>%
    mutate(type = paste(sort(unique(type)), collapse = ",")) %>%
    ungroup() %>%
    distinct() %>%
    select(Node, Source, ID, Timepoint, Sample_type, TcRb, V, type, colors)
}

summarize_communities = function(network, communities) {
  community_nodes = split(V(network)$name, communities$membership)
  community_ids = as.integer(names(community_nodes))

  summary_rows = lapply(seq_along(community_nodes), function(i) {
    nodes = community_nodes[[i]]
    comm_id = community_ids[[i]]

    community_specificity = collapse_unique_or_tbd(safe_token_extract(nodes, "[:_-]", 3))

    anbl_nodes = nodes[grepl("^anbl", nodes)]
    patient_nodes = if (length(anbl_nodes) > 0) {
      paste(sort(unique(anbl_nodes)), collapse = ",")
    } else {
      "TBD"
    }

    n_patients = if (length(anbl_nodes) > 0) {
      length(unique(safe_token_extract(anbl_nodes, "[:_]", 2)))
    } else {
      0L
    }

    known_source_regex = "mdavis|anbl|chp|hamm|wherry|mitchel|emerson|profyle"
    external_nodes = grep(known_source_regex, nodes, value = TRUE, invert = TRUE)
    known_external = if (length(external_nodes) > 0) {
      paste(sort(unique(external_nodes)), collapse = ",")
    } else {
      "TBD"
    }

    sample_type_range = if (length(anbl_nodes) > 0) {
      collapse_unique_or_tbd(safe_token_extract(anbl_nodes, "[:_-]", 5))
    } else {
      "TBD"
    }

    timepoint_range = if (length(anbl_nodes) > 0) {
      collapse_unique_or_tbd(safe_token_extract(anbl_nodes, "[:_-]", 4))
    } else {
      "TBD"
    }

    n_timepoints = if (identical(timepoint_range, "TBD")) {
      0L
    } else {
      length(unique(strsplit(timepoint_range, ",", fixed = TRUE)[[1]]))
    }

    tibble(
      Community_id = comm_id,
      Community_size = length(nodes),
      Community_specificity = community_specificity,
      PatientDerived_Nodes = patient_nodes,
      NumberOf_PatientDerived_Nodes = length(unique(anbl_nodes)),
      Number_of_Patients = n_patients,
      KnownExternalTCRs = known_external,
      Sample_Type_Range = sample_type_range,
      Timepoint_Range = timepoint_range,
      Number_of_Timepoints = n_timepoints
    )
  })

  bind_rows(summary_rows) %>%
    arrange(Community_id)
}

add_supernode_attribute = function(graph, community_stats, attr_name, value_col, default = NA) {
  values = rep(default, vcount(graph))
  community_idx = as.integer(community_stats$Community_id)
  in_range = !is.na(community_idx) & community_idx >= 1 & community_idx <= length(values)
  values[community_idx[in_range]] = community_stats[[value_col]][in_range]
  set_vertex_attr(graph, name = attr_name, value = values)
}

# LOAD INPUT ----
gliph_data = readr::read_csv(file.path(gliph_path, gliph_fname), show_col_types = FALSE) %>%
  filter(pattern != "single") %>%
  select(
    type, TcRb, V, pattern, Sample, Fisher_score, number_unique_cdr3, number_subject,
    final_score, vb_score, expansion_score, length_score, cluster_size_score, J, Freq
  ) %>%
  group_by(Sample, TcRb) %>%
  mutate(V = str_replace_all(paste(sort(unique(V)), collapse = "-or"), "(?!^)TRBV", "")) %>%
  ungroup() %>%
  distinct() %>%
  group_by(type) %>%
  filter(any(str_detect(Sample, "anbl"))) %>%
  ungroup() %>%
  mutate(
    ID = safe_token_extract(Sample, ":", 1),
    Source = case_when(
      grepl("anbl", Sample) ~ "anbl",
      grepl("pro", Sample) ~ "profyle_nb",
      grepl("tiger", Sample) ~ "tiger",
      grepl("chp", Sample) ~ "intercept_pc",
      grepl("mdavis", Sample) ~ "mdavid_ac",
      grepl("HammCancer", Sample) ~ "hamm_ac",
      grepl("mitchel", Sample) ~ "mitchel_ph",
      grepl("wherry", Sample) ~ "wherry_pc",
      grepl("EmersonChild", Sample) ~ "emerson_ph",
      grepl("EmersonAdult", Sample) ~ "emerson_ah",
      grepl("HammHealthy", Sample) ~ "hamm_ah",
      TRUE ~ "unknown"
    ),
    Timepoint = ifelse(
      grepl("anbl|chp|pro", Sample),
      safe_token_extract(safe_token_extract(Sample, "_", 2), "-", 1),
      "TBD"
    ),
    Sample_type = ifelse(
      grepl("anbl|chp|pro", Sample),
      safe_token_extract(safe_token_extract(Sample, "_", 2), "-", 2),
      "TBD"
    ),
    Node = str_c(TcRb, "_", Sample)
  ) %>%
  mutate(
    Sample = normalize_tokens(Sample),
    Node = normalize_tokens(Node),
    ID = normalize_tokens(ID)
  ) %>%
  mutate(
    colors = case_when(
      grepl("CMV|DENV|EBV|HCV|HPV|Influenza|MCPyV|HTLV-1|CEF|YFV", Sample) ~ group_col[["Viral"]],
      grepl("S-pneumoniae|M.tuberculosis", Sample) ~ group_col[["Bacterial"]],
      grepl("Human", Sample) ~ group_col[["Human"]],
      Source == "hamm_ah" ~ group_col[["Adult healthy"]],
      Source == "emerson_ah" ~ group_col[["Adult healthy"]],
      Source == "hamm_ac" ~ group_col[["Adult cancer"]],
      Source == "mdavid_ac" ~ group_col[["Adult cancer"]],
      Source == "mitchel_ph" ~ group_col[["Pediatric healthy"]],
      Source == "emerson_ph" ~ group_col[["Pediatric healthy"]],
      Source == "intercept_pc" ~ group_col[["Pediatric cancer"]],
      Source == "wherry_pc" ~ group_col[["Neuroblastoma"]],
      Source == "profyle_nb" ~ group_col[["Neuroblastoma"]],
      Source == "anbl" ~ group_col[["ANBL"]],
      grepl("NBSolid", Sample) ~ group_col[["Neuroblastoma"]],
      TRUE ~ NA_character_
    )
  )

if (run_diagnostics) {
  print(table(gliph_data$Source, useNA = "ifany"))
  print(table(gliph_data$Timepoint, useNA = "ifany"))
  print(table(gliph_data$Sample_type, useNA = "ifany"))
}

# RAW NETWORK GENERATION ----
vertex_data = build_vertex_table(gliph_data)
if (run_diagnostics) {
  janitor::get_dupes(vertex_data, Node)
}

EDGES = build_edges_from_signatures(gliph_data)
if (nrow(EDGES) == 0) {
  stop("No edges were generated from specificity signatures; cannot build graph.")
}

network = graph_from_data_frame(
  d = EDGES,
  vertices = vertex_data,
  directed = FALSE
)
network = simplify(network, remove.multiple = TRUE, remove.loops = TRUE)

# NETWORK TRIMMING VIA CLIQUE IDENTIFICATION ----
MaxClique_collection = max_cliques(network, min = 4, max = NULL)
if (length(MaxClique_collection) == 0) {
  stop("No cliques of size >= 4 were found; cannot continue with trimmed network.")
}
network = induced_subgraph(network, vids = sort(unique(unlist(MaxClique_collection))))
network = simplify(network, remove.multiple = TRUE, remove.loops = TRUE)

# POST-FILTRATION NETWORK STATS ----
components_out = igraph::components(network)

component_sizes = tibble(
  component_id = names(components_out$csize),
  component_size = as.integer(components_out$csize)
)

component_membership = tibble(
  Node = names(components_out$membership),
  component_id = as.integer(components_out$membership)
)

if (run_diagnostics) {
  print(igraph::count_components(network))
  print(range(components_out$csize))
  print(summary(components_out$csize))
  print(length(V(network)[V(network)$Source == "anbl"]))
  print(any_loop(network))
  print(any_multiple(network))
}

# COMMUNITY DETECTION ----
communities = cluster_leiden(
  network,
  objective_function = "CPM",
  resolution = 0.9,
  n_iterations = 1000
)

network = set_vertex_attr(
  network,
  name = "EdgeBetwCommunity",
  value = communities$membership
)

GLIPHII_Community_stats = summarize_communities(network, communities) %>%
  mutate(
    Community_specificity = str_replace_all(
      Community_specificity,
      c(",V\\d" = ",Wherry", "V\\d" = "Wherry")
    )
  ) %>%
  mutate(
    Community_specificity = str_replace_all(
      Community_specificity,
      c("(Wherry)(?=.*\\1)" = "")
    )
  )

# SPECIFICITY ANNOTATION ----
viral_terms = c("CEF", "Influenza", "MCPyV", "CMV", "DENV", "EBV", "HCV", "HPV", "HTLV-1", "YFV")
bacterial_terms = c("S.pneumoniae", "M.tuberculosis")
pediatric_healthy_terms = c("Normal", "MitchelChild", "EmersonChild", "Wherry")
pediatric_cancer_terms = c("Leukemia", "Lymphoma", "T.cellmalignancies", "Solidtumors")
pediatric_terms = c(pediatric_healthy_terms, pediatric_cancer_terms)
adult_healthy_terms = c("EmersonAdult", "HammHealthy")
adult_cancer_terms = c("MDavis", "HammCancer")
adult_terms = c(adult_healthy_terms, adult_cancer_terms)
neuroblastoma_terms = c("NB", "NBSolidtumors")

esc = function(x) str_replace_all(x, "([.\\^$|()?*+\\[\\]{}\\\\-])", "\\\\\\1")
pat = function(terms) paste0("(^|,)(?:", paste0(esc(terms), collapse = "|"), ")(?=,|$)")

specificity_annotation = tibble(
  Community_specificity = unique(GLIPHII_Community_stats$Community_specificity)
) %>%
  mutate(
    has_human = str_detect(Community_specificity, "(^|,)Human(?=,|$)"),
    has_viral = str_detect(Community_specificity, pat(viral_terms)),
    has_bact = str_detect(Community_specificity, pat(bacterial_terms)),
    has_ped_healthy = str_detect(Community_specificity, pat(pediatric_healthy_terms)),
    has_ped_cancer = str_detect(Community_specificity, pat(pediatric_cancer_terms)),
    has_pediatric = has_ped_healthy | has_ped_cancer,
    has_adult_healthy = str_detect(Community_specificity, pat(adult_healthy_terms)),
    has_adult_cancer = str_detect(Community_specificity, pat(adult_cancer_terms)),
    has_adult = has_adult_healthy | has_adult_cancer,
    has_nb = str_detect(Community_specificity, pat(neuroblastoma_terms)),
    none_pathogen_human = !has_viral & !has_bact & !has_human
  ) %>%
  mutate(
    Abstract_Annotation = case_when(
      (has_bact & has_human) | (has_viral & has_human) ~ "Cross-reactive",
      has_human & !has_viral & !has_bact ~ "Human",
      has_bact & has_viral & !has_human ~ "Viral and bacterial",
      has_viral & !has_bact & !has_human ~ "Viral",
      has_bact & !has_viral & !has_human ~ "Bacterial",
      has_ped_cancer & !has_ped_healthy & none_pathogen_human & !has_adult ~ "Pediatric cancer",
      has_adult_cancer & !has_adult_healthy & none_pathogen_human & !has_pediatric ~ "Adult cancer",
      has_adult_cancer & has_ped_cancer & !has_adult_healthy &
        none_pathogen_human & !has_ped_healthy ~ "Adult and pediatric cancer",
      has_pediatric & none_pathogen_human & has_adult ~ "Adult and pediatric",
      has_pediatric & none_pathogen_human & !has_adult ~ "Pediatric",
      !has_pediatric & none_pathogen_human & has_adult ~ "Adult",
      has_nb & none_pathogen_human & !has_adult & !has_pediatric ~ "Neuroblastoma",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    Color = case_when(
      Abstract_Annotation == "Viral" ~ "#001959",
      Abstract_Annotation == "Viral and bacterial" ~ "#0C325E",
      Abstract_Annotation == "Bacterial" ~ "#124560",
      Abstract_Annotation == "Cross-reactive" ~ "#185461",
      Abstract_Annotation == "Human" ~ "#26635F",
      Abstract_Annotation == "Adult healthy" ~ "#3B6C55",
      Abstract_Annotation == "Adult and pediatric healthy" ~ "#577646",
      Abstract_Annotation == "Pediatric healthy" ~ "#727E38",
      Abstract_Annotation == "Adult" ~ "#91852C",
      Abstract_Annotation == "Adult and pediatric" ~ "#B28D2E",
      Abstract_Annotation == "Pediatric" ~ "#D29243",
      Abstract_Annotation == "Adult cancer" ~ "#EC9961",
      Abstract_Annotation == "Adult and pediatric cancer" ~ "#FAA588",
      Abstract_Annotation == "Pediatric cancer" ~ "#FDB2AE",
      Abstract_Annotation == "Neuroblastoma" ~ "#FCBFD2",
      Abstract_Annotation == "ANBL" ~ "#F9CCF9",
      TRUE ~ NA_character_
    )
  )

GLIPHII_Community_stats = left_join(
  GLIPHII_Community_stats,
  specificity_annotation,
  by = "Community_specificity"
)

# NODE FEATURES ----
anbl_nodes = V(network)$name[grepl("^anbl", V(network)$name)]
node_degree = degree(network)
node_transitivity = igraph::transitivity(network, type = "local")

Node_Features = tibble(
  Node = anbl_nodes,
  Node_degree = as.numeric(node_degree[anbl_nodes]),
  Node_transitivity = as.numeric(node_transitivity[anbl_nodes])
)

# SUPER-NODE NETWORK ----
Abstract_Network = igraph::contract.vertices(
  graph = network,
  mapping = V(network)$EdgeBetwCommunity
)
Abstract_Network = igraph::simplify(Abstract_Network, remove.multiple = TRUE, remove.loops = TRUE)

Abstract_Network = add_supernode_attribute(
  Abstract_Network,
  GLIPHII_Community_stats,
  "EdgeBetwCommunity_Size",
  "Community_size",
  default = NA_real_
)
Abstract_Network = add_supernode_attribute(
  Abstract_Network,
  GLIPHII_Community_stats,
  "EdgeBetwCommunity_Specificity",
  "Community_specificity",
  default = NA_character_
)
Abstract_Network = add_supernode_attribute(
  Abstract_Network,
  GLIPHII_Community_stats,
  "EdgeBetwCommunity_AbstractSpecificityAnnotation",
  "Abstract_Annotation",
  default = NA_character_
)
Abstract_Network = add_supernode_attribute(
  Abstract_Network,
  GLIPHII_Community_stats,
  "SuperNode_Color",
  "Color",
  default = NA_character_
)
Abstract_Network = add_supernode_attribute(
  Abstract_Network,
  GLIPHII_Community_stats,
  "Number_of_Patients",
  "Number_of_Patients",
  default = NA_real_
)
Abstract_Network = add_supernode_attribute(
  Abstract_Network,
  GLIPHII_Community_stats,
  "Number_of_PatientDerivedTCRs",
  "NumberOf_PatientDerived_Nodes",
  default = NA_real_
)
Abstract_Network = add_supernode_attribute(
  Abstract_Network,
  GLIPHII_Community_stats,
  "Community_Sample_Type",
  "Sample_Type_Range",
  default = NA_character_
)
Abstract_Network = add_supernode_attribute(
  Abstract_Network,
  GLIPHII_Community_stats,
  "Community_Timepoints",
  "Timepoint_Range",
  default = NA_character_
)
Abstract_Network = add_supernode_attribute(
  Abstract_Network,
  GLIPHII_Community_stats,
  "NumberOfTimepoints",
  "Number_of_Timepoints",
  default = NA_real_
)

# SAVE FILES ----
output_dir = file.path(gliph_path, output_path)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(Abstract_Network, file = file.path(output_dir, "Abstract_Network.rds"))
saveRDS(communities, file = file.path(output_dir, "communities.rds"))
saveRDS(component_sizes, file = file.path(output_dir, "component_sizes.rds"))
saveRDS(component_membership, file = file.path(output_dir, "component_membership.rds"))
saveRDS(EDGES, file = file.path(output_dir, "EDGES.rds"))
saveRDS(gliph_data, file = file.path(output_dir, "gliph_data.rds"))
saveRDS(GLIPHII_Community_stats, file = file.path(output_dir, "GLIPHII_Community_stats.rds"))
saveRDS(network, file = file.path(output_dir, "network.rds"))
saveRDS(Node_Features, file = file.path(output_dir, "Node_Features.rds"))

# OPTIONAL: EXAMPLE COMMUNITY PLOT + SEQUENCE LOGO ----
if (run_example_plots) {
  if (example_community_id %in% V(network)$EdgeBetwCommunity) {
    nt = igraph::subgraph(network, which(V(network)$EdgeBetwCommunity %in% c(example_community_id)))

    pdf(
      file.path(output_dir, paste0("community_", example_community_id, ".pdf")),
      width = 6, height = 6, family = "Helvetica"
    )
    plot.igraph(
      nt,
      vertex.size = 25,
      vertex.color = case_when(
        V(nt)$Source == "anbl" ~ "#D4B9DA",
        V(nt)$Source == "profyle_nb" ~ "#DF65B0",
        V(nt)$Source == "intercept_pc" ~ "#DF65B0",
        TRUE ~ "#CE1256"
      ),
      vertex.frame.color = "transparent",
      vertex.shape = "circle",
      vertex.cex = 0.75,
      vertex.label = NA,
      edge.width = 1,
      edge.color = "#B8B8B8",
      rescale = TRUE
    )
    dev.off()

    svg(
      filename = file.path(output_dir, paste0("seqlogo_", example_community_id, ".svg")),
      family = "Helvetica",
      width = 2.0, height = 1.0,
      onefile = FALSE,
      bg = "transparent"
    )
    ggplot() +
      geom_logo(
        V(nt)$TcRb,
        method = "prob",
        seq_type = "aa",
        font = "akrobat_regular",
        col_scheme = make_col_scheme(chars = LETTERS, cols = rep("#000000", length(LETTERS)))
      ) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
      )
    dev.off()
  } else {
    warning("example_community_id not found in this run; no example plots generated.")
  }
}
