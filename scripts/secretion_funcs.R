#! /usr/bin/env R

# From dataframe with Whisson method, separates into HMM and regex
sep_whisson <- function(df) {
  if(sum(grepl("Whisson", df$method)) == 0) {
    print("Whisson method not used to find RXLRs")
    return(df)
  }
  df_sepwhisson <- df %>%
    filter(grepl("Whisson", method)) %>%
    mutate(is_rxlr = str_replace(is_rxlr, "neither", "N")) %>%
    mutate(Whisson2007_hmm = str_replace(is_rxlr, "hmm", "Y")) %>%
    mutate(Whisson2007_hmm = str_replace(Whisson2007_hmm, "re", "N")) %>%
    mutate(Whisson2007_re = str_replace(is_rxlr, "re", "Y")) %>%
    mutate(Whisson2007_re = str_replace(Whisson2007_re, "hmm", "N")) %>%
    select(-c(is_rxlr, method)) %>%
    pivot_longer(-ID, names_to = "method", values_to = "is_rxlr") %>%
    bind_rows(df) %>%
    arrange(method, ID) %>%
    filter(method != "Whisson2007")
  return(df_sepwhisson)
}

prune_tm_by_sp <- function(x, check_signalp = FALSE) {
  # instead of using filter(any()), checking all TM-helices per prot,
  # just check the stop pos of the most C-terminal helix.
  no_tm_prots <- x %>%
    # To be more explicit here, change NA in start,stop to -inf 
    group_by(Protein) %>%
    # This does keep the NAs but likely unintended behavior,
    # other solutions I found converted lots of values to NA
    slice_max(Stop, with_ties = FALSE) %>%
    # this is why originally I just chose 70 as a hard-limit
    # max(pr102v3_secreted$NN_Ymax_pos)
    # For definitely secreted proteins ymax goes to 66
    # Only removing if past middle of helix looks perfect for first 10 with Pred_hel > 0 I see
    # filter(pr102v3_secreted, HMM_Sprob_score > 0.9)
    # %>% arrange(desc(NN_Ymax_pos)) %>% filter(Pred_hel > 0) %>%
    # mutate(mid_helix = mean(Start, Stop))
    # Remove if the signal peptide is earlier than the middle of the helix
    mutate(mean_helix_pos = mean(c(Start, Stop), na.rm = TRUE)) %>%
    filter(NN_Ymax_pos > mean_helix_pos | is.na(mean_helix_pos))
    # mutate(Cterm_TM = max(Stop, na.rm = TRUE)) %>%
    # filter(Cterm_TM < 70 | is.na(Cterm_TM)) %>%
    # distinct(Protein, .keep_all = TRUE) %>%
    # mutate(Cterm_TM = na_if(Cterm_TM, "-Inf"))
    
    if(check_signalp){
      no_tm_prots <- filter(no_tm_prots, HMM_Sprob_score >= 0.9)
    }
    
    return(no_tm_prots)
}

# Examples - Not run
#sp3_cols_all <- c("ID","NN_Cmax_score","NN_Cmax_pos","NN_Cmax_pred","NN_Ymax_score","NN_Ymax_pos","NN_Ymax_pred","NN_Smax_score","NN_Smax_pos","NN_Smax_pred","NN_Smean_score","NN_Smean_pred","NN_D_score","NN_D_pred","HMM_type","HMM_Cmax_score","HMM_Cmax_pos","HMM_Cmax_pred","HMM_Sprob_score","HMM_Sprob_pred")
#pr102v3_sp3 <- readr::read_tsv(sp3_f, col_names = sp3_cols_all, skip = 1)
#pr102v3_sp3 <- select(pr102v3_sp3, "Protein" = ID, HMM_Sprob_score, NN_Ymax_pos)
#pr102v3_tm_raw <- read_tmhmm(tmhmm_f)
## Trust the filtering from earlier signalp runs, just join the dfs
## In the end I want a list of proteins that pass the TM test,
## and I can remove any proteins from rxlr hit list that aren't present
## This means ALL rows are important, especially if they have NO tms pred
#pr102v3_sp3_tm <- iSecrete::join_sp_to_tmhmm(pr102v3_tm_raw, pr102v3_sp3, is.secreted = FALSE)
#pr102v3_notm <- prune_tm_by_sp(pr102v3_sp3_tm)

# Pasted from Pram_effector_analysis.Rmd
# Where you filter out any prot where the signal pep is succeeded by
# at least one average tm position (between start and stop)
read_prune_tms <- function(sp3_f, tmhmm_f, isolate = NULL, quiet = FALSE, check_sigp = FALSE) {
  sp3_cols_all <- c("ID","NN_Cmax_score","NN_Cmax_pos","NN_Cmax_pred","NN_Ymax_score","NN_Ymax_pos","NN_Ymax_pred","NN_Smax_score","NN_Smax_pos","NN_Smax_pred","NN_Smean_score","NN_Smean_pred","NN_D_score","NN_D_pred","HMM_type","HMM_Cmax_score","HMM_Cmax_pos","HMM_Cmax_pred","HMM_Sprob_score","HMM_Sprob_pred")
  
  sp3_cands <- readr::read_tsv(sp3_f, col_names = sp3_cols_all, skip = 1,
                               show_col_types = FALSE) %>%
    # For signal peptide length use same as rxlr_motifs.py,
    # @line 250 comment defines it as NN_Ymax_pos-1 
    # pr102v3_sp3 <- select(pr102v3_sp3, "Protein" = ID, HMM_Sprob_score, NN_Ymax_pos)
    select("Protein" = ID, HMM_Sprob_score, NN_Ymax_pos)
  
  tm_raw <- iSecrete::read_tmhmm(tmhmm_f, show_col_types = FALSE)
  sp3_tm_cands <- iSecrete::join_sp_to_tmhmm(tm_raw, sp3_cands, is.secreted = FALSE)
  notm_cands <- prune_tm_by_sp(sp3_tm_cands, check_signalp = check_sigp)
  if(quiet == FALSE) {
    if(!is.null(isolate)) {
      print(isolate)
      notm_cands$isolate <- isolate
    }
    print(paste("Candidates tested in TMHMM:", nrow(tm_raw)))
    print(paste("Candidates remaining after filtering by TM (and SigP if check_sigp used):", nrow(notm_cands)))
  }
  
  return(notm_cands)
}
# Example (not run):
# notm_orfs1 <- read_prune_tms(
#  sp3_f = "annotation/effectors/output_data/PR-102_v3.1.orfs-min70long.start2stop.rxlr_effp3.outsp3_tabular.tmp",
#  tmhmm_f = "annotation/effectors/output_data/PR-102_v3.1.orfs-min70long.start2stop.rxlr_whisson_hmm_win_y.tmhmm.out")


# Read results in from rxlr_motif.py
read_rxlrmotifpy_out <- function(file) {
  x <- read_tsv(file) %>%
    mutate("method" = colnames(.)[2]) %>%
    select("ID" = `#ID`, "is_rxlr" = 2, 3)
  return(x)
}
# Examples (not run): 
# read_rxlrmotifpy_out(methods_out_files[1])
# rxlr_cands_df <- methods_out_files %>%
#   map_dfr(read_rxlrmotifpy_out)

# From dataframe with effectorp3 method, separate into cyto and apoplastic
sep_effp3 <- function(df, tm_df = NULL) {
  if(sum(grepl("effectorp3", df$method)) == 0) {
    print("EffectorP3 not used to find RXLRs")
    return(df)
  }
  df_sepeffp <- df %>%
    filter(grepl("effectorp3", method)) %>%
    # filter(!grepl("apo|N", is_rxlr)) %>%
    # Y to effp_symp and both to effp_symp_apo
    mutate(EffectorP_apo = str_replace(is_rxlr, "Y", "N")) %>%
    mutate(EffectorP_apo = str_replace(EffectorP_apo, "apo", "Y")) %>%
    mutate(EffectorP_apo = str_replace(EffectorP_apo, "both", "Y")) %>%
    mutate(EffectorP_sym = str_replace(is_rxlr, "apo", "N")) %>%
    mutate(EffectorP_sym = str_replace(EffectorP_sym, "both", "Y")) %>%
    select(-c(is_rxlr, method)) %>%
    pivot_longer(-ID, names_to = "method", values_to = "is_rxlr")
  if(!is.null(tm_df)) {
    df_sepeffp <- filter(df_sepeffp, ID %in% tm_df$Protein)
  }
  # {if !is.null(tm_df) filter(., ID %in% tm_df$Protein)}
  df_sepeffp <- df_sepeffp %>%
    bind_rows(df) %>%
    arrange(method, ID) %>%
    filter(method != "effectorp3")
  return(df_sepeffp)
}
# rxlr_sepeffp <- sep_effp3(rxlr_cands_df)
# rxlr_sepeffp
# rxlr_df2venn(rxlr_sepeffp, circle_order = c(1,2,3))

# Go from wide/pivoted df to venn diagram
rxlr_df2venn <- function(df, draw_venn = TRUE, circle_order = c(2,3,1,4),
                         abbr_methods = FALSE) {
  # https://stackoverflow.com/a/58837832/9120324
  df_venn <- df %>%
    filter(is_rxlr == "Y") %>%
    select(-is_rxlr) %>%
    group_by(method) %>%
    mutate(row = row_number())
  if(abbr_methods) {
    df_venn <- df_venn %>%
      mutate(method = str_replace(method, "200", "0")) %>%
      mutate(method = str_replace(method, "Bhattacharjee", "Bhat")) %>%
      mutate(method = str_replace(method, "Whisson", "Whis")) %>%
      mutate(method = str_replace(method, "EffectorP", "EffP"))
  }
  df_venn <- df_venn %>%
    pivot_wider(values_from = ID, names_from = method) %>%
    select(-row)
  if(draw_venn) {
    n_circ_diff = length(circle_order) - ncol(df_venn)
    if(n_circ_diff > 0) {
      print("Too many circles specified, remove this many from order:")
      return(abs(n_circ_diff))
    } else if(n_circ_diff < 0) {
      print("Not enough circles specified, please add this many:")
      return(n_circ_diff)
    } else {
      # Convert to list to remove NAs which inflate #s with false overlaps
      list_venn <- as.list(df_venn)
      list_venn <- map(list_venn, function(x) x[!is.na(x)])
      list_venn[circle_order] %>%
        # matches circle order in supp fig s1
        ggVennDiagram(.) +
        scale_fill_viridis_c(direction=-1) +
        theme(legend.position = "none")
    }}
}
