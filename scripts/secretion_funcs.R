#! /usr/bin/env R

# Read list of WY candidates - for Section 3.1, Step 3
read_wy_list <- function(x) {
  read_tsv(x, show_col_types = FALSE,
           col_names = c("ID"), id = "fpath") %>%
    mutate("isolate" = str_extract(fpath, "(?<=//).*(?=\\.orfs)")) %>%
  select(-fpath)
}

# Runs and parses regex search for RXLR-EER domains - Section 3.1, Step 4
effectr_eer_sp3 <- function(seqs_file, sp3_file, isolate_name = NULL) {
  # Read in seqs to scan 
  orfs <- read.fasta(seqs_file)
  # Read signalp3 output from all possible candidates
  sp3 <- read_tsv(sp3_file, show_col_types = FALSE, id = "fpath") %>%
    mutate("isolate" = str_extract(fpath, ".*(?=\\.orfs)")) %>%
    mutate(isolate = str_remove(isolate, ".*data/*")) %>%
    filter(HMM_Sprob_score >= 0.9) %>%
    rename("ID" = `#ID`)
  if(is.null(isolate_name)) {
    isolate_name <- pull(sp3, isolate)[1]
  }
  
  # Search new regexes
  re_win <- regex.search(orfs, motif = "Win2007")
  re_whis <- regex.search(orfs, motif = "Whisson2007")
  re_whis_rxlr <- regex.search(orfs, motif = "Whisson2007_rxlr")
  re_whis_eer <- regex.search(orfs, motif = "Whisson2007_eer")
  
  # Get all names passing regex
  re_win_sp_names <- data.frame("method" = "Win2007",
                                "ID" = names(re_win)[names(re_win) %in% sp3$ID])
  re_whis_sp_names <- data.frame("method" = "Whis2007",
                                 "ID" = names(re_whis)[names(re_whis) %in% sp3$ID])
  re_whis_rxlr_sp_names <- data.frame("method" = "Whis_rxlr",
                                      "ID" = names(re_whis_rxlr)[names(re_whis_rxlr) %in% sp3$ID])
  re_whis_eer_sp_names <- data.frame("method" = "Whis_eer",
                                 "ID" = names(re_whis_eer)[names(re_whis_eer) %in% sp3$ID])
  
  # Compile numbers of hits before and after filtering for signalp3 threshold
  n_hits_df <- tibble("regex_method" = c("secreted_ORFs", "Win2007", "Whisson2007", "Whisson_RXLR", "Whisson_EER"),
                      "n_hits_tot" = c(NA,
                                       length(re_win),
                                       length(re_whis),
                                       length(re_whis_rxlr),
                                       length(re_whis_eer)),
                      "n_hits_sp3" = c(nrow(sp3),
                                       nrow(re_win_sp_names),
                                       nrow(re_whis_sp_names),
                                       nrow(re_whis_rxlr_sp_names),
                                       nrow(re_whis_eer_sp_names)),
                      "isolate" = rep(isolate_name, 5)
  )
  
  print(n_hits_df)
  # Choose one list of names to return from function
  # Probably EER but I don't think my summary func is working for that motif yet
  all_names_hits <- bind_rows(
    re_whis_sp_names, re_win_sp_names, re_whis_rxlr_sp_names, re_whis_eer_sp_names) %>%
    mutate("isolate" = isolate_name)
  return(all_names_hits)
}

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
    # Only removing if past middle of helix looks perfect for first 10 with Pred_hel > 0 I see
    # Remove if the signal peptide is earlier than the middle of the helix
    mutate(mean_helix_pos = mean(c(Start, Stop), na.rm = TRUE)) %>%
    filter(NN_Ymax_pos > mean_helix_pos | is.na(mean_helix_pos))
    
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
## In the end I want a list of proteins that pass the TM test,
## and I can remove any proteins from rxlr hit list that aren't present
## This means ALL rows are important, especially if they have NO tms pred
#pr102v3_sp3_tm <- join_sp_to_tmhmm(pr102v3_tm_raw, pr102v3_sp3, is.secreted = FALSE)
#pr102v3_notm <- prune_tm_by_sp(pr102v3_sp3_tm)

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


#' Combine SignalP and TMHMM results
#'
#' Further, filter results on whether or not proteins are secreted.
#' This means right now, TMHMM results must be more comprehensive
#' than SignalP results.
#'
#' @param tmhmm Tibble from read_tmhmm()
#' @param sp Tibble from read_sp3(), sp5 not yet supported.
#' @param is.secreted Also filter the results down to only secretions.
#'
#' @seealso \code{/link{read_tmhmm}}
#'
#' @export
join_sp_to_tmhmm <- function(tmhmm, sp, is.secreted = TRUE) {
  tmhmm_sp <- tmhmm %>%
    pivot_tmhmm() %>%
    left_join(sp, by = "Protein")

  if(is.secreted) {
    tmhmm_sp <- filter(tmhmm_sp, Prob > 0.5)
  }

  return(tmhmm_sp)
}


#' Filter proteins on TMHMM results
#'
#' Removes proteins from joined TMHMM-SignalP data
#'
#' @param x Tibble from internally joined data
#'
#' @sealso \code{/link{join_sp_to_tmhmm}}
#'
#' @export
prune_tm <- function(x) {
  # instead of using filter(any()), checking all TM-helices per prot,
  # just check the stop pos of the most C-terminal helix.
  no_tm_prots <- x %>%
    group_by(Protein) %>%
    mutate(Cterm_TM = max(Stop, na.rm = TRUE)) %>%
    filter(Cterm_TM < 70 | is.na(Cterm_TM)) %>%
    distinct(Protein, .keep_all = TRUE) %>%
    mutate(Cterm_TM = na_if(Cterm_TM, "-Inf"))
}


#' Process TMHMM-2.0c output
#'
#' Read in a file from TMHMM-2.0c.
#' This function simply reads in the file and adds logical column names
#' @param file gff output file from SignalP-5
#' @param ... additional arguments passed to read_tsv
#' @keywords tmhmm
#' @return A tibble
#' @seealso \code{/link{pivot_tmhmm}} for conversion to "long" format.
#'
#' @export
read_tmhmm <- function(file, cols=c("Protein", "Length", "Exp_AA", "First_60", "Pred_hel", "Topology"), ...) {
  tmhmm_df <- readr::read_tsv(file, col_names = cols, ...)

  # clean data of text repeated in every row
  tmhmm_clean <- tmhmm_df %>%
    tidyr::separate(Length, c(NA, "Length"), convert = TRUE) %>%
    # merge extra - includes fractions of Exp_AA
    tidyr::separate(Exp_AA, c(NA, "Exp_AA"), convert = TRUE, extra = "merge") %>%
    # merge extra - includes fractions of First_60
    tidyr::separate(First_60, c(NA, "First_60"), convert = TRUE, extra = "merge") %>%
    tidyr::separate(Pred_hel, c(NA, "Pred_hel"), convert = TRUE) %>%
    # if you don't merge the extra pieces of topology, dplyr only keeps the first feat
    tidyr::separate(Topology, c(NA, "Topology"), convert = TRUE, extra = "merge")

  return(tmhmm_clean)
}


#' Tidy TMHMM-2.0c data
#'
#' After reading TMHMM-2.0c output into R, it's still in an unfriendly-format.
#' This tidies the "topology" string in the last column of the output
#' @param x tab-separated TMHMM-2.0c output
#' @keywords tmhmm
#' @return A tibble
#' @seealso \code{/link{read_tmhmm}} for reading input with proper names.
#'
#' @export
#' @examples
#' tmhmm_test <- tibble(Protein = c("A", "B", "C", "D", "E"),
#'                      Info = c(1, 2, 3, 4, 5),
#'                      Topology = c("o", "i", "i4-9",
#'                                   "o10-20i30-50",
#'                                   "i105-205o305-405i505-605"))
#'
#' pivot_tmhmm(tmhmm_test)
#'
#' # This can now be used in ggplot
#' # e.g. How much residue length is TM-helix folds?
#'
#' pivot_tmhmm(tmhmm_test) %>%
#'   mutate(local_in = Stop-Start+1) %>%
#'   group_by(Protein) %>%
#'   mutate(Total_in = sum(local_in, na.rm = TRUE)) %>%
#'   distinct(Protein, .keep_all = TRUE) %>%
#'   select(-Local_in) %>%
#'   ggplot2::ggplot(.) +
#'   ggplot2::geom_histogram(aes(Total_in))
pivot_tmhmm <- function(x) {
  # tidy raw tmhmm output (x) into long format (returns)
  # This was immensely helpful
  # https://www.r-bloggers.com/strsplit-but-keeping-the-delimiter/
  # I thought I needed to use regex for one tmhmm "term",
  # but that didnt end up working anyway
  # topo_regex <- "[io][[:digit:]]+-[[:digit:]]+"
  # solution aws a regex lookahead

  # first get each "term" on its own line
  # then separate each term into columns
  # c(Topology, c("Localization", "Start", "Stop"))
  # Topology is "o" or "i"
  # Start and stop are separated by hyphen

  # I want to split terms into their own line, but avoid the first character
  # so maybe check if i or o AND if previous char was a digit
  tidy_tmhmm <- x %>%
    tidytext::unnest_tokens(Topo_long, Topology,
                            token = stringr::str_split,
                            pattern = "(?=[[:alpha:]])" ) %>%
    dplyr::filter(Topo_long != "") %>%
    tidyr::separate(Topo_long, c("Localization", "pos"), "(?<=[[:alpha:]])") %>%
    tidyr::separate(pos, c("Start", "Stop"), sep = "-",
                    convert = TRUE, fill = "left")

  return(tidy_tmhmm)
}
