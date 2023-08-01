#!/usr/bin/env Rscript

# options(conflict.policy = list("can.mask" = c("base", "stats", "dplyr")))
require("optparse", quietly = TRUE)
require("dplyr", quietly = TRUE, warn.conflicts = FALSE)
# require("dplyr", quietly = TRUE)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Required. List of sequences from getorf (EMBOSS)"),
  make_option(c("-s", "--suffix"), type="character", default=".list",
              help="Suffix of filename passed to input [default= %default]")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# print(opt$input)
# print(opt$suffix)

# Get a GFF for sequences obtained from EMBOSS' getorf.
# This tool assumes genes have no introns.
# Also the syntax is somewhat specific to my braker gff,
# really just the naming convention for -exons and -cds though

if(is.null(opt$input)) {
  stop("Input list of contig names must be provided, see usage (--help)")
}

#seqs_list_f <- "annotation/effectors/output_data/PR-102_v3.1.orfs-min90long.start2stop.list"
seqs_list_f <- opt$input
#seqs_list_suffix <- ".list"
seqs_list_suffix <- opt$suffix

# get output filename
seqs_list_dir <- dirname(seqs_list_f)
seqs_list_bname <- stringr::str_remove(basename(seqs_list_f), seqs_list_suffix)
out_gff_f <- paste(file.path(seqs_list_dir, seqs_list_bname), "gff", sep = ".")
out_gff_f_orf <- paste(file.path(seqs_list_dir, seqs_list_bname), "orfs.gff", sep = ".")

# Use a prefix to rename orfs, breaks link between orf seqs and new gff
# prefix <- "PHRA102"

seqs_list_colnames <- c("orf", "start", "to", "stop", "is_rev1", "is_rev2")
seqs_list_raw <- readr::read_delim(seqs_list_f, delim = " ",
                               col_names = seqs_list_colnames) %>%
  # remove redundant cols
  select(-to) %>%
  select(-any_of(c("is_rev1", "is_rev2"))) %>%
  # remove non-numeric chars from seqnames
  mutate(across(c(start, stop), ~ stringr::str_remove(.x, "\\[|\\]"))) %>%
  mutate(across(c(start, stop), as.numeric)) %>%
  # commonly seq names lists start with ">" as artifact from fasta, remove it
  mutate(orf = stringr::str_remove(orf, "^>"))
  # handle fwd/reverse separately, then merge
  # select(-is_rev2)

# get reverse seqs, switch start and stop
seqs_list_rev <- seqs_list_raw %>%
  filter(start > stop) %>%
  mutate(end = start) %>%
  mutate(start = stop) %>%
  select(-stop) %>%
  mutate(strand = "-") %>%
  # adjust end pos to include stop codon
  # These are reverse, so START is actually the end.
  mutate(last_peptide = start) %>%
  mutate(start = last_peptide-3)

# fwd seqs and merge with rev seqs
seqs_list <- seqs_list_raw %>%
  filter(start < stop) %>%
  # adjust end pos to include stop codon
  mutate(end = stop+3) %>%
  select(-stop) %>%
  mutate(strand = "+") %>%
  bind_rows(seqs_list_rev)

# start with gff columns common for all
# The end pos of getorf is one codon before stop codon,
# need to adjust that so child stop_codon are within coords of parent gene.
# I visually confirmed in IGV that these windows dont include stop codons
# 2021/06/11 Apparently Reverse don't end by stop codon, -3 to their chrStart
# OHHH, I put them out of phase. The end of fwd is actually start of rev.
seqs_list_gff_common <- seqs_list %>%
  # getorf uses _[[:digit:]] to add IDs, find that 1+ times before end of string
  tidyr::separate(orf, c("seqid", "orf_num"), "_(?=[[:digit:]]+$)",
           remove = FALSE, convert = TRUE) %>%
  mutate(source = "getorf") %>%
  mutate(score = ".") #%>%
  # adjust end pos to include stop codon
  # mutate(last_peptide = end) %>%
  # mutate(end = last_peptide+3)

# add each feat with unique requirements
# ID follows the seqname from getorf, otherwise gff & seqs don't match up,
# so the cleaner naming scheme I tried to make wont work.
# It definitely doesnt work once you realize orf_num restarts every scaffold.
# order_ID will help sort into the right order after we combine
seqs_list_gff_gene <- seqs_list_gff_common %>%
  mutate(type = "gene") %>%
  mutate(phase = ".") %>%
  # mutate(ID_gene = paste("ID=", prefix, "orf", orf_num, sep = "")) %>%
  mutate(ID_gene = paste("ID=", orf, sep = "")) %>%
  mutate(attributes = ID_gene) %>%
  mutate(order_ID = seq(1, n()))
  # select(seqid, source, type, start, end, score, strand, phase, attributes)

# mRNA is child of gene
seqs_list_gff_mrna <- seqs_list_gff_gene %>%
  mutate(type = "mRNA") %>%
  # mutate(phase = ".") %>%
  
  # mutate(ID_mrna = paste("ID=", prefix, "orf", orf_num, ".1", sep = "")) %>%
  mutate(ID_mrna = paste("ID=", orf, ".1", sep = "")) %>%
  # mutate(parent = paste("Parent=", prefix, "orf", orf_num, sep = "")) %>%
  mutate(parent = gsub(pattern = "ID=", replacement = "Parent=", ID_gene)) %>%
  mutate(attributes = paste(ID_mrna, parent, sep = ";")) %>%
  mutate(order_ID = seq(1.1, n()+1))
  # select(seqid, source, type, start, end, score, strand, phase, attributes)

# exon and CDS are child of the mRNA
seqs_list_gff_exon <- seqs_list_gff_mrna %>%
  mutate(type = "exon") %>%
  # mutate(phase = ".") %>%
  # mutate(ID=paste("ID=", prefix, "orf", orf_num, ".1", "-exon1", sep = "")) %>%
  mutate(ID_exon=paste(ID_mrna, "-exon1", sep = "")) %>%
  # mutate(parent = paste("Parent=", prefix, "orf", orf_num, ".1", sep = "")) %>%
  mutate(parent = gsub(pattern = "ID=", replacement = "Parent=", ID_mrna)) %>%
  mutate(attributes = paste(ID_exon, parent, sep = ";")) %>%
  mutate(order_ID = seq(1.2, n()+1))
  # select(seqid, source, type, start, end, score, strand, phase, attributes)

seqs_list_gff_cds <- seqs_list_gff_mrna %>%
  mutate(type = "CDS") %>%
  mutate(phase = "0") %>%
  # mutate(ID = paste("ID=", prefix, "orf", orf_num, ".1", "-cds2", sep = "")) %>%
  mutate(ID_cds = paste(ID_mrna, "-cds2", sep = "")) %>%
  # mutate(parent = paste("Parent=", prefix, "orf", orf_num, ".1", sep = "")) %>%
  mutate(parent = gsub(pattern = "ID=", replacement = "Parent=", ID_mrna)) %>%
  mutate(attributes = paste(ID_cds, parent, sep = ";")) %>%
  mutate(order_ID = seq(1.3, n()+1))
  # select(seqid, source, type, start, end, score, strand, phase, attributes)

# Make sure start and end positions for
# start and stop codons (startc, stopc)
# are calculated/inferred correctly.
seqs_list_gff_startc <- seqs_list_gff_mrna %>%
  mutate(type = "start_codon") %>%
  mutate(end = start) %>%
  mutate(phase = ".") %>%
  mutate(ID_start = paste(ID_mrna, "-start_codon3", sep = "")) %>%
  mutate(parent = gsub(pattern = "ID=", replacement = "Parent=", ID_mrna)) %>%
  mutate(attributes = paste(ID_start, parent, sep = ";")) %>%
  mutate(order_ID = seq(1.4, n()+1))
seqs_list_gff_stopc <- seqs_list_gff_mrna %>%
  mutate(type = "stop_codon") %>%
  mutate(start = end) %>%
  mutate(phase = ".") %>%
  mutate(ID_start = paste(ID_mrna, "-stop_codon4", sep = "")) %>%
  mutate(parent = gsub(pattern = "ID=", replacement = "Parent=", ID_mrna)) %>%
  mutate(attributes = paste(ID_start, parent, sep = ";")) %>%
  mutate(order_ID = seq(1.5, n()+1))


seqs_list_gff <- bind_rows(seqs_list_gff_gene, seqs_list_gff_mrna,
                           seqs_list_gff_exon, seqs_list_gff_cds,
                           seqs_list_gff_startc, seqs_list_gff_stopc) %>%
  arrange(order_ID) %>%
  select(seqid, source, type, start, end,
         score, strand, phase, attributes) %>%
  filter(type != "exon") %>%
  # filter(type != "CDS") %>%
  filter(type != "start_codon") %>%
  filter(type != "stop_codon")
# XXX: Debugging 12/10
#head(seqs_list_gff, n = 5)
#quit(save="no")

readr::write_tsv(seqs_list_gff, file = out_gff_f, col_names = FALSE)


