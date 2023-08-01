#!/usr/bin/env python3

import sys
import argparse
import warnings
import copy
from natsort import natsort_keygen
import pandas as pd
import gffpandas.gffpandas as gffpd

def check_overlaps(ann_ref_pd, ann_sup_pd, overlaps_f = None, overlaps_ref_f = None):
	#print(ann_ref_pd.df.head())
	#print(ann_sup_pd.df.head())
	#ann_ref = ann_ref_pd.df
	#ann_sup = ann_sup_pd.df
	#ann_sup_genes = ann_sup_pd.filter_feature_of_type(['gene'])
	#print("Checking only genes for overlaps")
	#is_complem = ann_sup_genes.df.apply(lambda x: detect_overlaps(ref_pd=ann_ref_pd, seq_id=x['seq_id'], seq_type="gene", start=x['start'], end=x['end'], strand=x['strand']), axis=1)
	warnings.filterwarnings("ignore", category = UserWarning)
	is_complem = ann_sup_pd.df.apply(lambda x: detect_overlaps(ref_pd=ann_ref_pd, seq_id=x['seq_id'], start=x['start'], end=x['end'], strand=x['strand']), axis=1)
	if sum(is_complem) > 0:
		print("Checked %s features among %s in supp gff" % (len(is_complem), len(ann_sup_pd.df)))
		print("Found %s overlapping features" % (sum(is_complem)))
		if overlaps_f is not None:
			# Default shallow copy overwrites original, which I still need to return
			#overlaps_df = ann_sup_pd
			overlaps_df = copy.deepcopy(ann_sup_pd)
			overlaps_df.df = overlaps_df.df[is_complem]
			overlaps_df.to_gff3(overlaps_f)
		if overlaps_ref_f is not None:
			is_complem_ref = ann_ref_pd.df.apply(lambda x: detect_overlaps(ref_pd=ann_sup_pd, seq_id=x['seq_id'], start=x['start'], end=x['end'], strand=x['strand']), axis=1)
			print("Found %s overlapping features in reference" % (sum(is_complem_ref)))
			# Overwriting original here w/shallow copy OK because not returned from func. Save mem!
			overlaps_df_ref = ann_ref_pd
			overlaps_df_ref.df = overlaps_df_ref.df[is_complem_ref]
			overlaps_df_ref.to_gff3(overlaps_ref_f)
			
	else:
		print("No overlapping features found")
	# https://stackoverflow.com/a/10678448/9120324
	is_complem_keep = [not elem for elem in is_complem]
	# ann_sup_no_ovls = ann_sup_genes.df[is_complem_keep]
	# XXX Run on all feats, polyexonic genes in the supplemental gff could end up with orphaned exons
	# in the case where one exon overlaps but the other doesnt.
	# We are only adding from getorf which are single exonic so ignore this edge case for now
	ann_sup_pd.df = ann_sup_pd.df[is_complem_keep]
	ann_sup_genes = ann_sup_pd.filter_feature_of_type(['gene'])
	print("%s genes from supplemental gff were found to not overlap" % (len(ann_sup_genes.df)))
	print("Returning GFF with %s features not overlapping with ref" % (len(ann_sup_pd.df)))
	return(ann_sup_pd)

def detect_overlaps(ref_pd, seq_id, start, end, strand):
	'''Detects genes in ref GFF that overlap with each line of a supplementary GFF'''
	#ovl=ref_pd.overlaps_with(seq_id=seq_id, type=seq_type, start=start, end=end, strand=strand, complement = False)
	ovl=ref_pd.overlaps_with(seq_id=seq_id, start=start, end=end, strand=strand, complement = False)
	if len(ovl.df.index) > 0:
		#print("Gene overlaps ref:")
		is_ovl = True
	else:
		#print("No overlapping gene in ref")
		is_ovl = False
	return(is_ovl)

def merge_annots(annot_sup, annot_ref, sort_by_type):
	'''Supplement genes by concatenating dataframes and sorting'''
	print("Types of supplement and reference gff objects")
	# Bind rows
	gff_pds_concat_list = [annot_ref.df, annot_sup.df]
	gff_pd_concat = pd.concat(gff_pds_concat_list)
	# Sort
	# gene -> tRNA -> UTR -> exon -> UTR -> CDS -> mRNA -> UTR -> exon -> UTR -> CDS
	# chrom -> gene position -> type -> position
	# chroms and positions then attribute with parent attribute first
	# TODO: To make more flexible - create a key by finding first time each attribute appears
	# and noting the order it appears after 'gene'. Like, what would happen to start_codon?
	# type_order = ['gene', 'tRNA', 'mRNA', 'five_prime_UTR', 'exon', 'three_prime_UTR', 'CDS']
	# split attribute by = and - into att_value, gene, and att_list, don't remove original
	if sort_by_type is True:
		print("Sort by chr, start, end")
		gff_pd_concat[['field', 'gene', 'att_list']] = gff_pd_concat['attributes'].str.split('=|-|;', n=2, expand=True)
		# Sort by position to get genes in place,
		# then group by gene ID and sort with type order
		sort_order = ['seq_id', 'start', 'end']
		# https://stackoverflow.com/a/63890954/9120324	
		gff_pd_concat = gff_pd_concat.sort_values(sort_order, key = natsort_keygen())
		#print(gff_pd_concat)
		# https://stackoverflow.com/a/54301218/9120324
		type_order = {'gene':0, 'tRNA':1, 'mRNA':2, 'five_prime_UTR':3, 'exon':4, 'three_prime_UTR':5, 'CDS':6}
		print("Sorting GFF output...")
		gff_pd_concat_gene = gff_pd_concat.groupby(['gene'], sort = False, as_index = False)
		gff_pd_concat = gff_pd_concat_gene.apply(lambda x: x.sort_values(['seq_id', 'start', 'end', 'type'], key = lambda y: y.map(type_order)))
		#print(gff_pd_concat)
		# Newly added cols removed automatically during conversion to gff object
	if sort_by_type is False:
		print("New genes appended to end of GFF")
	annot_sup.df = gff_pd_concat
	merged_gff = annot_sup
	return(merged_gff)

def main(ref, supp, outfile, overlaps_f, overlaps_ref_f, sort_by_type):
	# Read gff files into pandas
	gff_ref_f = str(ref)
	gff_sup_f = str(supp)
	annot_ref = gffpd.read_gff3(gff_ref_f)
	annot_sup = gffpd.read_gff3(gff_sup_f)
	# Check for overlaps
	annot_sup_no_ovl = check_overlaps(annot_ref, annot_sup, overlaps_f, overlaps_ref_f)
	# Merge gffs
	merged_gff = merge_annots(annot_sup_no_ovl, annot_ref, sort_by_type)
	# Write file
	print("Writing final GFF to file")
	merged_gff.to_gff3(outfile)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Supplement GFF with additional annotations in another GFF.\n  Reads all GFFs into dataframes so may consume a lot of memory depending on GFF length.", epilog = "author: Nick Carleson (carleson@oregonstate.edu)", formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("-r", "--ref", metavar = "File", help = "Reference GFF", required = True)
	parser.add_argument("-s", "--sup", metavar = "File", help = "GFF file of features to supplement the reference with", required = True)
	parser.add_argument("-o", "--output", metavar = "File", help = "Output GFF file of ref & non-overlapping supplemental genes", required = True)
	parser.add_argument("-l", "--overlaps", metavar = "File", help = "Write overlapping features in supplemental GFF to file", required = False)
	parser.add_argument("-k", "--overlaps_ref", metavar = "File", help = "Write overlapping features in reference GFF to file", required = False)
	#parser.add_argument("--sort_by_type", help = "Sort output gff by type? If unsorted, new feats will be at end of GFF", required = False, dest='sort_by_type', default=True, action='store_true')
	parser.add_argument("--sorted", help = "Sort output gff by type", action='store_true')
	parser.add_argument("--unsorted", help = "Do not sort output. New feats will be at end of GFF. MUCH faster runtime.",\
dest='sorted', action='store_false')
	parser.set_defaults(sort=True)

	args = parser.parse_args()
	
	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)
	else:
		main(args.ref, args.sup, args.output, args.overlaps, args.overlaps_ref, args.sorted)

	
