import argparse
import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(parent_dir)
from lib import expression
import pandas as pd
from styleframe import StyleFrame, Styler, utils

def get_args():

	parser = argparse.ArgumentParser(
		formatter_class = argparse.ArgumentDefaultsHelpFormatter,
		description = "Calculate TPM and CPM"
	)

	parser.add_argument("countfiles", type = str, help = "Text file of all exon-exon junction files to be processed")
	parser.add_argument("output", type = str, help = "Output directory")

	args = parser.parse_args()

	return(args)

def merge_table():

	args = get_args()

	df = pd.read_csv(

		args.countfiles,
		sep = "\t",
		header = None

	)

	for index, row in df.iterrows():

		count_df = pd.read_csv(

			row[0],
			sep = "\t",
			dtype = {
				0: "str",
				1: "int32",
				2: "int32"
			}

		)

		sample = count_df.columns[2]

		if "count_all_df" in locals():

			count_all_df = pd.merge(

				count_all_df,
				count_df[["Geneid", sample]],
				on = "Geneid"

			)

		else:

			count_all_df = count_df

	count_all_df = count_all_df.rename(columns = {"Geneid": "gene_name"})
	count_all_df = count_all_df.sort_values("gene_name")
	return(count_all_df)

def main():

	args = get_args()

	print("Merge count tables...", file = sys.stderr)
	merged_table_df = merge_table()
	count_df = merged_table_df.drop(columns = ["Length"])

	count_df.to_csv(

		args.output + "counts.txt",
		sep = "\t",
		index = False

	)

	# copy count table
	count_a_df = merged_table_df.copy()
	count_b_df = merged_table_df.copy()
	del merged_table_df

	print("Calculate TPM...", file=sys.stderr)
	expression_processor = expression.ExpressionProcessor(count_a_df)
	tpm_df = expression_processor.TPM()
	del count_a_df

	tpm_df.to_csv(
		args.output + "TPM.txt",
		sep="\t",
		index=False
	)

	print("Calculate CPM...", file=sys.stderr)
	expression_processor = expression.ExpressionProcessor(count_b_df)
	cpm_df = expression_processor.CPM()
	del count_b_df

	cpm_df.to_csv(
		args.output + "CPM.txt",
		sep="\t",
		index=False
	)

	# Excel file
	print("Export to an excel file....", file = sys.stderr)

	# Style
	style = Styler(

		horizontal_alignment = utils.horizontal_alignments.left,
		border_type = utils.borders.default_grid,
		wrap_text = False

	)

	with StyleFrame.ExcelWriter(args.output + "TPM_CPM.xlsx") as writer:

		tpm_sf = StyleFrame(tpm_df)
		tpm_sf.set_column_width(columns = tpm_df.columns, width = 20)
		tpm_sf.apply_column_style(cols_to_style = tpm_df.columns, styler_obj = style, style_header = True)
		tpm_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "TPM")

		cpm_sf = StyleFrame(cpm_df)
		cpm_sf.set_column_width(columns = cpm_df.columns, width = 20)
		cpm_sf.apply_column_style(cols_to_style = cpm_df.columns, styler_obj = style, style_header = True)
		cpm_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "CPM")

if __name__ == '__main__':

    main()
