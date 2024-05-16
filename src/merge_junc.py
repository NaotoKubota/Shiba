import argparse
import sys
import pandas as pd

def get_args():

	parser = argparse.ArgumentParser(
		formatter_class = argparse.ArgumentDefaultsHelpFormatter,
		description = "Merge junction read counts"
	)

	parser.add_argument("juncfiles", type = str, help = "Text file of all junction files to be processed")
	parser.add_argument("output", type = str, help = "Output name")

	args = parser.parse_args()

	return(args)

def exonexon():

	print("Merge exon-exon junction count...", file = sys.stdout)

	args = get_args()

	df = pd.read_csv(

		args.juncfiles,
		sep = "\t",
		header = None

	)

	df = df[df[1] == "exon-exon"]

	for index, row in df.iterrows():

		print(row[0].split("/")[-1] + "...", file = sys.stdout)

		junc_df = pd.read_csv(

			row[0],
			sep = "\t",
			header = None,
			dtype = str

		)

		sample = row[0].split("/")[-1].rstrip("_exon-exon.junc")

		junc_df = junc_df.iloc[:, [0, 1, 2, 4, 10]]

		junc_df.columns = ["chr", "start", "end", "count", "block"]
		junc_df = junc_df.reset_index()

		junc_df.loc[(junc_df["chr"].str.isdecimal() == True) | (junc_df["chr"].str.len() <= 2), "chr"] = "chr" + junc_df["chr"]

		junc_df["count"] = junc_df["count"].astype("int32")
		junc_df["blockSize1"] = junc_df["block"].str.split(",", expand = True)[0].astype("int32")
		junc_df["start"] = junc_df["start"].astype("int32") + junc_df["blockSize1"]
		junc_df["blockSize2"] = junc_df["block"].str.split(",", expand = True)[1].astype("int32")
		junc_df["end"] = junc_df["end"].astype("int32") - junc_df["blockSize2"] + 1
		junc_df["ID"] = junc_df["chr"].astype(str) + ":" + junc_df["start"].astype(str) + "-" + junc_df["end"].astype(str)
		junc_df["sample"] = sample
		junc_df = junc_df[["ID", "sample", "count"]]
		# Group by ID and sample
		junc_df = junc_df.groupby(["ID", "sample"], as_index = False).sum()

		if "result_df" in locals():

			result_df = pd.concat(

				[result_df, junc_df],
				axis = 0

			)

		else:

			result_df = junc_df

	result_df["count"] = result_df["count"].astype("int32")
	# Check if there are duplicated junctions
	dup_df = result_df.groupby(["ID", "sample"], as_index = False).count()
	dup_df = dup_df[dup_df["count"] > 1]
	if dup_df.shape[0] > 0:

		print("Duplicated junctions:", file = sys.stdout)
		print(dup_df, file = sys.stdout)
		print("Please check the input files!", file = sys.stdout)
		sys.exit(1)

	result_df = result_df.pivot(

		index = "ID",
		columns = "sample",
		values = "count"

	).fillna(0).reset_index()

	result_df = result_df.rename(columns = {"index": "ID"})

	result_df["chr"] = result_df["ID"].str.split(":", expand = True)[0]
	result_df["start"] = result_df["ID"].str.split(":", expand = True)[1].str.split("-", expand = True)[0].astype("int32")
	result_df["end"] = result_df["ID"].str.split(":", expand = True)[1].str.split("-", expand = True)[1].astype("int32")
	result_df["chr-start"] = result_df["chr"] + "-" + result_df["start"].astype(str)
	result_df["chr-end"] = result_df["chr"] + "-" + result_df["end"].astype(str)

	col = [i for i in result_df.columns if i not in ["chr", "start", "end", "ID", "chr-start", "chr-end", "mean"]]
	for j in col:

		result_df = result_df.astype({j: "int32"})

	col = ["chr", "start", "end", "ID"] + col
	result_df = result_df[col]

	return(result_df)

def exonintron():

	print("Merge exon-intron junction count...", file = sys.stdout)

	args = get_args()

	df = pd.read_csv(

		args.juncfiles,
		sep = "\t",
		header = None

	)

	df = df[df[1] == "exon-intron"]

	for index, row in df.iterrows():

		print(row[0].split("/")[-1] + "...", file = sys.stdout)

		junc_df = pd.read_csv(

			row[0],
			sep = "\t",
			comment = "#",
			dtype = str

		)

		sample = row[0].split("/")[-1].rstrip("_exon-intron.junc")

		junc_df = junc_df.iloc[:, [1, 2, 3, 0, 6]]

		junc_df.columns = ["chr", "start", "end", "ID", "count"]
		junc_df["sample"] = sample

		junc_df.loc[(junc_df["chr"].str.isdecimal() == True) | (junc_df["chr"].str.len() <= 2), "chr"] = "chr" + junc_df["chr"]

		if "result_df" in locals():

			result_df = pd.concat(

				[result_df, junc_df],
				axis = 0

			)

		else:

			result_df = junc_df

	result_df["count"] = result_df["count"].astype("int32")
	result_df = result_df.pivot(

		index = ["chr", "start", "end", "ID"],
		columns = "sample",
		values = "count"

	).fillna(0).reset_index()

	col = [i for i in result_df.columns if i not in ["chr", "start", "end", "ID"]]
	for j in col:

		result_df = result_df.astype({j: "int32"})

	return(result_df)

def main():

	args = get_args()

	# exon-exon junctions
	exon_exon_junc_df = exonexon()

	# exon-intron junctions
	exon_intron_junc_df = exonintron()

	result_df = pd.concat(

		[exon_exon_junc_df, exon_intron_junc_df]

	)

	result_df = result_df.sort_values(["chr", "start"])

	result_df.to_csv(

		args.output,
		sep = "\t",
		index = False

	)

	junc_num = str(result_df.count()[0])
	print("Total number of junctions: " + junc_num, file = sys.stdout)

if __name__ == '__main__':

    main()
