import argparse
import sys
import os
import pandas as pd
import numpy as np
import collections
from collections import defaultdict
import multiprocessing as mp
import itertools
import time
import concurrent.futures

"""
This script converts a GTF file into a pandas DataFrame containing information about alternative splicing events.
"""

def get_args():
	"""
	Parses command line arguments.

	Returns:
		argparse.Namespace: An object containing the parsed arguments.
	"""

	parser = argparse.ArgumentParser(
		description = "Extract alternative splicing events from GTF file"
	)

	parser.add_argument("-i", "--gtf", type = str, help = "Input GTF file", required = True)
	parser.add_argument("-r", "--reference-gtf", type = str, help = "Reference GTF file", required = False)
	parser.add_argument("-o", "--output", type = str, help = "Output directory", required = True)
	parser.add_argument("-p", "--num-process", type = int, help = "Number of processors to use", default = 1)

	args = parser.parse_args()

	return(args)

def gtf(gtf, num_process) -> pd.DataFrame:
	"""
	Reads a GTF file and extracts exon information to create a pandas DataFrame.

	Args:
		gtf (str): The path to the GTF file.

	Returns:
		Dict: A dictionary containing information about the GTF file.
	"""

	gtf_df = pd.read_csv(

		gtf,
		sep = "\t",
		usecols = [
			0, 2, 3, 4, 6, 8
		],
		dtype = {
			0: "str",
			2: "str",
			3: "int32",
			4: "int32",
			6: "str",
			8: "str"
		},
		comment = "#",
		header = None

	)

	gtf_df = gtf_df[gtf_df[2] == "exon"]
	gtf_df = gtf_df.reset_index()
	gtf_df = gtf_df[[0, 3, 4, 6, 8]]
	gtf_df.columns = ["chr", "start", "end", "strand", "information"]

	gtf_info = gtf_df.information.values
	gene_id_dic = {}
	gene_name_dic = {}
	gene_id_list_dic = defaultdict(list)
	gene_name_list_dic = defaultdict(list)
	gene_id_col = []
	gene_name_col = []
	transcript_id_col = []

	for index in range(gtf_df.shape[0]):

		dic = {}

		l = gtf_info[index].split(";")[0:-1]

		for i in l:

			if '"' in i:

				key = i.split('"')[0].strip(" ")
				value = i.split('"')[1]

			else:

				key = i.strip(" ").split(" ")[0]
				value = i.strip(" ").split(" ")[1]

			dic[key] = value

		if "ref_gene_id" in dic:

			gene_id = dic["ref_gene_id"]

			gene_id_dic[dic["gene_id"]] = dic["ref_gene_id"]
			gene_id_list_dic[dic["gene_id"]] += [dic["ref_gene_id"]]

		else:

			gene_id = dic["gene_id"]

		if "gene_name" in dic:

			gene_name = dic["gene_name"]

			gene_name_dic[dic["gene_id"]] = dic["gene_name"]
			gene_name_list_dic[dic["gene_id"]] += [dic["gene_name"]]


		else:

			if "ref_gene_id" in dic:

				gene_name = dic["ref_gene_id"]

			else:

				gene_name = dic["gene_id"]

		transcript_id = dic["transcript_id"]

		gene_id_col += [gene_id]
		gene_name_col += [gene_name]
		transcript_id_col += [transcript_id]

	gtf_df["gene_id"] = gene_id_col
	gtf_df["gene_name"] = gene_name_col
	gtf_df["transcript_id"] = transcript_id_col

	gtf_gene_id = gtf_df.gene_id.values
	gtf_gene_name = gtf_df.gene_name.values
	gtf_chr = gtf_df.chr.values

	for index in range(gtf_df.shape[0]):

		if gtf_gene_id[index] in gene_id_list_dic:

			l = gene_id_list_dic[gtf_gene_id[index]]
			c = collections.Counter(l)
			most_common = c.most_common()[0][0]

			gtf_df.at[index, "gene_id"] = most_common

		if gtf_gene_name[index] in gene_name_list_dic:

			l = gene_name_list_dic[gtf_gene_name[index]]
			c = collections.Counter(l)
			most_common = c.most_common()[0][0]

			gtf_df.at[index, "gene_name"] = most_common

		if ~(gtf_chr[index].startswith("chr")) and (len(gtf_chr[index]) <= 2):

			gtf_df.at[index, "chr"] = "chr" + gtf_chr[index]

	gtf_df = gtf_df.sort_values(["gene_id", "transcript_id", "start"])
	gtf_df = gtf_df.reset_index()
	gtf_gene_id = gtf_df.gene_id.values
	gtf_gene_name = gtf_df.gene_name.values
	gtf_transcript_id = gtf_df.transcript_id.values
	gtf_chr = gtf_df.chr.values
	gtf_start = gtf_df.start.values
	gtf_end = gtf_df.end.values
	gtf_strand = gtf_df.strand.values
	gtf_dic = defaultdict(dict)
	transcript = ""
	gene = ""

	for index in range(gtf_df.shape[0]):

		gtf_dic[gtf_gene_id[index]]["gene_name"] = gtf_gene_name[index]
		gtf_dic[gtf_gene_id[index]]["chr"] = gtf_chr[index]
		gtf_dic[gtf_gene_id[index]]["strand"] = gtf_strand[index]

		try:
			gtf_dic[gtf_gene_id[index]]["start"] = np.append(gtf_dic[gtf_gene_id[index]]["start"], gtf_start[index])
		except:
			gtf_dic[gtf_gene_id[index]]["start"] = np.array([gtf_start[index]])

		try:
			gtf_dic[gtf_gene_id[index]]["end"] = np.append(gtf_dic[gtf_gene_id[index]]["end"], gtf_end[index])
		except:
			gtf_dic[gtf_gene_id[index]]["end"] = np.array([gtf_end[index]])

		try:
			gtf_dic[gtf_gene_id[index]]["exon_list"].add(gtf_chr[index] + ":" + str(gtf_start[index]) + "-" + str(gtf_end[index]))
		except:
			gtf_dic[gtf_gene_id[index]]["exon_list"] = {gtf_chr[index] + ":" + str(gtf_start[index]) + "-" + str(gtf_end[index])}

		if "start_dic" not in gtf_dic[gtf_gene_id[index]]:
			gtf_dic[gtf_gene_id[index]]["start_dic"] = {}
		try:
			gtf_dic[gtf_gene_id[index]]["start_dic"][str(gtf_start[index])].add(str(gtf_end[index]))
		except:
			gtf_dic[gtf_gene_id[index]]["start_dic"][str(gtf_start[index])] = {str(gtf_end[index])}

		if "end_dic" not in gtf_dic[gtf_gene_id[index]]:
			gtf_dic[gtf_gene_id[index]]["end_dic"] = {}
		try:
			gtf_dic[gtf_gene_id[index]]["end_dic"][str(gtf_end[index])].add(str(gtf_start[index]))
		except:
			gtf_dic[gtf_gene_id[index]]["end_dic"][str(gtf_end[index])] = {str(gtf_start[index])}

		if "transcript_exon_dic" not in gtf_dic[gtf_gene_id[index]]:
			gtf_dic[gtf_gene_id[index]]["transcript_exon_dic"] = {}
		if gtf_transcript_id[index] not in gtf_dic[gtf_gene_id[index]]["transcript_exon_dic"]:
			gtf_dic[gtf_gene_id[index]]["transcript_exon_dic"][gtf_transcript_id[index]] = set()
		try:
			gtf_dic[gtf_gene_id[index]]["transcript_exon_dic"][gtf_transcript_id[index]].add(gtf_chr[index] + ":" + str(gtf_start[index]) + "-" + str(gtf_end[index]))
		except:
			gtf_dic[gtf_gene_id[index]]["transcript_exon_dic"][gtf_transcript_id[index]] = {gtf_chr[index] + ":" + str(gtf_start[index]) + "-" + str(gtf_end[index])}

		if gtf_transcript_id[index] == transcript:

			if "intron_start_dic" not in gtf_dic[gtf_gene_id[index]]:
				gtf_dic[gtf_gene_id[index]]["intron_start_dic"] = {}
			try:
				gtf_dic[gtf_gene_id[index]]["intron_start_dic"][str(end)].add(str(gtf_start[index]))
			except:
				gtf_dic[gtf_gene_id[index]]["intron_start_dic"][str(end)] = {str(gtf_start[index])}

			if "intron_end_dic" not in gtf_dic[gtf_gene_id[index]]:
				gtf_dic[gtf_gene_id[index]]["intron_end_dic"] = {}
			try:
				gtf_dic[gtf_gene_id[index]]["intron_end_dic"][str(gtf_start[index])].add(str(end))
			except:
				gtf_dic[gtf_gene_id[index]]["intron_end_dic"][str(gtf_start[index])] = {str(end)}

			if "intron_list" not in gtf_dic[gtf_gene_id[index]]:
				gtf_dic[gtf_gene_id[index]]["intron_list"] = set()
			try:
				gtf_dic[gtf_gene_id[index]]["intron_list"].add(gtf_chr[index] + ":" + str(end) + "-" + str(gtf_start[index]))
			except:
				gtf_dic[gtf_gene_id[index]]["intron_list"] = {gtf_chr[index] + ":" + str(end) + "-" + str(gtf_start[index])}

			if "transcript_intron_dic" not in gtf_dic[gtf_gene_id[index]]:
				gtf_dic[gtf_gene_id[index]]["transcript_intron_dic"] = {}
			if gtf_transcript_id[index] not in gtf_dic[gtf_gene_id[index]]["transcript_intron_dic"]:
				gtf_dic[gtf_gene_id[index]]["transcript_intron_dic"][gtf_transcript_id[index]] = set()
			try:
				gtf_dic[gtf_gene_id[index]]["transcript_intron_dic"][gtf_transcript_id[index]].add(gtf_chr[index] + ":" + str(end) + "-" + str(gtf_start[index]))
			except:
				gtf_dic[gtf_gene_id[index]]["transcript_intron_dic"][gtf_transcript_id[index]] = {gtf_chr[index] + ":" + str(end) + "-" + str(gtf_start[index])}

		end = gtf_end[index]
		transcript = gtf_transcript_id[index]

	# Discard genes with only one transcript
	gtf_dic = {k: v for k, v in gtf_dic.items() if len(v["transcript_exon_dic"]) > 1}

	# Split gene list into number of processes
	gene_l_split = np.array_split(list(gtf_dic.keys()), num_process)
	# Split dictionry into number of processes
	gtf_dic_split = defaultdict(dict)
	for i in range(num_process):
		gtf_dic_split[i] = defaultdict(dict)
		for gene in gene_l_split[i]:
			gtf_dic_split[i][gene] = gtf_dic[gene]

	return(gtf_dic_split)

def gtf_exon_set(gtf_path) -> set:
	"""
	Reads a GTF file and extracts exon information to create a set of exon coordinates.

	Args:
		gtf (str): The path to the GTF file.

	Returns:
		set: A set of exon coordinates in the format "chr:start-end".
	"""

	gtf_df = pd.read_csv(

		gtf_path,
		sep = "\t",
		usecols = [
			0, 2, 3, 4
		],
		dtype = {
			0: "str",
			2: "str",
			3: "int32",
			4: "int32"
		},
		comment = "#",
		header = None

	)

	gtf_df = gtf_df[gtf_df[2] == "exon"][[0, 3, 4]]
	gtf_df.columns = ["chr", "start", "end"]
	gtf_df.loc[(~(gtf_df["chr"].str.startswith("chr")) & (gtf_df["chr"].str.len() <= 2)), "chr"] = "chr" + gtf_df["chr"]
	gtf_df["exon"] = gtf_df["chr"] + ":" + gtf_df["start"].astype(str) + "-" + gtf_df["end"].astype(str)
	gtf_exon_set = set(gtf_df["exon"])

	return(gtf_exon_set)


def se(gtf_dic) -> list:
	"""
	Make skipped exon list.

	Args:
		gtf_dic: A dictionary containing information about the GTF file.

	Returns:
		list: List of skipped exon events, where each event is represented as a list of the form
		[exon, inc1, inc2, exc, strand, gene, gene_name].

	"""

	event_l = []

	for gene in gtf_dic.keys():

		if "intron_list" not in gtf_dic[gene]:
			continue

		chr = gtf_dic[gene]["chr"]
		strand = gtf_dic[gene]["strand"]
		gene_name = gtf_dic[gene]["gene_name"]
		gene_start_values = gtf_dic[gene]["start"]
		gene_end_values = gtf_dic[gene]["end"]
		intron_list = gtf_dic[gene]["intron_list"]
		intron_start_dict = gtf_dic[gene]["intron_start_dic"]
		intron_end_dict = gtf_dic[gene]["intron_end_dic"]

		exon_list = gtf_dic[gene]["exon_list"]
		exon_list_unique = np.unique([[i.split(":")[1].split("-")[0], i.split(":")[1].split("-")[1]] for i in exon_list], axis = 0)
		exon_start = np.array([i[0] for i in exon_list_unique]).astype("int32")
		exon_end = np.array([i[1] for i in exon_list_unique]).astype("int32")

		for index in range(len(exon_start)):

			# (inc1)[exon](inc2)
			# (x1, y1)[y1, x2](x2, y2)
			# inc1: (x1, y1)
			# inc2: (x2, y2)
			# exc: (x1, y2)

			y1 = str(exon_start[index])
			x2 = str(exon_end[index])

			x1_list = intron_end_dict.get(y1, set())
			y2_list = intron_start_dict.get(x2, set())

			x1_y2_iter = itertools.product(x1_list, y2_list)
			for x1, y2 in x1_y2_iter:

				exon = chr + ":" + str(y1) + "-" + str(x2)
				inc1 = chr + ":" + str(x1) + "-" + str(y1)
				inc2 = chr + ":" + str(x2) + "-" + str(y2)
				exc = chr + ":" + str(x1) + "-" + str(y2)

				if exc in intron_list:

					event_l += [[exon, inc1, inc2, exc, strand, gene, gene_name]]

	return(event_l)


def mse(gtf_dic) -> list:
	'''
	Make multi-skipped exon list.

	Args:
		gtf_dic: A dictionary containing information about the GTF file.

	Returns:
		list: List of multi-skipped exon events, where each event is represented as a list of the form
		[exonlist, intronlist, mse_n, strand, gene, gene_name].
		exonlist is concatenated exon list with semi-colon (e.g. exon1;exon2;exon3).
		intronlist is concatenated intron list with semi-colon (e.g. intron1;intron2;intron3;intron4;exclusion_intron).
		mse_n is the number of exons skipped.
	'''

	event_l = []

	for gene in gtf_dic.keys():

		if "intron_list" not in gtf_dic[gene]:
			continue

		chr = gtf_dic[gene]["chr"]
		strand = gtf_dic[gene]["strand"]
		gene_name = gtf_dic[gene]["gene_name"]
		gene_start_values = gtf_dic[gene]["start"]
		gene_end_values = gtf_dic[gene]["end"]
		intron_list = gtf_dic[gene]["intron_list"]
		intron_start_dict = gtf_dic[gene]["intron_start_dic"]
		intron_end_dict = gtf_dic[gene]["intron_end_dic"]
		intron_dic = gtf_dic[gene]["transcript_intron_dic"]
		exon_dic = gtf_dic[gene]["transcript_exon_dic"]
		# Sort exons by start position, ascending order
		exon_dic = {k: sorted(list(v), key = lambda x: int(x.split(":")[1].split("-")[0])) for k, v in exon_dic.items()}
		# Transcript list sorted by exon number
		transcript_list = sorted(exon_dic, key = lambda x: len(exon_dic[x]))

		# Identify multi-skipped exon events until five-hundredth exon skipping
		for mse_n in range(2, 501):

			# Get transcript with at least mse_n+2 exons
			transcript_list = [transcript for transcript in exon_dic.keys() if len(exon_dic[transcript]) >= mse_n + 2]
			if len(transcript_list) == 0:
				break
			for transcript in transcript_list:

				exon_list_in_transcript = exon_dic[transcript]
				exon_start_in_transcript = np.array([i.split(":")[1].split("-")[0] for i in exon_list_in_transcript]).astype("int32")
				exon_end_in_transcript = np.array([i.split(":")[1].split("-")[1] for i in exon_list_in_transcript]).astype("int32")

				# Get combinations of n adjuscent index
				# e.g. (1, 2) when mse_n = 2 and exon number = 3, first and last exons are excluded
				# e.g. (1, 2), (2, 3) when mse_n = 2 and exon number = 4, first and last exons are excluded
				# e.g. (1, 2, 3), (2, 3, 4) when mse_n = 3 and exon number = 6, first and last exons are excluded
				# e.g. (1, 2, 3), (2, 3, 4), (3, 4, 5) when mse_n = 3 and exon number = 7, first and last exons are excluded
				idx_number_list = [i for i in range(len(exon_list_in_transcript) - mse_n)] # e.g. [0, 1] when mse_n = 2 and exon number = 4
				idx_list_list = [list(range(i + 1, i + 1 + mse_n)) for i in idx_number_list] # e.g. [[1, 2], [2, 3]] when mse_n = 2 and exon number = 4
				for idx_list in idx_list_list:

					# (inc_1)[exon_1](inc_2)[exon_2]...[exon_(mse_n-1)](inc_(mse_n))[exon_(mse_n)](inc_(mse_n+1))
					# (x1, y1)[y1, x2](x2, y2)[y2, x3]...[x(mse_n-1), y(mse_n)](x(mse_n), y(mse_n))[y(mse_n), x(mse_n+1)](x(mse_n+1), y(mse_n+1))
					# inc1: (x1, y1)
					# inc2: (x2, y2)
					# ...
					# inc(mse_n): (x(mse_n), y(mse_n))
					# inc(mse_n+1): (x(mse_n+1), y(mse_n+1))
					# exc: (x1, y(mse_n+1))

					all_exons = [chr + ":" + str(exon_start_in_transcript[i]) + "-" + str(exon_end_in_transcript[i]) for i in idx_list]
					exonlist = ";".join(all_exons)

					x1_list = intron_end_dict.get(str(exon_start_in_transcript[idx_list[0]]), set())
					y_mse_n_1_list = intron_start_dict.get(str(exon_end_in_transcript[idx_list[mse_n - 1]]), set())
					for x1, y_mse_n_1 in itertools.product(x1_list, y_mse_n_1_list):

						all_inclusion_introns = [chr + ":" + str(x1) + "-" + str(exon_start_in_transcript[idx_list[0]])] # inc1 (first intron)
						for i in range(mse_n - 2): # inc2 to inc(mse_n)
							all_inclusion_introns += [chr + ":"
								+ str(exon_end_in_transcript[idx_list[i + 1]])
								+ "-"
								+ str(exon_start_in_transcript[idx_list[i + 2]])
							]
						all_inclusion_introns += [chr + ":" + str(exon_end_in_transcript[idx_list[mse_n - 1]]) + "-" + str(y_mse_n_1)] # inc(mse_n+1)
						exc = chr + ":" + str(x1) + "-" + str(y_mse_n_1)
						all_introns = all_inclusion_introns + [exc]
						intronlist = ";".join(all_introns)

						# Check if all inclusion introns are and exclusion introns are NOT present in the same transcript
						if (set(all_inclusion_introns) <= intron_dic[transcript]) and (exc not in intron_dic[transcript]) and (exc in intron_list):

							event_l += [[exonlist, intronlist, mse_n, strand, gene, gene_name]]

	return(event_l)


def five(gtf_dic) -> list:
	"""
	Make alternative five prime ss list.

	Args:
		gtf_dic: A dictionary containing information about the GTF file.

	Returns:
		list: List of alternative five prime ss events, where each event is represented as a list of the form
		[exon_a, exon_b, intron_a, intron_b, strand, gene, gene_name].

	"""

	event_l = []

	for gene in gtf_dic.keys():

		if "intron_list" not in gtf_dic[gene]:
			continue

		chr = gtf_dic[gene]["chr"]
		strand = gtf_dic[gene]["strand"]
		gene_name = gtf_dic[gene]["gene_name"]
		gene_start_values = gtf_dic[gene]["start"]
		gene_end_values = gtf_dic[gene]["end"]
		intron_list = gtf_dic[gene]["intron_list"]
		intron_start_dict = gtf_dic[gene]["intron_start_dic"]
		intron_end_dict = gtf_dic[gene]["intron_end_dic"]
		intron_start_set = set(intron_start_dict.keys())
		intron_end_set = set(intron_end_dict.keys())

		if strand == "+":

			exon_dic = gtf_dic[gene]["start_dic"]

		else:

			exon_dic = gtf_dic[gene]["end_dic"]

		five_dic = {}

		for key in exon_dic.keys():

			if len(exon_dic[key]) != 1:

				five_dic[key] = exon_dic[key]

		for con in five_dic.keys():

			alt_list = five_dic[con]

			i_j_iter = itertools.combinations(alt_list, 2)
			for i, j in i_j_iter:

				set_ij = {str(i), str(j)}

				if strand == "+":

					if (set_ij & intron_start_set) == set_ij:

						intron_i = intron_start_dict[str(i)]
						intron_j = intron_start_dict[str(j)]

						intron_common = list(set(intron_i) & set(intron_j))

						for s in intron_common:

							if int(j) < int(i):
								exon_a_end = i
								exon_b_end = j
								intron_a_start = i
								intron_b_start = j
							else:
								exon_a_end = j
								exon_b_end = i
								intron_a_start = j
								intron_b_start = i

							exon_a = chr + ":" + str(con) + "-" + str(exon_a_end)
							exon_b = chr + ":" + str(con) + "-" + str(exon_b_end)

							intron_a = chr + ":" + str(intron_a_start) + "-" + str(s)
							intron_b = chr + ":" + str(intron_b_start) + "-" + str(s)

							event_l += [[exon_a, exon_b, intron_a, intron_b, strand, gene, gene_name]]

				else:

					if (set_ij & intron_end_set) == set_ij:

						intron_i = intron_end_dict[str(i)]
						intron_j = intron_end_dict[str(j)]

						intron_common = list(set(intron_i) & set(intron_j))

						for s in intron_common:

							if int(j) < int(i):
								exon_a_start = j
								exon_b_start = i
								intron_a_end = j
								intron_b_end = i
							else:
								exon_a_start = i
								exon_b_start = j
								intron_a_end = i
								intron_b_end = j

							exon_a = chr + ":" + str(exon_a_start) + "-" + str(con)
							exon_b = chr + ":" + str(exon_b_start) + "-" + str(con)

							intron_a = chr + ":" + str(s) + "-" + str(intron_a_end)
							intron_b = chr + ":" + str(s) + "-" + str(intron_b_end)

							event_l += [[exon_a, exon_b, intron_a, intron_b, strand, gene, gene_name]]

	return(event_l)


def three(gtf_dic) -> list:
	"""
	Make alternative three prime ss list.

	Args:
		gtf_dic: A dictionary containing information about the GTF file.

	Returns:
		list: List of alternative three prime ss events, where each event is represented as a list of the form
		[exon_a, exon_b, intron_a, intron_b, strand, gene, gene_name].

	"""

	event_l = []

	for gene in gtf_dic.keys():

		if "intron_list" not in gtf_dic[gene]:
			continue

		chr = gtf_dic[gene]["chr"]
		strand = gtf_dic[gene]["strand"]
		gene_name = gtf_dic[gene]["gene_name"]
		gene_start_values = gtf_dic[gene]["start"]
		gene_end_values = gtf_dic[gene]["end"]
		intron_list = gtf_dic[gene]["intron_list"]
		intron_start_dict = gtf_dic[gene]["intron_start_dic"]
		intron_end_dict = gtf_dic[gene]["intron_end_dic"]
		intron_start_set = set(intron_start_dict.keys())
		intron_end_set = set(intron_end_dict.keys())

		if strand == "+":

			exon_dic = gtf_dic[gene]["end_dic"]

		else:

			exon_dic = gtf_dic[gene]["start_dic"]

		three_dic = {}

		for key in exon_dic.keys():

			if len(exon_dic[key]) != 1:

				three_dic[key] = exon_dic[key]

		for con in three_dic.keys():

			alt_list = three_dic[con]

			i_j_iter = itertools.combinations(alt_list, 2)
			for i, j in i_j_iter:

				set_ij = {str(i), str(j)}

				if strand == "+":

					if (set_ij & intron_end_set) == set_ij:

						intron_i = intron_end_dict[str(i)]
						intron_j = intron_end_dict[str(j)]

						intron_common = list(set(intron_i) & set(intron_j))

						for s in intron_common:

							if int(j) < int(i):
								exon_a_start = j
								exon_b_start = i
								intron_a_end = j
								intron_b_end = i
							else:
								exon_a_start = i
								exon_b_start = j
								intron_a_end = i
								intron_b_end = j

							exon_a = chr + ":" + str(exon_a_start) + "-" + str(con)
							exon_b = chr + ":" + str(exon_b_start) + "-" + str(con)

							intron_a = chr + ":" + str(s) + "-" + str(intron_a_end)
							intron_b = chr + ":" + str(s) + "-" + str(intron_b_end)

							event_l += [[exon_a, exon_b, intron_a, intron_b, strand, gene, gene_name]]

				else:

					if (set_ij & intron_start_set) == set_ij:

						intron_i = intron_start_dict[str(i)]
						intron_j = intron_start_dict[str(j)]

						intron_common = list(set(intron_i) & set(intron_j))

						for s in intron_common:

							if int(j) < int(i):
								exon_a_end = i
								exon_b_end = j
								intron_a_start = i
								intron_b_start = j
							else:
								exon_a_end = j
								exon_b_end = i
								intron_a_start = j
								intron_b_start = i

							exon_a = chr + ":" + str(con) + "-" + str(exon_a_end)
							exon_b = chr + ":" + str(con) + "-" + str(exon_b_end)

							intron_a = chr + ":" + str(intron_a_start) + "-" + str(s)
							intron_b = chr + ":" + str(intron_b_start) + "-" + str(s)

							event_l += [[exon_a, exon_b, intron_a, intron_b, strand, gene, gene_name]]

	return(event_l)


def afe(gtf_dic) -> list:
	'''
	Make alternative first exon list.

	Args:
		gtf_dic: A dictionary containing information about the GTF file.

	Returns:
		list: List of alternative first exon events, where each event is represented as a list of the form
		[exon_a, exon_b, intron_a, intron_b, strand, gene, gene_name].
	'''

	event_l = []

	for gene in gtf_dic.keys():

		if "intron_list" not in gtf_dic[gene]:
			continue

		chr = gtf_dic[gene]["chr"]
		strand = gtf_dic[gene]["strand"]
		gene_name = gtf_dic[gene]["gene_name"]
		gene_start_values = gtf_dic[gene]["start"]
		gene_end_values = gtf_dic[gene]["end"]
		intron_list = gtf_dic[gene]["intron_list"]
		intron_start_dict = gtf_dic[gene]["intron_start_dic"]
		intron_end_dict = gtf_dic[gene]["intron_end_dic"]
		intron_dic = gtf_dic[gene]["transcript_intron_dic"]
		exon_dic = gtf_dic[gene]["transcript_exon_dic"]
		# Sort exons by start position, ascending order
		exon_dic = {k: sorted(list(v), key = lambda x: int(x.split(":")[1].split("-")[0])) for k, v in exon_dic.items()}
		# Transcript list sorted by exon number
		transcript_list = sorted(exon_dic, key = lambda x: len(exon_dic[x]))
		# Keep transcripts with at least two exons
		transcript_list = [transcript for transcript in transcript_list if len(exon_dic[transcript]) >= 2]
		if len(transcript_list) < 2:
			continue

		for transcript1, transcript2 in itertools.combinations(transcript_list, 2):

			# Get first exons
			first_exon_list = [exon_dic[transcript][0] if strand == "+" else exon_dic[transcript][-1] for transcript in [transcript1, transcript2]]
			first_exon_transcript1 = first_exon_list[0]
			first_exon_transcript1_start = first_exon_transcript1.split(":")[1].split("-")[0]
			first_exon_transcript1_end = first_exon_transcript1.split(":")[1].split("-")[1]
			first_exon_transcript2 = first_exon_list[1]
			first_exon_transcript2_start = first_exon_transcript2.split(":")[1].split("-")[0]
			first_exon_transcript2_end = first_exon_transcript2.split(":")[1].split("-")[1]
			# Get second exons
			second_exon_list = [exon_dic[transcript][1] if strand == "+" else exon_dic[transcript][-2] for transcript in [transcript1, transcript2]]
			second_exon_transcript1 = second_exon_list[0]
			second_exon_transcript1_start = second_exon_transcript1.split(":")[1].split("-")[0]
			second_exon_transcript1_end = second_exon_transcript1.split(":")[1].split("-")[1]
			second_exon_transcript2 = second_exon_list[1]
			second_exon_transcript2_start = second_exon_transcript2.split(":")[1].split("-")[0]
			second_exon_transcript2_end = second_exon_transcript2.split(":")[1].split("-")[1]

			if strand == "+":

				# [exc1](inc1), [exc2](inc2)
				# [first_exon_distal_start, first_exon_distal_end](x1, y1), [first_exon_proximal_start, first_exon_proximal_end](x2, y2)
				# exc1: [first_exon_distal_start, first_exon_distal_end]
				# inc1: (x1 = first_exon_distal_end, y1 = second_exon_start)
				# exc2: [first_exon_proximal_start, first_exon_proximal_end]
				# inc2: (x2 = first_exon_proximal_end, y2 = second_exon_start)

				# Check if the second exons have the same start position
				if second_exon_transcript1_start != second_exon_transcript2_start:
					continue

				# Set distal and proximal first exons
				if (first_exon_transcript1_start < first_exon_transcript2_start) and (first_exon_transcript1_end < first_exon_transcript2_end):

					first_exon_distal = first_exon_transcript1
					first_exon_distal_start = first_exon_transcript1_start
					first_exon_distal_end = first_exon_transcript1_end
					first_exon_proximal = first_exon_transcript2
					first_exon_proximal_start = first_exon_transcript2_start
					first_exon_proximal_end = first_exon_transcript2_end

				elif (first_exon_transcript1_start > first_exon_transcript2_start) and (first_exon_transcript1_end > first_exon_transcript2_end):

					first_exon_distal = first_exon_transcript2
					first_exon_distal_start = first_exon_transcript2_start
					first_exon_distal_end = first_exon_transcript2_end
					first_exon_proximal = first_exon_transcript1
					first_exon_proximal_start = first_exon_transcript1_start
					first_exon_proximal_end = first_exon_transcript1_end

				else:

					continue

				second_exon_start = second_exon_transcript1_start
				exon_a = chr + ":" + first_exon_distal_start + "-" + first_exon_distal_end
				exon_b = chr + ":" + first_exon_proximal_start + "-" + first_exon_proximal_end
				intron_a = chr + ":" + first_exon_distal_end + "-" + second_exon_start
				intron_b = chr + ":" + first_exon_proximal_end + "-" + second_exon_start
				intron_c = chr + ":" + first_exon_distal_end + "-" + first_exon_proximal_start # Intron connecting the distal and proximal first exons

			else: # strand == "-"

				# (inc1)[exc1], (inc2)[exc2]
				# (x1, y1)[first_exon_proximal_start, first_exon_proximal_end], (x2, y2)[first_exon_distal_start, first_exon_distal_end]
				# inc1: (x1 = second_exon_end, y1 = first_exon_proximal_start)
				# exc1: [first_exon_proximal_start, first_exon_proximal_end]
				# inc2: (x2 = second_exon_end, y2 = first_exon_distal_start)
				# exc2: [first_exon_distal_start, first_exon_distal_end]

				# Check if the second exons have the same end position
				if second_exon_transcript1_end != second_exon_transcript2_end:
					continue

				# Set distal and proximal first exons
				if (first_exon_transcript1_start < first_exon_transcript2_start) and (first_exon_transcript1_end < first_exon_transcript2_end):

					first_exon_distal = first_exon_transcript2
					first_exon_distal_start = first_exon_transcript2_start
					first_exon_distal_end = first_exon_transcript2_end
					first_exon_proximal = first_exon_transcript1
					first_exon_proximal_start = first_exon_transcript1_start
					first_exon_proximal_end = first_exon_transcript1_end

				elif (first_exon_transcript1_start > first_exon_transcript2_start) and (first_exon_transcript1_end > first_exon_transcript2_end):

					first_exon_distal = first_exon_transcript1
					first_exon_distal_start = first_exon_transcript1_start
					first_exon_distal_end = first_exon_transcript1_end
					first_exon_proximal = first_exon_transcript2
					first_exon_proximal_start = first_exon_transcript2_start
					first_exon_proximal_end = first_exon_transcript2_end

				else:

					continue

				second_exon_end = second_exon_transcript1_end
				exon_a = chr + ":" + first_exon_distal_start + "-" + first_exon_distal_end
				exon_b = chr + ":" + first_exon_proximal_start + "-" + first_exon_proximal_end
				intron_a = chr + ":" + second_exon_end + "-" + first_exon_distal_start
				intron_b = chr + ":" + second_exon_end + "-" + first_exon_proximal_start
				intron_c = chr + ":" + first_exon_proximal_end + "-" + first_exon_distal_start # Intron connecting the distal and proximal first exons

			# Check if no intron connecting the distal exon and the proximal exon
			if intron_c in intron_list:
				continue

			event_l += [[exon_a, exon_b, intron_a, intron_b, strand, gene, gene_name]]

	return(event_l)


def ale(gtf_dic) -> list:
	'''
	Make alternative last exon list.

	Args:
		gtf_dic: A dictionary containing information about the GTF file.

	Returns:
		list: List of alternative last exon events, where each event is represented as a list of the form
		[exon_a, exon_b, intron_a, intron_b, strand, gene, gene_name].
	'''

	event_l = []

	for gene in gtf_dic.keys():

		if "intron_list" not in gtf_dic[gene]:
			continue

		chr = gtf_dic[gene]["chr"]
		strand = gtf_dic[gene]["strand"]
		gene_name = gtf_dic[gene]["gene_name"]
		gene_start_values = gtf_dic[gene]["start"]
		gene_end_values = gtf_dic[gene]["end"]
		intron_list = gtf_dic[gene]["intron_list"]
		intron_start_dict = gtf_dic[gene]["intron_start_dic"]
		intron_end_dict = gtf_dic[gene]["intron_end_dic"]
		intron_dic = gtf_dic[gene]["transcript_intron_dic"]
		exon_dic = gtf_dic[gene]["transcript_exon_dic"]
		# Sort exons by start position, ascending order
		exon_dic = {k: sorted(list(v), key = lambda x: int(x.split(":")[1].split("-")[0])) for k, v in exon_dic.items()}
		# Transcript list sorted by exon number
		transcript_list = sorted(exon_dic, key = lambda x: len(exon_dic[x]))
		# Keep transcripts with at least two exons
		transcript_list = [transcript for transcript in transcript_list if len(exon_dic[transcript]) >= 2]
		if len(transcript_list) < 2:
			continue

		for transcript1, transcript2 in itertools.combinations(transcript_list, 2):

			# Get last exons
			last_exon_list = [exon_dic[transcript][-1] if strand == "+" else exon_dic[transcript][0] for transcript in [transcript1, transcript2]]
			last_exon_transcript1 = last_exon_list[0]
			last_exon_transcript1_start = last_exon_transcript1.split(":")[1].split("-")[0]
			last_exon_transcript1_end = last_exon_transcript1.split(":")[1].split("-")[1]
			last_exon_transcript2 = last_exon_list[1]
			last_exon_transcript2_start = last_exon_transcript2.split(":")[1].split("-")[0]
			last_exon_transcript2_end = last_exon_transcript2.split(":")[1].split("-")[1]
			# Get penultimate exons
			penultimate_exon_list = [exon_dic[transcript][-2] if strand == "+" else exon_dic[transcript][1] for transcript in [transcript1, transcript2]]
			penultimate_exon_transcript1 = penultimate_exon_list[0]
			penultimate_exon_transcript1_start = penultimate_exon_transcript1.split(":")[1].split("-")[0]
			penultimate_exon_transcript1_end = penultimate_exon_transcript1.split(":")[1].split("-")[1]
			penultimate_exon_transcript2 = penultimate_exon_list[1]
			penultimate_exon_transcript2_start = penultimate_exon_transcript2.split(":")[1].split("-")[0]
			penultimate_exon_transcript2_end = penultimate_exon_transcript2.split(":")[1].split("-")[1]

			if strand == "+":

				# (inc1)[exc1], (inc2)[exc2]
				# (x1, y1)[last_exon_proximal_start, last_exon_proximal_end], (x2, y2)[last_exon_distal_start, last_exon_distal_end]
				# inc1: (x1 = penultimate_exon_end, y1 = last_exon_proximal_start)
				# exc1: [last_exon_proximal_start, last_exon_proximal_end]
				# inc2: (x2 = penultimate_exon_end, y2 = last_exon_distal_start)
				# exc2: [last_exon_distal_start, last_exon_distal_end]

				# Check if the penultimate exons have the same end position
				if penultimate_exon_transcript1_end != penultimate_exon_transcript2_end:
					continue

				# Set proximal and distal last exons
				if (last_exon_transcript1_start < last_exon_transcript2_start) and (last_exon_transcript1_end < last_exon_transcript2_end):

					last_exon_proximal = last_exon_transcript1
					last_exon_proximal_start = last_exon_transcript1_start
					last_exon_proximal_end = last_exon_transcript1_end
					last_exon_distal = last_exon_transcript2
					last_exon_distal_start = last_exon_transcript2_start
					last_exon_distal_end = last_exon_transcript2_end

				elif (last_exon_transcript1_start > last_exon_transcript2_start) and (last_exon_transcript1_end > last_exon_transcript2_end):

					last_exon_proximal = last_exon_transcript2
					last_exon_proximal_start = last_exon_transcript2_start
					last_exon_proximal_end = last_exon_transcript2_end
					last_exon_distal = last_exon_transcript1
					last_exon_distal_start = last_exon_transcript1_start
					last_exon_distal_end = last_exon_transcript1_end

				else:

					continue

				penultimate_exon_end = penultimate_exon_transcript1_end
				exon_a = chr + ":" + last_exon_distal_start + "-" + last_exon_distal_end
				exon_b = chr + ":" + last_exon_proximal_start + "-" + last_exon_proximal_end
				intron_a = chr + ":" + penultimate_exon_end + "-" + last_exon_distal_start
				intron_b = chr + ":" + penultimate_exon_end + "-" + last_exon_proximal_start
				intron_c = chr + ":" + last_exon_proximal_end + "-" + last_exon_distal_start # Intron connecting the distal and proximal last exons

			else: # strand == "-"

				# [exc1](inc1), [exc2](inc2)
				# [last_exon_distal_start, last_exon_distal_end](x1, y1), [last_exon_proximal_start, last_exon_proximal_end](x2, y2)
				# exc1: [last_exon_distal_start, last_exon_distal_end]
				# inc1: (x1 = last_exon_distal_end, y1 = penultimate_exon_start)
				# exc2: [last_exon_proximal_start, last_exon_proximal_end]
				# inc2: (x2 = last_exon_proximal_end, y2 = penultimate_exon_start)

				# Check if the penultimate exons have the same start position
				if penultimate_exon_transcript1_start != penultimate_exon_transcript2_start:
					continue

				# Set distal and proximal last exons
				if (last_exon_transcript1_start < last_exon_transcript2_start) and (last_exon_transcript1_end < last_exon_transcript2_end):

					last_exon_distal = last_exon_transcript1
					last_exon_distal_start = last_exon_transcript1_start
					last_exon_distal_end = last_exon_transcript1_end
					last_exon_proximal = last_exon_transcript2
					last_exon_proximal_start = last_exon_transcript2_start
					last_exon_proximal_end = last_exon_transcript2_end

				elif (last_exon_transcript1_start > last_exon_transcript2_start) and (last_exon_transcript1_end > last_exon_transcript2_end):

					last_exon_distal = last_exon_transcript2
					last_exon_distal_start = last_exon_transcript2_start
					last_exon_distal_end = last_exon_transcript2_end
					last_exon_proximal = last_exon_transcript1
					last_exon_proximal_start = last_exon_transcript1_start
					last_exon_proximal_end = last_exon_transcript1_end

				else:

					continue

				penultimate_exon_start = penultimate_exon_transcript1_start
				exon_a = chr + ":" + last_exon_distal_start + "-" + last_exon_distal_end
				exon_b = chr + ":" + last_exon_proximal_start + "-" + last_exon_proximal_end
				intron_a = chr + ":" + last_exon_distal_end + "-" + penultimate_exon_start
				intron_b = chr + ":" + last_exon_proximal_end + "-" + penultimate_exon_start
				intron_c = chr + ":" + last_exon_distal_end + "-" + last_exon_proximal_start # Intron connecting the distal and proximal last exons

			# Check if no intron connecting the distal exon and the proximal exon
			if intron_c in intron_list:
				continue

			event_l += [[exon_a, exon_b, intron_a, intron_b, strand, gene, gene_name]]

	return(event_l)


def mxe(gtf_dic) -> list:
	"""
	Make mutually exclusive exons list.

	Args:
		gtf_dic: A dictionary containing information about the GTF file.

	Returns:
		list: List of mutually exclusive exons events, where each event is represented as a list of the form
		[exon_a, exon_b, intron_a1, intron_a2, intron_b1, intron_b2, strand, gene, gene_name].

	"""

	event_l = []

	for gene in gtf_dic.keys():

		if "intron_list" not in gtf_dic[gene]:
			continue

		chr = gtf_dic[gene]["chr"]
		strand = gtf_dic[gene]["strand"]
		gene_name = gtf_dic[gene]["gene_name"]
		gene_start_values = gtf_dic[gene]["start"]
		gene_end_values = gtf_dic[gene]["end"]
		intron_list = gtf_dic[gene]["intron_list"]
		intron_start_dic = gtf_dic[gene]["intron_start_dic"]
		intron_end_dic = gtf_dic[gene]["intron_end_dic"]
		intron_dic = gtf_dic[gene]["transcript_intron_dic"]
		exon_dic = gtf_dic[gene]["transcript_exon_dic"]

		exon_list = gtf_dic[gene]["exon_list"]
		exon_list_unique = np.unique([[i.split(":")[1].split("-")[0], i.split(":")[1].split("-")[1]] for i in exon_list], axis = 0)
		exon_start = np.array([i[0] for i in exon_list_unique]).astype("int32")
		exon_end = np.array([i[1] for i in exon_list_unique]).astype("int32")

		idx1_idx2_iter = itertools.combinations(range(len(exon_start)), 2)
		for idx1, idx2 in idx1_idx2_iter:

			retained_intron = chr + ":" + str(exon_start[idx1]) + "-" + str(exon_end[idx2])

			# exon_a is upstream of exon_b
			# Not retained intron
			if (exon_end[idx1] < exon_start[idx2]) and (retained_intron not in exon_list) and (str(exon_start[idx1]) in intron_end_dic) and (str(exon_end[idx1]) in intron_start_dic) and (str(exon_start[idx2]) in intron_end_dic) and (str(exon_end[idx2]) in intron_start_dic):

				intron_a1_start_list = intron_end_dic[str(exon_start[idx1])]
				intron_a1_end = exon_start[idx1]
				intron_a2_start = exon_end[idx1]
				intron_a2_end_list = intron_start_dic[str(exon_end[idx1])]
				intron_b1_start_list = intron_end_dic[str(exon_start[idx2])]
				intron_b1_end = exon_start[idx2]
				intron_b2_start = exon_end[idx2]
				intron_b2_end_list = intron_start_dic[str(exon_end[idx2])]

				intron_iter = itertools.product(intron_a1_start_list, intron_a2_end_list, intron_b1_start_list, intron_b2_end_list)
				for intron_a1_start, intron_a2_end, intron_b1_start, intron_b2_end in intron_iter:

					exon_a = chr + ":" + str(exon_start[idx1]) + "-" + str(exon_end[idx1])
					exon_b = chr + ":" + str(exon_start[idx2]) + "-" + str(exon_end[idx2])
					intron_a1 = chr + ":" + str(intron_a1_start) + "-" + str(intron_a1_end)
					intron_a2 = chr + ":" + str(intron_a2_start) + "-" + str(intron_a2_end)
					intron_b1 = chr + ":" + str(intron_b1_start) + "-" + str(intron_b1_end)
					intron_b2 = chr + ":" + str(intron_b2_start) + "-" + str(intron_b2_end)
					intron_c = chr + ":" + str(intron_a2_start) + "-" + str(intron_b1_end)
					intron_d = chr + ":" + str(intron_a1_start) + "-" + str(intron_b2_end)

					if (int(intron_a1_start) == int(intron_b1_start)) and (int(intron_a1_end) != int(intron_b1_end)) and (int(intron_a2_start) != int(intron_b2_start)) and (int(intron_a2_end) == int(intron_b2_end)) and (intron_c not in intron_list) and (intron_d not in intron_list):

						key1_key2_iter = itertools.permutations(intron_dic.keys(), 2)
						for key1, key2 in key1_key2_iter:

							if (intron_a1 in intron_dic[key1]) and (intron_a2 in intron_dic[key1]) and (intron_b1 in intron_dic[key2]) and (intron_b2 in intron_dic[key2]):

								# exons not present in the same transcript
								flag = False
								for key3 in exon_dic.keys():

									if (exon_a in exon_dic[key3]) and (exon_b in exon_dic[key3]):

										flag = True
										break

								if flag == False:

									event_l += [[exon_a, exon_b, intron_a1, intron_a2, intron_b1, intron_b2, strand, gene, gene_name]]

	return(event_l)


def ri(gtf_dic) -> list:
	"""
	Make retained introns list.

	Args:
		gtf_dic: A dictionary containing information about the GTF file.

	Returns:
		list: List of retained introns events, where each event is represented as a list of the form
		[exon_a, exon_b, exon_c, intron_a, strand, gene, gene_name].

	"""

	event_l = []

	for gene in gtf_dic.keys():

		if "intron_list" not in gtf_dic[gene]:
			continue

		chr = gtf_dic[gene]["chr"]
		strand = gtf_dic[gene]["strand"]
		gene_name = gtf_dic[gene]["gene_name"]
		intron_list = gtf_dic[gene]["intron_list"]
		exon_dic = gtf_dic[gene]["transcript_exon_dic"]
		exon_list = gtf_dic[gene]["exon_list"]
		exon_list_unique = np.unique([[i.split(":")[1].split("-")[0], i.split(":")[1].split("-")[1]] for i in exon_list], axis = 0)
		exon_start = np.array([i[0] for i in exon_list_unique]).astype("int32")
		exon_end = np.array([i[1] for i in exon_list_unique]).astype("int32")

		idx1_idx2_iter = itertools.combinations(range(len(exon_start)), 2)
		for idx1, idx2 in idx1_idx2_iter:

			# Retained intron
			retained_intron = chr + ":" + str(exon_end[idx1]) + "-" + str(exon_start[idx2])
			retained_exon = chr + ":" + str(exon_start[idx1]) + "-" + str(exon_end[idx2])
			exon_a = chr + ":" + str(exon_start[idx1]) + "-" + str(exon_end[idx1])
			exon_b = chr + ":" + str(exon_start[idx2]) + "-" + str(exon_end[idx2])
			exon_c = retained_exon
			intron_a = chr + ":" + str(exon_end[idx1]) + "-" + str(exon_start[idx2])

			# exon_a is upstream of exon_b
			if (exon_end[idx1] < exon_start[idx2]) and (retained_intron in intron_list) and (retained_exon in exon_list):

				# exons present in the same transcript
				for exon_key in exon_dic.keys():

					if (exon_a in exon_dic[exon_key]) and (exon_b in exon_dic[exon_key]):

						event_l += [[exon_a, exon_b, exon_c, intron_a, strand, gene, gene_name]]

						break

	return(event_l)


def main():
	## Main

	args = get_args()

	gtf_path = args.gtf
	reference_gtf_path = args.reference_gtf
	num_process = args.num_process
	output_dir = args.output

	print("Loading " + str(gtf_path) + "....", file = sys.stdout)
	# start_time = time.time()
	gtf_dic_split = gtf(gtf_path, num_process)
	# end_time = time.time()
	# print("Loading time: " + str(end_time - start_time) + " seconds", file = sys.stdout)

	if reference_gtf_path:

		print("Loading " + str(reference_gtf_path) + "....", file = sys.stdout)
		gtf_ref_exon_set = gtf_exon_set(reference_gtf_path)
		gtf_ref_dic = gtf(reference_gtf_path, 1)
		# Only intron_list is needed
		gtf_ref_intron_set_dict = {k: v["intron_list"] for k, v in gtf_ref_dic[0].items() if "intron_list" in v}
		gtf_ref_intron_set = set()
		for k in gtf_ref_intron_set_dict:
			gtf_ref_intron_set |= gtf_ref_intron_set_dict[k]

	# Skipped exon
	# start_time = time.time()
	print("Searching skipped exon....", file = sys.stdout)
	with concurrent.futures.ProcessPoolExecutor(max_workers=num_process) as executor:

		futures = [executor.submit(se, gtf_dic_split[i]) for i in range(num_process)]

	output_l = []

	for future in concurrent.futures.as_completed(futures):

		output_l += future.result()

	output_df = pd.DataFrame(

		output_l,
		columns = ["exon", "intron_a", "intron_b", "intron_c", "strand", "gene_id", "gene_name"]

	)

	output_df["pos_id"] = \
		output_df["exon"].str.split(":", expand = True)[0].astype(str) + "@" + \
		output_df["exon"].str.split(":", expand = True)[1].str.split("-", expand = True)[0].astype(str) + "-" + output_df["exon"].str.split(":", expand = True)[1].str.split("-", expand = True)[1].astype(str) + "@" + \
		output_df["intron_c"].str.split(":", expand = True)[1].str.split("-", expand = True)[0].astype(str) + "-" + output_df["intron_c"].str.split(":", expand = True)[1].str.split("-", expand = True)[1].astype(str)
	output_df = output_df.sort_values("exon")
	output_df = output_df.drop_duplicates(subset = "pos_id", keep = "first")
	output_df = output_df.reset_index()
	output_df["event_id_num"] = output_df.index + 1
	output_df["event_id"] = "SE_" + output_df["event_id_num"].astype(str)
	output_df = output_df[["event_id", "pos_id", "exon", "intron_a", "intron_b", "intron_c", "strand", "gene_id", "gene_name"]]

	if reference_gtf_path:

		output_df["label"] = output_df.apply(lambda x: "annotated" if (x["intron_a"] in gtf_ref_intron_set) and (x["intron_b"] in gtf_ref_intron_set) and (x["intron_c"] in gtf_ref_intron_set) else "unannotated", axis = 1)

	else:

		output_df["label"] = "annotated"

	SE_output_df = output_df.copy()
	del output_df

	# end_time = time.time()
	# print("Skipped exon search time: " + str(end_time - start_time) + " seconds", file = sys.stdout)

	# Alternative five prime ss
	# start_time = time.time()
	print("Searching alternative five prime ss....", file = sys.stdout)
	with concurrent.futures.ProcessPoolExecutor(max_workers=num_process) as executor:

		futures = [executor.submit(five, gtf_dic_split[i]) for i in range(num_process)]

	output_l = []

	for future in concurrent.futures.as_completed(futures):

		output_l += future.result()

	output_df = pd.DataFrame(

		output_l,
		columns = ["exon_a", "exon_b", "intron_a", "intron_b", "strand", "gene_id", "gene_name"]

	)

	output_df["pos_id"] = \
		output_df["intron_a"].str.split(":", expand = True)[0].astype(str) + "@" + \
		output_df["intron_a"].str.split(":", expand = True)[1].str.split("-", expand = True)[0].astype(str) + "-" + output_df["intron_a"].str.split(":", expand = True)[1].str.split("-", expand = True)[1].astype(str) + "@" + \
		output_df["intron_b"].str.split(":", expand = True)[1].str.split("-", expand = True)[0].astype(str) + "-" + output_df["intron_b"].str.split(":", expand = True)[1].str.split("-", expand = True)[1].astype(str)
	output_df = output_df.sort_values("exon_a")
	output_df = output_df.drop_duplicates(subset = "pos_id", keep = "first")
	output_df = output_df.reset_index()
	output_df["event_id_num"] = output_df.index + 1
	output_df["event_id"] = "FIVE_" + output_df["event_id_num"].astype(str)
	output_df = output_df[["event_id", "pos_id", "exon_a", "exon_b", "intron_a", "intron_b", "strand", "gene_id", "gene_name"]]

	if reference_gtf_path:

		output_df["label"] = output_df.apply(lambda x: "annotated" if (x["intron_a"] in gtf_ref_intron_set) and (x["intron_b"] in gtf_ref_intron_set) else "unannotated", axis = 1)

	else:

		output_df["label"] = "annotated"

	FIVE_output_df = output_df.copy()
	del output_df

	# end_time = time.time()
	# print("Alternative five prime ss search time: " + str(end_time - start_time) + " seconds", file = sys.stdout)

	# Alternative three prime ss
	# start_time = time.time()
	print("Searching alternative three prime ss....", file = sys.stdout)
	with concurrent.futures.ProcessPoolExecutor(max_workers=num_process) as executor:

		futures = [executor.submit(three, gtf_dic_split[i]) for i in range(num_process)]

	output_l = []

	for future in concurrent.futures.as_completed(futures):

		output_l += future.result()

	output_df = pd.DataFrame(

		output_l,
		columns = ["exon_a", "exon_b", "intron_a", "intron_b", "strand", "gene_id", "gene_name"]

	)

	output_df["pos_id"] = \
		output_df["intron_a"].str.split(":", expand = True)[0].astype(str) + "@" + \
		output_df["intron_a"].str.split(":", expand = True)[1].str.split("-", expand = True)[0].astype(str) + "-" + output_df["intron_a"].str.split(":", expand = True)[1].str.split("-", expand = True)[1].astype(str) + "@" + \
		output_df["intron_b"].str.split(":", expand = True)[1].str.split("-", expand = True)[0].astype(str) + "-" + output_df["intron_b"].str.split(":", expand = True)[1].str.split("-", expand = True)[1].astype(str)
	output_df = output_df.sort_values("exon_a")
	output_df = output_df.drop_duplicates(subset = "pos_id", keep = "first")
	output_df = output_df.reset_index()
	output_df["event_id_num"] = output_df.index + 1
	output_df["event_id"] = "THREE_" + output_df["event_id_num"].astype(str)
	output_df = output_df[["event_id", "pos_id", "exon_a", "exon_b", "intron_a", "intron_b", "strand", "gene_id", "gene_name"]]

	if reference_gtf_path:

		output_df["label"] = output_df.apply(lambda x: "annotated" if (x["intron_a"] in gtf_ref_intron_set) and (x["intron_b"] in gtf_ref_intron_set) else "unannotated", axis = 1)

	else:

		output_df["label"] = "annotated"

	THREE_output_df = output_df.copy()
	del output_df

	# end_time = time.time()
	# print("Alternative three prime ss search time: " + str(end_time - start_time) + " seconds", file = sys.stdout)

	# Mutually exclusive exon
	# start_time = time.time()
	print("Searching mutually exclusive exons....", file = sys.stdout)
	with concurrent.futures.ProcessPoolExecutor(max_workers=num_process) as executor:

		futures = [executor.submit(mxe, gtf_dic_split[i]) for i in range(num_process)]

	output_l = []

	for future in concurrent.futures.as_completed(futures):

		output_l += future.result()

	output_df = pd.DataFrame(

		output_l,
		columns = ["exon_a", "exon_b", "intron_a1", "intron_a2", "intron_b1", "intron_b2", "strand", "gene_id", "gene_name"]

	)

	output_df["pos_id"] = \
		output_df["intron_a1"].str.split(":", expand = True)[0].astype(str) + "@" + \
		output_df["intron_a1"].str.split(":", expand = True)[1].str.split("-", expand = True)[0].astype(str) + "@" + \
		output_df["exon_a"].str.split(":", expand = True)[1].str.split("-", expand = True)[0].astype(str) + "-" + output_df["exon_a"].str.split(":", expand = True)[1].str.split("-", expand = True)[1].astype(str) + "@" + \
		output_df["exon_b"].str.split(":", expand = True)[1].str.split("-", expand = True)[0].astype(str) + "-" + output_df["exon_b"].str.split(":", expand = True)[1].str.split("-", expand = True)[1].astype(str) + "@" + \
		output_df["intron_b2"].str.split(":", expand = True)[1].str.split("-", expand = True)[1].astype(str)
	output_df = output_df.sort_values("exon_a")
	output_df = output_df.drop_duplicates(subset = "pos_id", keep = "first")
	output_df = output_df.reset_index()
	output_df["event_id_num"] = output_df.index + 1
	output_df["event_id"] = "MXE_" + output_df["event_id_num"].astype(str)
	output_df = output_df[["event_id", "pos_id", "exon_a", "exon_b", "intron_a1", "intron_a2", "intron_b1", "intron_b2", "strand", "gene_id", "gene_name"]]

	if reference_gtf_path:

		output_df["label"] = output_df.apply(lambda x: "annotated" if (x["intron_a1"] in gtf_ref_intron_set) and (x["intron_a2"] in gtf_ref_intron_set) and (x["intron_b1"] in gtf_ref_intron_set) and (x["intron_b2"] in gtf_ref_intron_set) else "unannotated", axis = 1)

	else:

		output_df["label"] = "annotated"

	MXE_output_df = output_df.copy()
	del output_df

	# end_time = time.time()
	# print("Mutually exclusive exon search time: " + str(end_time - start_time) + " seconds", file = sys.stdout)

	# Retained intron
	# start_time = time.time()
	print("Searching retained intron....", file = sys.stdout)
	with concurrent.futures.ProcessPoolExecutor(max_workers=num_process) as executor:

		futures = [executor.submit(ri, gtf_dic_split[i]) for i in range(num_process)]

	output_l = []

	for future in concurrent.futures.as_completed(futures):

		output_l += future.result()

	output_df = pd.DataFrame(

		output_l,
		columns = ["exon_a", "exon_b", "exon_c", "intron_a", "strand", "gene_id", "gene_name"]

	)

	output_df["pos_id"] = \
		output_df["intron_a"].str.replace(":", "@")
	output_df = output_df.sort_values("exon_a")
	output_df = output_df.drop_duplicates(subset = "pos_id", keep = "first")
	output_df = output_df.reset_index()
	output_df["event_id_num"] = output_df.index + 1
	output_df["event_id"] = "RI_" + output_df["event_id_num"].astype(str)
	output_df = output_df[["event_id", "pos_id", "exon_a", "exon_b", "exon_c", "intron_a", "strand", "gene_id", "gene_name"]]

	if reference_gtf_path:

		output_df["label"] = output_df.apply(lambda x: "annotated" if (x["intron_a"] in gtf_ref_intron_set) and (x["exon_c"] in gtf_ref_exon_set) else "unannotated", axis = 1)

	else:

		output_df["label"] = "annotated"

	RI_output_df = output_df.copy()
	del output_df

	# end_time = time.time()
	# print("Retained intron search time: " + str(end_time - start_time) + " seconds", file = sys.stdout)

	# Multiple skipped exons
	# start_time = time.time()
	print("Searching multiple skipped exons....", file = sys.stdout)
	with concurrent.futures.ProcessPoolExecutor(max_workers=num_process) as executor:

		futures = [executor.submit(mse, gtf_dic_split[i]) for i in range(num_process)]

	output_l = []

	for future in concurrent.futures.as_completed(futures):

		output_l += future.result()

	output_df = pd.DataFrame(

		output_l,
		columns = ["exon", "intron", "mse_n", "strand", "gene_id", "gene_name"]

	)

	# pos_id = chromosome@exon_start-exon_end;exon_start-exon_end@exclusionintron_start-exclusionintron_end
	output_df["chr"] = output_df["exon"].str.split(":", expand = True)[0]
	output_df["exon_for_posid"] = output_df.apply(lambda x: x["exon"].replace(x["chr"] + ":", ""), axis = 1)
	output_df["exc"] = output_df["intron"].apply(lambda x: x.split(";")[-1])
	output_df["pos_id"] = \
		output_df["chr"] + "@" + \
		output_df["exon_for_posid"] + "@" + \
		output_df["exc"].str.split(":", expand = True)[1].str.split("-", expand = True)[0].astype(str) + "-" + output_df["exc"].str.split(":", expand = True)[1].str.split("-", expand = True)[1].astype(str)
	output_df = output_df.sort_values("exon")
	output_df = output_df.drop_duplicates(subset = "pos_id", keep = "first")
	output_df = output_df.reset_index()
	output_df["event_id_num"] = output_df.index + 1
	output_df["event_id"] = "MSE_" + output_df["event_id_num"].astype(str)
	output_df = output_df[["event_id", "pos_id", "mse_n", "exon", "intron", "strand", "gene_id", "gene_name"]]

	# Check if the intron is annotated
	if reference_gtf_path:

		output_df["label"] = output_df["intron"].apply(lambda x: "annotated" if set(x.split(";")) <= gtf_ref_intron_set else "unannotated")

	else:

		output_df["label"] = "annotated"

	MSE_output_df = output_df.copy()
	del output_df

	# end_time = time.time()
	# print("Multiple skipped exons search time: " + str(end_time - start_time) + " seconds", file = sys.stdout)

	# Alternative first exons
	# start_time = time.time()
	print("Searching alternative first exons....", file = sys.stdout)

	with concurrent.futures.ProcessPoolExecutor(max_workers=num_process) as executor:

		futures = [executor.submit(afe, gtf_dic_split[i]) for i in range(num_process)]

	output_l = []

	for future in concurrent.futures.as_completed(futures):

		output_l += future.result()

	output_df = pd.DataFrame(

		output_l,
		columns = ["exon_a", "exon_b", "intron_a", "intron_b", "strand", "gene_id", "gene_name"]

	)

	# pos_id = chromosome@exon_start-exon_end;exon_start-exon_end@intron_start-intron_end;intron_start-intron_end
	output_df["chr"] = output_df["exon_a"].str.split(":", expand = True)[0]
	output_df["exon_for_posid"] = output_df.apply(lambda x: x["exon_a"].replace(x["chr"] + ":", "") + ";" + x["exon_b"].replace(x["chr"] + ":", ""), axis = 1)
	output_df["intron_for_posid"] = output_df.apply(lambda x: x["intron_a"].replace(x["chr"] + ":", "") + ";" + x["intron_b"].replace(x["chr"] + ":", ""), axis = 1)
	output_df["pos_id"] = \
		output_df["chr"] + "@" + \
		output_df["exon_for_posid"] + "@" + \
		output_df["intron_for_posid"]
	output_df = output_df.sort_values(["exon_a", "exon_b"], ascending = [True, True])
	output_df = output_df.drop_duplicates(subset = "intron_for_posid", keep = "first")
	output_df = output_df.reset_index()
	output_df["event_id_num"] = output_df.index + 1
	output_df["event_id"] = "AFE_" + output_df["event_id_num"].astype(str)
	output_df = output_df[["event_id", "pos_id", "exon_a", "exon_b", "intron_a", "intron_b", "strand", "gene_id", "gene_name"]]

	# Check if the intron is annotated
	if reference_gtf_path:

		output_df["label"] = output_df.apply(lambda x: "annotated" if (x["intron_a"] in gtf_ref_intron_set) and (x["intron_b"] in gtf_ref_intron_set) else "unannotated", axis = 1)

	else:

		output_df["label"] = "annotated"

	AFE_output_df = output_df.copy()
	del output_df

	# end_time = time.time()
	# print("Alternative first exons search time: " + str(end_time - start_time) + " seconds", file = sys.stdout)

	# Alternative last exons
	# start_time = time.time()
	print("Searching alternative last exons....", file = sys.stdout)

	with concurrent.futures.ProcessPoolExecutor(max_workers=num_process) as executor:

		futures = [executor.submit(ale, gtf_dic_split[i]) for i in range(num_process)]

	output_l = []

	for future in concurrent.futures.as_completed(futures):

		output_l += future.result()

	output_df = pd.DataFrame(

		output_l,
		columns = ["exon_a", "exon_b", "intron_a", "intron_b", "strand", "gene_id", "gene_name"]

	)

	# pos_id = chromosome@exon_start-exon_end;exon_start-exon_end@intron_start-intron_end;intron_start-intron_end
	output_df["chr"] = output_df["exon_a"].str.split(":", expand = True)[0]
	output_df["exon_for_posid"] = output_df.apply(lambda x: x["exon_a"].replace(x["chr"] + ":", "") + ";" + x["exon_b"].replace(x["chr"] + ":", ""), axis = 1)
	output_df["intron_for_posid"] = output_df.apply(lambda x: x["intron_a"].replace(x["chr"] + ":", "") + ";" + x["intron_b"].replace(x["chr"] + ":", ""), axis = 1)
	output_df["pos_id"] = \
		output_df["chr"] + "@" + \
		output_df["exon_for_posid"] + "@" + \
		output_df["intron_for_posid"]
	output_df = output_df.sort_values(["exon_a", "exon_b"], ascending = [True, True])
	output_df = output_df.drop_duplicates(subset = "intron_for_posid", keep = "first")
	output_df = output_df.reset_index()
	output_df["event_id_num"] = output_df.index + 1
	output_df["event_id"] = "ALE_" + output_df["event_id_num"].astype(str)
	output_df = output_df[["event_id", "pos_id", "exon_a", "exon_b", "intron_a", "intron_b", "strand", "gene_id", "gene_name"]]

	# Check if the intron is annotated
	if reference_gtf_path:

		output_df["label"] = output_df.apply(lambda x: "annotated" if (x["intron_a"] in gtf_ref_intron_set) and (x["intron_b"] in gtf_ref_intron_set) else "unannotated", axis = 1)

	else:

		output_df["label"] = "annotated"

	ALE_output_df = output_df.copy()
	del output_df

	# end_time = time.time()
	# print("Alternative last exons search time: " + str(end_time - start_time) + " seconds", file = sys.stdout)

	# Export
	os.makedirs(output_dir, exist_ok = True)
	SE_output_df.to_csv(

		output_dir + "/EVENT_SE.txt",
		sep = "\t",
		index = False

	)

	FIVE_output_df.to_csv(

		output_dir + "/EVENT_FIVE.txt",
		sep = "\t",
		index = False

	)

	THREE_output_df.to_csv(

		output_dir + "/EVENT_THREE.txt",
		sep = "\t",
		index = False

	)

	MXE_output_df.to_csv(

		output_dir + "/EVENT_MXE.txt",
		sep = "\t",
		index = False

	)

	RI_output_df.to_csv(

		output_dir + "/EVENT_RI.txt",
		sep = "\t",
		index = False

	)

	MSE_output_df.to_csv(

		output_dir + "/EVENT_MSE.txt",
		sep = "\t",
		index = False

	)

	AFE_output_df.to_csv(

		output_dir + "/EVENT_AFE.txt",
		sep = "\t",
		index = False

	)

	ALE_output_df.to_csv(

		output_dir + "/EVENT_ALE.txt",
		sep = "\t",
		index = False

	)

	print("Done!", file = sys.stdout)

if __name__ == '__main__':

	main()

