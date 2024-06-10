import warnings
warnings.simplefilter('ignore')
import argparse
import sys
import os
import pandas as pd
import numpy as np
import plotly.express as px

def get_args():

	parser = argparse.ArgumentParser(
		formatter_class = argparse.ArgumentDefaultsHelpFormatter,
		description = "Make plots for alternative splicing events"
	)

	parser.add_argument("-i", "--input", type = str, help = "Directory that contains result files")
	parser.add_argument("-e", "--experiment-table", type = str, help = "Experiment table file")
	parser.add_argument("-o", "--output", type = str, help = "Directory for output files")

	args = parser.parse_args()

	return(args)

def load_experiment_table(experiment_table: str):

	# Load experiment table
	experiment_table_df = pd.read_csv(

		experiment_table,
		sep = "\t",
		usecols = ["sample", "group"]

	)

	return experiment_table_df

def load_tpm_pca_table(input_dir: str, experiment_table_df: pd.DataFrame, output_dir: str):

	# Load PCA matrix for TPM
	pca_tpm_df = pd.read_csv(

		input_dir + "/pca/tpm_pca.tsv",
		sep = "\t"

	)
	pca_tpm_df = pca_tpm_df.rename(columns = {"Unnamed: 0": "sample"})
	pca_tpm_df = pd.merge(pca_tpm_df, experiment_table_df, on = "sample")
	pca_tpm_df = pca_tpm_df.sort_values("group")
	# Load contribution table for TPM
	contribution_tpm_df = pd.read_csv(

		input_dir + "/pca/tpm_contribution.tsv",
		sep = "\t",
		names = ["PC", "contribution"]

	)
	contribution_tpm_PC1 = str((contribution_tpm_df.iloc[0][1]*100).round(2))
	contribution_tpm_PC2 = str((contribution_tpm_df.iloc[1][1]*100).round(2))

	return pca_tpm_df, contribution_tpm_PC1, contribution_tpm_PC2

def load_psi_pca_table(input_dir: str, experiment_table_df: pd.DataFrame, output_dir: str):

	# Load PCA matrix for PSI
	pca_psi_df = pd.read_csv(

		input_dir + "/pca/psi_pca.tsv",
		sep = "\t"

	)
	pca_psi_df = pca_psi_df.rename(columns = {"Unnamed: 0": "sample"})
	pca_psi_df = pd.merge(pca_psi_df, experiment_table_df, on = "sample")
	pca_psi_df = pca_psi_df.sort_values("group")
	# Load contribution table for PSI
	contribution_psi_df = pd.read_csv(

		input_dir + "/pca/psi_contribution.tsv",
		sep = "\t",
		names = ["PC", "contribution"]

	)
	contribution_psi_PC1 = str((contribution_psi_df.iloc[0][1]*100).round(2))
	contribution_psi_PC2 = str((contribution_psi_df.iloc[1][1]*100).round(2))

	return pca_psi_df, contribution_psi_PC1, contribution_psi_PC2

def plots_pca(name: str, pca_df: pd.DataFrame, contribution_PC1: str, contribution_PC2: str, output_dir: str):

	fig = px.scatter(

		pca_df,
		x = "PC1",
		y = "PC2",
		color = "group",
		opacity = 0.5,
		hover_data = ["sample"]

	)

	fig.update_traces(
		marker = dict(
			size = 8,
			line = dict(width = 0, color = 'DarkSlateGrey')),
			selector = dict(mode = 'markers')
	)

	fig.update_layout(
		title = dict(
			text = "PCA for " + name,
			font = dict(size = 26, color = 'black'),
			xref = 'paper',
			x = 0.5,
			y = 0.95,
			xanchor = 'center',
		)
	)

	fig.update_layout(
		width = 550,
		height = 400,
		font_family = "Arial",
		xaxis_title = "PC1 ({}%)".format(contribution_PC1),
		yaxis_title = "PC2 ({}%)".format(contribution_PC2),
		legend_title = "Group",
	)

	fig.write_html(output_dir + "/data/pca_" + name + ".html")

def plots(AS: str, input_dir: str, output_dir: str):

	# load data
	df = pd.read_csv(

		input_dir + "/splicing/PSI_" + AS + ".txt",
		sep = "\t"

	)

	if AS == "FIVE" or AS == "THREE":

		Ref_group = list(df.columns)[12]
		Exp_group = list(df.columns)[15]

	elif AS == "MXE":

		Ref_group = list(df.columns)[16]
		Exp_group = list(df.columns)[21]

	else:

		Ref_group = list(df.columns)[13]
		Exp_group = list(df.columns)[17]

	# Volcano plot
	if not df.empty:

		df['-log10(q)'] = -np.log10(df["q"])
		df.loc[(df["Diff events"] == "Yes") & (df["dPSI"] > 0.1), "group"] = "up"
		df.loc[(df["Diff events"] == "Yes") & (df["dPSI"] < -0.1), "group"] = "down"
		df = df.fillna({"group": "others"})

		fig = px.scatter(

			df,
			x = "dPSI",
			y = "-log10(q)",
			color = "group",
			symbol="label",
			opacity = 0.5,
			category_orders = {"group": ["up", "down", "others"], "label": ["annotated", "unannotated"]},
			color_discrete_sequence = ["salmon", "steelblue", "lightgrey"],
			hover_data = ["gene_name", "event_id", Ref_group, Exp_group, "q"]

		)

		fig.update_traces(
			marker = dict(
				size = 8,
				line = dict(width = 0, color = 'DarkSlateGrey')),
				selector = dict(mode = 'markers')
		)

		fig.update_layout(
			title = dict(
				text = "Volcano plot",
				font = dict(size = 26, color = 'black'),
				xref = 'paper',
				x = 0.5,
				y = 0.95,
				xanchor = 'center'
			)
		)

	else:

		# If there is no data, make empty plot
		fig = px.scatter(

			x = [None],
			y = [None]

		)

		# Remove xlabel and ylabel
		fig.update_xaxes(title = None)
		fig.update_yaxes(title = None)

		# Add "No data" text to the empty plot
		fig.add_annotation(

			text = "No data",
			xref = "paper",
			yref = "paper",
			x = 0.5,
			y = 0.5,
			showarrow = False,
			font = dict(size = 26, color = 'black')

		)

	fig.update_layout(
		width = 550,
		height = 400,
		font_family = "Arial",
	)

	fig.write_html(output_dir + "/data/volcano_" + AS + ".html")

	# Scatter
	if not df.empty:

		fig = px.scatter(

			df,
			x = Ref_group,
			y = Exp_group,
			color = "group",
			symbol="label",
			opacity = 0.5,
			category_orders = {"group": ["up", "down", "others"], "label": ["annotated", "unannotated"]},
			color_discrete_sequence = ["salmon", "steelblue", "lightgrey"],
			hover_data = ["gene_name", "event_id", "dPSI", "q"]

		)

		fig.update_traces(
			marker = dict(size=8,
							line=dict(width=0, color='DarkSlateGrey')),
							selector=dict(mode='markers')
		)

		fig.update_layout(title=dict(text = "Scatter plot",
										font=dict(size=26,
												color='black'),
										xref='paper',
										x=0.5,
										y=0.95,
										xanchor='center'
									)
							)

	else:

		# If there is no data, make empty plot
		fig = px.scatter(

			x = [None],
			y = [None]

		)

		# Remove xlabel and ylabel
		fig.update_xaxes(title = None)
		fig.update_yaxes(title = None)

		# Add "No data" text to the empty plot
		fig.add_annotation(

			text = "No data",
			xref = "paper",
			yref = "paper",
			x = 0.5,
			y = 0.5,
			showarrow = False,
			font = dict(size = 26, color = 'black')

		)

	fig.update_layout(
		width = 550,
		height= 400,
		font_family = "Arial",
	)

	fig.write_html(output_dir + "/data/scatter_" + AS + ".html")

	# Bar
	if not df.empty:

		count_df = df[df["group"] != "others"].groupby(["group", "label"], as_index = False).count()[["group", "label", "event_id"]]
		fig = px.bar(
			count_df,
			x = "group",
			y = "event_id",
			color = "label",
			labels = {"event_id": "Count"},
			category_orders = {"group": ["up", "down"], "label": ["annotated", "unannotated"]},
			color_discrete_sequence = ["#66c2a5", "#fc8d62"],
			barmode = "relative"
		)

		fig.update_layout(title=dict(text = "Number of events",
										font=dict(size=26,
												color='black'),
										xref='paper',
										x=0.5,
										y=0.95,
										xanchor='center'
									)
							)

	else:

		# if there is no data, make empty plot
		fig = px.bar(

			x = [0],
			y = [0]

		)

		# Remove xlabel and ylabel
		fig.update_xaxes(title = None)
		fig.update_yaxes(title = None)

		# Add "No data" text to the empty plot
		fig.add_annotation(

			text = "No data",
			xref = "paper",
			yref = "paper",
			x = 0.5,
			y = 0.5,
			showarrow = False,
			font = dict(size = 26, color = 'black')

		)

	fig.update_layout(
		width = 550,
		height = 400,
		font_family = "Arial",
	)

	fig.write_html(output_dir + "/data/bar_" + AS + ".html")

def write_summary_html(output_dir: str):

	# PCA
	lines_strip_pca_dict = {}
	for pca in ["TPM", "PSI"]:

		with open(output_dir + "/data/pca_" + pca + ".html", 'r') as f:

			lines = f.readlines()

		lines_strip = [line.strip() for line in lines[3:12]]
		lines_strip_pca_dict[pca] = lines_strip

	# Splicing events
	events = ["SE", "FIVE", "THREE", "MXE", "RI"]
	plottypes = ["volcano", "scatter", "bar"]
	lines_strip_dict = {}

	for event in events:

		lines_strip_dict[event] = {}

		for plottype in plottypes:

			with open(output_dir + "/data/" + plottype + "_" + event + ".html", 'r') as f:

				lines = f.readlines()

			lines_strip = [line.strip() for line in lines[3:12]]
			lines_strip_dict[event][plottype] = lines_strip

	summary_html = '''
	<!DOCTYPE html>
	<html lang="en">
	<head>
	<title>Shiba results summary</title>
	<meta http-equiv="Cache-Control" content="no-store">
	</head>
	<body>

		<font size="7">Shiba results summary</font>
		<p><font size="5">Principal component analysis</font></p>
		<TABLE WIDTH="1200" CELLSPACING="0" CELLPADDING="0">

			<TR>
			<TD WIDTH="600">

				{pca_tpm_l1}
				{pca_tpm_l2}
				{pca_tpm_l3}
				{pca_tpm_l4}
				{pca_tpm_l5}
				{pca_tpm_l6}
				{pca_tpm_l7}
				{pca_tpm_l8}
				{pca_tpm_l9}

			</TD>
			<TD WIDTH="600">

				{pca_psi_l1}
				{pca_psi_l2}
				{pca_psi_l3}
				{pca_psi_l4}
				{pca_psi_l5}
				{pca_psi_l6}
				{pca_psi_l7}
				{pca_psi_l8}
				{pca_psi_l9}

			</TD>
			</TR>

		</TABLE>

		<p><font size="5">Skipped exon</font></p>
		<TABLE WIDTH="1800" CELLSPACING="0" CELLPADDING="0">

			<TR>
			<TD WIDTH="600">

				{volcano_se_l1}
				{volcano_se_l2}
				{volcano_se_l3}
				{volcano_se_l4}
				{volcano_se_l5}
				{volcano_se_l6}
				{volcano_se_l7}
				{volcano_se_l8}
				{volcano_se_l9}

			</TD>
			<TD WIDTH="600">

				{scatter_se_l1}
				{scatter_se_l2}
				{scatter_se_l3}
				{scatter_se_l4}
				{scatter_se_l5}
				{scatter_se_l6}
				{scatter_se_l7}
				{scatter_se_l8}
				{scatter_se_l9}

			</TD>
			<TD WIDTH="600">

				{bar_se_l1}
				{bar_se_l2}
				{bar_se_l3}
				{bar_se_l4}
				{bar_se_l5}
				{bar_se_l6}
				{bar_se_l7}
				{bar_se_l8}
				{bar_se_l9}

			</TD>
			</TR>

		</TABLE>

		<p><font size="5">Alternative 5'ss</font></p>
		<TABLE WIDTH="1800" CELLSPACING="0" CELLPADDING="0">

			<TR>
			<TD WIDTH="600">

				{volcano_five_l1}
				{volcano_five_l2}
				{volcano_five_l3}
				{volcano_five_l4}
				{volcano_five_l5}
				{volcano_five_l6}
				{volcano_five_l7}
				{volcano_five_l8}
				{volcano_five_l9}

			</TD>
			<TD WIDTH="600">

				{scatter_five_l1}
				{scatter_five_l2}
				{scatter_five_l3}
				{scatter_five_l4}
				{scatter_five_l5}
				{scatter_five_l6}
				{scatter_five_l7}
				{scatter_five_l8}
				{scatter_five_l9}

			</TD>
			<TD WIDTH="600">

				{bar_five_l1}
				{bar_five_l2}
				{bar_five_l3}
				{bar_five_l4}
				{bar_five_l5}
				{bar_five_l6}
				{bar_five_l7}
				{bar_five_l8}
				{bar_five_l9}

			</TD>
			</TR>

		</TABLE>

		<p><font size="5">Alternative 3'ss</font></p>
		<TABLE WIDTH="1800" CELLSPACING="0" CELLPADDING="0">

			<TR>
			<TD WIDTH="600">

				{volcano_three_l1}
				{volcano_three_l2}
				{volcano_three_l3}
				{volcano_three_l4}
				{volcano_three_l5}
				{volcano_three_l6}
				{volcano_three_l7}
				{volcano_three_l8}
				{volcano_three_l9}

			</TD>
			<TD WIDTH="600">

				{scatter_three_l1}
				{scatter_three_l2}
				{scatter_three_l3}
				{scatter_three_l4}
				{scatter_three_l5}
				{scatter_three_l6}
				{scatter_three_l7}
				{scatter_three_l8}
				{scatter_three_l9}

			</TD>
			<TD WIDTH="600">

				{bar_three_l1}
				{bar_three_l2}
				{bar_three_l3}
				{bar_three_l4}
				{bar_three_l5}
				{bar_three_l6}
				{bar_three_l7}
				{bar_three_l8}
				{bar_three_l9}

			</TD>
			</TR>

		</TABLE>

		<p><font size="5">Mutually exclusive exons</font></p>
		<TABLE WIDTH="1800" CELLSPACING="0" CELLPADDING="0">

			<TR>
			<TD WIDTH="600">

				{volcano_mxe_l1}
				{volcano_mxe_l2}
				{volcano_mxe_l3}
				{volcano_mxe_l4}
				{volcano_mxe_l5}
				{volcano_mxe_l6}
				{volcano_mxe_l7}
				{volcano_mxe_l8}
				{volcano_mxe_l9}

			</TD>
			<TD WIDTH="600">

				{scatter_mxe_l1}
				{scatter_mxe_l2}
				{scatter_mxe_l3}
				{scatter_mxe_l4}
				{scatter_mxe_l5}
				{scatter_mxe_l6}
				{scatter_mxe_l7}
				{scatter_mxe_l8}
				{scatter_mxe_l9}

			</TD>
			<TD WIDTH="600">

				{bar_mxe_l1}
				{bar_mxe_l2}
				{bar_mxe_l3}
				{bar_mxe_l4}
				{bar_mxe_l5}
				{bar_mxe_l6}
				{bar_mxe_l7}
				{bar_mxe_l8}
				{bar_mxe_l9}

			</TD>
			</TR>

		</TABLE>

		<p><font size="5">Retained intron</font></p>
		<TABLE WIDTH="1800" CELLSPACING="0" CELLPADDING="0">

			<TR>
			<TD WIDTH="600">

				{volcano_ri_l1}
				{volcano_ri_l2}
				{volcano_ri_l3}
				{volcano_ri_l4}
				{volcano_ri_l5}
				{volcano_ri_l6}
				{volcano_ri_l7}
				{volcano_ri_l8}
				{volcano_ri_l9}

			</TD>
			<TD WIDTH="600">

				{scatter_ri_l1}
				{scatter_ri_l2}
				{scatter_ri_l3}
				{scatter_ri_l4}
				{scatter_ri_l5}
				{scatter_ri_l6}
				{scatter_ri_l7}
				{scatter_ri_l8}
				{scatter_ri_l9}

			</TD>
			<TD WIDTH="600">

				{bar_ri_l1}
				{bar_ri_l2}
				{bar_ri_l3}
				{bar_ri_l4}
				{bar_ri_l5}
				{bar_ri_l6}
				{bar_ri_l7}
				{bar_ri_l8}
				{bar_ri_l9}

			</TD>
			</TR>

		</TABLE>

	</body>
	</html>
	'''.format(
		pca_tpm_l1 = lines_strip_pca_dict["TPM"][0],
		pca_tpm_l2 = lines_strip_pca_dict["TPM"][1],
		pca_tpm_l3 = lines_strip_pca_dict["TPM"][2],
		pca_tpm_l4 = lines_strip_pca_dict["TPM"][3],
		pca_tpm_l5 = lines_strip_pca_dict["TPM"][4],
		pca_tpm_l6 = lines_strip_pca_dict["TPM"][5],
		pca_tpm_l7 = lines_strip_pca_dict["TPM"][6],
		pca_tpm_l8 = lines_strip_pca_dict["TPM"][7],
		pca_tpm_l9 = lines_strip_pca_dict["TPM"][8],
		pca_psi_l1 = lines_strip_pca_dict["PSI"][0],
		pca_psi_l2 = lines_strip_pca_dict["PSI"][1],
		pca_psi_l3 = lines_strip_pca_dict["PSI"][2],
		pca_psi_l4 = lines_strip_pca_dict["PSI"][3],
		pca_psi_l5 = lines_strip_pca_dict["PSI"][4],
		pca_psi_l6 = lines_strip_pca_dict["PSI"][5],
		pca_psi_l7 = lines_strip_pca_dict["PSI"][6],
		pca_psi_l8 = lines_strip_pca_dict["PSI"][7],
		pca_psi_l9 = lines_strip_pca_dict["PSI"][8],
		volcano_se_l1 = lines_strip_dict["SE"]["volcano"][0],
		volcano_se_l2 = lines_strip_dict["SE"]["volcano"][1],
		volcano_se_l3 = lines_strip_dict["SE"]["volcano"][2],
		volcano_se_l4 = lines_strip_dict["SE"]["volcano"][3],
		volcano_se_l5 = lines_strip_dict["SE"]["volcano"][4],
		volcano_se_l6 = lines_strip_dict["SE"]["volcano"][5],
		volcano_se_l7 = lines_strip_dict["SE"]["volcano"][6],
		volcano_se_l8 = lines_strip_dict["SE"]["volcano"][7],
		volcano_se_l9 = lines_strip_dict["SE"]["volcano"][8],
		scatter_se_l1 = lines_strip_dict["SE"]["scatter"][0],
		scatter_se_l2 = lines_strip_dict["SE"]["scatter"][1],
		scatter_se_l3 = lines_strip_dict["SE"]["scatter"][2],
		scatter_se_l4 = lines_strip_dict["SE"]["scatter"][3],
		scatter_se_l5 = lines_strip_dict["SE"]["scatter"][4],
		scatter_se_l6 = lines_strip_dict["SE"]["scatter"][5],
		scatter_se_l7 = lines_strip_dict["SE"]["scatter"][6],
		scatter_se_l8 = lines_strip_dict["SE"]["scatter"][7],
		scatter_se_l9 = lines_strip_dict["SE"]["scatter"][8],
		bar_se_l1 = lines_strip_dict["SE"]["bar"][0],
		bar_se_l2 = lines_strip_dict["SE"]["bar"][1],
		bar_se_l3 = lines_strip_dict["SE"]["bar"][2],
		bar_se_l4 = lines_strip_dict["SE"]["bar"][3],
		bar_se_l5 = lines_strip_dict["SE"]["bar"][4],
		bar_se_l6 = lines_strip_dict["SE"]["bar"][5],
		bar_se_l7 = lines_strip_dict["SE"]["bar"][6],
		bar_se_l8 = lines_strip_dict["SE"]["bar"][7],
		bar_se_l9 = lines_strip_dict["SE"]["bar"][8],
		volcano_five_l1 = lines_strip_dict["FIVE"]["volcano"][0],
		volcano_five_l2 = lines_strip_dict["FIVE"]["volcano"][1],
		volcano_five_l3 = lines_strip_dict["FIVE"]["volcano"][2],
		volcano_five_l4 = lines_strip_dict["FIVE"]["volcano"][3],
		volcano_five_l5 = lines_strip_dict["FIVE"]["volcano"][4],
		volcano_five_l6 = lines_strip_dict["FIVE"]["volcano"][5],
		volcano_five_l7 = lines_strip_dict["FIVE"]["volcano"][6],
		volcano_five_l8 = lines_strip_dict["FIVE"]["volcano"][7],
		volcano_five_l9 = lines_strip_dict["FIVE"]["volcano"][8],
		scatter_five_l1 = lines_strip_dict["FIVE"]["scatter"][0],
		scatter_five_l2 = lines_strip_dict["FIVE"]["scatter"][1],
		scatter_five_l3 = lines_strip_dict["FIVE"]["scatter"][2],
		scatter_five_l4 = lines_strip_dict["FIVE"]["scatter"][3],
		scatter_five_l5 = lines_strip_dict["FIVE"]["scatter"][4],
		scatter_five_l6 = lines_strip_dict["FIVE"]["scatter"][5],
		scatter_five_l7 = lines_strip_dict["FIVE"]["scatter"][6],
		scatter_five_l8 = lines_strip_dict["FIVE"]["scatter"][7],
		scatter_five_l9 = lines_strip_dict["FIVE"]["scatter"][8],
		bar_five_l1 = lines_strip_dict["FIVE"]["bar"][0],
		bar_five_l2 = lines_strip_dict["FIVE"]["bar"][1],
		bar_five_l3 = lines_strip_dict["FIVE"]["bar"][2],
		bar_five_l4 = lines_strip_dict["FIVE"]["bar"][3],
		bar_five_l5 = lines_strip_dict["FIVE"]["bar"][4],
		bar_five_l6 = lines_strip_dict["FIVE"]["bar"][5],
		bar_five_l7 = lines_strip_dict["FIVE"]["bar"][6],
		bar_five_l8 = lines_strip_dict["FIVE"]["bar"][7],
		bar_five_l9 = lines_strip_dict["FIVE"]["bar"][8],
		volcano_three_l1 = lines_strip_dict["THREE"]["volcano"][0],
		volcano_three_l2 = lines_strip_dict["THREE"]["volcano"][1],
		volcano_three_l3 = lines_strip_dict["THREE"]["volcano"][2],
		volcano_three_l4 = lines_strip_dict["THREE"]["volcano"][3],
		volcano_three_l5 = lines_strip_dict["THREE"]["volcano"][4],
		volcano_three_l6 = lines_strip_dict["THREE"]["volcano"][5],
		volcano_three_l7 = lines_strip_dict["THREE"]["volcano"][6],
		volcano_three_l8 = lines_strip_dict["THREE"]["volcano"][7],
		volcano_three_l9 = lines_strip_dict["THREE"]["volcano"][8],
		scatter_three_l1 = lines_strip_dict["THREE"]["scatter"][0],
		scatter_three_l2 = lines_strip_dict["THREE"]["scatter"][1],
		scatter_three_l3 = lines_strip_dict["THREE"]["scatter"][2],
		scatter_three_l4 = lines_strip_dict["THREE"]["scatter"][3],
		scatter_three_l5 = lines_strip_dict["THREE"]["scatter"][4],
		scatter_three_l6 = lines_strip_dict["THREE"]["scatter"][5],
		scatter_three_l7 = lines_strip_dict["THREE"]["scatter"][6],
		scatter_three_l8 = lines_strip_dict["THREE"]["scatter"][7],
		scatter_three_l9 = lines_strip_dict["THREE"]["scatter"][8],
		bar_three_l1 = lines_strip_dict["THREE"]["bar"][0],
		bar_three_l2 = lines_strip_dict["THREE"]["bar"][1],
		bar_three_l3 = lines_strip_dict["THREE"]["bar"][2],
		bar_three_l4 = lines_strip_dict["THREE"]["bar"][3],
		bar_three_l5 = lines_strip_dict["THREE"]["bar"][4],
		bar_three_l6 = lines_strip_dict["THREE"]["bar"][5],
		bar_three_l7 = lines_strip_dict["THREE"]["bar"][6],
		bar_three_l8 = lines_strip_dict["THREE"]["bar"][7],
		bar_three_l9 = lines_strip_dict["THREE"]["bar"][8],
		volcano_mxe_l1 = lines_strip_dict["MXE"]["volcano"][0],
		volcano_mxe_l2 = lines_strip_dict["MXE"]["volcano"][1],
		volcano_mxe_l3 = lines_strip_dict["MXE"]["volcano"][2],
		volcano_mxe_l4 = lines_strip_dict["MXE"]["volcano"][3],
		volcano_mxe_l5 = lines_strip_dict["MXE"]["volcano"][4],
		volcano_mxe_l6 = lines_strip_dict["MXE"]["volcano"][5],
		volcano_mxe_l7 = lines_strip_dict["MXE"]["volcano"][6],
		volcano_mxe_l8 = lines_strip_dict["MXE"]["volcano"][7],
		volcano_mxe_l9 = lines_strip_dict["MXE"]["volcano"][8],
		scatter_mxe_l1 = lines_strip_dict["MXE"]["scatter"][0],
		scatter_mxe_l2 = lines_strip_dict["MXE"]["scatter"][1],
		scatter_mxe_l3 = lines_strip_dict["MXE"]["scatter"][2],
		scatter_mxe_l4 = lines_strip_dict["MXE"]["scatter"][3],
		scatter_mxe_l5 = lines_strip_dict["MXE"]["scatter"][4],
		scatter_mxe_l6 = lines_strip_dict["MXE"]["scatter"][5],
		scatter_mxe_l7 = lines_strip_dict["MXE"]["scatter"][6],
		scatter_mxe_l8 = lines_strip_dict["MXE"]["scatter"][7],
		scatter_mxe_l9 = lines_strip_dict["MXE"]["scatter"][8],
		bar_mxe_l1 = lines_strip_dict["MXE"]["bar"][0],
		bar_mxe_l2 = lines_strip_dict["MXE"]["bar"][1],
		bar_mxe_l3 = lines_strip_dict["MXE"]["bar"][2],
		bar_mxe_l4 = lines_strip_dict["MXE"]["bar"][3],
		bar_mxe_l5 = lines_strip_dict["MXE"]["bar"][4],
		bar_mxe_l6 = lines_strip_dict["MXE"]["bar"][5],
		bar_mxe_l7 = lines_strip_dict["MXE"]["bar"][6],
		bar_mxe_l8 = lines_strip_dict["MXE"]["bar"][7],
		bar_mxe_l9 = lines_strip_dict["MXE"]["bar"][8],
		volcano_ri_l1 = lines_strip_dict["RI"]["volcano"][0],
		volcano_ri_l2 = lines_strip_dict["RI"]["volcano"][1],
		volcano_ri_l3 = lines_strip_dict["RI"]["volcano"][2],
		volcano_ri_l4 = lines_strip_dict["RI"]["volcano"][3],
		volcano_ri_l5 = lines_strip_dict["RI"]["volcano"][4],
		volcano_ri_l6 = lines_strip_dict["RI"]["volcano"][5],
		volcano_ri_l7 = lines_strip_dict["RI"]["volcano"][6],
		volcano_ri_l8 = lines_strip_dict["RI"]["volcano"][7],
		volcano_ri_l9 = lines_strip_dict["RI"]["volcano"][8],
		scatter_ri_l1 = lines_strip_dict["RI"]["scatter"][0],
		scatter_ri_l2 = lines_strip_dict["RI"]["scatter"][1],
		scatter_ri_l3 = lines_strip_dict["RI"]["scatter"][2],
		scatter_ri_l4 = lines_strip_dict["RI"]["scatter"][3],
		scatter_ri_l5 = lines_strip_dict["RI"]["scatter"][4],
		scatter_ri_l6 = lines_strip_dict["RI"]["scatter"][5],
		scatter_ri_l7 = lines_strip_dict["RI"]["scatter"][6],
		scatter_ri_l8 = lines_strip_dict["RI"]["scatter"][7],
		scatter_ri_l9 = lines_strip_dict["RI"]["scatter"][8],
		bar_ri_l1 = lines_strip_dict["RI"]["bar"][0],
		bar_ri_l2 = lines_strip_dict["RI"]["bar"][1],
		bar_ri_l3 = lines_strip_dict["RI"]["bar"][2],
		bar_ri_l4 = lines_strip_dict["RI"]["bar"][3],
		bar_ri_l5 = lines_strip_dict["RI"]["bar"][4],
		bar_ri_l6 = lines_strip_dict["RI"]["bar"][5],
		bar_ri_l7 = lines_strip_dict["RI"]["bar"][6],
		bar_ri_l8 = lines_strip_dict["RI"]["bar"][7],
		bar_ri_l9 = lines_strip_dict["RI"]["bar"][8]

	)

	with open(output_dir + "/summary.html", 'w') as f:
		f.write(summary_html)

	return 0

def main():

	args = get_args()
	input_dir = args.input
	output_dir = args.output
	# Make directory
	os.makedirs(output_dir + "/data", exist_ok = True)

	print("Making plots....", file = sys.stdout)
	# PCA
	experiment_table_df = load_experiment_table(args.experiment_table)
	# TPM
	pca_tpm_df, contribution_tpm_PC1, contribution_tpm_PC2 = load_tpm_pca_table(input_dir, experiment_table_df, output_dir)
	plots_pca("TPM", pca_tpm_df, contribution_tpm_PC1, contribution_tpm_PC2, output_dir)
	# PSI
	pca_psi_df, contribution_psi_PC1, contribution_psi_PC2 = load_psi_pca_table(input_dir, experiment_table_df, output_dir)
	plots_pca("PSI", pca_psi_df, contribution_psi_PC1, contribution_psi_PC2, output_dir)

	AS_list = ["SE", "FIVE", "THREE", "MXE", "RI"]

	for AS in AS_list:

		plots(AS, input_dir, output_dir)

	write_summary_html(output_dir)

	print("Plots: " + output_dir + "/summary.html", file = sys.stdout)

if __name__ == '__main__':

    main()
