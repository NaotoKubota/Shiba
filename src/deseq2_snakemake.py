import argparse
import subprocess
import pandas as pd
import os

def str2bool(v):
    if v == "True":
        return True
    elif v == "False":
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def get_args():

	parser = argparse.ArgumentParser(
		formatter_class = argparse.ArgumentDefaultsHelpFormatter,
		description = "DESeq2"
	)

	parser.add_argument("--count", type = str, help = "Read count files")
	parser.add_argument("--experiment-table", type = str, help = "Experiment table")
	parser.add_argument("--reference", type = str, help = "Reference")
	parser.add_argument("--alternative", type = str, help = "Alternative")
	parser.add_argument("--output", type = str, help = "Output file")

	args = parser.parse_args()

	return args

def deseq2(count, experiment_table, reference, alternative, output):

	df = pd.read_csv(experiment_table, sep = "\t")
	count_df = df.groupby("group").count()
	count_reference = count_df["sample"][reference]
	count_alternative = count_df["sample"][alternative]

	if count_reference >= 2 and count_alternative >= 2:

		# Get path of directory where this script is located
		workdir = os.path.dirname(os.path.abspath(__file__))
		subprocess.call(["/opt_shiba/R/4.1.3/bin/Rscript", workdir + "/deseq2.R", experiment_table, count, reference, alternative, output])

	else:

		print("Error: The number of samples is less than 2.")
		# Save the error message to a file
		with open(output, "w") as f:
			f.write("Error: The number of samples is less than 2.\n")

def main():

	args = get_args()
	deseq2(args.count, args.experiment_table, args.reference, args.alternative, args.output)

if __name__ == "__main__":

	main()
