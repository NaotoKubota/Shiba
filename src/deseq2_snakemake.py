import argparse
import pandas as pd
import os
import sys
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(parent_dir)
from lib import general
import logging
# Configure logging
logger = logging.getLogger(__name__)

def str2bool(v):
    if v == "True":
        return True
    elif v == "False":
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def parse_args():
	## Get arguments from command line

	parser = argparse.ArgumentParser(
		formatter_class = argparse.ArgumentDefaultsHelpFormatter,
		description = "DESeq2"
	)
	parser.add_argument("--count", type = str, help = "Read count files")
	parser.add_argument("--experiment-table", type = str, help = "Experiment table")
	parser.add_argument("--reference", type = str, help = "Reference")
	parser.add_argument("--alternative", type = str, help = "Alternative")
	parser.add_argument("--output", type = str, help = "Output file")
	parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
	args = parser.parse_args()
	return args

def deseq2(count, experiment_table, reference, alternative, output):

	# Load experiment table
	df = pd.read_csv(experiment_table, sep = "\t")
	# Count the number of samples in each group
	count_df = df.groupby("group").count()
	# Get the number of samples in the reference and alternative groups
	count_reference = count_df["sample"][reference]
	count_alternative = count_df["sample"][alternative]

	# Check if the number of samples is greater than or equal to 2
	if count_reference >= 2 and count_alternative >= 2:
		# Get path of directory where this script is located
		run_command = ["Rscript", os.path.join(current_dir, "deseq2.R"), experiment_table, count, reference, alternative, output]
		returncode = general.execute_command(run_command)
		if returncode != 0:
			logger.error(f"Error executing the command. Exiting...")
			sys.exit(1)
	else:
		logger.error("Error: The number of samples is less than 2.")
		# Save the error message to a file
		with open(output, "w") as f:
			f.write("Error: The number of samples is less than 2.\n")

def main():

	# Parse arguments
	args = parse_args()
	# Set up logging
	logging.basicConfig(
		format = "[%(asctime)s] %(levelname)7s %(message)s",
		level = logging.DEBUG if args.verbose else logging.INFO
	)
	logger.debug(args)

	# Run DESeq2
	logger.info("Running DESeq2...")
	deseq2(args.count, args.experiment_table, args.reference, args.alternative, args.output)

	# Finish
	logger.info("Done.")

if __name__ == "__main__":

	main()
