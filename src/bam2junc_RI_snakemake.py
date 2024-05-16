import argparse
import os
import subprocess
import pysam


def get_args():
	## Get arguments from command line

	parser = argparse.ArgumentParser(
		formatter_class = argparse.ArgumentDefaultsHelpFormatter,
		description = "bam2junc_RI_snakemake.py"
	)

	parser.add_argument('-b', '--bam', type=str, help='Input bam file')
	parser.add_argument('-r', '--RI', type=str, help='Input RI file')
	parser.add_argument('-o', '--junc', type=str, help='Output junction file')
	parser.add_argument('-t', '--threads', type=int, help='Number of threads')

	args = parser.parse_args()

	return args


def bam2junc(bam, RI, output, threads):

	shell = ["featureCounts", "-a", RI, "-o", output, "-F", "SAF", "--fracOverlapFeature", "1.0", "-T", str(threads), "-O"]

	bamfile = pysam.AlignmentFile(bam, "rb")
	for read in bamfile.fetch():
		if read.is_paired:
			subprocess.call(shell + ["-p", bam])
		else:
			subprocess.call(shell + [bam])
		break

def main():

	args = get_args()
	bam2junc(args.bam, args.RI, args.junc, args.threads)

if __name__ == "__main__":

	main()
