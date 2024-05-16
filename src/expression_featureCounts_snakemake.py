import argparse
import os
import subprocess
import pysam

def get_args():
	## Get arguments from command line

	parser = argparse.ArgumentParser(
		formatter_class = argparse.ArgumentDefaultsHelpFormatter,
		description = "expression_featureCounts_snakemake.py"
	)

	parser.add_argument('-b', '--bam', type=str, help='Input bam file')
	parser.add_argument('-g', '--gtf', type=str, help='Input gtf file')
	parser.add_argument('-o', '--output', type=str, help='Output count file')
	parser.add_argument('-t', '--threads', type=int, help='Number of threads')

	args = parser.parse_args()

	return args

def bam2junc(bam, gtf, output, threads):

	shell = ["featureCounts", "-a", gtf, "-o", output, "-T", str(threads), "-t", "exon", "-g", "gene_id"]

	bamfile = pysam.AlignmentFile(bam, "rb")
	for read in bamfile.fetch():
		if read.is_paired:
			subprocess.call(shell + ["-p", "-B", bam])
		else:
			subprocess.call(shell + [bam])
		break

def main():

	args = get_args()
	bam2junc(args.bam, args.gtf, args.output, args.threads)

if __name__ == "__main__":

	main()
