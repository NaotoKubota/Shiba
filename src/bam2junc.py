import argparse
import os
import subprocess
import logging
import pandas as pd
import pysam

# Configure logging
logger = logging.getLogger(__name__)

def get_args():
	parser = argparse.ArgumentParser(
		description="Pipeline for processing junction read counts."
	)
	parser.add_argument("-i", "--input", required=True, help="Experiment table")
	parser.add_argument("-r", "--ri_event", required=True, help="Intron retention event file")
	parser.add_argument("-o", "--output", required=True, help="Output junction read counts file")
	parser.add_argument("-p", "--processors", type=int, default=1, help="Number of processors to use (default: 1)")
	parser.add_argument("-a", "--anchor", type=int, default=8, help="Minimum anchor length (default: 8)")
	parser.add_argument("-m", "--min_intron", type=int, default=70, help="Minimum intron size (default: 70)")
	parser.add_argument("-M", "--max_intron", type=int, default=500000, help="Maximum intron size (default: 500000)")
	parser.add_argument("-s", "--strand", default="XS", help="Strand specificity (default: XS)")
	parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
	return parser.parse_args()

def check_file(file_path, description):
	if not os.path.isfile(file_path):
		logger.error(f"{description} does not exist: {file_path}")
		raise FileNotFoundError(f"{file_path} not found.")

def run_command(command, log_file=None):
	logger.info(f"Executing: {command}")
	with open(log_file, "a") if log_file else None as log:
		result = subprocess.run(command, shell=True, text=True, stdout=log, stderr=log)
	if result.returncode != 0:
		logger.error(f"Command failed: {command}")
		raise RuntimeError(f"Command failed with exit code {result.returncode}")

def prepare_output_dir(output_path):
	output_dir = os.path.dirname(output_path)
	logs_dir = os.path.join(output_dir, "logs")
	os.makedirs(logs_dir, exist_ok=True)
	return output_dir, logs_dir

def create_saf_file(ri_event, output_dir):
	saf_file = os.path.join(output_dir, "RI.saf")
	logger.info(f"Generating SAF file: {saf_file}")
	saf_data = []
	with open(ri_event, "r") as ri_file:
		next(ri_file)  # Skip header
		for line in ri_file:
			intron, strand = line.strip().split("\t")[5:7]
			chrom = intron.split(":")[0]
			start = intron.split(":")[1].split("-")[0]
			start_plus_1 = str(int(start) + 1)
			end = intron.split(":")[1].split("-")[1]
			end_minus_1 = str(int(end) - 1)
			saf_data.append([f"{chrom}:{start}-{start_plus_1}", chrom, start, start_plus_1, strand])
			saf_data.append([f"{chrom}:{end_minus_1}-{end}", chrom, end_minus_1, end, strand])
	saf_df = pd.DataFrame(saf_data, columns=["GeneID", "Chr", "Start", "End", "Strand"])
	saf_df.drop_duplicates().to_csv(saf_file, sep="\t", index=False)
	return saf_file

def is_paired_end(bam_file):
	"""
	Determine if a BAM file is paired-end.
	Returns True if paired-end, False otherwise.
	"""
	logger.debug(f"Checking if BAM file is paired-end: {bam_file}")
	with pysam.AlignmentFile(bam_file, "rb") as bam:
		for read in bam:
			if read.is_paired:
				return True
		return False

def process_samples(experiment_file, strand, anchor, min_intron, max_intron, output_dir, logs_dir, saf_file, processors):
	junc_files = []
	with open(experiment_file, "r") as experiment:
		for line in experiment:
			line = line.strip()
			if not line or line.startswith("sample"):
				continue
			sample, bam, _group = line.split(maxsplit=2)
			sample_dir = os.path.dirname(bam)
			bam_index = f"{bam}.bai"

			# Ensure BAM index exists
			if not os.path.isfile(bam_index):
				logger.info(f"Indexing BAM file for sample {sample}: {bam}")
				run_command(f"samtools index {bam}")

			# Extract exon-exon junctions
			logger.info(f"Counting exon-exon junctions for sample {sample}...")
			exon_junc_file = os.path.join(sample_dir, f"{sample}_exon-exon.junc")
			run_command(
				f"regtools junctions extract -s {strand} -a {anchor} -m {min_intron} -M {max_intron} "
				f"-o {exon_junc_file} {bam}",
				os.path.join(logs_dir, "regtools.log"),
			)
			junc_files.append((exon_junc_file, "exon-exon"))

			# Count exon-intron junctions
			logger.info(f"Counting exon-intron junctions for sample {sample}...")
			exon_intron_file = os.path.join(sample_dir, f"{sample}_exon-intron.junc")
			if is_paired_end(bam):
				feature_counts_cmd = (
					f"featureCounts -a {saf_file} -o {exon_intron_file} -F SAF --fracOverlapFeature 1.0 "
					f"-T {processors} -O -p {bam}"
				)
			else:
				feature_counts_cmd = (
					f"featureCounts -a {saf_file} -o {exon_intron_file} -F SAF --fracOverlapFeature 1.0 "
					f"-T {processors} -O {bam}"
				)
			run_command(feature_counts_cmd, os.path.join(logs_dir, "featureCounts.log"))
			junc_files.append((exon_intron_file, "exon-intron"))

	return junc_files

def merge_junction_files(junc_files, output_file):
	exon_exon_files = [j[0] for j in junc_files if j[1] == "exon-exon"]
	exon_intron_files = [j[0] for j in junc_files if j[1] == "exon-intron"]

	def process_junction_files(files, junction_type):
		logger.info(f"Merging {junction_type} junction files...")
		result_df = None
		for file in files:
			sample_name = os.path.basename(file).rstrip("_exon-exon.junc") if junction_type == "exon-exon" else os.path.basename(file).rstrip("_exon-intron.junc")
			logger.debug(f"Sample name: {sample_name}")
			df = pd.read_csv(file, sep="\t", comment="#" if junction_type == "exon-intron" else None, header = None, dtype = {0: str})
			if junction_type == "exon-exon":
				df = df.iloc[:, [0, 1, 2, 4, 10]]
				df.columns = ["chr", "start", "end", "count", "block"]
				df = df.reset_index()
				df.loc[(df["chr"].str.isdecimal() == True) | (df["chr"].str.len() <= 2), "chr"] = "chr" + df["chr"]
				df["blockSize1"] = df["block"].str.split(",", expand = True)[0].astype("int32")
				df["start"] = df["start"].astype("int32") + df["blockSize1"]
				df["blockSize2"] = df["block"].str.split(",", expand = True)[1].astype("int32")
				df["end"] = df["end"].astype("int32") - df["blockSize2"] + 1
				df["ID"] = df["chr"].astype(str) + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
			else:
				df = df.iloc[:, [1, 2, 3, 0, 6]]
				# Remove the first row
				df = df.iloc[1:]
				df.columns = ["chr", "start", "end", "ID", "count"]
				df.loc[(df["chr"].str.isdecimal() == True) | (df["chr"].str.len() <= 2), "chr"] = "chr" + df["chr"]
			df["sample"] = sample_name
			df = df[["ID", "sample", "count"]] if junction_type == "exon-exon" else df
			result_df = pd.concat([result_df, df]) if result_df is not None else df
		result_df["count"] = result_df["count"].astype("int32")
		result_df = result_df.drop_duplicates()
		return result_df

	exon_exon_df = process_junction_files(exon_exon_files, "exon-exon")
	logger.debug("Check if there are duplicated junctions")
	dup_df = exon_exon_df.groupby(["ID", "sample"], as_index=False).count()
	dup_df = dup_df[dup_df["count"] > 1]
	if dup_df.shape[0] > 0:
		logger.error("Duplicated junctions:")
		logger.error(dup_df)
		logger.error("Please check the input files!")
		raise ValueError("Duplicated junctions found.")
	exon_exon_df = exon_exon_df.pivot(
		index="ID",
		columns="sample",
		values="count"
	).fillna(0).reset_index()
	exon_exon_df["chr"] = exon_exon_df["ID"].str.split(":", expand=True)[0]
	exon_exon_df["start"] = exon_exon_df["ID"].str.split(":", expand=True)[1].str.split("-", expand=True)[0].astype("int32")
	exon_exon_df["end"] = exon_exon_df["ID"].str.split(":", expand=True)[1].str.split("-", expand=True)[1].astype("int32")
	col = [i for i in exon_exon_df.columns if i not in ["chr", "start", "end", "ID"]]
	for j in col:
		exon_exon_df = exon_exon_df.astype({j: "int32"})
	col = ["chr", "start", "end", "ID"] + col
	exon_exon_df = exon_exon_df[col]

	exon_intron_df = process_junction_files(exon_intron_files, "exon-intron")
	exon_intron_df = exon_intron_df.pivot(
		index = ["chr", "start", "end", "ID"],
		columns = "sample",
		values = "count"
	).fillna(0).reset_index()
	col = [i for i in exon_intron_df.columns if i not in ["chr", "start", "end", "ID"]]
	for j in col:
		exon_intron_df = exon_intron_df.astype({j: "int32"})

	final_df = pd.concat([exon_exon_df, exon_intron_df]).sort_values(["chr", "start"])
	final_df.to_csv(output_file, sep="\t", index = False)
	logger.info(f"Junction read counts merged into {output_file}")

def main():

	# Parse arguments
	args = get_args()
	# Set up logging
	logging.basicConfig(
		format = "[%(asctime)s] %(levelname)7s %(message)s",
		level = logging.DEBUG if args.verbose else logging.INFO
	)
	logger.info("Processing junction read counts...")
	logger.debug(args)

	check_file(args.input, "Experiment table")
	check_file(args.ri_event, "Intron retention event file")

	output_dir, logs_dir = prepare_output_dir(args.output)
	saf_file = create_saf_file(args.ri_event, output_dir)
	logger.info("Extracting junctions from BAM files...")
	junc_files = process_samples(
		args.input, args.strand, args.anchor, args.min_intron, args.max_intron, output_dir, logs_dir, saf_file, args.processors
	)
	logger.debug(junc_files)
	logger.info("Merging junction read counts...")
	merge_junction_files(junc_files, args.output)

	# Cleanup
	os.remove(saf_file)
	logger.info("Junction read counts processing completed!")

if __name__ == "__main__":
	main()
