import argparse
import os
import logging
import subprocess

# Configure logging
logger = logging.getLogger(__name__)

def parse_args():
	parser = argparse.ArgumentParser(
		description="Pipeline for transcript assembly using StringTie2"
	)
	parser.add_argument("-i", "--input", required=True, help="Experiment table")
	parser.add_argument("-r", "--reference", required=True, help="Reference GTF file")
	parser.add_argument("-o", "--output", required=True, help="Assembled GTF file")
	parser.add_argument("-p", "--processors", type=int, default=1, help="Number of processors to use (default: 1)")
	parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
	return parser.parse_args()

def check_file(file_path, description):
	if not os.path.isfile(file_path):
		logger.error(f"{description} does not exist: {file_path}")
		raise FileNotFoundError(f"{file_path} not found.")

def validate_output_file(output_file):
	if os.path.isdir(output_file):
		logger.error(f"Output path is a directory, not a file: {output_file}")
		raise IsADirectoryError(f"{output_file} is a directory, not a file.")
	if os.path.exists(output_file) and not os.access(output_file, os.W_OK):
		logger.error(f"Output file is not writable: {output_file}")
		raise PermissionError(f"{output_file} is not writable.")

def run_command(command):
	logger.info(f"Executing: {command}")
	result = subprocess.run(command, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if result.returncode != 0:
		logger.error(f"Command failed: {command}")
		logger.error(result.stderr)
		raise RuntimeError(f"Command failed with exit code {result.returncode}")

def create_directory(directory):
	os.makedirs(directory, exist_ok=True)

def main():

	# Parse arguments
	args = parse_args()
	# Set up logging
	logging.basicConfig(
		format = "[%(asctime)s] %(levelname)7s %(message)s",
		level = logging.DEBUG if args.verbose else logging.INFO
	)

	experiment_file = args.input
	reference_gtf = args.reference
	assembled_gtf = args.output
	num_processors = args.processors

	logger.info("Starting transcript assembly...")
	logger.debug(args)

	# Validate input files
	check_file(experiment_file, "Experiment table")
	check_file(reference_gtf, "Reference GTF")

	# Validate output file
	validate_output_file(assembled_gtf)

	# Prepare output directory
	output_dir = os.path.dirname(assembled_gtf)
	create_directory(output_dir)

	gtf_list_path = os.path.join(output_dir, "gtf_list.txt")
	with open(experiment_file, "r") as experiment, open(gtf_list_path, "w") as gtf_list:
		for line in experiment:
			line = line.strip()
			if not line or line.startswith("sample"):
				continue

			sample, bam_file, _group = line.split(maxsplit=2)
			sample_dir = os.path.dirname(bam_file)
			bam_index = f"{bam_file}.bai"

			logger.info(f"Processing sample: {sample}")
			logger.debug(f"BAM file: {bam_file}")

			# Check or create BAM index
			if not os.path.isfile(bam_index):
				logger.info(f"Creating BAM index for {bam_file}")
				run_command(f"samtools index {bam_file}")
			else:
				logger.debug(f"Found BAM index for {bam_file}")

			# Run StringTie2 for assembly
			sample_gtf = os.path.join(sample_dir, f"{sample}.assembled.gtf")
			stringtie_command = (
				f"stringtie -p {num_processors} -G {reference_gtf} -o {sample_gtf} {bam_file}"
			)
			run_command(stringtie_command)

			# Add to GTF list
			gtf_list.write(f"{sample_gtf}\n")

	# Merge GTF files
	logger.info("Merging GTF files...")
	with open(gtf_list_path, "r") as gtf_list_file:
		gtf_files = " ".join(gtf_list_file.read().splitlines())

	merge_command = f"stringtie --merge -p {num_processors} -G {reference_gtf} -o {assembled_gtf} {gtf_files}"
	run_command(merge_command)

	# Cleanup
	os.remove(gtf_list_path)
	logger.info("Transcript assembly completed successfully!")

if __name__ == "__main__":

	main()
