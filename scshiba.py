import argparse
import os
import sys
import logging
import yaml
import subprocess
# Configure logger
logger = logging.getLogger(__name__)
# Set version
VERSION = "v0.4.1"

def parse_args():
    parser = argparse.ArgumentParser(
                description=f"""scShiba {VERSION} - Pipeline for identification of differential RNA splicing in single-cell RNA-seq data

Step 1: gtf2event.py
    - Converts GTF files to event format.
Step 2: sc2junc.py
    - Counts junction reads from STARsolo output files.
Step 3: scpsi.py
    - Calculates PSI values and perform differential analysis.""",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("config", help="Config file in yaml format")
    parser.add_argument("-p", "--process", type=int, default=1, help="Number of processors to use (default: 1)")
    parser.add_argument("-s", "--start-step", type=int, default=0, help="Start the pipeline from the specified step (default: 0, run all steps)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode")
    return parser.parse_args()

def main():

	# Get arguments
	args = parse_args()

	# Add parent directory to path
	current_dir = os.path.dirname(os.path.abspath(__file__))
	parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
	sys.path.append(parent_dir)
	from lib import general

	# Set up logging
	logging.basicConfig(
		format = "[%(asctime)s] %(levelname)7s %(message)s",
		level = logging.DEBUG if args.verbose else logging.INFO
	)

	# Validate input and config
	logger.info("Running scShiba...")
	logger.debug(f"Arguments: {args}")
	# Get number of processors
	processors = str(args.process)

	# Load config
	logger.info("Loading configuration...")
	config_path = args.config
	config = general.load_config(config_path)

	# Check essential config keys
	missing_keys = general.check_config(config, ["workdir", "experiment_table", "gtf"])
	if missing_keys:
		logger.error(f"Missing required keys in configuration file: {', '.join(missing_keys)}")
		sys.exit(1)
	else:
		logger.info(f"workdir: {config['workdir']}")
		logger.info(f"experiment_table: {config['experiment_table']}")
		logger.info(f"gtf: {config['gtf']}")

	# Prepare output directory
	output_dir = config["workdir"]
	logger.debug("Making output directory...")
	os.makedirs(output_dir, exist_ok=True)
	log_file = os.path.join(output_dir, "scShiba.log")
	logger.debug(f"Log file: {log_file}")
	logging.FileHandler(log_file)

	# Get parent directory of this script
	logger.debug("Getting script directory...")
	script_dir = os.path.dirname(os.path.realpath(__file__))

	# Steps
	experiment_table = config["experiment_table"]
	gtf = config["gtf"]
	steps = [
		{
			"name": "Step 1: gtf2event.py",
			"command": [
				"python", os.path.join(script_dir, "src", "gtf2event.py"),
				"-i", gtf,
				"-o", os.path.join(output_dir, "events"),
				"-p", processors
			]
		},
		{
			"name": "Step 2: sc2junc.py",
			"command": [
				"python", os.path.join(script_dir, "src", "sc2junc.py"),
				"-i", experiment_table,
				"-o", os.path.join(output_dir, "junctions", "junctions.bed")
			]
		},
		{
			"name": "Step 3: scpsi.py",
			"command": [
				"python", os.path.join(script_dir, "src", "scpsi.py"),
				"-p", "1",
				"-r", config['reference_group'],
				"-a", config['alternative_group'],
				"-f", str(config['fdr']),
				"-d", str(config['delta_psi']),
				"-m", str(config['minimum_reads']),
				"--excel" if config['excel'] else "",
				os.path.join(output_dir, "junctions", "junctions.bed"),
				os.path.join(output_dir, "events"),
				os.path.join(output_dir, "results")
			]
		}
	]

	logger.info("Starting pipeline execution...")

	# Get start step
	start_step = args.start_step
	if start_step > 0:
		logger.info(f"Starting from step {start_step}...")
		steps = steps[start_step-1:]

	# Execute steps
	for step in steps:
		logger.info(f"Executing {step['name']}...")
		command_to_run = step["command"] + ["-v"] if args.verbose else step["command"]
		logger.debug(command_to_run)
		returncode = general.execute_command(command_to_run)
		if returncode != 0:
			logger.error(f"Error executing {step['name']}. Exiting...")
			sys.exit(1)

	# Finish
	logger.info(f"scShiba finished! Results saved in {output_dir}")

if __name__ == "__main__":
    main()
