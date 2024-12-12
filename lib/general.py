# Modules used in shiba.py and scshiba.py
import os
import subprocess
import sys
import yaml
import logging

def load_config(config_path):
    """
    Loads a YAML configuration file.

    Parameters:
    config_path (str): Path to the YAML configuration file.

    Returns:
    dict: The loaded configuration as a dictionary.
    """
    try:
        with open(config_path, 'r') as file:
            config = yaml.safe_load(file)
        return config
    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    except yaml.YAMLError as e:
        raise ValueError(f"Error parsing YAML configuration file: {e}")

def check_config(config, keys):
    missing_keys = [key for key in keys if key not in config or not config[key]]
    return missing_keys

def execute_command(command):
    result = subprocess.run(command, shell=False)
    return result.returncode
