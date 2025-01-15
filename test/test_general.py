import unittest
from unittest.mock import patch, mock_open, MagicMock
import subprocess
import os
import sys
import yaml

# Add parent directory to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from lib.general import load_config, check_config, execute_command

class TestGeneralModule(unittest.TestCase):
    def test_load_config_success(self):
        mock_yaml_content = """
        key1: value1
        key2: value2
        """
        with patch("builtins.open", mock_open(read_data=mock_yaml_content)):
            with patch("yaml.safe_load", return_value=yaml.safe_load(mock_yaml_content)) as mock_safe_load:
                config = load_config("config.yaml")
                self.assertEqual(config, {'key1': 'value1', 'key2': 'value2'})
                mock_safe_load.assert_called_once()

    def test_load_config_file_not_found(self):
        with patch("builtins.open", side_effect=FileNotFoundError):
            with self.assertRaises(FileNotFoundError):
                load_config("nonexistent_config.yaml")

    def test_load_config_yaml_error(self):
        invalid_yaml_content = "{invalid: yaml,"
        with patch("builtins.open", mock_open(read_data=invalid_yaml_content)):
            with patch("yaml.safe_load", side_effect=yaml.YAMLError):
                with self.assertRaises(ValueError):
                    load_config("invalid_config.yaml")

    def test_check_config_all_keys_present(self):
        config = {'key1': 'value1', 'key2': 'value2'}
        keys = ['key1', 'key2']
        missing_keys = check_config(config, keys)
        self.assertEqual(missing_keys, [])

    def test_check_config_missing_keys(self):
        config = {'key1': 'value1'}
        keys = ['key1', 'key2']
        missing_keys = check_config(config, keys)
        self.assertEqual(missing_keys, ['key2'])

    def test_check_config_empty_values(self):
        config = {'key1': 'value1', 'key2': None}
        keys = ['key1', 'key2']
        missing_keys = check_config(config, keys)
        self.assertEqual(missing_keys, ['key2'])

    @patch("subprocess.run")
    def test_execute_command_success(self, mock_subprocess):
        mock_subprocess.return_value = MagicMock(returncode=0)
        returncode = execute_command(["echo", "Hello"])
        self.assertEqual(returncode, 0)
        mock_subprocess.assert_called_once_with(["echo", "Hello"], shell=False)

    @patch("subprocess.run")
    def test_execute_command_failure(self, mock_subprocess):
        mock_subprocess.return_value = MagicMock(returncode=1)
        returncode = execute_command(["false"])
        self.assertEqual(returncode, 1)
        mock_subprocess.assert_called_once_with(["false"], shell=False)

if __name__ == "__main__":
    unittest.main()
