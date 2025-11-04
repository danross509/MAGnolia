#!/usr/bin/env python

## Originally written by David Ross for use within __
## See git repository (https://github.com/) for full license text.

# USAGE: ./dram_update_config.py -f $dram_config_file

import argparse

parser = argparse.ArgumentParser(
                    prog='dram_update_config',
                    description='Update the dram_config_loc parameter in the existing database_download_config to avoid redownloading',
                    epilog='')

parser.add_argument('-d', '--dram_config')      # dram config file path
parser.add_argument('-c', '--db_config')        # database_download.config path
parser.add_argument('-v', '--verbose',
                    action='store_true')        # on/off flag

args = parser.parse_args()

DRAM_CONFIG = args.dram_config
DATABASE_CONFIG = args.db_config

print(DRAM_CONFIG)
print(DATABASE_CONFIG)

# Read the existing database_download.config file
with open(DATABASE_CONFIG, "r", encoding="utf8") as file:
    text = file.read()

# Replace the dram config loc parameter 
text = text.replace("dram_config_loc = false", f"dram_config_loc = "{DRAM_CONFIG}"")

# Replace the database_download.config with the updated dram_config_loc
with open(DATABASE_CONFIG, "w") as file:
  file.write(text)