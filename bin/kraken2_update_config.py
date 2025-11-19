#!/usr/bin/env python

## Originally written by David Ross for use within __
## See git repository (https://github.com/) for full license text.

# USAGE: ./kraken2_update_config.py -f $kraken2_config_file

import argparse

parser = argparse.ArgumentParser(
                    prog='kraken2_update_config',
                    description='Update the kraken2_db parameter in the designated databases.config to avoid redownloading',
                    epilog='')

parser.add_argument('-d', '--kraken2_db')       # kraken2 db file path
parser.add_argument('-c', '--db_config')        # databases.config path
parser.add_argument('-v', '--verbose',
                    action='store_true')        # on/off flag

args = parser.parse_args()

KRAKEN2_DB = args.kraken2_db
DATABASE_CONFIG = args.db_config

print(KRAKEN2_DB)
print(DATABASE_CONFIG)

# Read the existing database_download.config file
with open(DATABASE_CONFIG, "r", encoding="utf8") as file:
    text = file.read()

# Replace the dram config loc parameter 
text = text.replace('kraken2_db = false', f'kraken2_db = "{KRAKEN2_DB}"')

# Replace the database_download.config with the updated dram_config_loc
with open(DATABASE_CONFIG, "w") as file:
  file.write(text)