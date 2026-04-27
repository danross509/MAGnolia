#!/usr/bin/env python3

## Originally written by David Ross for use within MAGnolia

# USAGE: ./eggnog_update_config.py -d $eggnog_db_dir -c $config_file

import argparse

parser = argparse.ArgumentParser(
                    prog='eggnog_update_config',
                    description='Update the eggnog_db parameter in the designated databases.config to avoid redownloading',
                    epilog='')

parser.add_argument('-d', '--eggnog_db')       # eggnog db file path
parser.add_argument('-c', '--db_config')        # databases.config path
parser.add_argument('-v', '--verbose',
                    action='store_true')        # on/off flag

args = parser.parse_args()

EGGNOG_DB = args.eggnog_db
DATABASE_CONFIG = args.db_config

print(EGGNOG_DB)
print(DATABASE_CONFIG)

# Read the existing databases.config file
with open(DATABASE_CONFIG, "r", encoding="utf8") as file:
    text = file.read()

# Replace the eggnog db file path 
text = text.replace('eggnog_db = false', f'eggnog_db = "{EGGNOG_DB}"')

# Replace the databases.config with the updated eggnog_db parameter
with open(DATABASE_CONFIG, "w") as file:
  file.write(text)