#!/usr/bin/env python3

## Originally written by David Ross for use within MAGnolia

# USAGE: ./bakta_update_config.py -d $bakta_db_dir -c $config_file

import argparse

parser = argparse.ArgumentParser(
                    prog='bakta_update_config',
                    description='Update the bakta_db parameter in the designated databases.config to avoid redownloading',
                    epilog='')

parser.add_argument('-d', '--bakta_db')         # bakta db file path
parser.add_argument('-c', '--db_config')        # databases.config path
parser.add_argument('-v', '--verbose',
                    action='store_true')        # on/off flag

args = parser.parse_args()

BAKTA_DB = args.bakta_db
DATABASE_CONFIG = args.db_config

#print(BAKTA_DB)
#print(DATABASE_CONFIG)

# Read the existing databases.config file
with open(DATABASE_CONFIG, "r", encoding="utf8") as file:
    text = file.read()

# Replace the bakta db file path 
text = text.replace('bakta_db = false', f'bakta_db = "{BAKTA_DB}"')

# Replace the databases.config with the updated bakta_db parameter
with open(DATABASE_CONFIG, "w") as file:
  file.write(text)