#!/usr/bin/env python3

## Originally written by David Ross for use within MAGnolia

# USAGE: ./checkm_update_config.py -d $checkm_db_dir -c $config_file

import argparse

parser = argparse.ArgumentParser(
                    prog='checkm_update_config',
                    description='Update the checkm_db parameter in the designated databases.config to avoid redownloading',
                    epilog='')

parser.add_argument('-d', '--checkm_db')       # checkm db file path
parser.add_argument('-c', '--db_config')        # databases.config path
parser.add_argument('-v', '--verbose',
                    action='store_true')        # on/off flag

args = parser.parse_args()

CHECKM_DB = args.checkm_db
DATABASE_CONFIG = args.db_config

print(CHECKM_DB)
print(DATABASE_CONFIG)

# Read the existing databases.config file
with open(DATABASE_CONFIG, "r", encoding="utf8") as file:
    text = file.read()

# Replace the checkm db file path 
text = text.replace('checkm_db = false', f'checkm_db = "{CHECKM_DB}"')

# Replace the databases.config with the updated checkm_db parameter
with open(DATABASE_CONFIG, "w") as file:
  file.write(text)