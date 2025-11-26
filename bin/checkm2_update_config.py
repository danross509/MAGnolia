#!/usr/bin/env python

## Originally written by David Ross for use within __
## See git repository (https://github.com/) for full license text.

# USAGE: ./checkm2_update_config.py -d $checkm2_db_dir -c $config_file

import argparse

parser = argparse.ArgumentParser(
                    prog='checkm2_update_config',
                    description='Update the checkm2_db parameter in the designated databases.config to avoid redownloading',
                    epilog='')

parser.add_argument('-d', '--checkm2_db')       # checkm2 db file path
parser.add_argument('-c', '--db_config')        # databases.config path
parser.add_argument('-v', '--verbose',
                    action='store_true')        # on/off flag

args = parser.parse_args()

CHECKM2_DB = args.checkm2_db
DATABASE_CONFIG = args.db_config

print(CHECKM2_DB)
print(DATABASE_CONFIG)

# Read the existing databases.config file
with open(DATABASE_CONFIG, "r", encoding="utf8") as file:
    text = file.read()

# Replace the checkm2 db file path 
text = text.replace('checkm2_db = false', f'checkm2_db = "{CHECKM2_DB}"')

# Replace the databases.config with the updated checkm2_db parameter
with open(DATABASE_CONFIG, "w") as file:
  file.write(text)