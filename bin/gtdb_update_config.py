#!/usr/bin/env python

## Originally written by David Ross for use within __
## See git repository (https://github.com/) for full license text.

# USAGE: ./gtdb_update_config.py -d $gtdb_db_dir -c $config_file

import argparse

parser = argparse.ArgumentParser(
                    prog='gtdb_update_config',
                    description='Update the gtdb_db parameter in the designated databases.config to avoid redownloading',
                    epilog='')

parser.add_argument('-d', '--gtdb_db')         # bakta db file path
parser.add_argument('-c', '--db_config')        # databases.config path
parser.add_argument('-v', '--verbose',
                    action='store_true')        # on/off flag

args = parser.parse_args()

GTDB_DB = args.gtdb_db
DATABASE_CONFIG = args.db_config

# Read the existing databases.config file
with open(DATABASE_CONFIG, "r", encoding="utf8") as file:
    text = file.read()

# Replace the bakta db file path 
text = text.replace('gtdb_db = false', f'gtdb_db = "{GTDB_DB}"')

# Replace the databases.config with the updated bakta_db parameter
with open(DATABASE_CONFIG, "w") as file:
  file.write(text)