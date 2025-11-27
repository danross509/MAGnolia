#!/usr/bin/env python

## Originally written by David Ross for use within __
## See git repository (https://github.com/) for full license text.

# USAGE: ./kraken2_update_config.py -d $kraken2_db_dir -c $config_file

import argparse

parser = argparse.ArgumentParser(
                    prog='kraken2_update_config',
                    description='Update the kraken2_db parameter in the designated databases.config to avoid redownloading',
                    epilog='')

parser.add_argument('-d', '--kraken2_db')       # kraken2 db file path
parser.add_argument('-c', '--db_config')        # databases.config path
parser.add_argument('-e', '--executables')      # path to kraken2 executables
parser.add_argument('-v', '--verbose',
                    action='store_true')        # on/off flag

args = parser.parse_args()

KRAKEN2_DB = args.kraken2_db
DATABASE_CONFIG = args.db_config
EXEC = args.executables

print(KRAKEN2_DB)
print(DATABASE_CONFIG)

# Read the existing databases.config file
with open(DATABASE_CONFIG, "r", encoding="utf8") as file:
    text = file.read()

# Replace the kraken2 db file path 
text = text.replace('kraken2_db = false', f'kraken2_db = "{KRAKEN2_DB}"')
text = text.replace('kraken2_exec = false', f'kraken2_exec = "{EXEC}"')

# Replace the databases.config with the updated kraken2_db parameter
with open(DATABASE_CONFIG, "w") as file:
  file.write(text)