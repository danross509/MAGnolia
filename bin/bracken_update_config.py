#!/usr/bin/env python

## Originally written by David Ross for use within __
## See git repository (https://github.com/) for full license text.

# USAGE: ./bracken_update_config.py -d $bracken_built -c $bracken_config_file

import argparse

parser = argparse.ArgumentParser(
                    prog='bracken_update_config',
                    description='Update the bracken_build_exists parameter in the designated databases.config to avoid rebuilding',
                    epilog='')

parser.add_argument('-d', '--bracken_built')    # bracken_built = true
parser.add_argument('-c', '--db_config')        # databases.config path
parser.add_argument('-v', '--verbose',
                    action='store_true')        # on/off flag

args = parser.parse_args()

BUILT = args.bracken_built
DATABASE_CONFIG = args.db_config

#print(BUILT)
#print(DATABASE_CONFIG)

# Read the existing databases.config file
with open(DATABASE_CONFIG, "r", encoding="utf8") as file:
    text = file.read()

# Replace the bracken_build_exists parameter 
text = text.replace('bracken_build_exists = false', f'bracken_build_exists = {BUILT}')

# Replace the databases.config with the updated bracken_build_exists parameter
with open(DATABASE_CONFIG, "w") as file:
  file.write(text)