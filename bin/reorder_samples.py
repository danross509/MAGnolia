#!/usr/bin/env python

## Originally written by David Ross for use within __
## See git repository (https://github.com/) for full license text.

# USAGE: ./reorder_samples.py
# Used in sample preparation step to reorder the samples in samples_tmp.csv for easier manipulation by the user

import pandas as pd

df = pd.read_csv("samples_tmp.csv", header=0, na_filter=False)

df.sort_values(by='sampleID', inplace=True)
        
df.to_csv("samples.csv", index=False)