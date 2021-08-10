#!/usr/bin/env python

import os
import pandas as pd

# reading .odf (libre office spread-sheets)
import ezodf
import odf


## ----------------
## Path Environment
## ----------------

def get_data_dir(path_file, path_dir):
    
    rel_path = os.path.relpath(path_file, path_dir)
    data_dir = rel_path.split(os.sep)[0]
    
    return data_dir


def get_sequencing_type(sample_prefix):

    if 'ChIP' in sample_prefix or 'INPUT' in sample_prefix:
        return 'ChIP'

    elif 'RIP' in sample_prefix:
        return 'RIP'

    elif 'totalRNA' in sample_prefix:
        return 'totalRNA'

    elif 'RNA' in sample_prefix or 'pA' in sample_prefix :
        return 'RNA'

    else:
        return 'unknown'


## -------------
## Spread-sheets
## -------------

def read_odf_doc(file):
        
    print("Importing odf file %s ... \n" % file)
    doc = ezodf.opendoc(file)

    print("Spreadsheet contains %d sheet(s)." % len(doc.sheets))
    for sheet in doc.sheets:
        print("-"*40)
        print("   Sheet name : '%s'" % sheet.name)
        print("Size of Sheet : (rows=%d, cols=%d)" % (sheet.nrows(), sheet.ncols()) )

    # convert the first sheet to a pandas.DataFrame
    sheet = doc.sheets[0]
    df_dict = {}
    n_cols = 0
    
    for i, row in enumerate(sheet.rows()):
        
        # row is a list of cells
        # assume the header is on the first row
        if i == 0:
            #import pdb; pdb.set_trace()
            # columns as lists in a dictionary
            df_dict = {cell.value:[] for cell in row}
            
            # remove empty columns
            try:
                del df_dict[None]
            except:
                pass
            
            # create index for the column headers
            col_index = {j:cell.value for j, cell in enumerate(row)}
            # remove empty colums
            col_index = {kk:vv for (kk, vv) in col_index.items() if not isinstance(vv, type(None))}
            
            # init number of columns
            n_cols = len(col_index)
            
            continue
            
        for j, cell in enumerate(row):
            
            ## only use non-empty cols!
            if j <= n_cols - 1:
                # use header instead of column index
                df_dict[col_index[j]].append(cell.value)
            
    # and convert to a DataFrame
    df = pd.DataFrame(df_dict)
    # drop empty columns
    df = df[~df.isnull().all(axis='columns')]
    
    print("\nDone.")

    return df
