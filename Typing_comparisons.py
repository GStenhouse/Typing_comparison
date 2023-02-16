# This is a script to take an excel spreadsheet of isolates and their types (by various typing schemes) and
# compare the granualrity of the grouping by the various typing methods. 

# The script will generate descriptive stats of the level of similarity between the grouping of isolates 
# between each typing scheme - based on the numbers of isolates in each group (ignoring population struture).
# These descriptive stats will be used to generate stats describing the skew of granularity between two typing 
# schemes.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Import packages ##

from collections import Counter
import numpy as np
import pandas as pd
import itertools

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###############
#  FUNCTIONS  #
###############

# To create a dictionary of pairwise pivot tables
def pairingtab(df_in, list):
    pivot_dict = {}

    for pair in list:
        x = pair[0]
        y = pair[1]
        pivot = pd.crosstab(df_in[x], df_in[y], margins=True)

        key = x + " - " + y

        pivot_dict[key] = pivot
    return pivot_dict

# Identify matching groups, with and without isolate numbers
def compare(pivot_dict):
    piv_comp_dict = {}            # create empty dicts to collect all  output data - pivot : row and matching columns
    piv_pair_dict = {}                # pivot : pair and number of matches

    for key, piv in pivot_dict.items():
        row_names = piv.index.to_series() # create list of row names
        col_cat = piv.columns.to_list()   # create list of column names

        pair_out = []              # create empty list to take the pairs and number objects
        row_match_dict = {}     # create empty dictionaries for collecting the data - list of matching columns by row
        row_pair_dict = {}          # group pairs and number of matchin isolates
                
        for i in range(0,len(piv.index)):
            row = piv.iloc[i]              # extract row
            r_name = row_names.iloc[i]         # collect row name
            col_matches = []                       # create empty list for names of columns with isolates in row
            row_count_dict = {}                     # create dict for collecting number of column groups which match row group

            column_no = -1      # set column number to -1 (so that it changes to 0 before cell is checked)

            if r_name != 'All':     # skip grand total rows   
                for cell in row.values.tolist():
                    column_no += 1              # add 1 to column number
                
                    if cell  > 0:
                        c_name = col_cat[column_no]     # collect corresponding column name to non-zero cell

                        if c_name != 'All':     # ignore grand total columns
                            col_matches.append(c_name)        # add matching column name to list
                            pair = str(r_name) + " " + str(c_name) + "\t" + str(cell)
                        
            row_match_dict[r_name] = col_matches        # link row name with matching columns
            pair_out.append(pair)

        piv_comp_dict[key] = row_match_dict
        piv_pair_dict[key] = pair_out


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##########
#  MAIN  #
##########

# Read in spreadsheet
df_in = pd.read_excel('~/OneDrive - The University of Liverpool (1)/Typing_test_data.xlsx')

# define number of columns as object
x = len(df_in.columns)

# create list of column names
columns = df_in.columns.to_list()

# create a series of pairwise comparison pivot tables.
# create list of column names as pairwise comparisons
pair_list = itertools.permutations(columns, r=2)

# use pairwise list to run pair wise pivot table creation
pivot_dict = pairingtab(df_in, pair_list)

# use output to identify matching groups
compare(pivot_dict)





# Create subset dataframes with all the non-intersecting categories removed, tables of just intersecting categories??