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

###################################################
#                                                 #
# From Excel to pairwise comparisons, pivot style #
#                                                 #
###################################################

# Read in spreadsheet
df_in = pd.read_excel('~/OneDrive - The University of Liverpool (1)/Typing_test_data.xlsx')

#print(df_in)

# define number of columns as object
x = len(df_in.columns)

#print(x)

# create list of column names
columns = df_in.columns.to_list()

#print(columns)

# create a series of pairwise comparison pivot tables.
# create list of column names as pairwise comparisons
pair_list = itertools.permutations(columns, r=2)

# use pairwise list to run loop and iteratively add pivots to dictionary, pair as a string = key : pivot
pivot_dict = {}

for pair in pair_list:
    x = pair[0]
    y = pair[1]
    pivot = pd.crosstab(df_in[x], df_in[y], margins=True)

    key = x + y

    pivot_dict[key] = pivot




key = "HC200250 SNP Threshold"
pivot = pivot_dict[key]
row_names = pivot.index.to_series() # create list of row names
col_cat = pivot.columns.to_list()   # create list of column names
pivot_comp_dict = {}              # create empty dictionaries for collecting the data
pivot_count_dict = {}
comp = []


for i in range(0,len(pivot.index)):
    row = pivot.iloc[i]              # extract row
    r_name = row_names.iloc[i]         # collect row name
    name_dict = {}                       # create empty dictionaries to collect data for future dataframes
    row_count_dict = {} 
    

    counter = 0         # set counters to 0
    column_no = -1
            
    for cell in row.values.tolist():
        column_no += 1
        if cell > 0:
            counter += 1            # add 1 every time the typing systems intersect
            c_name = col_cat[column_no]
            name_dict[r_name] = c_name         
                    
    comp.append(name_dict)
    row_count_dict[r_name] = counter
    row_count_df = pd.DataFrame.from_dict(row_count_dict, orient = "index")
    
pivot_comp_dict[key] = name_dict
pivot_count_dict[key] = row_count_df


pivot_comp_df = pd.DataFrame.from_dict(pivot_comp_dict)
pivot_count_df = pd.DataFrame.from_dict(pivot_count_dict)





# Create subset dataframes with all the non-intersecting categories removed, tables of just intersecting categories??