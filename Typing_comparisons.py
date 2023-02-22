# This is a script to take an excel spreadsheet of isolates and their types (by various typing schemes) and
# compare the granularity of the grouping by the various typing methods. 

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
def pairingtab(df_in, list, ignore=None):
    ct_dict = {}

    for pair in list:
        x = pair[0]
        y = pair[1]
        if x != ignore and y != ignore:
            pivot = pd.crosstab(df_in[x], df_in[y])
            key = x + " - " + y
            ct_dict[key] = pivot
    return ct_dict


# Identify matching groups, with and without isolate numbers
def compare(piv_dict):
    piv_comp_dict = {}            # create empty dicts to collect all  output data - pivot : row and matching columns
    piv_pair_dict = {}                # pivot : pair and number of matches
    out =[]
    for key, piv in piv_dict.items():
        row_names = piv.index.to_series() # create list of row names
        pair_out = []              # create empty list to take the pairs and number objects
        match_dict = {}     # create empty dictionaries for collecting the data - list of matching columns by row
           
        for i in range(0,len(row_names)):    
            matches = []                       # create empty list for names of matching columns
            index = -1                          # set index value for extracting the matching column name
            name = row_names.iloc[i]            # collect row name
            data = piv.iloc[i]                      # extract row
            for cell in data:
                index += 1              # add 1 to row number
                if cell  > 0 :
                    o_name = data.index[index]       # extract name of matching column
                    matches.append(o_name)        # add matching column name to list
                    pair = str(name) + " - " + str(o_name) + ": " + str(cell)   # collect row, column and no. of matching isolates
                    pair_out.append(pair)           # collect in to a list         
            match_dict[name] = matches      # add row name and all matching column names to dict
        
        piv_comp_dict[key] = match_dict     # create dict of matches dict - by pivot name
        piv_pair_dict[key] = pair_out       # create dict of match lists - by pivot name
    out.append(piv_comp_dict)           # create a 2 item list of the two dictionaries - to be returned
    out.append(piv_pair_dict)
    return out


# generate pairwise comparison stats of granular similarity
def comp_stats(piv_dict, in_list):
    out_dict = {}        # create empty dictionary to collect results

    for key, piv in piv_dict.items():
        col_names = piv.columns.to_series()     # extract column names
        row_names = piv.index.to_series()       # extract row names
        n_cols = len(col_names)   # determine the total number of isolates groups defined by typing scheme level - cols
        n_rows = len(row_names)   # determine the total number of isolates groups defined by typing scheme level - rows          
        output = []
        no_excess = 0       # set counters to zero
        vol_excess = 0

        for row in row_names.to_list():
            matches = len(in_list[0][key][row])
            if matches > 1:
                no_excess += 1
                vol_excess = vol_excess + matches    

        if n_rows > n_cols:
            diff = n_rows - n_cols
            perc_diff = diff/n_rows
        elif n_rows < n_cols:
            diff = n_cols - n_rows
            perc_diff = diff/n_cols
        else:
            diff = 0
            perc_diff = "na"
        
        if no_excess == 0:
            ave_excess = "na"
        
        else:
            ave_excess = vol_excess / no_excess

        ratio = str(n_rows) + ":" + str(n_cols)
        output.append(ratio)
        output.append(diff) 
        output.append(perc_diff)
        output.append(no_excess) 
        output.append(vol_excess)
        output.append(ave_excess)  
        out_dict[key] = output
    return(out_dict)

# identify pairs within a certain threshold similarity
def simi(piv_dict, stats, threshold=0.45):
    similar = []
    for key in piv_dict.keys():
        if stats[key][2] <= threshold:
            similar.append(key)

    return(similar)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##########
#  MAIN  #
##########

# Read in spreadsheet
df_in = pd.read_excel('~/OneDrive - The University of Liverpool (1)/Typing_test_data.xlsx')

# create list of column names
columns = df_in.columns.to_list()

# create a series of pairwise comparison pivot tables.
# create list of column names as pairwise comparisons
pair_list = itertools.permutations(columns, r=2)

# use pairwise list to run pair wise pivot table creation
pivot_dict = pairingtab(df_in, pair_list, ignore="Accession")

# use output to identify matching groups - don't need to do separate by columns comp because pairingtab create list with pair both ways round
gran_comp_list = compare(pivot_dict)

# generate comparison stats
# find the number of times, per pivot pair, that there are exact matches, as a proportion of the total number of matches
pair_list = itertools.permutations(columns, r=2)
comp_stats_out = comp_stats(pivot_dict, gran_comp_list)

# compile comparisons in to dataframe
comp_stats_df = pd.DataFrame.from_dict(comp_stats_out, orient = "index", columns = ("Rows:Columns", "Difference", "Pecentage difference", "Frequency excess columns", "Number of excess columns", "Average excess"))

# and then send to csv file
comp_stats_df.to_csv("test_gran_comp.csv", sep=",")

# identify similar typing levels
similars = simi(pivot_dict, comp_stats_out)


    



    





