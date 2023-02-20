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
def pairingtab(df_in, list):
    ct_dict = {}

    for pair in list:
        x = pair[0]
        y = pair[1]
        pivot = pd.crosstab(df_in[x], df_in[y])

        key = x + " - " + y

        ct_dict[key] = pivot
    return ct_dict


# Identify matching groups, with and without isolate numbers
def compare(piv_dict, d):
    piv_comp_dict = {}            # create empty dicts to collect all  output data - pivot : row and matching columns
    piv_pair_dict = {}                # pivot : pair and number of matches
    out =[]
    for key, piv in piv_dict.items():
        row_names = piv.index.to_series() # create list of row names
        col_names = piv.columns.to_list()   # create list of column names
        pair_out = []              # create empty list to take the pairs and number objects
        match_dict = {}     # create empty dictionaries for collecting the data - list of matching columns by row
           
        if d in ["rows", "Rows", "by row", "By rows"] :
            for i in range(0,len(row_names)):    
                matches = []                       # create empty list for names of matching columns
                index = -1                          # set index value for extracting the matching column name
                name = row_names.iloc[i]            # collect row name
                data = piv.iloc[i]                      # extract row
                for cell in data.values.tolist():
                    index += 1              # add 1 to column number
        
                    if cell  > 0:
                        o_name = col_names[index]       # extract name of matching column
                        matches.append(o_name)        # add matching column name to list
                    else:
                        break

                    pair = str(name) + " - " + str(o_name) + ": " + str(cell)   # collect row, column and no. of matching isolates
                    pair_out.append(pair)           # collect in to a list            
                match_dict[name] = matches      # add row name and all matching column names to dict

        elif d in ["columns", "Columns", "by column", "By column"] :
            for i in range(0,len(col_names)):
                matches = []                       # create empty list for names of matching rows
                index = -1                          # set index value for extracting the matching row name
                name = col_names[i]                 # collect column name
                data = piv[name]                        # extract column
                for cell in data.values.tolist():
                    index += 1              # add 1 to row number

                    if cell  > 0 :
                        o_name = row_names.iloc[index]       # extract name of matching row
                        matches.append(o_name)        # add matching row name to list
                    else:
                        break

                    pair = str(name) + " - " + str(o_name) + ": " + str(cell)
                    pair_out.append(pair)        
                match_dict[name] = matches

        else:
            print("Error: No direction specified")
            return
        
        piv_comp_dict[key] = match_dict     # create dict of matches dict - by pivot name
        piv_pair_dict[key] = pair_out       # create dict of match lists - by pivot name
    out.append(piv_comp_dict)           # create a 2 item list of the two dictionaries - to be returned
    out.append(piv_pair_dict)
    return out


# generate pairwise comparison stats of granular similarity
def comp_stats(names, comparisons):
    out_dict = {}
    out_dict["Pair"] = "Total matches, Exact matches", "Excess matches", "Exact proportion", "Excess porportion", "Skew", "skew proportion"        # create empty dictionary to collect results


    for duo in names:
        pair = duo[0] + " - " + duo[1]
        matches = comparisons[0][pair]     # extract matches by pair
    
        tot_count = 0   # set counters
        exact_count = 0
        exceed_count = 0
        skew = 0


        for group in matches.values():        # extract specific row/column
            tot_count = tot_count + len(group)      # count total matches

            if len(group) == 1:         # count exact matches
                exact_count += 1
            
            elif len(group) > 1:          # count number of instances where 1 group matches with multiple groups
                exceed_count += 1
                skew = skew + (len(group) - 1)
    
        perc_act = exact_count/tot_count    # define totals as proportion of total matches
        perc_eed = exceed_count/tot_count
        perc_skew = skew/tot_count

        output = repr(tot_count) + ", " + repr(exact_count) + ", " + repr(exceed_count) + ", " + repr(perc_act) + ", " + repr(perc_eed) + ", " + repr(skew) + ", " + repr(perc_skew)
        out_dict[pair] = output
    return(out_dict)


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
pivot_dict = pairingtab(df_in, pair_list)

# use output to identify matching groups
row_gran_comp_list = compare(pivot_dict, d="rows")

col_gran_comp_list = compare(pivot_dict, d="columns")

# generate comparison stats
# find the number of times, per pivot pair, that there are exact matches, as a proportion of the total number of matches
pair_list = itertools.permutations(columns, r=2)
row_comp_stats = comp_stats(pair_list, row_gran_comp_list)

pair_list = itertools.permutations(columns, r=2)
column_comp_stats = comp_stats(pair_list, col_gran_comp_list)



    





