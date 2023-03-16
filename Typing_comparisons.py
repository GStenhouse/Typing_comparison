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
import seaborn as sns
import matplotlib.pyplot as plt
from dython.nominal import associations

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###############
#  FUNCTIONS  #
###############

# To create a dictionary of pairwise pivot tables
def pairingtab(df_in, list, ignore=None):
    out = []
    ct_dict = {}
    ctm_dict = {}
    for pair in list:
        x = pair[0]
        y = pair[1]
        if x not in ignore and y not in ignore:
            pivot_m = pd.crosstab(df_in[x], df_in[y], margins=True)
            pivot = pd.crosstab(df_in[x], df_in[y])
            key = x + " - " + y
            ctm_dict[key] = pivot_m
            ct_dict[key] = pivot

    out.append(ctm_dict)
    out.append(ct_dict)
    return out


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
def simi(piv_dict, stats, threshold=0.15):
    similar = []
    for key in piv_dict.keys(): 
        if stats[key][2] != "na" and stats[key][2] <= threshold:
            similar.append(key)

    return(similar)


# Use crosstabs to create a new dataframe with the group size instead of group name for each isolate
# Then use this new ataframe to create density plots
def density_plt(df, piv_dict, id="Accession", xlab=None, ylab=None):
    key_set = []            # create an empty list to fill with a unique set of typing scheme names 
    g_size_dict = {}        # create an empty dictionary to collect group sizes by typing scheme 
    for pair, piv in piv_dict.items():
        keys = pair.split(" - ")
        i_key = keys[1]
        g_size = piv.loc["All"]     # collect group sizes
        g_size_dict[i_key] = g_size
        key_set = key_set + keys
    key_set = set(key_set)

    n_df = df[[id, "Year"]]        # create new series / dataframe of column with isolate IDs and year
    for key in key_set:
        data = []
        o_df = df[key]
        for cell in o_df:
            size = g_size_dict[key].loc[cell]        # extract group size based on which group each isolate belongs to
            data.append(size)
        ns = pd.Series(data,name=key)     # turn sizes in to a series, with typing level name as column name
        n_df = pd.concat([n_df,ns], axis=1)     # add new series to the new dataframe

    for key in key_set:
        subsets = []        # create empty list for collecting year based subsets
        df = n_df[[key,"Year"]]
        colours = ["navy", "mediumblue", "royalblue", "cornflowerblue", "deepskyblue","lightblue"]
        for i in range(0,6):  # create subsets by year
            subdf = df[(df["Year"].between(2016,(2016+i), inclusive="both"))]     # create year based subsets
            subsets.append(subdf)
        sns.set(style="darkgrid")
        for i in range(0,6):  # create subsets by year
            fig = sns.kdeplot(subsets[i][key], fill=True, color=colours[i], warn_singular=False)          # create density plots with each year having a different hade curve
        plt.legend(["2016","2017", "2018","2019","2020","2021"], loc='upper center')
        plt.title(key)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        filename = "density_plots/" + str(key) + ".png"
        plt.savefig(filename)     # save plots as image files
        plt.close()

        Fig, axs = plt.subplots(1)
        x_multi = []
        for i in range (0,6):
            x_multi.append(subsets[i][key])
        axs.hist(x_multi, density=True, color=colours)
        axs.legend(["2016","2017", "2018","2019","2020","2021"], loc='upper center')
        axs.set_title(key)
        axs.set_xlabel(xlab)
        axs.set_ylabel(ylab)
        filename_h = "density_plots/" + str(key) + "_hist.png"
        plt.savefig(filename_h)     # save plots as image files
        plt.close()
    
    return(n_df)        # collect group size dataframe for validation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##########
#  MAIN  #
##########

# Read in spreadsheet
df_in = pd.read_excel('~/OneDrive - The University of Liverpool (1)/Scripting/Typing_data.xlsx')

# create list of column names
columns = df_in.columns.to_list()

# create a series of pairwise comparison pivot tables.
# create list of column names as pairwise comparisons
pair_list = itertools.permutations(columns, r=2)

# use pairwise list to run pair wise pivot table creation
Ignore = ["Accession","Uberstrain","Year","Sex","Age Group", "Foreign Travel","Continent of Travel", "Genotype Name", "blaCTX-M-27 gene", "Name"]
pivot_dicts = pairingtab(df_in, pair_list, ignore=Ignore)

# use output to identify matching groups - don't need to do separate by columns comp because pairingtab create list with pair both ways round
gran_comp_list = compare(pivot_dicts[1])

# generate comparison stats
# find the number of times, per pivot pair, that there are exact matches, as a proportion of the total number of matches
pair_list = itertools.permutations(columns, r=2)
comp_stats_out = comp_stats(pivot_dicts[1], gran_comp_list)

# compile comparisons in to dataframe
comp_stats_df = pd.DataFrame.from_dict(comp_stats_out, orient = "index", columns = ("Rows:Columns", "Difference", "Pecentage difference", "Frequency excess columns", "Number of excess columns", "Average excess"))

# and then send to csv file
comp_stats_df.to_csv("test_gran_comp.csv", sep=",")

# identify similar typing levels
similars = simi(pivot_dicts[1], comp_stats_out)

# create df with a group size for each isolate and typing level and generate density plots - with group size along x-axis
g_size_out_df = density_plt(df_in, pivot_dicts[0], xlab = "Group size (number of isolates)", ylab = "Proportion on population")

# send group size dataframe to file for validataion
g_size_out_df.to_csv("g_size.csv", sep=",")

# Theil's U 
associations(df_in, nominal_columns='auto', numerical_columns=None, mark_columns=False, nom_nom_assoc='theil', num_num_assoc='pearson', nom_num_assoc='correlation_ratio', symmetric_nom_nom=True, symmetric_num_num=False, display_rows='all', display_columns='all', hide_rows=None, hide_columns=None, cramers_v_bias_correction=True, nan_strategy="replace", nan_replace_value="0", ax=None, figsize=None, annot=True, fmt='.2f', sv_color='silver', cbar=True, vmax=1.0, vmin=0.0, plot=True, compute_only=False, clustering=False, title=None, filename=None, multiprocessing=True, max_cpu_cores=None)

theil_df = associations(df_in, nominal_columns='auto', numerical_columns=None, mark_columns=False, nom_nom_assoc='theil', symmetric_nom_nom=False, hide_rows=Ignore, hide_columns=Ignore, multiprocessing=True, compute_only=True)

sns.heatmap(theil_df)
plt.show


    





