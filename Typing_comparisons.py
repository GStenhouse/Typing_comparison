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
def compair(piv_dict):
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
        if n_rows > n_cols or n_rows < n_cols:
            diff = n_cols - n_rows
            perc_diff = diff/(n_cols + n_rows)
        else:
            diff = 0
            perc_diff = "0"       
        if no_excess == 0:
            ave_excess = "0"
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
        if stats[key][2] != "0" and stats[key][2] <= threshold:
            similar.append(key)
    return(similar)


def year_sets(df, start_date = 2016, no_years = int(6)):
    subsets = {}
    for i in range(0,no_years):  # create subsets by year
        subdf = df[(df["Year"].between(start_date,(start_date+i), inclusive="both"))]     # create year based subsets
        subsets[i] = subdf
    return(subsets)


def identifiers(piv_dict):
    out = []
    pair_list = []
    key_set = []            # create an empty list to fill with a unique set of typing scheme names
    for pair, piv in piv_dict.items():
        pair_list.append(pair)
        keys = pair.split(" - ")
        key_set = key_set + keys
    key_set = set(key_set)
    out.append(pair_list)
    out.append(key_set)
    return(out)


def sizing(df, piv_dict, subsetting=False, start_date = 2016, no_years = int(6), id="Accession"):
    ident = identifiers(piv_dict)
    if subsetting:
        subsets = (year_sets(df, start_date, no_years))
        g_size_dict = {}        # create an empty dictionary to collect group sizes by typing scheme - for each subset
        for i in range(0,6):  # create subset dfs
            df = subsets[i]
            g_size = {}
            for column in df:
                data = df[column]
                size = data.value_counts(dropna=False)     # get df of row totals
                g_size[column] = size
            g_size_dict[i] = g_size
                
        gs_dfs = {}     # create empty dict to contain dfs of group size for year groups key = subset number 
        for i in range(0,no_years):       
            df = subsets[i]
            n_df = df[[id, "Year"]]        # create new series / dataframe of column with isolate IDs and year
            for key in ident[1]:
                data = []
                tmp_df = df[key]
                for cell in tmp_df:
                    size = g_size_dict[i][key].loc[cell]        # extract group size based on which group each isolate belongs to
                    data.append(size)
                ns = pd.Series(data,name=key)     # turn sizes in to a series, with typing level name as column name
                n_df = pd.concat([n_df,ns], axis=1)     # add new series to the new dataframe
            gs_dfs[i] = n_df
        out = gs_dfs
    
    else:
        g_size_dict = {}        # create an empty dictionary to collect group sizes by typing scheme 
        for pair in ident[0]:
            keys = pair.split(" - ")
            i_key = keys[1]
            g_size = piv_dict[pair].loc["All"]     # collect group sizes
            g_size_dict[i_key] = g_size

        n_df = df[[id, "Year"]]        # create new series / dataframe of column with isolate IDs and year
        for key in ident[1]:
            data = []
            o_df = df[key]
            for cell in o_df:
                size = g_size_dict[key].loc[cell]        # extract group size based on which group each isolate belongs to
                data.append(size)
            ns = pd.Series(data,name=key)     # turn sizes in to a series, with typing level name as column name
            n_df = pd.concat([n_df,ns], axis=1)     # add new series to the new dataframe
        out = n_df
    return(out)


# define plotting functions - first density plots and the histograms - to be used within other functions to generated the plots in the same way with various input data
def d_plot(input, legend=None, title=None, x_lab=None, y_lab=None, filename="d_plot.png", colours = "black", no_years = int(6)):
    sns.set(style="darkgrid")
    for i in range(0,no_years):  # create subsets by year
        fig = sns.kdeplot(input[i][title], fill=True, color=colours[i], warn_singular=False)          # create density plots with each year having a different hade curve
    plt.legend(legend, loc='upper center')
    plt.title(title)
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.savefig(filename)     # save plots as image files
    plt.close()

def h_plot(input, legend=None,title=None, x_lab=None, y_lab=None, filename="h_plot.png", colours = "black", no_years = int(6)):
    Fig, axs = plt.subplots(1)
    x_multi = []
    for i in range (0,no_years):
        x_multi.append(input[i][title])
    axs.hist(x_multi, density=True, color=colours)
    axs.legend(legend, loc='upper center')
    axs.set_title(title)
    axs.set_xlabel(x_lab)
    axs.set_ylabel(y_lab)
    plt.savefig(filename)     # save plots as image files
    plt.close()


# the previously defined functions of generating yearly subsets, and define groups sizes and plotting functions to generate plots showing yearly breakdown of group size frenquency for each typing level
def density_plt(df, piv_dict, xlab=None, ylab=None, subsetting=False, no_years = int(6), colours = list(("navy", "mediumblue", "royalblue", "cornflowerblue", "deepskyblue","lightblue")), leg = list(("2016","2017","2018","2019","2020","2021"))):
    ident = identifiers(piv_dict)
    size_dfs = sizing(df, piv_dict)
    if subsetting:    
        for key in ident[1]:
            f_name = "density_plots/" + str(key) + "_split.png"
            d_plot(size_dfs, legend=leg, title=key, x_lab=xlab, y_lab=ylab, filename=f_name, colours=colours, no_years = no_years)
            
            f_name_h = "density_plots/" + str(key) + "_split_hist.png"
            h_plot(size_dfs, legend=leg, title=key, x_lab=xlab, y_lab=ylab, filename=f_name_h, colours=colours, no_years = no_years)
    else:
        for key in ident[1]:
            subsets = []        # create empty list for collecting year based subsets
            df = size_dfs[[key,"Year"]]
            subsets = year_sets(df, start_date = 2016, no_years = int(6))
            
            f_name = "density_plots/" + str(key) + ".png"
            d_plot(subsets, legend=leg, title=key, x_lab=xlab, y_lab=ylab, filename=f_name, colours=colours, no_years = no_years)
            
            f_name_h = "density_plots/" + str(key) + "_hist.png"
            h_plot(subsets, legend=leg, title=key, x_lab=xlab, y_lab=ylab, filename=f_name_h, colours=colours, no_years = no_years)
  

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
gran_comp_list = compair(pivot_dicts[1])

# generate comparison stats
# find the number of times, per pivot pair, that there are exact matches, as a proportion of the total number of matches
pair_list = itertools.permutations(columns, r=2)
comp_stats_out = comp_stats(pivot_dicts[1], gran_comp_list)

# compile comparisons in to dataframe
comp_stats_df = pd.DataFrame.from_dict(comp_stats_out, orient = "index", columns = ("Rows:Columns", "Difference", "Pecentage difference", "Frequency excess columns", "Number of excess columns", "Average excess"))

# and then send to csv file
comp_stats_df.to_csv("gran_comp.csv", sep=",")

# identify similar typing levels
similars = simi(pivot_dicts[1], comp_stats_out)

# create df with a group size for each isolate and typing level and generate density plots - with group size along x-axis : sizing based on whole dataset
density_plt(df_in, pivot_dicts[0], subsetting=False, xlab = "Group size (number of isolates)", ylab = "Proportion on population")

# send group size dataframe to file for validataion
g_size_out_df = sizing(df_in, pivot_dicts[0], subsetting=False, start_date = 2016, no_years = int(6))
g_size_out_df.to_csv("g_size.csv", sep=",")

# Re-do density plots, this time with group size determined for each subset, rather than from total dataset
density_plt(df_in, pivot_dicts[0], subsetting=True, xlab = "Group size (number of isolates)", ylab = "Proportion on population")

# send group size dataframe to file for validataion
split_g_size_out_df = sizing(df_in, pivot_dicts[0], subsetting=True, start_date = 2016, no_years = int(6))
split_g_size_out_df.to_csv("g_size_split.csv", sep=",")

# Theil's U 
ignore_2 = ["Accession","Uberstrain", "HC1100 (cgST Cplx)", "HC2350 (subsp.)", "HC2000", "HC1500", "Name"]

associations(df_in, nominal_columns='auto', numerical_columns=None, mark_columns=False, nom_nom_assoc='theil', num_num_assoc='pearson', nom_num_assoc='correlation_ratio', symmetric_nom_nom=True, symmetric_num_num=False, hide_rows=ignore_2, hide_columns=ignore_2, cramers_v_bias_correction=True, nan_strategy="replace", nan_replace_value="0", ax=None, figsize=(25,25), annot=True, fmt='.2f', sv_color='silver', cbar=True, vmax=1.0, vmin=0.0, plot=True, compute_only=False, clustering=False, title=None, multiprocessing=True)

