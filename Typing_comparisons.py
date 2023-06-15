## This is a script to take an excel spreadsheet of isolates and their types (by various typing schemes) and
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
import sys
from xlsxwriter import Workbook

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###############
#  FUNCTIONS  #
###############

# To create a dictionary of pairwise pivot tables
def pairingtab(df_in, list, ignore=None): # ignore = list of columns to ignore in the pairwise comparisons
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
    return(out)


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
    return(out)


# generate pairwise comparison stats of granular similarity between typing schemes
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


# identify pairs within a certain threshold similarity - based on the number of isolates in each group
def simi(piv_dict, stats, threshold=0.15):
    similar = []
    for key in piv_dict.keys(): 
        if stats[key][2] != "0" and stats[key][2] <= threshold:
            similar.append(key)
    return(similar)


# function for subsetting dataframe by year
def year_sets(df, start_date = 2016, no_years = int(6), inclusive = "both"):
    subsets = {}
    for i in range(0,no_years):  # create subsets by year
        subdf = df[(df["Year"].between(start_date,(start_date+i), inclusive=inclusive))]     # create year based subsets
        subsets[i] = subdf
    return(subsets)


# dataframe for collecting column names and pairs from cross tabulation to drive for loops iterating through columns / typing scheme level names
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


# function for determining the group size for each isolates across all typing scheme levels (columns) and generating new dataframe(s) with the group size rather than name for each isolate - for determining frequency by group size
def sizing(df, piv_dict, subsetting=False, start_date = 2016, no_years = int(6), id="Accession", inclusive = "both"):
    ident = identifiers(piv_dict)
    if subsetting:
        subsets = (year_sets(df, start_date=start_date, no_years=no_years, inclusive=inclusive))
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


# define plotting functions - first density plots and then histograms - to be used within other functions to generated the plots in the same way with various input datasets
# density plots
def d_plot(input, legend=None, title=None, x_lab=None, y_lab=None, filename="d_plot.png", colours = "black", no_years = int(6)):
    sns.set(style="darkgrid")
    for i in range(0,no_years):  # create subsets by year
        data = input[i][title]
        data = data.dropna()
        if data.empty == True:
            break
        else:
            fig = sns.kdeplot(data, fill=True, color=colours[i], warn_singular=False)          # create density plots with each year having a different hade curve
    plt.legend(legend, loc='upper center')
    plt.title(title)
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.savefig(filename, bbox_inches = "tight")     # save plots as image files
    plt.close()

# histograms
def h_plot(input, legend=None,title=None, x_lab=None, y_lab=None, filename="h_plot.png", colours = "", no_years = int(6)):
    x_multi = []
    c = []
    for i in range (0,no_years):
        data = input[i][title]
        data = data.dropna()
        if data.empty == True:
            break
        else:
            x_multi.append(data)
            c.append(colours[i])
        Fig, axs = plt.subplots(1)
        axs.hist(x_multi, density=True, color=c)
        axs.legend(legend, loc='upper center')
        axs.set_title(title)
        axs.set_xlabel(x_lab)
        axs.set_ylabel(y_lab)
        plt.savefig(filename, bbox_inches = "tight")     # save plots as image files
        plt.close()


# the previously defined functions of generating yearly subsets, and define groups sizes and plotting functions to generate plots showing yearly breakdown of group size frenquency for each typing level
def density_plt(df, piv_dict, subsetting=False, xlab=None, ylab=None, no_years = int(6), colours = list(("navy", "mediumblue", "royalblue", "cornflowerblue", "deepskyblue","lightblue")), leg = list(("2016","2017","2018","2019","2020","2021"))):
    ident = identifiers(piv_dict)
    size_dfs = sizing(df, piv_dict, subsetting=subsetting)
    if subsetting:    
        for key in ident[1]:
            f_name = "density_plots/" + str(key) + "_split.png"
            d_plot(size_dfs, legend=leg, title=(key), x_lab=xlab, y_lab=ylab, filename=f_name, colours=colours, no_years=no_years)
            
            f_name_h = "density_plots/" + str(key) + "_split_hist.png"
            h_plot(size_dfs, legend=leg, title=(key), x_lab=xlab, y_lab=ylab, filename=f_name_h, colours=colours, no_years=no_years)
    else:
        for key in ident[1]:
            subsets = []        # create empty list for collecting year based subsets
            df = size_dfs[[key,"Year"]]
            subsets = year_sets(df, start_date = 2016, no_years = int(6))
            
            f_name = "density_plots/" + str(key) + ".png"
            d_plot(subsets, legend=leg, title=key, x_lab=xlab, y_lab=ylab, filename=f_name, colours=colours, no_years = no_years)
            
            f_name_h = "density_plots/" + str(key) + "_hist.png"
            h_plot(subsets, legend=leg, title=key, x_lab=xlab, y_lab=ylab, filename=f_name_h, colours=colours, no_years=no_years)
  

# function to assess the number of small groups and numbers of isolates witin them identified each year, with option to generate density plots / histograms
def size_stats(df, piv_dict, threshold = int(10), greater_than = True, plotting = False, subsetting = True, inclusive = "neither", start_date = int(2016), no_years = int(6), colours = list(("navy", "mediumblue", "royalblue", "cornflowerblue", "deepskyblue","lightblue")), leg = list(("2016","2017","2018","2019","2020","2021")), order = ""):
    ident = identifiers(piv_dict)
    size_dfs = sizing(df, piv_dict, subsetting=False)
    if not order: 
        ident = identifiers(piv_dict)
        order = ident[1] 
    # create df subsets using threshold to exclude those not of interest
    def value_sets(in_df, size_df, greater_than=greater_than,threshold=threshold):
        n_df = in_df["Year"] 
        for key in order:
            a = in_df[key]
            b = size_df[key]
            col = key + "_gsize"
            df = pd.concat([a,b], keys=[key,col], axis = 1)
            if greater_than:
                subset = df[(df[col] >= threshold)]
                n_df = pd.concat([n_df,subset], axis =1)
            else:
                subset = df[(df[col] <= threshold)]
                n_df = pd.concat([n_df,subset], axis =1)
        return(n_df)
    v_df = value_sets(df, size_dfs)
    
    # create yearly subset dfs of group sizes according to threshold and group name dataframe 
    if subsetting:
        if inclusive == "neither":
            subs = {}
            for i in range(0, no_years):
                year = start_date + i
                subset = v_df[(v_df["Year"] == year)]
                subs[i] = subset        
        elif inclusive == "both":
            subs = year_sets(v_df, start_date=start_date, no_years=no_years)
            for i in range(0, no_years):
                year = start_date + i
                subset = subs[i]
        else:
            sys.exit("Error inclusive can only be both or neither")
        # count the number of unique group names for each column for each subset
        s_list = []
        y_list = []
        for i in range(0, no_years):
            subset = subs[i]
            year = start_date + i
            y_list.append(year)
            s = pd.Series(dtype=int)
            for key in order:
                s[key] = (subset[key].nunique(dropna=True))
            s_list.append(s)
        n = pd.DataFrame(s_list, index = y_list)
        out = n
    else:
        n = v_df.nunique(axis=0, dropna=True)
        out = pd.DataFrame(n)

    # generate plots
    if plotting:     
        if greater_than:
            relative = "greater_than"
            rel = "greater than"
        else:
            relative = "less_than"
            rel = "less than"
        ylab = " ".join(("Number of groups with n", rel, str(threshold)))
        if subsetting:
            # create single bar chart of total n all years
            if inclusive == "both":
                inc = "_inc_"
            else:
                inc = "_exc_"
            N = n.sum(axis=0)      
            N.plot.bar()
            plt.ylabel(str(ylab))
            filenam = "density_plots/" + "_groupsize_" + relative + "_" + str(threshold) + inc + "subset.png"
            plt.savefig(filenam, bbox_inches = "tight")
            plt.close()
            # create bar plots coloured by year
            u = n.iloc[0:no_years].transpose()
            if inclusive == "neither":
                stack = True
            else:
                stack = False
            u.plot.bar(color=colours, stacked=stack)
            plt.ylabel(str(ylab))
            filename = "density_plots/" + "_groupsize_" + relative + "_" + str(threshold) + inc + "year_stack.png"
            plt.savefig(filename, bbox_inches = "tight")
            plt.close()
            # level specific plots
            for key in order:
                D = n[key].iloc[0:no_years]
                D.plot.bar(color=colours)
                plt.ylabel(str(ylab))
                plt.title(key)
                if inclusive == "both":
                    f_name = "density_plots/" + str(key) + "_groupsize_" + relative + "_" + str(threshold) + inc + "_inclu_years_split.png"
                else:    
                    f_name = "density_plots/" + str(key) + "_groupsize_" + relative + "_" + str(threshold) + inc + "sing_years_split.png"
                plt.savefig(f_name, bbox_inches = "tight")
                plt.close()
        else:
            filenam = "density_plots/" + "_groupsize_" + relative + "_" + str(threshold) + "not_subset.png"     
            n.plot.bar()
            plt.ylabel(str(ylab))
            plt.savefig(filenam, bbox_inches = "tight")
            plt.close() 
    return(out)            
           

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##########
#  MAIN  #
##########

# Read in spreadsheet
df_in = pd.read_excel('~/OneDrive - The University of Liverpool (1)/Scripting/S_sonnei_spreadsheet.xlsx')

# set up necessary variables
columns = df_in.columns.to_list()

Ignore = ["Accession","Uberstrain","Year","Sex","Age Group", "Foreign Travel","Continent of Travel", "Genotype Name", "blaCTX-M-27 gene", "Name", "Epi cluster"]

ignore_2 = ["Accession","Uberstrain", "HC1100 (cgST Cplx)", "HC2350 (subsp.)", "HC2000", "HC1500", "Name"]

hues = list(("navy", "mediumblue", "royalblue", "cornflowerblue", "deepskyblue","lightblue"))

Led = list(("2016","2017","2018","2019","2020","2021"))

order = list(("Genotype", "Lineage", "Clade", "Sub-Clade", "Sub-sub-Clade", "Sub-sub-sub-Clade", "250 SNP Threshold", "100 SNP Threshold", "50 SNP Threshold",  "25 SNP Threshold", "10 SNP Threshold", "5 SNP Threshold", "0 SNP Threshold", "SNP Address", "ST", "HC0 (indistinguishable)", "HC2", "HC5", "HC10", "HC20", "HC50", "HC100", "HC200", "HC400", "HC1100 (cgST Cplx)", "HC1500", "HC2000", "HC2350 (subsp.)"))


#############################################
## Analyses of full dataset

# create a series of pairwise comparison pivot tables.
# create list of column names as pairwise comparisons
pair_list = itertools.permutations(columns, r=2)

# use pairwise list to run pair wise pivot table creation
pivot_dicts = pairingtab(df_in, pair_list, ignore=Ignore)

# use output to identify matching groups - don't need to do separate by columns comp because pairingtab create list with pair both ways round
gran_comp_list = compair(pivot_dicts[1])

# generate comparison stats
# find the number of times, per pivot pair, that there are exact matches, as a proportion of the total number of matches
pair_list = itertools.permutations(columns, r=2)
comp_stats_out = comp_stats(pivot_dicts[1], gran_comp_list)

# compile comparisons in to dataframe
comp_stats_df = pd.DataFrame.from_dict(comp_stats_out, orient = "index", columns = ["Rows:Columns", "Difference", "Pecentage difference", "Frequency excess columns", "Number of excess columns", "Average excess"])

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

writer = pd.ExcelWriter("g_size_split.xlsx", engine='xlsxwriter')
for sheet, frame in split_g_size_out_df.items():
    year = int(2016) + int(sheet)
    frame.to_excel(writer, sheet_name = str(year))
writer.close()

# Theil's U 
theil = associations(df_in, nominal_columns='auto', numerical_columns=None, mark_columns=False, nom_nom_assoc='theil', num_num_assoc='pearson', nom_num_assoc='correlation_ratio', symmetric_nom_nom=True, symmetric_num_num=False, hide_rows=ignore_2, hide_columns=ignore_2, cramers_v_bias_correction=True, nan_strategy="replace", nan_replace_value=0, ax=None, figsize=(25,25), annot=True, fmt='.2f', sv_color='silver', cbar=True, vmax=1.0, vmin=0.0, plot=True, compute_only=False, clustering=False, title=None, multiprocessing=True)

pd.DataFrame(theil.items()).to_csv("theils_U.csv", sep=",")

# Group size stats by threshold
thresh_dict = {}

thresh_dict["less_10_inc"] = size_stats(df_in, pivot_dicts[0], threshold = int(10), greater_than = False, plotting = True, subsetting = True, inclusive = "both", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict["less_50_inc"] = size_stats(df_in, pivot_dicts[0], threshold = int(50), greater_than = False, plotting = True, subsetting = True, inclusive = "both", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict["great_10_inc"] = size_stats(df_in, pivot_dicts[0], threshold = int(10), greater_than = True, plotting = True, subsetting = True, inclusive = "both", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict["great_50_inc"] = size_stats(df_in, pivot_dicts[0], threshold = int(50), greater_than = True, plotting = True, subsetting = True, inclusive = "both", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

# single year subsets
thresh_dict["less_10_exc"] = size_stats(df_in, pivot_dicts[0], threshold = int(10), greater_than = False, plotting = True, subsetting = True, inclusive = "neither", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict["less_50_exc"] = size_stats(df_in, pivot_dicts[0], threshold = int(50), greater_than = False, plotting = True, subsetting = True, inclusive = "neither", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict["great_10_exc"] = size_stats(df_in, pivot_dicts[0], threshold = int(10), greater_than = True, plotting = True, subsetting = True, inclusive = "neither", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict["great_50_exc"] = size_stats(df_in, pivot_dicts[0], threshold = int(50), greater_than = True, plotting = True, subsetting = True, inclusive = "neither", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

# save group size dataframes into excel spreadsheets
writer = pd.ExcelWriter('thresholdstats.xlsx', engine='xlsxwriter')

for sheet, frame in thresh_dict.items():
    frame.to_excel(writer, sheet_name = sheet)

#critical last step
writer.close()


#############################################
## Analyses minus the identified likely-MSM clades

# Subset out MSM clades
df_less_MSM = df_in[df_in["Epi cluster"].isnull()]

# create a series of pairwise comparison pivot tables.
# create list of column names as pairwise comparisons
pair_list = itertools.permutations(columns, r=2)

# use pairwise list to run pair wise pivot table creation
pivot_dicts_2 = pairingtab(df_less_MSM, pair_list, ignore=Ignore)

# use output to identify matching groups - don't need to do separate by columns comp because pairingtab create list with pair both ways round
gran_comp_list_2 = compair(pivot_dicts_2[1])

# generate comparison stats
# find the number of times, per pivot pair, that there are exact matches, as a proportion of the total number of matches
pair_list = itertools.permutations(columns, r=2)
comp_stats_out_2 = comp_stats(pivot_dicts_2[1], gran_comp_list_2)

# compile comparisons in to dataframe
comp_stats_df_2 = pd.DataFrame.from_dict(comp_stats_out_2, orient = "index", columns = ["Rows:Columns", "Difference", "Pecentage difference", "Frequency excess columns", "Number of excess columns", "Average excess"])

# and then send to csv file
comp_stats_df_2.to_csv("gran_comp_wo_MSM.csv", sep=",")

# identify similar typing levels
similars_2 = simi(pivot_dicts_2[1], comp_stats_out_2)

# create df with a group size for each isolate and typing level and generate density plots - with group size along x-axis : sizing based on whole dataset
density_plt(df_less_MSM, pivot_dicts_2[0], subsetting=False, xlab = "Group size (number of isolates)", ylab = "Proportion on population")

# send group size dataframe to file for validataion
g_size_out_df_2 = sizing(df_less_MSM, pivot_dicts_2[0], subsetting=False, start_date = 2016, no_years = int(6))
g_size_out_df_2.to_csv("g_size_wo_MSM.csv", sep=",")

# Re-do density plots, this time with group size determined for each subset, rather than from total dataset
density_plt(df_less_MSM, pivot_dicts_2[0], subsetting=True, xlab = "Group size (number of isolates)", ylab = "Proportion on population")

# send group size dataframe to file for validataion
split_g_size_out_df_2 = sizing(df_less_MSM, pivot_dicts_2[0], subsetting=True, start_date = 2016, no_years = int(6))

writer = pd.ExcelWriter("g_size_split_wo_MSM.xlsx", engine='xlsxwriter')
for sheet, frame in split_g_size_out_df_2.items():
    year = int(2016) + int(sheet)
    frame.to_excel(writer, sheet_name = str(year))
writer.close()

# Theil's U 
theil_2 = associations(df_less_MSM, nominal_columns='auto', numerical_columns=None, mark_columns=False, nom_nom_assoc='theil', num_num_assoc='pearson', nom_num_assoc='correlation_ratio', symmetric_nom_nom=True, symmetric_num_num=False, hide_rows=ignore_2, hide_columns=ignore_2, cramers_v_bias_correction=True, nan_strategy="replace", nan_replace_value=0, ax=None, figsize=(25,25), annot=True, fmt='.2f', sv_color='silver', cbar=True, vmax=1.0, vmin=0.0, plot=True, compute_only=False, clustering=False, title=None, multiprocessing=True)

pd.DataFrame(theil_2.items()).to_csv("theils_U_wo_MSM.csv", sep=",")

# Group size stats by threshold
thresh_dict_2 = {}

thresh_dict_2["less_10_inc"] = size_stats(df_less_MSM, pivot_dicts_2[0], threshold = int(10), greater_than = False, plotting = True, subsetting = True, inclusive = "both", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict_2["less_50_inc"] = size_stats(df_less_MSM, pivot_dicts_2[0], threshold = int(50), greater_than = False, plotting = True, subsetting = True, inclusive = "both", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict_2["great_10_inc"] = size_stats(df_less_MSM, pivot_dicts_2[0], threshold = int(10), greater_than = True, plotting = True, subsetting = True, inclusive = "both", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict_2["great_50_inc"] = size_stats(df_less_MSM, pivot_dicts_2[0], threshold = int(50), greater_than = True, plotting = True, subsetting = True, inclusive = "both", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

# single year subsets
thresh_dict_2["less_10_exc"] = size_stats(df_less_MSM, pivot_dicts_2[0], threshold = int(10), greater_than = False, plotting = True, subsetting = True, inclusive = "neither", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict_2["less_50_exc"] = size_stats(df_less_MSM, pivot_dicts_2[0], threshold = int(50), greater_than = False, plotting = True, subsetting = True, inclusive = "neither", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict_2["great_10_exc"] = size_stats(df_less_MSM, pivot_dicts_2[0], threshold = int(10), greater_than = True, plotting = True, subsetting = True, inclusive = "neither", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict_2["great_50_exc"] = size_stats(df_less_MSM, pivot_dicts_2[0], threshold = int(50), greater_than = True, plotting = True, subsetting = True, inclusive = "neither", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

# save group size dataframes into excel spreadsheets
writer = pd.ExcelWriter('thresholdstats_wo_MSM.xlsx', engine='xlsxwriter')

for sheet, frame in thresh_dict_2.items():
    frame.to_excel(writer, sheet_name = sheet)

#critical last step
writer.close()


#############################################
## Analyses of the identified likely-MSM clades

# Subset out MSM clades
df_MSM = df_in[df_in["Epi cluster"].notnull()]

# create a series of pairwise comparison pivot tables.
# create list of column names as pairwise comparisons
pair_list = itertools.permutations(columns, r=2)

# use pairwise list to run pair wise pivot table creation
pivot_dicts_3 = pairingtab(df_MSM, pair_list, ignore=Ignore)

# use output to identify matching groups - don't need to do separate by columns comp because pairingtab create list with pair both ways round
gran_comp_list_3 = compair(pivot_dicts_3[1])

# generate comparison stats
# find the number of times, per pivot pair, that there are exact matches, as a proportion of the total number of matches
pair_list = itertools.permutations(columns, r=2)
comp_stats_out_3 = comp_stats(pivot_dicts_3[1], gran_comp_list_3)

# compile comparisons in to dataframe
comp_stats_df_3 = pd.DataFrame.from_dict(comp_stats_out_3, orient = "index", columns = ["Rows:Columns", "Difference", "Pecentage difference", "Frequency excess columns", "Number of excess columns", "Average excess"])

# and then send to csv file
comp_stats_df_3.to_csv("gran_comp_wo_MSM.csv", sep=",")

# identify similar typing levels
similars_3 = simi(pivot_dicts_3[1], comp_stats_out_3)

# create df with a group size for each isolate and typing level and generate density plots - with group size along x-axis : sizing based on whole dataset
density_plt(df_MSM, pivot_dicts_3[0], subsetting=False, xlab = "Group size (number of isolates)", ylab = "Proportion on population")

# send group size dataframe to file for validataion
g_size_out_df_3 = sizing(df_MSM, pivot_dicts_3[0], subsetting=False, start_date = 2016, no_years = int(6))
g_size_out_df_3.to_csv("g_size_wo_MSM.csv", sep=",")

# Re-do density plots, this time with group size determined for each subset, rather than from total dataset
density_plt(df_MSM, pivot_dicts_3[0], subsetting=True, xlab = "Group size (number of isolates)", ylab = "Proportion on population")

# send group size dataframe to file for validataion
split_g_size_out_df_3 = sizing(df_MSM, pivot_dicts_3[0], subsetting=True, start_date = 2016, no_years = int(6))

writer = pd.ExcelWriter("g_size_split_MSM.xlsx", engine='xlsxwriter')
for sheet, frame in split_g_size_out_df_3.items():
    year = int(2016) + int(sheet)
    frame.to_excel(writer, sheet_name = str(year))
writer.close()

# Theil's U 
theil_3 = associations(df_MSM, nominal_columns='auto', numerical_columns=None, mark_columns=False, nom_nom_assoc='theil', num_num_assoc='pearson', nom_num_assoc='correlation_ratio', symmetric_nom_nom=True, symmetric_num_num=False, hide_rows=ignore_2, hide_columns=ignore_2, cramers_v_bias_correction=True, nan_strategy="replace", nan_replace_value=0, ax=None, figsize=(25,25), annot=True, fmt='.2f', sv_color='silver', cbar=True, vmax=1.0, vmin=0.0, plot=True, compute_only=False, clustering=False, title=None, multiprocessing=True)

pd.DataFrame(theil_3.items()).to_csv("theils_U_MSM.csv", sep=",")

# Group size stats by threshold
thresh_dict_3 = {}

thresh_dict_3["less_10_inc"] = size_stats(df_MSM, pivot_dicts_3[0], threshold = int(10), greater_than = False, plotting = True, subsetting = True, inclusive = "both", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict_3["less_50_inc"] = size_stats(df_MSM, pivot_dicts_3[0], threshold = int(50), greater_than = False, plotting = True, subsetting = True, inclusive = "both", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict_3["great_10_inc"] = size_stats(df_MSM, pivot_dicts_3[0], threshold = int(10), greater_than = True, plotting = True, subsetting = True, inclusive = "both", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict_3["great_50_inc"] = size_stats(df_MSM, pivot_dicts_3[0], threshold = int(50), greater_than = True, plotting = True, subsetting = True, inclusive = "both", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

# single year subsets
thresh_dict_3["less_10_exc"] = size_stats(df_MSM, pivot_dicts_3[0], threshold = int(10), greater_than = False, plotting = True, subsetting = True, inclusive = "neither", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict_3["less_50_exc"] = size_stats(df_MSM, pivot_dicts_3[0], threshold = int(50), greater_than = False, plotting = True, subsetting = True, inclusive = "neither", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict_3["great_10_exc"] = size_stats(df_MSM, pivot_dicts_3[0], threshold = int(10), greater_than = True, plotting = True, subsetting = True, inclusive = "neither", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

thresh_dict_3["great_50_exc"] = size_stats(df_MSM, pivot_dicts_3[0], threshold = int(50), greater_than = True, plotting = True, subsetting = True, inclusive = "neither", start_date = int(2016), no_years = int(6), colours = hues, leg = Led, order = order)

# save group size dataframes into excel spreadsheets
writer = pd.ExcelWriter('thresholdstats_MSM.xlsx', engine='xlsxwriter')

for sheet, frame in thresh_dict_3.items():
    frame.to_excel(writer, sheet_name = sheet)

#critical last step
writer.close()
