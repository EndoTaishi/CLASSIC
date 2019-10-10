#!/home/crp102/.conda/envs/python_dev/bin/python3.7
'''
This script is used to generate comparative plots between fluxnet observation data
and simulated data.

Ensure the outputs from the fluxnet run are referenced properly. The simulatedData variable
should point directly to the 'outputs' folder generated by the fluxnet run using submission
script with ID 93c72a01fe4a0351 (check the ppp $RUNPATH if you can't find it)

Ensure the observational files are also referenced properly.The observationalData variable
should point directly to the 'outputs' folder generated using the fluxnet conversion script
with ID 8f08daf22e65f5d5

Title: comparative_plot_generator.py
Use: $ python3 comparative_plot_generator.py
Author: Matthew Fortier
Script ID: d1f7cac6a0096afe
'''

import os
import sys
import math
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import ticker
import statistics
import scipy
from datetime import datetime
from sklearn.metrics import mean_squared_error, r2_score

#######################################################
output_dir = ""
observationalData = "/home/acrp001/mattfortier2019/code/CLASSIC/inputFiles/observationalDataFLUXNET"
simulatedRuns = []

#######################################################

'''
Variables

This array contains a dictionary for each variable you wish to compare. Each dictionary
has the following attributes:

- observational_attribute should be the name of the observational attribute being compared, and
should correspond to a {observational_attribute}_monthly_fluxnet.csv file. It also dictates how
the final graph will be labeled

- observational_multiplier should be a scalar value used to multiply with the observatinoal
data if needed. It this step is unnecessary, it may be set to 1.

- units indicates the string that will be displayed beside the attribute on the final graph. If
the substring "month" appears within {units}, the input data will be multiplied by 1 / 8.64E7.

- seconds_to_months should be True if you need to convert from a value per second, to a value
per month. In general, fluxes should mark this as True, whereas pools should mark this as False

- input_attributes is a list that will contain one or more simulated attributes to be used
in linear combination to compare to the observational attribute in question. Each string
in the input_attributes list should correspond to a {input_attribute}_monthly.csv file
in each fluxnet site's simulation output.

- input_multipliers is a list that will contain the scalar values which will be multiplied
with the corresponding values in the input_attributes list. If no multiplication is necessary,
it can be set to an array of 1's the same length as the input_attributes list.

- show is a string that can be set to either "model", "obs", or "both". This determines what
plots will be generated, and what will be included in those plots. For example, if you only want
to see what the observational value of a variable is, select "obs". This will disable scatter plots
for that variable and just show the observational value. If you want to see everything, set it to
"both"
'''

#######################################################
variables = [
    {
        "observational_attribute": "gpp",
        "observational_multiplier": 1,
        "units": "kg C/m$^2$/month",
        "seconds_to_months": True,
        "input_attributes": ["gpp"],
        "input_multipliers": [1],
        "show": "both"

    },
    {
        "observational_attribute": "reco",
        "observational_multiplier": 1,
        "units": "kg C/m$^2$/month",
        "seconds_to_months": True,
        "input_attributes": ["ra", "rh"], # reco is 1*ra + 1*rh
        "input_multipliers": [1, 1],
        "show": "both"

    },
    {
        "observational_attribute": "nee",
        "observational_multiplier": -1, # "nee" needs to be inverted
        "units": "kg C/m$^2$/month",
        "seconds_to_months": True,
        "input_attributes": ["nep"],
        "input_multipliers": [1],
        "show": "both"

    },
    {
        "observational_attribute": "hfss",
        "observational_multiplier": 1,
        "units": "W/m$^2$",
        "seconds_to_months": False,
        "input_attributes": ["hfss"],
        "input_multipliers": [1],
        "show": "both"

    },
    {
        "observational_attribute": "hfls",
        "observational_multiplier": 1,
        "units": "W/m$^2$",
        "seconds_to_months": False,
        "input_attributes": ["hfls"],
        "input_multipliers": [1],
        "show": "both"

    },
    {
        "observational_attribute": "rns",
        "observational_multiplier": 1,
        "units": "W/m$^2$",
        "seconds_to_months": False,
        "input_attributes": ["rss", "rls"],
        "input_multipliers": [1, 1],
        "show": "both"

    },
    {
        "observational_attribute": "hfg",
        "observational_multiplier": 1,
        "units": "W/m$^2$",
        "seconds_to_months": False,
        "input_attributes": ["hfg"],
        "input_multipliers": [1],
        "show": "both"
    },
    {
        "observational_attribute": "imbalance",
        "observational_multiplier": 1,
        "units": "W/m$^2$",
        "seconds_to_months": False,
        "input_attributes": ["rss", "rls", "hfg", "hfss", "hfls"],
        "input_multipliers": [1, 1, -1, -1, -1],
        "show": "obs"
    }
]
#######################################################


def main():
    avail_colours = ['k', 'm', 'c', 'r', 'g', 'b']

    colours = ['k']
    for run in simulatedRuns:
        colours.append(run["colour"])

    for colour in colours:
        avail_colours.remove(colour)


    # Check if the output directory is specified in the args
    global output_dir
    for i, arg in enumerate(sys.argv):
        if arg == "-o":
            output_dir = sys.argv[i+1]
            sys.argv.remove(arg)
            sys.argv.remove(output_dir)

    if output_dir == "":
        print("Error parsing output_dir")
        sys.exit(1)

    # Make sure the top-level output directory exists. Otherwise, create it.
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Add any command-line output folders to the list of simulated runs
    if len(sys.argv) > 1:
        for item in sys.argv[1:]:
            sim = {}
            sim["outputs"] = item
            name = re.match(r".*\/([^\/]*)\/outputs", item)
            if not name:
                sim["name"] = "CLASSIC"
            else:
                sim["name"] = name.group(1)
            sim["colour"] = avail_colours.pop(0)
            simulatedRuns.append(sim)
            colours.append(sim["colour"])

    # Iterate through the entire script for each variable
    for var in variables:
        print("Processing " + var["observational_attribute"].upper() + " site data...")

        # Confirm that the obs input file exists, then open it
        obsFile = observationalData + "/" + var["observational_attribute"] + "_monthly_fluxnet.csv"
        if not os.path.isfile(obsFile):
            print("Could not find file: " + obsFile)
            sys.exit(1)
            df_obs_full = pd.read_csv(obsFile, usecols=["time", var["observational_attribute"], "sitename", "IGBPcode"])
        df_obs_full['time'] = pd.to_datetime(df_obs_full['time'])

        # Keep all adjusted dataframes in a single dictionary
        dataframes = {}
        # Keep track of every site's IGBP classification
        IGBP = {}
	# Dictionary to store statistics
        stats = {}
        # Obtain a list of all the sites in the simulated data
        sites = os.listdir(simulatedRuns[0]["outputs"])
        for site in sites:
            if len(site) > 6: # to exclude the DE-Tha. Remove later...
                continue
            df_obs = df_obs_full[df_obs_full["sitename"] == site]
            df_obs = df_obs.reset_index(drop=True)
            site_code = df_obs['IGBPcode'][0]
            # If we're able to obtain meaningful data from the comparison:
            if generate_comparative_frame(df_obs, site, dataframes, var):
                # Add the site into the collection of processed sites (IGBP)
                if site_code in IGBP.keys():
                    IGBP[site_code].append(site)
                else:
                    IGBP[site_code] = [site]

        if var["show"] == "both":
            print("Running statistics...")
            run_statistics(dataframes, IGBP, var, stats)
        print("Generating timeseries and scatter matrices...")
        generate_plot_matrices(dataframes, colours, IGBP, var, stats)

        #if var["observational_attribute"] == "gpp":
        #    generate_functional_relationship_plots(dataframes, colours, IGBP)


        print("Generating seasonal matrices...")
        seasonal_plots(dataframes, colours, IGBP, var)
        print("\n")



###########################################################################
# Process the observational and simulated site data into a single DataFrame
###########################################################################

def generate_comparative_frame(df_obs, site, dataframes, var):
    # Presently, we're skipping any site with a long name (needs to be renamed)
    if len(site) > 6:
        print("Skipping " + site + "...")
        return False

    simulatedFrames = []
    for run in simulatedRuns:

        # Check that the sim input file exists
        simFile = run["outputs"] + "/" + site + "/csv/" + var["input_attributes"][0] + "_monthly.csv"
        if not os.path.isfile(simFile):
            print("Could not find file: " + simFile)
            return False

        # Read simulation data into dataframes
        df_sim = pd.read_csv(simFile, usecols=["time", var["input_attributes"][0]])
        df_sim['time'] = pd.to_datetime(df_sim['time'])
        df_sim = df_sim.replace(-9999, np.nan)
        df_sim[var["input_attributes"][0]] = df_sim[var["input_attributes"][0]].apply(lambda x: x*var["input_multipliers"][0])
        df_sim = df_sim.rename(columns={var["input_attributes"][0]: var["observational_attribute"]})

        # Add other simulated files as needed
        if len(var["input_attributes"]) > 1:
            for i, att in enumerate(var["input_attributes"][1:]):
                simFile = run["outputs"] + "/" + site + "/csv/" + att + "_monthly.csv"
                if not os.path.isfile(simFile):
                    print("Could not find file: " + simFile)
                    return False
                df_att = pd.read_csv(simFile, usecols=["time", att])
                df_att['time'] = pd.to_datetime(df_att['time'])
                df_att = df_att.replace(-9999, np.nan)
                df_att[att] = df_att[att].apply(lambda x: x*var["input_multipliers"][i+1])
                df_sim[var["observational_attribute"]] = df_sim[var["observational_attribute"]] + df_att[att]

        # Add this dataFrame to the list of simulatedFrames
        df_sim = df_sim.rename(columns={var["observational_attribute"]: run["name"]})
        simulatedFrames.append(df_sim)

    # Only include dates where we have overlapping data
        # Drop all unnecessary columns
        df_obs = df_obs.drop(columns=['IGBPcode', 'sitename'])

        # Remove -9999 values (csv sentinel value for no observation), apply multipliers
        df_obs = df_obs[df_obs[var["observational_attribute"]] != -9999]
        df_obs = df_obs.reset_index(drop=True) # go back to indexing from 0
        #df_obs = df_obs.replace(-9999, np.nan)
        df_obs[var["observational_attribute"]] = df_obs[var["observational_attribute"]].apply(lambda x: x*var["observational_multiplier"])

        if min(simulatedFrames[0][simulatedRuns[0]["name"]].count(), df_obs[var["observational_attribute"]].count()) == 0:
            print("Missing data for " + site + ". Skipping...")
            return False
        df_obs.sort_values(by=['time'])
        for df_sim in simulatedFrames:
            df_sim.sort_values(by=['time'])
        mintime = df_obs['time'][0] if simulatedFrames[0]['time'][0] < df_obs['time'][0] else simulatedFrames[0]['time'][0]
        maxtime = simulatedFrames[0]['time'][simulatedFrames[0]['time'].count()-1] if simulatedFrames[0]['time'][simulatedFrames[0]['time'].count()-1] \
            < df_obs['time'][df_obs['time'].count()-1] else df_obs['time'][df_obs['time'].count()-1]
        df_obs = df_obs[df_obs['time'] >= mintime]
        df_obs = df_obs[df_obs['time'] <= maxtime]
        df_obs = df_obs.reset_index(drop=True) # go back to indexing from 0
        for df_sim in simulatedFrames:
            df_sim.drop(df_sim[df_sim['time'] < mintime].index, inplace=True)
            df_sim.drop(df_sim[df_sim['time'] > maxtime].index, inplace=True)
            df_sim.reset_index(drop=True, inplace=True) # go back to indexing from 0
        if min(simulatedFrames[0][simulatedRuns[0]["name"]].count(), df_obs[var["observational_attribute"]].count()) == 0:
            print("No overlapping data for " + site + ". Skipping...")
            return False

    # Rename column
    df_obs = df_obs.rename(columns={var["observational_attribute"]: 'observed'})

    # Index by time, then concatenate the datasets by time
        simulatedFrames.insert(0, df_obs)
    for df_sim in simulatedFrames:
        if var["seconds_to_months"]:
            seconds_to_months(df_sim)
        df_sim.set_index('time', inplace=True)

    df_combined = pd.concat(simulatedFrames, axis=1, sort=False)
    df_combined = df_combined.replace(pd.NaT, np.nan)
    f = df_combined.apply(lambda col: col.first_valid_index()).max()
    l = df_combined.apply(lambda col: col.last_valid_index()).min()
    df_combined = df_combined.loc[f:l]
    dataframes[site] = df_combined
    return True


########################
# Generate plot matrices
########################

# Generates matrix of plots for the dataframes, as well as a matrix of scatter plots if desired.
def generate_plot_matrices(dataframes, colours, IGBP, var, stats):
    # Generate plot matrix
    num_plots = len(dataframes)+1
    num_cols = 5
    num_rows = math.ceil(num_plots/num_cols)
    # Create the figure, make the size based on the number of plots in the matrix.
    fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, squeeze=False)
    fig.suptitle(f"{var['observational_attribute'].upper()}\n({var['units']})", fontsize=20, x=0.125, y=1-(1/(4*num_rows)))
    fig.set_figheight(3*num_rows)
    fig.set_figwidth(3*num_cols)
    i = 0 # this will keep track of the row we're on
    j = 1 # this will keep track of the column we're on (skip column 0 on the first row, the legend goes there)
    plt.delaxes(axes[0,0])
    # Iterate through each IGBP code (eg. ENF, OSH, etc.)
    for code in IGBP.keys():
        for site in IGBP[code]:
            for c, run in enumerate(list(dataframes[site])):
                if var["show"] == "model" and c == 0:
                    continue
                if var["show"] == "obs" and c > 0:
                    continue
                dataframes[site][run].plot(title=f"{site}({code})", ax=axes[i][j], style=".-", lw=0.5, ms=3, legend=False, color=colours[c])
            j += 1
            if j == num_cols:
                j = 0
                i += 1
    for ii, ax1 in enumerate(axes):
        for jj, ax in enumerate(ax1):
            ax.set_position([0.075 + jj*0.19, 1-((ii+1)/num_rows)+(0.2/num_rows), 0.15, (1/num_rows)*(0.6)])
            ax.xaxis.set_label_text("")
            if len(ax.xaxis.get_ticklabels()) > 4:
                for n, label in enumerate(ax.xaxis.get_ticklabels()):
                    if n % 2 != 0:
                        label.set_visible(False)
    for ii in range(num_plots, num_cols*num_rows):
        plt.delaxes(axes[math.floor(ii/num_cols), ii%num_cols])
    handles, labels = axes[0][1].get_legend_handles_labels()
    fig.legend(handles, labels, loc=(0.08, 1-(3/(4*num_rows))))
    plt.savefig(output_dir + "/" + var["observational_attribute"] + "_time_series.pdf")
    plt.close()
    # Create the scatter plot matrix, if we're using observational data
    if var["show"] == "both":
        handles = []
        labels = []
        fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, squeeze=False)
        fig.suptitle(f"{var['observational_attribute'].upper()}\n({var['units']})", fontsize=20, x=0.125, y=1-(1/(4*num_rows)))
        fig.set_figheight(3*num_rows)
        fig.set_figwidth(3*num_cols)
        i = 0 # this will keep track of the row we're on
        j = 1 # this will keep track of the column we're on (skip column 0 on the first row, the legend goes there)
        plt.delaxes(axes[0,0])
        for code in IGBP.keys():
            for site in IGBP[code]:
                if i == 0 and j == 1:
                    for k, run in enumerate(list(dataframes[site])[1:]):
                        handles.append(Line2D([0], [0], color=colours[k+1]))
                        labels.append(run)
                floor = df_floor(dataframes[site])
                ceiling = df_ceil(dataframes[site])
                axes[i][j].set_xlim([floor, ceiling])
                axes[i][j].set_ylim([floor, ceiling])
                for k, run in enumerate(list(dataframes[site])[1:]):
                    dataframes[site].plot(kind="scatter", x="observed", y=run, color=colours[k+1], s=3, ax=axes[i][j], title=f"{site}({code})", label=run, zorder=3, alpha=0.7, legend=False)
                fits = find_linear_fits(dataframes[site])
                for k, run in enumerate(list(dataframes[site])[1:]):
                    axes[i][j].text(x=0.02, y=0.9-0.1*k, s="y={:5.3f}x+{:4.1f} / {} / {}".format(fits[k][0], fits[k][1], stats[site][k]["r2"][:5], stats[site][k]["mse"][:5]), color=colours[k+1], transform=axes[i][j].transAxes)
                for k, m in enumerate(fits):
                    axes[i][j].plot([floor, ceiling], [floor*m[0]+m[1], ceiling*m[0]+m[1]], color=colours[k+1], zorder=4, alpha=0.7)
                axes[i][j].plot([floor, ceiling], [floor, ceiling], color=colours[0], zorder=1, alpha=0.7)
                axes[i][j].set_position([0.075 + j*0.19, 1-((i+1)/num_rows)+(0.2/num_rows), 0.15, (1/num_rows)*(0.6)])
                axes[i][j].xaxis.set_label_text("Observed")
                if j == 0 or (j == 1 and i == 0):
                    axes[i][j].yaxis.set_label_text("Simulated")
                else:
                    axes[i][j].yaxis.set_label_text("")
                j += 1
                if j == num_cols:
                    j = 0
                    i += 1
        for ii in range(num_plots, num_cols*num_rows):
            plt.delaxes(axes[math.floor(ii/num_cols), ii%num_cols])
        fig.text(0.072, 1-(6/(7*num_rows)), "line eqn / r$^2$ / rmse")
        fig.legend(handles, labels, loc=(0.083, 1-(3/(4*num_rows))))
        plt.savefig(output_dir + "/" + var["observational_attribute"] + "_scatter" + ".pdf")
        plt.close()


def seasonal_plots(dataframes, colours, IGBP, var):
    # Generate plot matrix
    num_plots = len(dataframes)+1
    num_cols = 5
    num_rows = math.ceil(num_plots/num_cols)
    # Create the figure, make the size based on the number of plots in the matrix.
    fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, squeeze=False)
    fig.suptitle(f"Seasonal {var['observational_attribute'].upper()}\n({var['units']})", fontsize=20, x=0.125, y=1-(1/(4*num_rows)))
    fig.set_figheight(3*num_rows)
    fig.set_figwidth(3*num_cols)
    i = 0 # this will keep track of the row we're on
    j = 1 # this will keep track of the column we're on (skip column 0 on the first row, the legend goes there)
    plt.delaxes(axes[0,0])
    custom_legend = []
    custom_labels = []
    # Iterate through each IGBP code (eg. ENF, OSH, etc.)
    for code in IGBP.keys():
        for site in IGBP[code]:
            df = dataframes[site]
            strt = df.index[0].year
            end = df.index[-1].year
            # Make the legend (but only once!)
            if len(custom_legend) == 0:
                runs = list(df)
                for c, run in enumerate(list(df)):
                    if var["show"] == "model" and c == 0:
                        continue
                    if var["show"] == "obs" and c > 0:
                        continue
                    custom_legend.append(Line2D([0], [0], color=colours[c]))
                    custom_labels.append(run)
            df = df.reset_index()
            df['Month'] = df['time'].dt.month
            df = df.sort_values(by=['Month'])
            for c, run in enumerate(runs):
                if var["show"] == "model" and c == 0:
                    continue
                if var["show"] == "obs" and c > 0:
                    continue
                sns.lineplot(data=df, x='Month', y=run, ax=axes[i][j], legend=False, ci="sd", color=colours[c], alpha=(0.6 if c == 0 else 1)).set_title(f"{site}({code})")
            # set all the axis properties
            axes[i][j].xaxis.set_major_locator(ticker.MaxNLocator(6))
            axes[i][j].set_xticklabels(["", "Feb","Apr","Jun","Aug","Oct","Dec"])
            axes[i][j].set_xlim([1,12])
            axes[i][j].set(xlabel="", ylabel="")
            axes[i][j].text(x=1, y=0.95, s=f"{str(strt)} - {str(end)} (n={str(end-strt)})", verticalalignment='top', horizontalalignment='right', transform=axes[i][j].transAxes)
            axes[i][j].set_position([0.075 + j*0.19, 1-((i+1)/num_rows)+(0.2/num_rows), 0.15, (1/num_rows)*(0.6)])
            axes[i][j].xaxis.set_label_text("Month")
            j += 1
            if j == num_cols:
                j = 0
                i += 1
    # delete the empty axes
    for ii in range(num_plots, num_cols*num_rows):
        plt.delaxes(axes[math.floor(ii/num_cols), ii%num_cols])
    handles, labels = axes[0][1].get_legend_handles_labels()
    fig.legend(custom_legend, custom_labels, loc=(0.08, 1-(3/(4*num_rows))))
    plt.savefig(output_dir + "/" + var["observational_attribute"] + "_seasonal.pdf")
    plt.close()


def generate_functional_relationship_plots(dataframes, colours, IGBP):
    df = pd.DataFrame()
    sns.set(style="whitegrid")
    for code in IGBP.keys():
        for site in IGBP[code]:

            df_combined = dataframes[site]

            rain_file = simulatedRuns[0]["outputs"] + "/" + site + "/csv/pr_monthly.csv"
            if not os.path.isfile(rain_file):
                print("Could not find monthly precipitation file for " + site)
                return False
            temp_file = simulatedRuns[0]["outputs"] + "/" + site + "/csv/tas_monthly.csv"
            if not os.path.isfile(temp_file):
                print("Could not find monthly temperature file for " + site)
                return False

            df_rain = pd.read_csv(rain_file, usecols=["time", "pr"])
            df_temp = pd.read_csv(temp_file, usecols=["time", "tas"])
            df_rain['time'] = pd.to_datetime(df_rain['time'])
            df_temp['time'] = pd.to_datetime(df_temp['time'])
            mintime = df_combined.index[0]
            maxtime = df_combined.index[-1]
            df_rain = df_rain[df_rain['time'] >= mintime]
            df_rain = df_rain[df_rain['time'] <= maxtime]
            df_temp = df_temp[df_temp['time'] >= mintime]
            df_temp = df_temp[df_temp['time'] <= maxtime]
            df_rain.set_index('time', inplace=True)
            df_temp.set_index('time', inplace=True)
            df_rain.loc[:,'pr'] *= 86400
            df_combined = df_combined.reset_index()
            df_combined = months_to_days(df_combined)
            df_combined.set_index('time', inplace=True)
            df_combined = pd.concat([df_combined, df_rain, df_temp], axis=1, sort=False)
            if df.shape[0] == 0:
                df = df_combined
            else:
                df = pd.concat([df, df_combined])
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 9))
    fig.suptitle("Functional Relationship - GPP and Precipitation")
    for c, run in enumerate(list(df)[:-2]):
        sns.scatterplot(data=df, x="pr", y=run, legend=False, ax=ax, size=0.1, color=colours[c], linewidth=0, alpha=0.3)
    ax.xaxis.set_label_text("Precip (mL/m$^2$/day)")
    ax.yaxis.set_label_text("GPP (g C/m$^2$/day)")
    ax.set_xscale("symlog")
    plt.savefig(output_dir + "/functional_precip.pdf")

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 9))
    fig.suptitle("Functional Relationship - GPP and Temperature")
    for c, run in enumerate(list(df)[:-2]):
        sns.scatterplot(data=df, x="tas", y=run, legend=False, ax=ax, size=0.1, color=colours[c], linewidth=0, alpha=0.3)
    ax.xaxis.set_label_text("Temp (C)")
    ax.yaxis.set_label_text("GPP (g C/m$^2$/day)")
    plt.savefig(output_dir + "/functional_temp.pdf")



##################
# Helper functions
##################

# Returns a list of linear fits for the dataframe (comparing the first column to every other column)
def find_linear_fits(df):
    points = df.values
    points = points[~np.isnan(points).any(axis=1)]
    fits = []
    def best_fit(X, Y):
        xbar = sum(X)/len(X)
        ybar = sum(Y)/len(Y)
        n = len(X) # or len(Y)

        numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
        denum = sum([xi**2 for xi in X]) - n * xbar**2

        b = numer / denum
        a = ybar - b * xbar
        return a, b
    x=points[:,0]
    y=points[:,1]
    a, b = best_fit(list(x),list(y))
    z = [b, a]
    fits.append(z)
    #fits.append(np.polyfit(x, y, 1))
    return fits


# Returns an appropriate floor for axis range from a dataframe
def df_floor(df):
    floor = np.nanmin(df.iloc[:, 0].values)
    for x in range(len(list(df)[1:])):
        p = np.nanmin(df.iloc[:, x+1].values)
        if p < floor:
            floor = p
    return 0 if floor > 0 else floor*1.1


# Returns an appropriate ceiling for axis range from a dataframe
def df_ceil(df):
    ceil = np.nanmax(df.iloc[:, 0].values)
    for x in range(len(list(df)[1:])):
        p = np.nanmax(df.iloc[:, x+1].values)
        if p > ceil:
            ceil = p
    return ceil*1.1 if ceil > 0 else 0


# If comparing to observational data, run statistics on the model accuracy
def run_statistics(dataframes, IGBP, var, stats):
    with open(output_dir + '/fluxnet_stats_' + var["observational_attribute"] + '.csv', 'w') as f:
        run_names = [x["name"] for x in simulatedRuns]
        # Write headers for the csv file
        f.write('site,biome,')
        for stat in ['RMSE', 'R2', 'NSE', 'Pval', 'Pconf']:
            for run in run_names:
                f.write(f'{stat}({run}),')
        f.write('\n')
        for code in IGBP.keys():
            for site in IGBP[code]:
                stats[site] = []
                f.write(site + "," + code + ",")
                df = dataframes[site]
                # We have to transform this data from a dataframe into a set of lists,
                # one for the observational data and one for each simulated run.
                raw_vals = df.values.tolist()
                formatted_vals = []
                # if n = months and k = number of sim runs + 1 (obs), we need to go
                # from "n x k" to "k x n"
                for item in raw_vals[0]:
                    formatted_vals.append([])
                for i, row in enumerate(raw_vals):
                    for val in row:
                        if math.isnan(val):
                            break
                    else:
                        for j, val in enumerate(row):
                            formatted_vals[j].append(val)
                # formatted_vals now contains 1 list for each set of data (obs or sim)
                # Write the stats for this site
                for run in formatted_vals[1:]:
                    stats[site].append({})
                for i, run in enumerate(formatted_vals[1:]):
                    mse = str(math.sqrt(mean_squared_error(formatted_vals[0], run)))
                    f.write(mse + ",")
                    stats[site][i]["mse"] = mse
                for i, run in enumerate(formatted_vals[1:]):
                    r2 = str(r2_score(formatted_vals[0], run))
                    f.write(r2 + ",")
                    stats[site][i]["r2"] = r2
                for i, run in enumerate(formatted_vals[1:]):
                    ns = str(nash_sutcliffe(formatted_vals[0], run))
                    f.write(ns + ",")
                    stats[site][i]["ns"] = ns
                for i, run in enumerate(formatted_vals[1:]):
                    pr1 = str(scipy.stats.pearsonr(formatted_vals[0], run)[0])
                    f.write(pr1 + ",")
                    stats[site][i]["pr1"] = pr1
                for i, run in enumerate(formatted_vals[1:]):
                    pr2 = str(scipy.stats.pearsonr(formatted_vals[0], run)[1])
                    f.write(pr2 + ",")
                    stats[site][i]["pr2"] = pr2
                f.write("\n")


# Quick implementation of Nash-Sutcliffe statistic
def nash_sutcliffe(obs, sim):
    mean_obs = sum(obs)/len(obs)
    numerator = sum((i[0]-i[1])**2 for i in zip(sim, obs))
    denominator = sum((i-mean_obs)**2 for i in obs)
    return numerator/denominator


def cumulative_monthly_average(odf):
    df = odf.groupby([lambda x: x.month]).mean().cumsum()
    return df


# Determine the number of seconds in a given month (for data transformation)
def seconds_to_months(df):
    month_days = {'January': 31, 'February': 28, 'March': 31, 'April': 30, 'May': 31, 'June': 30, 'July': 31, \
                  'August': 31, 'September': 30, 'October': 31, 'November': 30, 'December': 31}
    arr = df.values
    for row in arr:
        row[1] = row[1]*86400*month_days[row[0].month_name()]
    df.iloc[:,:] = arr
    return df

def months_to_days(df):
    month_days = {'January': 31, 'February': 28, 'March': 31, 'April': 30, 'May': 31, 'June': 30, 'July': 31, \
                  'August': 31, 'September': 30, 'October': 31, 'November': 30, 'December': 31}
    arr = df.values
    for row in arr:
        for i in range(1, len(row)):
            row[i] = row[i]*1000/month_days[row[0].month_name()]
    df.iloc[:,:] = arr
    return df

# Call main
if __name__ == "__main__":
  main()
