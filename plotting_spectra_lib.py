import os
import numpy as np
import pandas as pd
# Plotlyc
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
# baseline removal function
from baseline_removal import baseline_als
from scipy.signal import savgol_filter
# matplotlib 
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
####
# The Goal of this file is to make standard figues
####
def plot_spectra(specimens, labels, plot_range=[[],[]],baseline_remove=True, baseline_params=[10e3,0.001],smooth=False,smooth_params=[9,3], norm=True, save=False, folder=""):
    """
    :params specimens should be a list of string paths
    :params labels should be a list of string
    :params title is the title of the figure
    :params norm if true the spectra are normalized
    """
    fields = ['Energy','Intensity']
    spectra = []
    baseline = ""
    # importing the spectra
    for specimen in specimens:
        # reading data from files
        df = pd.read_csv(str(specimen), skiprows=28, skipfooter=1,names=fields, engine='python', encoding = "ISO-8859-1")
        # !!!!! removing the first point that are uncorrect
        df = df[(df.Energy > 38)]
        # selecting only the range considered
        if plot_range[0] != []:
            df = df[(df.Energy < plot_range[0][1]) & (df.Energy > plot_range[0][0])]
            #print(plot_range[0])
        # baseline removal or/and smoothing with savgol filter
        if baseline_remove and smooth:
            baseline = baseline_als(savgol_filter(df.Intensity ,smooth_params[0],smooth_params[1]), lam=baseline_params[0], p=baseline_params[1])
            df.Intensity = df.Intensity - baseline
        elif baseline_remove:
            baseline = baseline_als(df.Intensity, lam=baseline_params[0], p=baseline_params[1])
            df.Intensity = df.Intensity - baseline
        else:
            # if baseline is not active it returns just the intensity
            #print("[+] baseline removal is not active")
            baseline = df.Intensity
        # Normalizing 
        if norm:
            df.Intensity = df.Intensity / max(df.Intensity)
        # adding to list of spectra
        spectra.append(df)
    ## Plotly figures
    fig = make_subplots(rows=1, cols=1,  subplot_titles=(''))
    # plotting the spectra
    for idx, spectrum in enumerate(spectra):

        fig.add_trace(
            go.Scatter(x=spectrum['Energy'], y= spectrum['Intensity'], name=labels[idx]),
            row=1, col=1
        )
        fig.add_trace(
            go.Scatter(x=spectrum['Energy'], y= baseline, name="baseline"),
            row=1, col=1
        )
        #print(spectrum)
    ## checking the inputs
    if plot_range[0] == []:
        # Update xaxis properties
        fig.update_xaxes(title_text="Energy (eV)", row=1, col=1,linecolor='black',gridcolor='lightgrey', ticks="outside")
    else:
        fig.update_xaxes(title_text="Energy (eV)", row=1, col=1,linecolor='black',gridcolor='lightgrey', range=plot_range[0], ticks="outside")
    if plot_range[1] == []:
        # Update yaxis properties
        fig.update_yaxes(title_text="Number of counts", row=1, col=1, linecolor='black',gridcolor='lightgrey', ticks="outside")
    else:
        fig.update_yaxes(title_text="Number of counts", row=1, col=1, linecolor='black',gridcolor='lightgrey', range=plot_range[1], ticks="outside")
    ## end check

    fig.update_layout(height=700, width=1400, title_text="", font=dict(size=25), 
        yaxis = dict(
        showexponent = 'all',
        exponentformat = 'power'
    ), plot_bgcolor="#fff", paper_bgcolor="#fff")
    # Save spectra in files in side a directory to be used with Origin
    if save==True:
        if not os.path.exists(folder):
            # Create the directory
            os.makedirs(folder)
        for idx, spectrum in enumerate(spectra):
            spectrum.to_csv(folder + "\\"+ labels[idx] +".csv", index=False)
    return fig, spectra

# plotting using matplotlib
def plot_spectra_matplotlib(specimens, labels, plot_range=[[],[]], range_step=2, fig_size=(), baseline_remove=True, baseline_params=[10e3,0.001],smooth=False,smooth_params=[9,3], offset=False, offset_value=[], norm=True, save=False, folder="", save_figure=False, folder_figure="", figure_name="figure"):
    fields = ['Energy','Intensity']
    spectra = []
    baseline = ""
    for specimen in specimens:
        # reading data from files
        df = pd.read_csv(str(specimen), skiprows=28, skipfooter=1,names=fields, engine='python', encoding = "ISO-8859-1")
        # !!!!! removing the first point that are uncorrect
        df = df[(df.Energy > 38)]
        # selecting only the range considered
        if plot_range[0] != []:
            df = df[(df.Energy < plot_range[0][1]) & (df.Energy > plot_range[0][0])]
            #print(plot_range[0])
        # baseline removal or/and smoothing with savgol filter
        if baseline_remove and smooth:
            baseline = baseline_als(savgol_filter(df.Intensity ,smooth_params[0],smooth_params[1]), lam=baseline_params[0], p=baseline_params[1])
            df.Intensity = df.Intensity - baseline
        elif baseline_remove:
            baseline = baseline_als(df.Intensity, lam=baseline_params[0], p=baseline_params[1])
            df.Intensity = df.Intensity - baseline
        else:
            # if baseline is not active it returns just the intensity
            #print("[+] baseline removal is not active")
            baseline = df.Intensity
        # Normalizing 
        if norm:
            df.Intensity = df.Intensity / max(df.Intensity)
        # adding to list of spectra
        spectra.append(df)
    ## Plotly figures
    fig = make_subplots(rows=1, cols=1,  subplot_titles=(''))
    # plotting the spectra
    font = {'size'   : 20, 'weight':'light'}
    matplotlib.rc('font', **font)
    if fig_size == ():
        fig, ax = plt.subplots(figsize=(16,12))
    else:
        fig, ax = plt.subplots(figsize=fig_size)
    for idx, spectrum in enumerate(spectra):
        # Create the x and y data
        # Create the plot
        if offset == True:
            plt.plot(spectrum.Energy + offset_value[idx], spectrum.Intensity, linewidth=3,label=labels[idx])
        else:
            plt.plot(spectrum.Energy, spectrum.Intensity, linewidth=3,label=labels[idx])
        
        #plt.plot(li_metal_5kv_100na_100s.Energy[li_metal_5kv_100na_100s.Energy < 70], baseline)
        # Set the axis labels and title
        plt.xlabel('Energy (eV)')
        plt.ylabel('Number of counts')

        # Set the ticks on both sides for x and y without numbers for the seconds
        ax.tick_params(which='both', direction='in', right=True, top=True)
        # For the minor ticks, use no labels; default NullFormatter.
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        # ticks size
        ax.tick_params(which='both', width=1)
        ax.tick_params(which='major', length=7)
        ax.tick_params(which='minor', length=4)

        
    ## checking the inputs
    if plot_range[0] == []:
        # Update xaxis properties
        pass
    else:
       plt.xlim(plot_range[0][0], plot_range[0][1])
       # Set the tick parameters
       plt.xticks(np.arange(plot_range[0][0], plot_range[0][1], range_step))
    if plot_range[1] == []:
        # Update yaxis properties
        pass
    else:
        plt.ylim(plot_range[1][0], plot_range[1][1])
    ## end check
    # Save spectra in files in side a directory to be used with Origin
    if save==True:
        if not os.path.exists(folder):
            # Create the directory
            os.makedirs(folder)
        for idx, spectrum in enumerate(spectra):
            spectrum.to_csv(folder + "\\"+ labels[idx] +".csv", index=False)
    plt.legend()
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
    # Save figure
    if save_figure == True:
        if not os.path.exists(folder_figure):
            os.makedirs(folder_figure)
        fig.savefig(folder_figure+"/"+figure_name+".jpg",dpi=600,bbox_inches = 'tight')
    return ax, fig