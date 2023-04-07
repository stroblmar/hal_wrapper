# ====================================================================================
# Functions to help with loading/formatting data for the CA simulations and for 
#  plotting the results.
# ====================================================================================
import numpy as np
import pandas as pd
import sys
import os
from math import isnan
import statsmodels.api as sm
from scipy import stats
if 'matplotlib' not in sys.modules:
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")
import myUtils as utils
from OnLatticeModel import OnLatticeModel

# ======================================== Visualisation ==========================================
def plot_simulation(dataDf, timeColumn="Time", feature='CA153', treatmentColumn="DrugConcentration",
              hue=None, style=None, legend=False, palette=None,
              plotDrug=True, plotDrugAsBar=True, drugBarPosition=0.85, treatmentScheduleList=None,
              drugColorMap={"Encorafenib": "blue", "Binimetinib": "green", "Nivolumab": sns.xkcd_rgb["goldenrod"]},
              xlim=None, ylim=None, y2lim=1,
              markInitialSize=False, markPositiveCutOff=False, plotHorizontalLine=False, lineYPos=1, despine=False,
              titleStr="", decorateX=True, decorateY=True, decorateY2=True,
              markersize=12, linestyle="None", linecolor='black',
              ax=None, figsize=(10, 8), outName=None, **kwargs):
    '''
    Plot longitudinal treatment data, together with annotations of drug administration and events responsible for
    changes in treatment dosing (e.g. toxicity).
    :param dataDf: Pandas data frame with longitudinal data to be plotted.
    :param timeColumn: Name (str) of the column with the time information.
    :param feature: Name (str) of the column with the metric to be plotted on the y-axis (e.g. PSA, CA125, etc).
    :param treatmentColumn: Name (str) of the column with the information about the dose administered.
    :param plotDrug: Boolean; whether or not to plot the treatment schedule.
    :param plotDrugAsBar: Boolean, whether to plot drug as bar across the top, or as shading underneath plot.
    :param drugBarPosition: Position of the drug bar when plotted across the top.
    :param drugColorMap: Color map for colouring the shading when using different drugs.
    :param lw_events: Line width for vertical event lines.
    :param xlim: x-axis limit.
    :param ylim: y-axis limit.
    :param y2lim: y2-axis limit.
    :param markInitialSize: Boolean, whether or not to draw horizontal line at height of fist data point.
    :param plotHorizontalLine: Boolean, whether or not to draw horizontal line at position specified at lineYPos.
    :param lineYPos: y-position at which to plot horizontal line.
    :param despine: Boolean, whether or not to despine the plot.
    :param titleStr: Title to put on the figure.
    :param decorateX: Boolean, whether or not to add labels and ticks to x-axis.
    :param decorateY: Boolean, whether or not to add labels and ticks to y-axis.
    :param decorateY2: Boolean, whether or not to add labels and ticks to y2-axis.
    :param markersize: Size of markers for feature variable.
    :param linestyle: Feature variable line style.
    :param linecolor: Feature variable line color.
    :param ax: matplotlib axis to plot on. If none provided creates a new figure.
    :param figsize: Tuple, figure dimensions when creating new figure.
    :param outName: Name under which to save figure.
    :param kwargs: Other kwargs to pass to plotting functions.
    :return:
    '''
    if ax is None: fig, ax = plt.subplots(1, 1, figsize=figsize)
    # Plot the data
    sns.lineplot(x=timeColumn, y=feature, hue=hue, style=style,
                 color=linecolor, legend=legend, palette=palette,
                 marker='o', markersize=markersize, markeredgewidth=2,
                 ax=ax, data=dataDf)

    # Plot the drug concentration
    if plotDrug:
        ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
        if plotDrugAsBar:
            if treatmentScheduleList is None:
                tVec = dataDf[timeColumn]
                drugConcentrationVec = utils.TreatmentListToTS(
                    treatmentList=utils.ExtractTreatmentFromDf(dataDf, timeColumn=timeColumn,
                                                               treatmentColumn=treatmentColumn),
                    tVec=tVec)
            else:
                tVec = dataDf[timeColumn].unique()
                drugConcentrationVec = utils.TreatmentListToTS(treatmentScheduleList,
                    tVec=tVec)
            drugConcentrationVec[drugConcentrationVec < 0] = 0
            drugConcentrationVec = np.array([x / (np.max(drugConcentrationVec) + 1e-12) for x in drugConcentrationVec])
            drugConcentrationVec = drugConcentrationVec / (1 - drugBarPosition) + drugBarPosition
            ax2.fill_between(tVec, drugBarPosition, drugConcentrationVec,
                             step="post", color="black", alpha=1., label="Drug Concentration")
            ax2.axis("off")
        else:
            currDrugBarPosition = drugBarPosition
            drugBarHeight = (1-drugBarPosition)/len(drugColorMap.keys())
            for drug in drugColorMap.keys():
                drugConcentrationVec = utils.TreatmentListToTS(
                    treatmentList=utils.ExtractTreatmentFromDf(dataDf, timeColumn=timeColumn,
                                                               treatmentColumn="%s Dose (mg)"%drug),
                    tVec=dataDf[timeColumn])
                drugConcentrationVec[drugConcentrationVec < 0] = 0
                # Normalise drug concentration to 0-1 (1=max dose(=initial dose))
                drugConcentrationVec = np.array([x / (np.max(drugConcentrationVec) + 1e-12) for x in drugConcentrationVec])
                # Rescale to make it fit within the bar at the top of the plot
                drugConcentrationVec = drugConcentrationVec * drugBarHeight + currDrugBarPosition
                ax2.fill_between(dataDf[timeColumn], currDrugBarPosition, drugConcentrationVec, step="post",
                                 color=drugColorMap[drug], alpha=0.5, label="Drug Concentration")
                ax2.hlines(xmin=dataDf[timeColumn].min(), xmax=dataDf[timeColumn].max(), 
                          y=currDrugBarPosition, linewidth=3, color="black")
                currDrugBarPosition += drugBarHeight
            # Line at the top of the drug bars
            ax2.hlines(xmin=dataDf[timeColumn].min(), xmax=dataDf[timeColumn].max(), 
                          y=currDrugBarPosition, linewidth=3, color="black")
        # Format y2 axis
        if y2lim is not None: ax2.set_ylim([0, y2lim])
        ax2.tick_params(labelsize=28)
        if not decorateY2:
            ax2.set_yticklabels("")

    # Format the plot
    if xlim is not None: ax.set_xlim(0, xlim)
    if ylim is not None: ax.set_ylim(0, ylim)
    if despine: sns.despine(ax=ax, trim=True, offset=50)

    # Draw horizontal lines (e.g. initial size)
    if plotHorizontalLine or markInitialSize or markPositiveCutOff:
        xlim = ax.get_xlim()[1]
        if markInitialSize: lineYPos = dataDf.loc[dataDf[timeColumn] == 0, feature]
        if markPositiveCutOff: lineYPos = 0.5 # cut-off value for positive is 0.5 copies/uL
        ax.hlines(xmin=0, xmax=xlim, y=lineYPos, linestyles=':', linewidth=4)

    # Decorate the plot
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(titleStr)
    if not decorateX:
        ax.set_xticklabels("")
    if not decorateY:
        ax.set_yticklabels("")
    plt.tight_layout()
    if outName is not None: plt.savefig(outName)

# =========================================================================================================
# Functions to help with game assay simulations
# =========================================================================================================
def simulate_well(rFrac, wellId, seed=None, outDir=None, jarFileName="../code.jar", modelConfigDic={}, **args):
    '''
    Function to simulate an individual well. Outputs a results file with detailed longitudinal measurements.
    '''
    # Setup model object
    model = OnLatticeModel(jarFileName = jarFileName)
    modelConfigDic = modelConfigDic.copy() # make local copy so don't edit parent
    modelConfigDic['rFrac'] = rFrac
    if seed is not None: modelConfigDic['seed'] = seed
    if outDir is not None: modelConfigDic['outDir'] = os.path.join(outDir, wellId)
    if outDir is not None: modelConfigDic['imageOutDir'] = os.path.join(outDir, wellId)
    modelConfigDic = {**modelConfigDic, **args}
    model.SetParams(**modelConfigDic)
    # Run
    model.RunSimulation()
    return pd.read_csv(os.path.join(modelConfigDic['outDir'],"RepId_%d.csv"%modelConfigDic['seed']))