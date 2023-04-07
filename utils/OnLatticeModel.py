# ====================================================================================
# Class to simulate spheroid growth using one of 6 ODE models
# ====================================================================================
import numpy as np
import pandas as pd
import math
import os
import sys
if 'matplotlib' not in sys.modules:
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")
# ====================================================================================
class OnLatticeModel():
    def __init__(self, **kwargs):
        # Initialise parameters
        self.modelConfigDic = {
                               # -------------------- Experimental Setup --------------------
                               "xDim": 100, "yDim": 100, # Dimensions of the domain (in lattice sites)
                               "nReplicates": 10, # Number of replicates
                               "seed": -1, # Random number seed
                               "initialSize": 250, # Total initial cell number (absolute)
                               "initialSizeProp": 0, # Initial cell density relative to (physical) carrying capacity; takes precedence over the absolute number if it is specified.
                               "rFrac": 0.05, # Initial resistance fraction in [0,1]
                               "tEnd": 5, # End time in days
                               "dt": 1, # Time step in days
                               # TODO: Test whether the schedule gets passed correctly
                               "treatmentScheduleList": None, # Treatment schedule in format [[tStart, tEnd, drugConcentration]]
                               # -------------------- Cell Properties --------------------
                               "divisionRate_S": 0.027, "divisionRate_R": 0.027, # Proliferation rate of cells in d^-1.
                               "movementRate_S": 0.0, "movementRate_R": 0.0, # Movement rate of cells in d^-1.
                               "deathRate_S": 0.027, "deathRate_R": 0.027, # Natural death rate of cells in d^-1.
                               "drugKillProportion": 0.75, # Probability that dividing sensitive cell is killed by drug.
                               "turnover": None, "cost": None, # Shorthand way of specifying cost and turnover
                               # -------------------- Output - Text --------------------
                               "verboseLevel": 0, "printFrequency": 100, # Command line verbosity. This goes to stdout.
                               "logCellCountFrequency": 1, # Frequency (in time units) at which the cell counts are written to file. <0 indicates no logging.
                               "outDir": "./tmp/", # Directory which to save output files to.
                               # -------------------- Output - Visualisation --------------------
                               # TODO: Clean up how visualisation are managed (and headless vs visualiseB stuff)
                            #    "headless": True, # Whether to show graphics window. Want this off to run on cluster etc
                            #    "visualiseB": False, # Whether to generate visualisation. Turn off for speed increase.
                               "scaleFactor": 2, "pause": 0, # Visualsation configs
                               "imageOutDir": "./tmp/", "imageFrequency": -1, # Saving of simulation images. Negative freq turns it off.") 
                               # -------------------- Output - Model File --------------------
                               # TODO: Get this implemented
                               "saveModelState": False, # Whether or not to save the model object at the end of the simulation.
                               "savedModelFileName": None # Name of model file to load when continuing a previous run.
                               }
        self.jarFileName = kwargs.get('jarFileName', './onLatticeModel.jar')
        self.resultsDf = None

        # Set the parameters
        self.SetParams(**kwargs)

        # Configure the solver
        self.suppressOutputB = kwargs.get('suppressOutputB',
                                          False)  # If true, suppress output of ODE solver (including warning messages)

    # =========================================================================================
    # Function to set the parameters
    def SetParams(self, **kwargs):
        for key in self.modelConfigDic.keys():
            self.modelConfigDic[key] = kwargs.get(key, self.modelConfigDic[key])
        # The 'cost' and 'turnover' parameters provide shorthands to set the global deathrate and resistant cell proliferation rate
        if self.modelConfigDic['turnover'] != None:
            self.modelConfigDic['deathRate_S'] = self.modelConfigDic['turnover']*self.modelConfigDic['divisionRate_S']
            self.modelConfigDic['deathRate_R'] = self.modelConfigDic['turnover']*self.modelConfigDic['divisionRate_S']
        if self.modelConfigDic['cost'] != None:
            self.modelConfigDic['divisionRate_R'] = (1-self.modelConfigDic['cost'])*self.modelConfigDic['divisionRate_S']
        # Ensure that the paths for the output directories have a file separator at the end so that
        # java then creates the correct file names.
        self.modelConfigDic['outDir'] = os.path.join(self.modelConfigDic['outDir'], "")
        if 'imageOutDir' in self.modelConfigDic: self.modelConfigDic['imageOutDir'] = os.path.join(self.modelConfigDic['imageOutDir'], "")

    def ConvertTreatmentScheduleToStr(self, treatmentScheduleList):
        treatmentScheduleStr = "["
        for interval in treatmentScheduleList:
            treatmentScheduleStr += "[%1.2f,%1.2f,%1.2f]," % tuple(interval)
        treatmentScheduleStr = treatmentScheduleStr[:-1]
        treatmentScheduleStr += "]"
        return treatmentScheduleStr

    def RunSimulation(self, printCommand=False):
        '''
        Parses arguments into command line command and run this command in the terminal
        '''
        argStr = " "
        for var in self.modelConfigDic.keys():
            if self.modelConfigDic[var] is None: continue # For vars set to None use default vals from java side
            if var in ["cost", "turnover"]: continue # Certain parameters are only used in the wrapper and not passed to java
            if isinstance(self.modelConfigDic[var], bool): 
                if self.modelConfigDic[var]: argStr += "--%s " % (var)
            else:
                argStr += "--%s %s " % (var, self.modelConfigDic[var])
        if printCommand: print("java -jar %s" % self.jarFileName + argStr)
        os.system("java -jar %s" % self.jarFileName + argStr)

    def LoadSimulations(self, normalise=True):
        tmpList = []
        replicateIdList = range(self.modelConfigDic['nReplicates']) if self.modelConfigDic['nReplicates']>1 else [self.modelConfigDic['seed']]
        for replicateId in replicateIdList:
            currDfName = os.path.join(self.modelConfigDic['outDir'], "RepId_%d.csv" % (replicateId))
            tmpDf = pd.read_csv(currDfName)
            tmpDf = tmpDf[['Time', 'NCells_S', 'NCells_R', 'NCells', 'DrugConcentration']]
            tmpDf['ReplicateId'] = replicateId
            tmpList.append(tmpDf)
        resultsDf = pd.concat(tmpList)
        resultsDf.rename(columns={"NCells": "TumourSize",
                                  "NCells_S": "S",
                                  "NCells_R": "R"}, inplace=True)
        if normalise:
            resultsDf['TumourSize'] /= 1e4
            resultsDf['S'] /= 1e4
            resultsDf['R'] /= 1e4
        return resultsDf

    def NormaliseToInitialSize(self, dataDf):
        dataDf['S'] /= dataDf.TumourSize.iloc[0]
        dataDf['R'] /= dataDf.TumourSize.iloc[0]
        dataDf['TumourSize'] /= dataDf.TumourSize.iloc[0]

    # =========================================================================================
    # Function to simulate the model
    def Simulate(self, treatmentScheduleList=None, scaleTumourVolume=False, printCommand=False, **kwargs):
        # Allow configuring the solver at this point as well
        self.jarFileName = kwargs.get('jarFileName', self.jarFileName)
        if treatmentScheduleList is not None: self.modelConfigDic["treatmentScheduleList"] = self.ConvertTreatmentScheduleToStr(treatmentScheduleList)
        # if (self.modelConfigDic['imageFreq']==-1) and ('imageOutDir' in self.modelConfigDic): self.modelConfigDic.pop('imageOutDir')
        self.SetParams() # Calling this to make sure parameters and file paths are properly checked in case they were reset at any point

        # Run the simulations
        self.RunSimulation(printCommand=printCommand)
        # self.SetParams(fromScratch=False) # TODO: get continuation to work again

        # Load data
        self.resultsDf = self.LoadSimulations()
        self.resultsDf = self.resultsDf.groupby(by="Time").mean()
        self.resultsDf.reset_index(inplace=True)
        self.resultsDf.drop(columns="ReplicateId", inplace=True)
        if scaleTumourVolume: self.NormaliseToInitialSize(self.resultsDf)

    # =========================================================================================
    # Function to plot the model predictions
    def Plot(self, scaleTumourVolume=True, aggregateData=True, progressBar=False, drugBarPosition=0.85,
             xlim=None, ylim=1.3, y2lim=1, decorateX=True, decorateY=True, axisLabels=False, markersize=10,
             labelsize=28,
             titleStr="", ax=None, figsize=(10, 8), outName=None, **kwargs):
        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)

        # Load data
        predictionDf = self.LoadSimulations()
        if scaleTumourVolume: self.NormaliseToInitialSize(predictionDf)

        # Plot the size the we will see on the images
        sns.lineplot(x="Time", y="TumourSize", ci='sd',
                     lw=kwargs.get('linewidthA', kwargs.get('linewidth', 7)),
                     color=kwargs.get('colorA', '#094486'),
                     estimator='mean' if aggregateData else None, legend=False,
                     data=predictionDf, ax=ax)

        # Plot the individual populations
        sns.lineplot(x="Time", y="S", ci='sd',
                     lw=kwargs.get('linewidth', 7), color=kwargs.get('colorS', "#0F4C13"),
                     estimator='mean' if aggregateData else None, legend=False,
                     data=predictionDf, ax=ax)
        ax.lines[1].set_linestyle(kwargs.get('linestyleS', '--'))
        sns.lineplot(x="Time", y="R", ci='sd',
                     lw=kwargs.get('linewidth', 7), color=kwargs.get('colorR', '#710303'),
                     estimator='mean' if aggregateData else None, legend=False,
                     data=predictionDf, ax=ax)
        ax.lines[2].set_linestyle(kwargs.get('linestyleR', '-.'))

        # Plot the drug concentration
        ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
        exampleReplicateId = predictionDf.ReplicateId.unique()[0]
        timeVec = predictionDf.Time[predictionDf.ReplicateId == exampleReplicateId]
        drugConcentrationVec = predictionDf.DrugConcentration[predictionDf.ReplicateId == exampleReplicateId]
        drugConcentrationVec = drugConcentrationVec / (1 - drugBarPosition) + drugBarPosition
        ax2.fill_between(timeVec,
                         drugBarPosition, drugConcentrationVec, color="black",
                         alpha=1., label="Drug Concentration")
        ax2.axis("off")

        # Format the plot
        if xlim is not None: ax.set_xlim(0, xlim)
        ax.set_ylim(0, ylim)
        ax2.set_ylim([0, y2lim])
        ax.set_xlabel("Time in Days" if axisLabels else "", fontdict={'fontsize': 28})
        ax.set_ylabel("PSA (Normalised)" if axisLabels else "", fontdict={'fontsize': 28})
        ax.set_title(titleStr)
        ax.tick_params(labelsize=labelsize)
        ax2.tick_params(labelsize=labelsize)
        if not decorateX:
            ax.set_xticklabels("")
        if not decorateY:
            ax.set_yticklabels("")
        plt.tight_layout()
        if outName is not None: plt.savefig(outName)

    # =========================================================================================
    # # Function to plot the model predictions
    # def Make_Gif(self, treatmentScheduleList=None, seed=1, imageFreq=10, printCommand=False,
    #          titleStr="", ax=None, figsize=(10, 8), outName=None, **kwargs):
    #     if ax is None:
    #         fig = plt.figure(figsize=figsize)
    #         ax = fig.add_subplot(111)
    #
    #     # Create images for Gif
    #     configDic_gifMaker = self.modelConfigDic.copy()
    #     if treatmentScheduleList is None:
    #         if self.modelConfigDic["treatmentScheduleList"] is None:
    #             raise ValueError("No treatment schedule defined!")
    #     else:
    #         configDic_gifMaker["treatmentScheduleList"] = self.ConvertTreatmentScheduleToStr(treatmentScheduleList)
    #     configDic_gifMaker['seed'] = seed
    #     configDic_gifMaker['nReplicates'] = 1
    #     configDic_gifMaker['imageFreq'] = imageFreq
    #     RunSimulation(configDic_gifMaker, jarFileName=self.jarFileName, printCommand=printCommand)
    #
    #     # Format the plot
    #     ax.set_title(titleStr)
    #     plt.tight_layout()
    #     if outName is not None: plt.savefig(outName)