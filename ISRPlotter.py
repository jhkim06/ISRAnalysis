import gc
import time
import os
import json
import pandas as pd
import numpy as np
import ROOT as rt
from root_numpy import hist2array, array2hist
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, LogLocator, NullFormatter, LogFormatter)
                               
from matplotlib.pyplot import cm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.ticker as mtick

from scipy.odr import *
from scipy import stats

import matplotlib as mpl
mpl.rcParams['hatch.linewidth'] = 0.3
from matplotlib.lines import Line2D

from scipy import optimize
import mplhep as hep

def root_to_numpy(hist) :
            
    values = []
    bins = []
    error = []  
                
    nbinsx = hist.GetNbinsX()
    for ibin in range(nbinsx):
        values.append(hist.GetBinContent(ibin+1))
        error.append(hist.GetBinError(ibin+1))
        bins.append(hist.GetXaxis().GetBinLowEdge(ibin+1))
                    
    bins.append(bins[-1]+hist.GetBinWidth(nbinsx))
                    
    values = np.array(values)
    bins = np.array(bins)
    error = np.array(error)
                    
    return (values, bins), error

def divide(nominator, denominator, output_type='root'):# 'root', 'numpy'

    ratio = nominator.Clone("ratio")
    ratio.Divide(denominator)

    if output_type == 'numpy' :
        # (values, bins, error)
        return root_to_numpy(ratio)

    else :
        # return root histogram
        return ratio

dashes = "--------------------------------------------------------------------------------"
# root TH1 to DataFrame and plot using matplotlib
class ISRPlotter :
    
    def __init__(self, inputHistFilePath, jasonConfigFilePath, doSystematic=False, verbose=True, setQuantile=False) :

        print(dashes)
        print("CREATE ISRPlotter......")

        self.setQuantile = setQuantile

        self.labels_list = []
        self.plots_list = []

        # Set using self.binDef
        self.binDef={} # dictionary of TUnfoldBinning object
        self.massBins=[] # list of tuple, ex) [(40.,64.), (64., 81.), (81., 101.), (101., 200.), (200., 320.)]

        self.histTypes=["Data", "DataBkgMCSubtracted", "SigMC", "BkgMC", "MCTotal", "Histogram"]

        # Dictionary of raw histograms according DataFrames
        self.rawHistsDict = {}

        for histType in self.histTypes :
            self.rawHistsDict[histType] = dict()
    
        # Dictionary of DataFrames containing nominal and systematic histogram content
        self.dfs = {} # self.dfs["Data"] self.dfs["DataBkgMCSubtracted"] self.dfs["SigMC"] self.dfs["BkgMC"] self.dfs["MCTotal"]

        for histType in self.histTypes :
            self.dfs[histType] = dict()

        self.normalisation = 1.
        self.bkgUsed = False

        self.nSystematics = 0 # Number of systematics
        self.systematic_marker={
            "IsoSF":"o", "IdSF":"o", "recoSF":"o", "trgSF":"o", "trgDZSF":"o", "PU":"o", "bveto":"v", "L1Prefire":"o", "Unfold": "^", "LepMomScale": "<", "LepMomRes": ">", 
            "Scale": "*", "AlphaS": "+", "PDF": "X", "UnfoldingMatrix": "o", "fsr" : "o"
        }

        self.lumi={"2016": "35.9/fb", "2017": "41.5/fb", "2018": "59.7/fb"}

        # FIXME make a function to read json configuration file
        # Read json file to get information of histograms
        with open(jasonConfigFilePath, 'r') as jsonFile :

            self.config = json.load(jsonFile)

            self.plotPrefix=self.config['plotPrefix'] # prefix for output plot file
            self.analysis=self.config['Analysis']
            self.year=self.config['Year']
            self.channel=self.config["Channel"]
            self.topDirName=self.config['topDirName']
            self.histPrefix=self.config["HistPrefix"]
            self.useTUnfoldBin=self.config["useTUnfoldBin"]
            self.tunfoldBinNames=self.config["TUnfoldBinNames"]
            self.variables=self.config["Variables"]
            self.variablePostfix=self.config["VariablePostfix"]
            self.steeringTUnfold=self.config['Steering']

            if doSystematic : 
                self.systematics=self.config['Systematics']

                if "UnfoldingInput" in self.systematics["measurement"] :
                    n_max=int(((self.systematics["measurement"]["UnfoldingInput"][1]).split("_"))[1])
                    prefix=((self.systematics["measurement"]["UnfoldingInput"][1]).split("_"))[0]

                    self.systematics["measurement"]["UnfoldingInput"].pop()
                    self.systematics["measurement"]["UnfoldingInput"]=[prefix+"_{}".format(str(i)) for i in range(0,n_max+1)]

                if "UnfoldingMatrix" in self.systematics["measurement"] :
                    n_max=int(((self.systematics["measurement"]["UnfoldingMatrix"][1]).split("_"))[1])
                    prefix=((self.systematics["measurement"]["UnfoldingMatrix"][1]).split("_"))[0]

                    self.systematics["measurement"]["UnfoldingMatrix"].pop()
                    self.systematics["measurement"]["UnfoldingMatrix"]=[prefix+"_{}".format(str(i)) for i in range(0,n_max+1)]

                if "PDF" in self.systematics["theory"] :
                    n_max=int(((self.systematics["theory"]["PDF"][1]).split("_"))[1])
                    pdf_prefix=((self.systematics["theory"]["PDF"][1]).split("_"))[0]

                    self.systematics["theory"]["PDF"].pop()
                    self.systematics["theory"]["PDF"]=[pdf_prefix+"{:0>3}".format(str(i)) for i in range(1,n_max+1)]
            else :
                self.systematics=dict()

            self.samples=self.config['Samples']
            self.stackOrder=self.config['StackOrder']

            # Out directory
            self.outDirPath="output/"+self.year+"/"+self.channel+"/"
            if "ChannelPostfix" in self.config :
                self.channelPostfix = self.config['ChannelPostfix']
                self.outDirPath="output/"+self.year+"/"+self.channel+"_"+self.channelPostfix+"/"
            if not os.path.exists(self.outDirPath) :
                os.makedirs(self.outDirPath)

        if verbose==True :
            print('This is {} {} data of {} analysis...'.format(self.year, self.channel, self.analysis))
            print("Systematics saved in the input root file...")
            print(self.systematics)
            print("Samples saved in the input root file...")
            print(self.samples)
            if self.useTUnfoldBin :
                print("TUnfoldBinning is used")

        # Open input histogram root file
        self.inRootFile=rt.TFile.Open(inputHistFilePath, 'READ')

        if self.useTUnfoldBin :
            
            for var in ["2D_dimass_dipt", "2D_dipt_dimass"] :
                self.binDef[var]=self.inRootFile.Get(self.topDirName+"/"+var+"/truth") 

            # todo if the bin has two axis then use the second axis to draw each distribution of the first axis
            temp_tvecd=self.binDef["2D_dipt_dimass"].GetDistributionBinning(1)
            temp_mass_bin_edges=temp_tvecd.GetMatrixArray()
            self.massBins.extend([ (temp_mass_bin_edges[index], temp_mass_bin_edges[index+1]) for index in range(temp_tvecd.GetNrows()-1)]) # list comprehension

        else :
            if self.channel == "electron" :
                self.massBins.extend([(50, 64), (64, 81), (81, 101), (101, 200), (200, 320), (320,830)])
            else :
                self.massBins.extend([(40, 64), (64, 81), (81, 101), (101, 200), (200, 320), (320,830)])

        print("Read histograms from input file...")
        start_time = time.time()

        count_nSystematics = True
        
        # Get raw TH1 histograms
        # Convert them into DataFrame
        for variable in self.variables : # Mass 
              
            # dict[variable]
            varDir=variable
            if self.useTUnfoldBin :
                varDir=variable.split("__")[0] # In case, TUnfoldBinning used
                if len(self.variablePostfix)>0 :

                    if varDir == "Mass" :
                        varDir=varDir+"_"+self.variablePostfix[1]
                    if varDir == "Pt" :
                        varDir=varDir+"_"+self.variablePostfix[0]

            # dict[variable][sample]
            for combinedName in self.samples :
                # initialize 
                if variable not in self.rawHistsDict[combinedName.split("_")[1]] : # rawHistDict["Data"]["Pt__0"]
                    self.rawHistsDict[combinedName.split("_")[1]][variable] = dict() 
                    temp_dict=self.rawHistsDict[combinedName.split("_")[1]][variable]
                    if "BkgMC" in combinedName : 
                        self.bkgUsed = True
                else :
                    temp_dict=self.rawHistsDict[combinedName.split("_")[1]][variable]
               
                if variable not in self.rawHistsDict["MCTotal"] :
                    self.rawHistsDict["MCTotal"][variable]=dict()

                isFirstFileinCombinedName=True # first sample in the combinedName
                isLastFileinCombinedName=False

                for fileIndex, mcFileName in enumerate(self.samples[combinedName]) :

                    if fileIndex == len(self.samples[combinedName])-1 :
                        isLastFileinCombinedName=True

                    temp_dict[mcFileName]=dict()

                    # Set nominal histograms
                    sysName = "Nominal"
                    postfix = "Nominal"

                    temp_dict[mcFileName][sysName]=dict()
                    temp_dict[mcFileName][sysName][postfix]=dict()

                    self.setRawHistsDict(varDir, combinedName, mcFileName, sysName, postfix, variable, temp_dict, 
                    isFirstFileinCombinedName, isLastFileinCombinedName)

                    # Set systemaitc histograms
                    for sysCategory in self.systematics.keys() :
                        for sysName, postfixs in self.systematics[sysCategory].items() :
                            # just for varifying
                            if count_nSystematics :
                                self.nSystematics += 1

                            temp_dict[mcFileName][sysName]=dict()
                            if isFirstFileinCombinedName : 
                                temp_dict[combinedName][sysName]=dict()
                            for postfix in postfixs :
                                temp_dict[mcFileName][sysName][postfix]=dict()
                                if isFirstFileinCombinedName : 
                                    temp_dict[combinedName][sysName][postfix]=dict()

                                self.setRawHistsDict(varDir, combinedName, mcFileName, sysName, postfix, variable, temp_dict, 
                                isFirstFileinCombinedName, isLastFileinCombinedName)

                    count_nSystematics = False
                    isFirstFileinCombinedName=False

        gc.collect() 
        end_time = time.time()
        print("DONE")
        print(dashes)
        print("time elapsed: {:.2f} s".format(end_time-start_time))
        print(dashes)
    
        print("Combine histograms.....")
        start_time = time.time()

        self.setTotalHists(self.rawHistsDict["Data"])
        self.setTotalHists(self.rawHistsDict["SigMC"])
        if self.bkgUsed : 
            self.setTotalHists(self.rawHistsDict["BkgMC"])
            self.setBkgSubtractedDataHis()
        self.setMCTotalHists() # signal + background MC
        if self.rawHistsDict["Histogram"] != 0 : 
            self.setTotalHists(self.rawHistsDict["Histogram"])    

        end_time = time.time()
        print("DONE")
        print(dashes)
        print("time elapsed: {:.2f} s".format(end_time-start_time))
        print(dashes)

        print("Convert histogram to DataFrames.....")
        start_time = time.time()

        # DataFrame
        for histType in self.histTypes :
        
            if "DataBkgMCSubtracted" == histType or "BkgMC" == histType : 
                if self.bkgUsed == False : continue
            
            self.setDataFrameWithUnc(histType)

        end_time = time.time()
        print("DONE")
        print(dashes)
        print("time elapsed: {:.2f} s".format(end_time-start_time))
        print(dashes)

        print("Calculate up/down systematics.....")
        start_time = time.time()

        for histType in self.histTypes :
        
            if "DataBkgMCSubtracted" == histType or "BkgMC" == histType : 
                if self.bkgUsed == False : continue

            self.calculateCombinedUnc(histType, "total")
            self.calculateCombinedUnc(histType, "theory")
            self.calculateCombinedUnc(histType, "measurement")

            if self.useTUnfoldBin :
                
                self.calculateCombinedUnc(histType, "total", "upDownUnc_meanValue")
                self.calculateCombinedUnc(histType, "theory", "upDownUnc_meanValue")
                self.calculateCombinedUnc(histType, "measurement", "upDownUnc_meanValue")

                if self.setQuantile :
                    self.calculateCombinedUnc(histType, "total", "upDownUnc_qValue", True)
                    self.calculateCombinedUnc(histType, "theory", "upDownUnc_qValue", True)
                    self.calculateCombinedUnc(histType, "measurement", "upDownUnc_qValue", True)

        end_time = time.time()
        print("DONE")
        print(dashes)
        print("time elapsed: {:.2f} s".format(end_time-start_time))
        print(dashes)

    def setRawHistsDict(self, varDir, combinedName, inputHistFileName, sysName, sysPostfix, variable, targetDict, isFirstFile, isLastFile) :

        if sysName == "Nominal" and sysPostfix == "Nominal" : 
            histToRead = self.topDirName + "/"+varDir + "/" + self.histPrefix + inputHistFileName 
        else :
            histToRead = self.topDirName + "/"+varDir + "/" + self.histPrefix + inputHistFileName + '_' + sysName + sysPostfix 

        # read histogram
        temp_TH1=self.inRootFile.Get(histToRead)
        if type(temp_TH1) != rt.TH1D :
            histToRead = self.topDirName + "/"+varDir + "/" + self.histPrefix + inputHistFileName 
            temp_TH1=self.inRootFile.Get(histToRead)
        
        if self.useTUnfoldBin :
            temp_TH1 = self.binDef[varDir.split("__")[0]].ExtractHistogram(inputHistFileName + variable + sysName + sysPostfix, temp_TH1, 
                                                                            0, True, self.steeringTUnfold[variable]) 

        targetDict[inputHistFileName][sysName][sysPostfix]["TH1"] = temp_TH1
        targetDict[inputHistFileName][sysName][sysPostfix]["DataFrame"] = self.convertTH1toDataFrame(temp_TH1)

        del temp_TH1

    def getOutBaseDir(self) :
        return self.outDirPath

    def loglinear_func(self, p, x):
        return 2.*p[0]*np.log(x)+p[1]

    def setBkgSubtractedDataHis(self) :

        for variable in self.variables :
            self.rawHistsDict["DataBkgMCSubtracted"][variable]=dict()
            self.rawHistsDict["DataBkgMCSubtracted"][variable]["total"]=dict()

            sysName = "Nominal" 
            postfix = "Nominal"
    
            self.rawHistsDict["DataBkgMCSubtracted"][variable]["total"][sysName]=dict()
            self.rawHistsDict["DataBkgMCSubtracted"][variable]["total"][sysName][postfix]=dict() 

            temp_dict=self.rawHistsDict["DataBkgMCSubtracted"][variable]["total"][sysName][postfix]

            data_dict=self.rawHistsDict["Data"][variable]["total"][sysName][postfix]
            mc_bkg_dict=self.rawHistsDict["BkgMC"][variable]["total"][sysName][postfix]

            temp_dict["TH1"]=data_dict["TH1"].Clone("data_bkg_subtracted")
            temp_dict["TH1"].Add(mc_bkg_dict["TH1"], -1)

            temp_dict["DataFrame"]=data_dict["DataFrame"].copy()
            temp_dict["DataFrame"].content=data_dict["DataFrame"].content-mc_bkg_dict["DataFrame"].content
            

            for sysCategory in self.systematics.keys() :
                for sysName, postfixs in self.systematics[sysCategory].items() :
                    self.rawHistsDict["DataBkgMCSubtracted"][variable]["total"][sysName]=dict()
                    for postfix in postfixs :
                        self.rawHistsDict["DataBkgMCSubtracted"][variable]["total"][sysName][postfix]=dict()
                        temp_dict=self.rawHistsDict["DataBkgMCSubtracted"][variable]["total"][sysName][postfix]

                        data_dict=self.rawHistsDict["Data"][variable]["total"][sysName][postfix]
                        mc_bkg_dict=self.rawHistsDict["BkgMC"][variable]["total"][sysName][postfix]

                        temp_dict["TH1"]=data_dict["TH1"].Clone("data_bkg_subtracted")
                        temp_dict["TH1"].Add(mc_bkg_dict["TH1"], -1)

                        temp_dict["DataFrame"]=data_dict["DataFrame"].copy()
                        temp_dict["DataFrame"].content=data_dict["DataFrame"].content-mc_bkg_dict["DataFrame"].content

    def setMCTotalHists(self) :

        temp_dict = self.rawHistsDict["MCTotal"]
        for variable in self.variables :
            # dict[variable]
            #temp_dict[variable]=dict()
            temp_dict[variable]["total"]=dict()
            # 
            sysName = "Nominal"
            postfix = "Nominal"

            temp_dict[variable]["total"][sysName]=dict()
            temp_dict[variable]["total"][sysName][postfix]=dict()

            # Lets combine MC histograms
            temp_dict[variable]["total"][sysName][postfix]["TH1"] = \
            self.rawHistsDict["SigMC"][variable]["total"][sysName][postfix]["TH1"].Clone("Clone_"+variable+sysName+postfix)
            temp_dict[variable]["total"][sysName][postfix]["DataFrame"] = \
            self.rawHistsDict["SigMC"][variable]["total"][sysName][postfix]["DataFrame"].copy()

            if self.bkgUsed :
                temp_dict[variable]["total"][sysName][postfix]["TH1"].Add(self.rawHistsDict["BkgMC"][variable]["total"][sysName][postfix]["TH1"])
                temp_dict[variable]["total"][sysName][postfix]["DataFrame"].content= \
                self.rawHistsDict["BkgMC"][variable]["total"][sysName][postfix]["DataFrame"].content+temp_dict[variable]["total"][sysName][postfix]["DataFrame"].content

            for sysCategory in self.systematics.keys() :
                for sysName, postfixs in self.systematics[sysCategory].items() :
                    temp_dict[variable]["total"][sysName]=dict()
                    for postfix in postfixs :
                        temp_dict[variable]["total"][sysName][postfix]=dict()

                        # Lets combine MC histograms
                        temp_dict[variable]["total"][sysName][postfix]["TH1"] = \
                        self.rawHistsDict["SigMC"][variable]["total"][sysName][postfix]["TH1"].Clone("Clone_"+variable+sysName+postfix)
                        temp_dict[variable]["total"][sysName][postfix]["DataFrame"] = \
                        self.rawHistsDict["SigMC"][variable]["total"][sysName][postfix]["DataFrame"].copy()

                        if self.bkgUsed :
                            temp_dict[variable]["total"][sysName][postfix]["TH1"].Add(self.rawHistsDict["BkgMC"][variable]["total"][sysName][postfix]["TH1"])
                            temp_dict[variable]["total"][sysName][postfix]["DataFrame"].content= \
                            self.rawHistsDict["BkgMC"][variable]["total"][sysName][postfix]["DataFrame"].content+temp_dict[variable]["total"][sysName][postfix]["DataFrame"].content

    def setDataFrameWithUnc(self, dictName="Data") :

        in_dict=self.rawHistsDict[dictName]
        out_dict=self.dfs[dictName]

        for variable in in_dict.keys() :
            out_dict[variable]=dict()
            print(variable)

            for sample in in_dict[variable].keys() : # Note loop over samples defined in the "dictionary!"
                
                out_dict[variable][sample]=dict()
                out_dict[variable][sample]["rawUnc"]=in_dict[variable][sample]["Nominal"]["Nominal"]["DataFrame"].copy()
                out_dict[variable][sample]["upDownUnc"]=in_dict[variable][sample]["Nominal"]["Nominal"]["DataFrame"].copy()
                # TODO Create column for statistical uncertainty

                # DataFrame for mean values
                # low mass cut, high mass cut, mean value
                # For Mass, make a DataFrame for all the mass bins
                # For pT, make a DataFrame for a mass bin

                q = 0.7
                if self.useTUnfoldBin :
                    out_dict[variable][sample]["rawUnc_meanValue"]    = self.createMeanDataFrame(variable, in_dict[variable][sample]["Nominal"]["Nominal"]["TH1"])
                    out_dict[variable][sample]["upDownUnc_meanValue"] = self.createMeanDataFrame(variable, in_dict[variable][sample]["Nominal"]["Nominal"]["TH1"])

                    # quantile
                    if "Pt" in variable and self.setQuantile: 
                        #print("call createQuantileDataFrame")
                        out_dict[variable][sample]["rawUnc_qValue"]    = self.createQuantileDataFrame(variable, in_dict[variable][sample]["Nominal"]["Nominal"]["TH1"], q, True)
                        out_dict[variable][sample]["upDownUnc_qValue"] = self.createQuantileDataFrame(variable, in_dict[variable][sample]["Nominal"]["Nominal"]["TH1"], q, True)

                for sysCategory in self.systematics.keys() :
                    for sysName, postfixs in self.systematics[sysCategory].items() :

                        if sysName == "Nominal" : continue
                        #print(" variable {} sample {}".format(variable, sample))
                        #print(" systematic name {} number of variations {} ".format(sysName, len(postfixs)))
                        #print(" Set differences....")

                        # temp dictionary
                        temp_rawUnc_dict = {}
                        temp_rawUnc_meanValue_dict = {}
                        if self.setQuantile :
                            temp_rawUnc_qValue_dict = {}

                        for postfix in postfixs :

                            if "fsr" not in sysName :

                                temp_rawUnc_dict[sysName+"_"+postfix] = \
                                in_dict[variable][sample][sysName][postfix]["DataFrame"].content-in_dict[variable][sample]["Nominal"]["Nominal"]["DataFrame"].content

                                if self.useTUnfoldBin :
                                    temp_df=self.createMeanDataFrame(variable, in_dict[variable][sample][sysName][postfix]["TH1"])

                                    temp_rawUnc_meanValue_dict[sysName+"_"+postfix] = \
                                    temp_df["mean"]-out_dict[variable][sample]["rawUnc_meanValue"]["mean"]

                                    if "Pt" in variable and self.setQuantile:
                                        temp_quantile_df=self.createQuantileDataFrame(variable, in_dict[variable][sample][sysName][postfix]["TH1"], q)

                                        temp_rawUnc_qValue_dict[sysName+"_"+postfix] = \
                                        temp_quantile_df["q_"+str(q)]-out_dict[variable][sample]["rawUnc_qValue"]["q_"+str(q)]

                            else :
                                current_index=postfixs.index(postfix)
                                the_other_postfix=None

                                if current_index == 0 : the_other_postfix=postfixs[1]
                                else : the_other_postfix=postfixs[0]

                                #out_dict[variable][sample]["rawUnc"][sysName+"_"+postfix]= \
                                #in_dict[variable][sample][sysName][postfix]["DataFrame"].content-in_dict[variable][sample][sysName][the_other_postfix]["DataFrame"].content

                                temp_rawUnc_dict[sysName+"_"+postfix] = \
                                in_dict[variable][sample][sysName][postfix]["DataFrame"].content-in_dict[variable][sample][sysName][the_other_postfix]["DataFrame"].content

                                if self.useTUnfoldBin :
                                    temp_df1=self.createMeanDataFrame(variable, in_dict[variable][sample][sysName][postfix]["TH1"])
                                    temp_df2=self.createMeanDataFrame(variable, in_dict[variable][sample][sysName][the_other_postfix]["TH1"])

                                    temp_rawUnc_meanValue_dict[sysName+"_"+postfix] = \
                                    temp_df1["mean"]-temp_df2["mean"]

                                    if "Pt" in variable and self.setQuantile:
                                        temp_quantile_df1=self.createQuantileDataFrame(variable, in_dict[variable][sample][sysName][postfix]["TH1"], q)
                                        temp_quantile_df2=self.createQuantileDataFrame(variable, in_dict[variable][sample][sysName][the_other_postfix]["TH1"], q)

                                        temp_rawUnc_qValue_dict[sysName+"_"+postfix] = \
                                        temp_quantile_df1["q_"+str(q)]-temp_quantile_df2["q_"+str(q)]

                        temp_rawUnc_df = pd.DataFrame(temp_rawUnc_dict)
                        temp_rawUnc_meanValue_df = pd.DataFrame(temp_rawUnc_meanValue_dict)
                        if self.setQuantile :
                            temp_rawUnc_qValue_df = pd.DataFrame(temp_rawUnc_qValue_dict)

                        out_dict[variable][sample]["rawUnc"] = pd.concat([out_dict[variable][sample]["rawUnc"], temp_rawUnc_df], axis = 1)
                        out_dict[variable][sample]["rawUnc_meanValue"] = pd.concat([out_dict[variable][sample]["rawUnc_meanValue"], temp_rawUnc_meanValue_df], axis = 1)

                        if "Pt" in variable and self.setQuantile:
                            out_dict[variable][sample]["rawUnc_qValue"] = pd.concat([out_dict[variable][sample]["rawUnc_qValue"], temp_rawUnc_qValue_df], axis = 1)

                        #print(" Set up/down....")

                        # For systematics taking r.m.s of variations as systematic up/down
                        if sysName == "PDF" or sysName == "UnfoldingInput" or sysName == "UnfoldingMatrix":

                            out_dict[variable][sample]["upDownUnc"][sysName+'_Up']  =np.sqrt(out_dict[variable][sample]["rawUnc"].filter(like=sysName).var(axis=1))
                            out_dict[variable][sample]["upDownUnc"][sysName+'_Down']= -1. * np.sqrt(out_dict[variable][sample]["rawUnc"].filter(like=sysName).var(axis=1))

                            out_dict[variable][sample]["upDownUnc_meanValue"][sysName+'_Up']  =\
                            np.sqrt(out_dict[variable][sample]["rawUnc_meanValue"].filter(like=sysName).var(axis=1))
                            out_dict[variable][sample]["upDownUnc_meanValue"][sysName+'_Down']=\
                            -1. * np.sqrt(out_dict[variable][sample]["rawUnc_meanValue"].filter(like=sysName).var(axis=1))

                            if "Pt" in variable and self.setQuantile:
                                out_dict[variable][sample]["upDownUnc_qValue"][sysName+'_Up']  =\
                                np.sqrt(out_dict[variable][sample]["rawUnc_qValue"].filter(like=sysName).var(axis=1))
                                out_dict[variable][sample]["upDownUnc_qValue"][sysName+'_Down']=\
                                -1. * np.sqrt(out_dict[variable][sample]["rawUnc_qValue"].filter(like=sysName).var(axis=1))

                        else :
                            # Note sysName MUST BE UNIQUE
                            regex_ = "(^" + sysName + "_)"
                            temp_max = out_dict[variable][sample]["rawUnc"].filter(regex=regex_).max(axis=1)
                            temp_min = out_dict[variable][sample]["rawUnc"].filter(regex=regex_).min(axis=1)

                            if temp_max.equals(temp_min) :
                                temp_max = temp_max.abs()
                                temp_min = -1. * temp_min.abs()

                            temp_max=temp_max.fillna(0) 
                            temp_min=temp_min.fillna(0)

                            temp_symmetric_max = (temp_max.abs() + temp_min.abs()) / 2.
                            temp_symmetric_min = (temp_max.abs() + temp_min.abs()) / -2.

                            out_dict[variable][sample]["upDownUnc"][sysName+'_Up'] = temp_symmetric_max
                            out_dict[variable][sample]["upDownUnc"][sysName+'_Down'] = temp_symmetric_min

                            if self.useTUnfoldBin : 
                                # for mean value
                                temp_mean_max = out_dict[variable][sample]["rawUnc_meanValue"].filter(regex=regex_).max(axis=1)
                                temp_mean_min = out_dict[variable][sample]["rawUnc_meanValue"].filter(regex=regex_).min(axis=1) 

                                if temp_mean_max.equals(temp_mean_min) :
                                    temp_mean_max = temp_mean_max.abs()
                                    temp_mean_min = -1. * temp_mean_min.abs() 

                                temp_mean_max = temp_mean_max.fillna(0)
                                temp_mean_min = temp_mean_min.fillna(0)

                                temp_mean_symmetric_max = (temp_mean_max.abs() + temp_mean_min.abs()) / 2.
                                temp_mean_symmetric_min = (temp_mean_max.abs() + temp_mean_min.abs()) / -2.

                                out_dict[variable][sample]["upDownUnc_meanValue"][sysName+'_Up'] = temp_mean_symmetric_max
                                out_dict[variable][sample]["upDownUnc_meanValue"][sysName+'_Down'] = temp_mean_symmetric_min

                                # for quantile
                                if "Pt" in variable and self.setQuantile:
                                    temp_quantile_max = out_dict[variable][sample]["rawUnc_qValue"].filter(regex=regex_).max(axis=1)
                                    temp_quantile_min = out_dict[variable][sample]["rawUnc_qValue"].filter(regex=regex_).min(axis=1) 

                                    if temp_quantile_max.equals(temp_quantile_min) :
                                        temp_quantile_max = temp_quantile_max.abs()
                                        temp_quantile_min = -1. * temp_quantile_min.abs() 

                                    temp_quantile_max=temp_quantile_max.fillna(0)
                                    temp_quantile_min=temp_quantile_min.fillna(0)

                                    temp_quantile_symmetric_max = (temp_quantile_max.abs() + temp_quantile_min.abs()) / 2.
                                    temp_quantile_symmetric_min = (temp_quantile_max.abs() + temp_quantile_min.abs()) / -2.

                                    out_dict[variable][sample]["upDownUnc_qValue"][sysName+'_Up'] = temp_quantile_symmetric_max
                                    out_dict[variable][sample]["upDownUnc_qValue"][sysName+'_Down'] = temp_quantile_symmetric_min


    def createMeanDataFrame(self, variable, TH1_hist) :

        if "dipt_dimass" in variable and self.useTUnfoldBin:
            num_str      = variable.split("__")[1]
            nth_mass_bin = int(num_str)
            temp_mean    = TH1_hist.GetMean()
            temp_dict    = {"low mass cut": self.massBins[nth_mass_bin][0], "high mass cut": self.massBins[nth_mass_bin][1], "mean": temp_mean, "stat_error": TH1_hist.GetMeanError()}
            temp_df      = pd.DataFrame([temp_dict], columns=['low mass cut','high mass cut','mean', 'stat_error'])

            return temp_df

        else :
            if "dimass_dipt" in variable :
                # Get mean mass for all mass bins
                #self.massBins
                # Create DataFrame
                row_list=[]

                for massBin in self.massBins :
                    TH1_hist.GetXaxis().SetRangeUser(massBin[0], massBin[1])
                    temp_mean = TH1_hist.GetMean()
                    temp_dict = {"low mass cut": massBin[0], "high mass cut": massBin[1], "mean": temp_mean, "stat_error": TH1_hist.GetMeanError()}
                    row_list.append(temp_dict)

                temp_df = pd.DataFrame(row_list, columns=['low mass cut','high mass cut','mean', 'stat_error'])

                return temp_df

            else :

                row_list=[]
                temp_mean=TH1_hist.GetMean()
                #print("low mass cut: {} high mass cut: {} mean : {}".format(massBin[0], massBin[1], temp_mean))
                temp_dict={"low mass cut": variable.split("_")[-2].split("mll")[1], "high mass cut": variable.split("_")[-1], "mean": temp_mean, "stat_error": TH1_hist.GetMeanError()}
                row_list.append(temp_dict)
                temp_df = pd.DataFrame(row_list, columns=['low mass cut','high mass cut','mean', 'stat_error'])

                return temp_df

    def createQuantileDataFrame(self, variable, TH1_hist, prob, set_stat = False) :

        #if "Pt" in variable :
        if not variable.isalpha() and self.useTUnfoldBin:
            num_str = variable.split("__")[1]
            nth_mass_bin = int(num_str)

            temp_TH1_hist = TH1_hist.Clone("temp_TH1_hist") 
            temp_TH1_hist.Scale(1., "width") # divide by binwidth
        
            length = 1
            q = np.zeros(length)
            prob = np.array([prob])
            temp_TH1_hist.GetQuantiles(length, q, prob)

            # stat
            # 1000 toys
            if set_stat : 
                n_toys = 10
                toy_quantiles = []
                for i in range(n_toys) :
                    #print(str(i) + " th toy")

                    temp_toy = TH1_hist.Clone("temp_toy")
                    temp_toy.Reset("ICESM")
                    temp_toy.FillRandom(temp_TH1_hist, 100000)
                    temp_toy.Scale(1., "width")

                    temp_q = np.zeros(1)
                    temp_toy.GetQuantiles(length, temp_q, prob)
                    toy_quantiles.append(temp_q[0])
                
                stat_error = np.std(toy_quantiles)

            else :
                stat_error = 0.

            temp_dict={"low mass cut": self.massBins[nth_mass_bin][0], "high mass cut": self.massBins[nth_mass_bin][1], "q_" + str(prob[0]): q[0], "stat_error": stat_error}
            temp_df=pd.DataFrame([temp_dict], columns=['low mass cut','high mass cut','q_' + str(prob[0]), 'stat_error'])

            del temp_TH1_hist

            return temp_df

    def calculateCombinedUnc(self, dictName="Data", sys_to_combine="total", column_name="upDownUnc", is_quantile = False) :

        in_dict=self.dfs[dictName]

        # total uncertainty (stat+sys), total theory, total measurement
        for variable in in_dict.keys() :

            # for quantile, mass is not considered
            if is_quantile : 
                if "Pt" not in variable :
                    continue

            for sample in in_dict[variable].keys() :

                combinedSys_Up = None
                combinedSys_Down = None
                first_sys = True

                for index, sysCategory in enumerate(self.systematics.keys()) :

                    if sys_to_combine != "total" :
                        if sysCategory != sys_to_combine : continue

                    else :
                        if first_sys == True :

                            combinedSys_Up=np.square(in_dict[variable][sample][column_name]['stat_error']).copy()
                            combinedSys_Down=np.square(in_dict[variable][sample][column_name]['stat_error']).copy()
                            first_sys=False

                    for sysName, postfixs in self.systematics[sysCategory].items() :
                        if sysName == "Nominal" : continue # Nominal don't have Up, Down...
                        if first_sys :
                            combinedSys_Up=np.square(in_dict[variable][sample][column_name][sysName+'_Up']).copy()
                            combinedSys_Down=np.square(in_dict[variable][sample][column_name][sysName+'_Down']).copy()
                            first_sys=False
                        else :
                            combinedSys_Up=combinedSys_Up+np.square(in_dict[variable][sample][column_name][sysName+'_Up'])
                            combinedSys_Down=combinedSys_Down+np.square(in_dict[variable][sample][column_name][sysName+'_Down'])

                if combinedSys_Up is not None :
                    in_dict[variable][sample][column_name][sys_to_combine+"_Up"]=np.sqrt(combinedSys_Up)
                else :
                    in_dict[variable][sample][column_name][sys_to_combine+"_Up"]=0.

                if combinedSys_Down is not None :
                    in_dict[variable][sample][column_name][sys_to_combine+"_Down"]=np.sqrt(combinedSys_Down)
                else :
                    in_dict[variable][sample][column_name][sys_to_combine+"_Down"]=0.

    def setTotalHists(self, mc_dict) :

        for variable in mc_dict.keys() :
            mc_dict[variable]["total"]=dict()

            sysName = "Nominal"
            postfix = "Nominal"

            mc_dict[variable]["total"][sysName]=dict()
            mc_dict[variable]["total"][sysName][postfix]=dict()

            first_mc=True 
            for sample in mc_dict[variable].keys() : 

                if sample == "total" : continue

                if first_mc :
                    first_mc=False
                    mc_dict[variable]["total"][sysName][postfix]["TH1"]=\
                    mc_dict[variable][sample][sysName][postfix]["TH1"].Clone("Clone_"+variable+sample+sysName+postfix)

                    mc_dict[variable]["total"][sysName][postfix]["DataFrame"]=\
                    mc_dict[variable][sample][sysName][postfix]["DataFrame"].copy()
                else :
                    mc_dict[variable]["total"][sysName][postfix]["TH1"].Add(mc_dict[variable][sample][sysName][postfix]["TH1"])

                    mc_dict[variable]["total"][sysName][postfix]["DataFrame"].content= \
                    mc_dict[variable]["total"][sysName][postfix]["DataFrame"].content+mc_dict[variable][sample][sysName][postfix]["DataFrame"].content


            for sysCategory in self.systematics.keys() :

                for sysName, postfixs in self.systematics[sysCategory].items() :

                    mc_dict[variable]["total"][sysName]=dict()

                    for postfix in postfixs :
                        mc_dict[variable]["total"][sysName][postfix]=dict()

                        # Lets combine MC histograms
                        first_mc=True

                        for sample in mc_dict[variable].keys() :

                            if sample in self.samples.keys(): continue # self.samples.keys() is a list of combined sample
                            if sample == "total" : continue

                            if first_mc :
                                first_mc=False
                                mc_dict[variable]["total"][sysName][postfix]["TH1"]=\
                                mc_dict[variable][sample][sysName][postfix]["TH1"].Clone("Clone_"+variable+sample+sysName+postfix)

                                mc_dict[variable]["total"][sysName][postfix]["DataFrame"]=\
                                mc_dict[variable][sample][sysName][postfix]["DataFrame"].copy()
                            else :
                                mc_dict[variable]["total"][sysName][postfix]["TH1"].Add(mc_dict[variable][sample][sysName][postfix]["TH1"])

                                mc_dict[variable]["total"][sysName][postfix]["DataFrame"].content= \
                                mc_dict[variable]["total"][sysName][postfix]["DataFrame"].content+mc_dict[variable][sample][sysName][postfix]["DataFrame"].content


    def convertTH1toDataFrame(self, temp_TH1, existing_df = None) :

        temp_content, temp_binEdge=hist2array(temp_TH1, return_edges=True)
        binWidth = (temp_binEdge[0][1:] - temp_binEdge[0][:-1])
        nBin=len(temp_binEdge[0])-1

        temp_stat_error = []
        for ibin in range(1, len(temp_content)+1) :
            temp_stat_error.append(temp_TH1.GetBinError(ibin))

        temp_stat_error_array = np.array(temp_stat_error)

        pd_series_binIndex    = pd.Series(range(1,nBin+1), range(1, nBin+1), name="bin_index")
        pd_series_binWidth    = pd.Series(binWidth, range(1, nBin+1), name="bin_width")
        pd_series_lowBinEdge  = pd.Series(temp_binEdge[0][0:-1], range(1, nBin+1), name="low_bin_edge")
        pd_series_highBinEdge = pd.Series(temp_binEdge[0][1:], range(1, nBin+1), name="high_bin_edge")

        dict_temp={
            'bin_width':     pd_series_binWidth,
            'low_bin_edge':   pd_series_lowBinEdge,
            'high_bin_edge':  pd_series_highBinEdge,
            'content':       temp_content,
            'stat_error':    temp_stat_error_array,
        }

        pd_temp=pd.DataFrame(dict_temp, index=pd_series_binIndex)

        if existing_df is None :
            return pd_temp
        else :
            existing_df.update(pd_temp)
            del pd_temp
            return None

    def setXaxisLabel(self, variable) :

        channelName = "e^{+}e^{-}"
        if self.channel == "muon" :
            channelName = "\mu^{+}\mu^{-}"

        if "Lep" in variable :
            channelName = "e"
            if self.channel == "muon" :
                channelName = "\mu"

        varName = "$p_{T}^{\mathit{"+channelName+"}}$"
        unit = "[GeV]"

        if variable == "Mass" :
            varName = "Mass$^{\mathit{"+channelName+"}}$"
        if "Eta" in variable :
            varName = "$\eta^{\mathit{"+channelName+"}}$"
            unit = ""
        if "MET" in variable :
            varName = "MET"
        if "nPV" in variable :
            varName = "Number of Vertex"
            unit = ""

        return varName, unit

    def rebinDataFrame(self, source_df, target_df) :

        output_source=pd.DataFrame(index=target_df.index, columns=target_df.columns)
        target_nbins = len(target_df.index)
        source_nbins = len(source_df.index)

        for target_index in range(target_nbins) :
            low = target_df.iloc[target_index].low_bin_edge
            high = target_df.iloc[target_index].high_bin_edge
            
            low_source_index = -1
            high_source_index = -1
            # find index matched to the target index
            for source_index in range(source_nbins) :
                current_source_low=source_df.iloc[source_index].low_bin_edge
                current_source_high=source_df.iloc[source_index].high_bin_edge
                if low > current_source_low : continue
                if low == current_source_low :
                    low_source_index = source_index
                
                    if high == current_source_high :
                        high_source_index = source_index
                        break
                if high == current_source_high :
                    high_source_index = source_index
                    
            #print("target low {} high {} index {}, source index {} {}".format(low, high, target_index+1, low_source_index+1, high_source_index+1))
            #print(source_df.iloc[low_source_index:high_source_index+1]["content"].sum(), np.sqrt(np.square(source_df.iloc[low_source_index:high_source_index+1]["stat_error"]).sum()))
            for column_name in target_df.columns :
                if "bin" in column_name :
                    output_source.iloc[target_index][column_name] = target_df.iloc[target_index][column_name]
                elif "content" == column_name :
                    output_source.iloc[target_index][column_name] = source_df.iloc[low_source_index:high_source_index+1]["content"].sum()
                else :
                    output_source.iloc[target_index][column_name] = np.sqrt(np.square(source_df.iloc[low_source_index:high_source_index+1][column_name]).sum())

        output_source=output_source.astype('float')

        return output_source

    # Draw data and all MC and the ratio plot between data and total MC
    def drawSubPlot(self, axis, variable, divde_by_bin_width = False, setLogy=False,
                    write_xaxis_title=False, write_yaxis_title=False, setLogx = False, showLegend = False, show_additional_legend=False,
                    min=1e-1, max=1e9,
                    ratio_max = 1.35, ratio_min = 0.65, optimzeXrange=False, minimum_content=1, show_ratio=False, ext_objects=None, ext_names=None, 
                    denominator="total_Data", internal_names=None,
                    setRatioLogy = False, showNEvents=False, showChi2=False, ratioName=None, normNominator=False,
                    oneColumn=False, setMarker="o", setColor='red', xFactor=1., colors=None) :

        if len(axis) == 2 :
            top_axis = axis[0]
            bottom_axis = axis[1]
        else :
            top_axis = axis[0] 

        denominator_hist_type = denominator.split("_")[1]
        denominator_hist_name = denominator.split("_")[0]
        
        if len(denominator.split("_")) > 2 :
            denominator_label_name = denominator.split("_")[2:]

        denominator_print_name = denominator.split("_")[1]
        denominator_df=None
        denominator_hist=None

        nominator_hist_type="MCTotal"
        nominator_hist_name="total"
        nominator_hist=None

        df_filter = self.dfs[denominator_hist_type][variable][denominator_hist_name]["upDownUnc"]['content'] > minimum_content

        # set denominator
        denominator_df   = self.dfs[denominator_hist_type][variable][denominator_hist_name]["upDownUnc"][df_filter].copy()
        denominator_hist = self.rawHistsDict[denominator_hist_type][variable][denominator_hist_name]["Nominal"]["Nominal"]["TH1"]

        # set nominators
        nominator_df      = self.dfs[nominator_hist_type][variable][nominator_hist_name]["upDownUnc"][df_filter].copy()
        nominator_hist    = self.rawHistsDict[nominator_hist_type][variable][nominator_hist_name]["Nominal"]["Nominal"]["TH1"]
            
        channelName = "e^{+}e^{-}"
        if self.channel == "muon" :
            channelName = "\mu^{+}\mu^{-}"

        if show_ratio : 
            ratio, ratio_error = divide(nominator_hist, denominator_hist, output_type='numpy')

        # bins
        bins = denominator_df.low_bin_edge.values
        last_edge = denominator_df.high_bin_edge.values[-1]
        bins = np.append(bins, last_edge)

        hep.cms.label("Preliminary", data=True, lumi=50, year=2016, ax=top_axis, fontsize=20)
        hep.histplot((denominator_df.content, bins), ax=top_axis, binwnorm=1., yerr=True, histtype='errorbar', label="test", color='black')
        hep.histplot((nominator_df.content, bins), ax=top_axis, binwnorm=1., yerr=True, histtype='step', label="test", color='red')
        top_axis.set_xlim(bins[0], bins[-1])
        
        print(ratio)
        hep.histplot((ratio[0], bins), ax=bottom_axis, yerr=ratio_error, histtype='step', label="ratio", color='red')
        bottom_axis.set_xlim(bins[0], bins[-1])
        bottom_axis.set_ylim(ratio_min, ratio_max)

        gc.collect()

    def drawHistPlot(self, *variables, divde_by_bin_width = False, setLogy=False, setLogx=False,
                     min = 1e-2, max = 1e9, ratio_max=1.35, ratio_min=0.65,
                     optimzeXrange=False, minimum_content=1,
                     figSize=(10,6), show_ratio=True, ext_object_list=None, ext_names_list=None,
                     denominator="Data_Data", internal_names=None, setRatioLogy=False,
                     showNEvents=False, showChi2=False, outPdfPostfix=None, ratioName=None,
                     normNominator=False, oneColumn=False, colors=None, top=0.9, bottom=0.1) :

        hep.style.use("CMS")
        # FIXME make a function for setting
        if show_ratio == True :
            num_rows = 2
            
            if oneColumn :
        
                heights = []
                for i in range(len(variables) + 1) :
                    if i == 0 :
                        heights.append(1.)
                    else :
                        heights.append(0.3)
                fig, axes = plt.subplots(len(variables) + 1, 1, sharex=False, figsize=figSize, gridspec_kw={'height_ratios': heights})
                
            else :
            
                fig, axes = plt.subplots(num_rows, len(variables), sharex=False, figsize=figSize, gridspec_kw={'height_ratios':[1, 0.3]})

        else :
            num_rows = 1
            if oneColumn : 
                fig, axes = plt.subplots(num_rows, 1, sharex=False, figsize=figSize)
            else :
                fig, axes = plt.subplots(num_rows, len(variables), sharex=False, figsize=figSize)

        plt.tight_layout()

        write_xaxis_title = True
        write_yaxis_title = True

        show_legend = False

        color=['red', 'blue', 'green', 'cyan', 'magenta', 'yellow', 'black', 'white']
        factor = 1.
        write_xaxis_title=False
        show_additional_legend=False
        
        for index, variable in enumerate(variables) :

            if len(variables) == 1:
                show_legend = True
                write_xaxis_title=True
                show_additional_legend=True

                if show_ratio : 
                    axes_tuple = (axes[0], axes[1])

                else :
                    axes_tuple = (axes, )

            # draw multiple histograms(variables) at once
            else :
                if index > 0 :
                    write_yaxis_title=False

                if index == len(variables)-1 :
                    write_xaxis_title=True
                    show_additional_legend=True

                if oneColumn : 
                    show_legend = True
                    factor *= 10.
                else :
                    if len(variables) == index + 1:
                        show_legend = True

                if show_ratio : 
                    if oneColumn :
                        axes_tuple = (axes[0], axes[index+1]) 
                    else :
                        axes_tuple = (axes[0][index], axes[1][index])

                else :
                    if oneColumn :
                        axes_tuple = (axes,) 
                    else :
                        axes_tuple = (axes[index], )

            
            self.drawSubPlot(axes_tuple, variable=variable, divde_by_bin_width=divde_by_bin_width, setLogy=setLogy,
                             write_xaxis_title=write_xaxis_title, write_yaxis_title=write_yaxis_title, setLogx=setLogx, 
                             showLegend=show_legend, show_additional_legend=show_additional_legend, min = min, max = max, ratio_max=ratio_max, ratio_min=ratio_min, optimzeXrange=optimzeXrange,
                             minimum_content=minimum_content, show_ratio=show_ratio, ext_objects=ext_object_list, ext_names=ext_names_list, 
                             denominator=denominator, internal_names=internal_names, setRatioLogy=setRatioLogy, 
                             showNEvents=showNEvents, showChi2=showChi2, ratioName=ratioName, normNominator=normNominator, oneColumn=oneColumn, setMarker=Line2D.filled_markers[index], setColor=color[index%len(color)],
                             xFactor=factor, colors=colors)

        outPdfName = self.outDirPath+self.plotPrefix+"_"+variable+"_"+self.channel+"_"+self.year+".pdf" 
        if outPdfPostfix is not None :
            outPdfName = self.outDirPath+self.plotPrefix+"_"+variable+"_"+outPdfPostfix+".pdf" 
        
        plt.savefig(outPdfName, format="pdf", dpi=300)
        plt.close(fig)

        self.labels_list.clear()
        self.plots_list.clear()

    def combinedPtDataFrame(self, name="Data") :

        in_df=self.dfs[name]

        combined_pt_df = None
        for nth_mass_bin in range(len(self.massBins)) :
            if nth_mass_bin == 0 :
                combined_pt_df=in_df["2D_dipt_dimass__"+str(nth_mass_bin)]["total"]["upDownUnc_meanValue"].copy()
            else :
                #combined_pt_df=combined_pt_df.append(in_df["Pt__"+str(nth_mass_bin)]["total"]["upDownUnc_meanValue"].copy(), ignore_index=True)
                combined_pt_df=pd.concat([combined_pt_df, in_df["2D_dipt_dimass__"+str(nth_mass_bin)]["total"]["upDownUnc_meanValue"].copy()], ignore_index=True)
        return combined_pt_df

    def combinedQuantilePtDataFrame(self, name="Data") :

        in_df=self.dfs[name]

        combined_pt_df = None
        for nth_mass_bin in range(len(self.massBins)) :
            if nth_mass_bin == 0 :
                combined_pt_df=in_df["Pt__"+str(nth_mass_bin)]["total"]["upDownUnc_qValue"].copy()
            else :
                combined_pt_df=combined_pt_df.append(in_df["Pt__"+str(nth_mass_bin)]["total"]["upDownUnc_qValue"].copy(), ignore_index=True)
        return combined_pt_df


    def combinedMassPtDataFrame(self, name="Data") :

        in_df=self.dfs[name]

        combined_pt_df = None
        combined_mass_df = None

        for key in in_df.keys() :

            if "mll_mll" in key :
                if combined_mass_df is None :
                    combined_mass_df = in_df[key]["total"]["upDownUnc_meanValue"].copy() 
                else :
                    combined_mass_df = combined_mass_df.append(in_df[key]["total"]["upDownUnc_meanValue"].copy(), ignore_index=True)
            else : 
                if combined_pt_df is None :
                    combined_pt_df = in_df[key]["total"]["upDownUnc_meanValue"].copy() 
                else :
                    combined_pt_df = combined_pt_df.append(in_df[key]["total"]["upDownUnc_meanValue"].copy(), ignore_index=True)

        return combined_mass_df, combined_pt_df


    def doLogLinearFit(self, ax, mass_df, pt_df, mass_unc, pt_unc, line_color, useDF=True, printPar=True) :

        print("do log-linear fit!")
        x_err=[]
        y_err=[]
        for i in range(len(mass_unc.T)) :
            x_err.append(max(mass_unc.T[i]))
            y_err.append(max(pt_unc.T[i]))

        loglinear=Model(self.loglinear_func)
        if useDF :
            data=RealData(mass_df["mean"], pt_df["mean"], sx=x_err, sy=y_err)
        else :
            data=RealData(mass_df, pt_df, sx=x_err, sy=y_err)
        odr=ODR(data, loglinear, beta0=[1.0, 0.0])
        out=odr.run()
        out.pprint()

        xn = np.linspace(40., 700, 1000)
        yn = self.loglinear_func(out.beta, xn)
        ax.plot(xn, yn, color=line_color, linewidth=0.5)
        # prepare parameters for confidence interval curves
        nstd = 1. # to draw 1-sigma intervals
        popt_up = out.beta + nstd * out.sd_beta
        popt_dw = out.beta - nstd * out.sd_beta
        print(out.sd_beta)
        print(type(out.sd_beta))

        # calculate y values for 1 sigma
        fit_up = self.loglinear_func(popt_up, xn)
        fit_dw = self.loglinear_func(popt_dw, xn)

        if printPar :
            ax.text(0.05, .05, "$<p_{T}^{DY}>$" + "=({:.2f}$\pm${:.2f})+({:.2f}$\mp${:.2f})".format(out.beta[1], out.sd_beta[1], out.beta[0], out.sd_beta[0]) + "x $log <m_{DY}>^{2}$"
, fontsize='xx-large', transform=ax.transAxes)

        # plot 1 sigma
        #ax.plot(xn, fit_up, '#0076D4', dashes=[9, 4.5], label='1 Sigma uncertainty', linewidth=0.8)
        #ax.plot(xn, fit_dw, '#0076D4', dashes=[9, 4.5], linewidth=0.8)

    def drawISRUncertaintyPlot(self, variable = "Pt", ymin=0.95, ymax=1.05) :

        print("draw isr uncertainty plot")
        print(self.nSystematics)

        fig, ax = plt.subplots(figsize=(8, 6))
        plt.subplots_adjust(left=0.12, right=0.97, bottom=0.15, top=0.9)
        ax.text(0., 1.05, "CMS Work in progress", fontsize='xx-large', transform=ax.transAxes)
        ax.text(1., 1.05, "(13 TeV, " + self.year + ")", fontsize=20, transform=ax.transAxes, ha='right')

        ax.tick_params(bottom=True, left=True, right=True, which="both", direction='in')  

        mass_bin_index = [index for index in range(1, len(self.massBins)+1)]
        ax.set_xlim(0,len(self.massBins)+1)
        ax.set_ylim(ymin, ymax)
        ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))

        channelName = "e^{+}e^{-}"
        if self.channel == "muon" :
            channelName = "\mu^{+}\mu^{-}"

        ax.set_xlabel("Mass bins", fontsize=20, ha='right', x=1.0)
        temp_df = None
        if variable == "Pt" :
            ax.set_ylabel("Variation/ Nominal of Mean $p_{T}^{" + channelName + "}$", fontsize=20, ha='right', y=1.0)
            temp_df=self.combinedPtDataFrame("Data")
        if variable == "Mass" :
            ax.set_ylabel("Variable/ Nominal of Mean $M^{" + channelName + "}$", fontsize=20, ha='right', y=1.0)
            temp_df=self.Data["Mass"]["total"]["upDownUnc_meanValue"]

        ax.errorbar(mass_bin_index, (temp_df.total_Up+temp_df["mean"])/temp_df["mean"], fmt='--', ms = 4., color="black", label="Total uncertainty", linewidth=0.7)
        ax.errorbar(mass_bin_index, (temp_df["mean"]-temp_df.total_Down)/temp_df["mean"], fmt='--', ms = 4., color="black", linewidth=0.7)

        ax.errorbar(mass_bin_index, (temp_df.stat_error+temp_df["mean"])/temp_df["mean"], fmt='--', ms = 4., color="black", label="Stat. uncertainty", linewidth=0.5)
        ax.errorbar(mass_bin_index, (temp_df["mean"]-temp_df.stat_error)/temp_df["mean"], fmt='--', ms = 4., color="black", linewidth=0.5)

        legend_handle = []
        label_handle = []

        color=iter(cm.rainbow(np.linspace(0,1,self.nSystematics))) 
        for sysCategory in self.systematics.keys() :
            for sysName in self.systematics[sysCategory] :

                color_=next(color)
                
                temp_handle = ax.errorbar(mass_bin_index, (temp_df[sysName+"_Up"].abs()+temp_df["mean"])/temp_df["mean"], fmt=self.systematic_marker[sysName], ms = 4., color=color_)
                legend_handle.append(temp_handle)
                label_handle.append(sysName + " Up")

        color=iter(cm.rainbow(np.linspace(0,1,self.nSystematics))) 
        for sysCategory in self.systematics.keys() :
            for sysName in self.systematics[sysCategory] :

                color_=next(color)
                
                temp_handle = ax.errorbar(mass_bin_index, (temp_df["mean"]-temp_df[sysName+"_Down"].abs())/temp_df["mean"], fmt=self.systematic_marker[sysName], fillstyle = "none", ms = 4., color=color_)

                legend_handle.append(temp_handle)
                label_handle.append(sysName + " Down")

        ax.legend(tuple(legend_handle), tuple(label_handle), loc='upper left', fancybox=False, framealpha=0.0, ncol=2, fontsize=5)

        plt.tight_layout()
        plt.savefig(self.outDirPath+self.plotPrefix+"_ISRUncertainty.pdf", format="pdf", dpi=300)
        plt.close(fig)

    def drawISRPlot(self, *list_to_plot, do_linear_fit=False, ymin=13, ymax=30, xmin=30, xmax=4e2, outPdfPostfix=None) :
        print("draw isr plot, do_linear_fit {}".format(do_linear_fit))
        color=iter(cm.rainbow(np.linspace(0,1,len(list_to_plot))))
        isData=False

        fig, ax = plt.subplots(figsize=(8, 8))
        plt.subplots_adjust(left=0.12, right=0.97, bottom=0.15, top=0.9)
        ax.text(0., 1.05, "CMS Work in progress", fontsize='xx-large', transform=ax.transAxes)
        ax.text(1., 1.05, "(13 TeV, " + self.year + ")", fontsize=20, transform=ax.transAxes, ha='right')

        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.set_xscale("log")
        ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))

        channelName = "e^{+}e^{-}"
        if self.channel == "muon" :
            channelName = "\mu^{+}\mu^{-}"

        ax.set_xlabel("Mean $M^{" + channelName + "}$ [GeV]", fontsize=20, ha='right', x=1.0)
        ax.set_ylabel("Mean $p_{T}^{" + channelName + "}$ [GeV]", fontsize=20, ha='right', y=1.0)

        for index, name in enumerate(list_to_plot) :
            if "Data" in name :
                isData=True

            temp_mass_df=self.dfs[name]["2D_dimass_dipt"]["total"]["upDownUnc_meanValue"]
            temp_pt_df=self.combinedPtDataFrame(name)

            temp_mass_total_up =   np.sqrt(np.square(temp_mass_df["total_Up"]) + np.square(temp_mass_df["stat_error"]))
            temp_mass_total_down = np.sqrt(np.square(temp_mass_df["total_Down"]) + np.square(temp_mass_df["stat_error"]))

            temp_pt_total_up =   np.sqrt(np.square(temp_pt_df["total_Up"]) + np.square(temp_pt_df["stat_error"]))
            temp_pt_total_down = np.sqrt(np.square(temp_pt_df["total_Down"]) + np.square(temp_pt_df["stat_error"]))

            mass_systematic=self.makeErrorNumpy(temp_mass_total_up, temp_mass_total_down)
            pt_systematic=self.makeErrorNumpy(temp_pt_total_up, temp_pt_total_down)

            color_=next(color)

            if "Data" in name :
                color_ = 'black'

            ax.errorbar(temp_mass_df["mean"], temp_pt_df["mean"], xerr=mass_systematic, yerr=pt_systematic, fmt='o', ms = 4., color=color_, label=name)

            if do_linear_fit :
                if isData :
                    self.doLogLinearFit(ax, temp_mass_df, temp_pt_df, mass_systematic, pt_systematic, color_)
                    print("mass_systematic: ", mass_systematic)
                    print("pt_systematic: ", pt_systematic)
            isData=False

        ax.grid(True, which='both', axis='x', color='black', linewidth=0.3, linestyle="--")
        ax.grid(True, axis='y', color='black', linewidth=0.3, linestyle="--")

        ax.tick_params(bottom=True, top=True, left=True, right=True, which='both', direction='in')
        ax.tick_params(length=10, which='major')
        ax.tick_params(length=5, which='minor')
        ax.yaxis.set_minor_locator(AutoMinorLocator())

        ax.legend(loc='best', fontsize=15, fancybox=False, framealpha=0.0)

        outPdfName = self.outDirPath+self.plotPrefix+"_ISR.pdf"
        if outPdfPostfix is not None :
            outPdfName = self.outDirPath+self.plotPrefix+"_"+outPdfPostfix+"_ISR.pdf"
            
        plt.savefig(outPdfName, format="pdf", dpi=300)
        plt.close(fig)

    def drawISRPlots(self, *objects_to_plot, names_in_objects, do_linear_fit=None, labels=None, markers=None, colors=None, facecolor_list = None, 
                     ymin=13, ymax=30, xmin=30, xmax=4e2, outPdfPostfix=None, years = None, both_lepton = False, ratios = False, nominators = None, 
                     showRatio=False, showCDF=False, combined_result=None, showCMS=False) :

        allMeanMassDF = None
        allErrorMassList = None

        allMeanPtDF = None
        allErrorPtList = None

        print("draw isr plot, do_linear_fit {}".format(do_linear_fit))
        color=iter(cm.rainbow(np.linspace(0,1,len(names_in_objects))))
        isData=False

        if showRatio == False :
            fig, ax = plt.subplots(figsize=(8, 8))
        else :
            fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=False, gridspec_kw={'height_ratios':[1, 0.3]})
            ax = axes[0]
            bottom = axes[1] 

        plt.tight_layout()
        plt.subplots_adjust(left=0.12, right=0.97, bottom=0.15, top=0.9, hspace=0.05)
        ax.text(0., 1.05, "CMS Work in progress", fontsize=20, transform=ax.transAxes)
        if years is not None :
            ax.text(1., 1.05, "(13 TeV " + years + ")", fontsize=20, transform=ax.transAxes, ha='right')
        else :
            ax.text(1., 1.05, "(13 TeV, " + self.year + ")", fontsize=20, transform=ax.transAxes, ha='right')

        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.set_xscale("log")

        channelName = "e^{+}e^{-}"
        if self.channel == "muon" :
            channelName = "\mu^{+}\mu^{-}"
        if showCMS:
            channelName="ll";

        ax.set_xlabel("Mean $M^{" + channelName + "}$ [GeV]", fontsize=20, ha='right', x=1.0)
        #ax.set_ylabel("Mean $p_{T}^{" + channelName + "}$ [GeV]", fontsize=20, ha='right', y=1.0)
        ax.set_ylabel("$p_{T}^{" + channelName + "}$ [GeV]", fontsize=20, ha='right', y=1.0)

        if both_lepton :
            ax.set_xlabel("Mean $M^{ll}$ [GeV]", fontsize=20, ha='right', x=1.0)
            ax.set_ylabel("Mean $p_{T}^{ll}$ [GeV]", fontsize=20, ha='right', y=1.0)

        ax.grid(True, which='both', axis='x', color='black', linewidth=0.3, linestyle="--")
        ax.grid(True, axis='y', color='black', linewidth=0.3, linestyle="--")

        ax.tick_params(bottom=True, top=True, left=True, right=True, which='both', direction='in')
        ax.tick_params(length=10, which='major')
        ax.tick_params(length=5, which='minor')
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        if showRatio :
            ax.set_xticklabels([]) 
        else :
            ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))

        #if ratios :
        #    ax.set_ylabel("R", fontsize=20, ha='right', y=1.0) 
        if showCMS == False: 
            for index, name in enumerate(names_in_objects) :

                isData=False
                if name == "Data" :
                    isData=True

                label_name=name
                if labels is not None :
                    label_name=labels[index]

                if markers is None :
                    current_marker='o'
                else :
                    current_marker=markers[index]

                if colors is None :
                    current_color="black"
                else :
                    current_color=colors[index]

                if facecolor_list is not None :
                    current_facecolor=facecolor_list[index]

                if index==0 :
                    temp_mass_df =self.getDict(name)["2D_dimass_dipt"]["total"]["upDownUnc_meanValue"]
                    temp_pt_df   =self.combinedPtDataFrame(name)
                    if self.setQuantile : 
                        temp_quantile_pt_df = self.combinedQuantilePtDataFrame(name)

                else :
                    if objects_to_plot[index-1].useTUnfoldBin :
                        temp_mass_df=objects_to_plot[index-1].getDict(name)["2D_dimass_dipt"]["total"]["upDownUnc_meanValue"]
                        temp_pt_df=objects_to_plot[index-1].combinedPtDataFrame(name)
                        if self.setQuantile :
                            temp_quantile_pt_df=objects_to_plot[index-1].combinedQuantilePtDataFrame(name)
                    else :
                        temp_mass_df, temp_pt_df = objects_to_plot[index-1].combinedMassPtDataFrame(name)

                temp_mass_total_up =   np.sqrt(np.square(temp_mass_df["total_Up"]) + np.square(temp_mass_df["stat_error"]))
                temp_mass_total_down = np.sqrt(np.square(temp_mass_df["total_Down"]) + np.square(temp_mass_df["stat_error"]))

                temp_pt_total_up =   np.sqrt(np.square(temp_pt_df["total_Up"]) + np.square(temp_pt_df["stat_error"]))
                temp_pt_total_down = np.sqrt(np.square(temp_pt_df["total_Down"]) + np.square(temp_pt_df["stat_error"]))

                mass_systematic=self.makeErrorNumpy(temp_mass_total_up, temp_mass_total_down)
                pt_systematic=self.makeErrorNumpy(temp_pt_total_up, temp_pt_total_down)

                if self.setQuantile :

                    temp_quantile_pt_total_up =   np.sqrt(np.square(temp_quantile_pt_df["total_Up"]) + np.square(temp_quantile_pt_df["stat_error"]))
                    temp_quantile_pt_total_down = np.sqrt(np.square(temp_quantile_pt_df["total_Down"]) + np.square(temp_quantile_pt_df["stat_error"]))
                    pt_quantile_systematic=self.makeErrorNumpy(temp_quantile_pt_total_up, temp_quantile_pt_total_down)

                mass_systematic=self.makeErrorNumpy(temp_mass_total_up, temp_mass_total_down)
                pt_systematic=self.makeErrorNumpy(temp_pt_total_up, temp_pt_total_down)

                if index==0 or objects_to_plot[index-1].useTUnfoldBin :

                    if ratios == False :
                        ax.errorbar(temp_mass_df["mean"], temp_pt_df["mean"], xerr=mass_systematic, yerr=pt_systematic, fmt=current_marker, color=current_color, mfc="none", label=label_name, linewidth=0.5, ms = 4)
                        if self.setQuantile :
                            ax.errorbar(temp_mass_df["mean"], temp_quantile_pt_df["q_0.7"], xerr=mass_systematic, yerr=pt_systematic, fmt=current_marker, color=current_color, mfc="none", label="q=0.7", linewidth=0.5, ms = 4)

                        if index == 0 :
                            allMeanMassDF = temp_mass_df.copy(deep=True)
                            allErrorMassList = mass_systematic.copy()

                            allMeanPtDF = temp_pt_df.copy(deep=True)
                            allErrorPtList = pt_systematic.copy()
                        else :
                            #allMeanMassDF.append(temp_mass_df)
                            allMeanMassDF=pd.concat([allMeanMassDF, temp_mass_df])
                            np.append(allErrorMassList[0], mass_systematic[0])
                            np.append(allErrorMassList[1], mass_systematic[1])

                            #allMeanPtDF.append(temp_pt_df)
                            allMeanPtDF=pd.concat([allMeanPtDF, temp_pt_df])
                            np.append(allErrorPtList[0], pt_systematic[0])
                            np.append(allErrorPtList[1], pt_systematic[1])

                        if showRatio :
                            if index == 0 :
                                bottom.set_ylim(-5, 5)
                                bottom.set_xlim(xmin,xmax)
                                bottom.set_xlabel("Mean $M^{" + channelName + "}$ [GeV]", fontsize=20, ha='right', x=1.0)
                                bottom.set_ylabel("$\\delta \\barp_{T}$", fontsize=20, ha='center', y=0.5)
                                bottom.set_xscale("log")
                                bottom.tick_params(bottom=True, top=True, left=True, right=True, which='both', direction='in')
                                bottom.tick_params(length=10, which='major')
                                bottom.tick_params(length=5, which='minor')
                                bottom.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
                                bottom.axhline(0., color='black', linewidth=1, linestyle="--")

                                denominator_mass = temp_mass_df["mean"]
                                denominator_pt = temp_pt_df["mean"]
                                denominator_mass_err = mass_systematic
                                denominator_pt_err = pt_systematic

                            else :

                                print("draw bottom")
                                temp_ratio = temp_pt_df["mean"] - denominator_pt
                                bottom.errorbar(temp_mass_df["mean"], temp_ratio, xerr=mass_systematic, yerr=pt_systematic, fmt=current_marker, color=current_color, mfc=current_facecolor, label=label_name, linewidth=1., ms = 4)
                    else :
                        # comparison to ratio method
                        temp_pt_df_nominator = nominators[index].combinedPtDataFrame("SigMC")
                        temp_mass_df_nominator =  nominators[index].getDict("SigMC")["Mass"]["total"]["upDownUnc_meanValue"]

                        if index==0 :
                            temp_pt_df_denominator = self.combinedPtDataFrame("SigMC")
                            temp_mass_df_denominator = self.getDict("SigMC")["Mass"]["total"]["upDownUnc_meanValue"]
                        else :
                            temp_pt_df_denominator = objects_to_plot[index-1].combinedPtDataFrame("SigMC")
                            temp_mass_df_denominator =  objects_to_plot[index-1].getDict("SigMC")["Mass"]["total"]["upDownUnc_meanValue"] 

                        temp_pt_df_ratio = temp_pt_df_nominator["mean"]/ temp_pt_df_denominator["mean"]
                        temp_mass_df_ratio = temp_mass_df_nominator["mean"]/ temp_mass_df_denominator["mean"]
                        #ax.errorbar(temp_mass_df["mean"], temp_pt_df_ratio, fmt=current_marker, color=current_color, mfc=current_facecolor, label=label_name, linewidth=1., ms = 4)
                        ax.errorbar(temp_mass_df_ratio * temp_mass_df["mean"], temp_pt_df_ratio * temp_pt_df["mean"], 
                                    fmt=current_marker, color=current_color, mfc=current_facecolor, label=label_name, linewidth=1., ms = 4)

                        # results from unfolding method
                        temp_pt_df_unfolding = nominators[index].combinedPtDataFrame("Data")
                        temp_mass_df_unfolding = nominators[index].getDict("Data")["Mass"]["total"]["upDownUnc_meanValue"]

                        temp_mass_total_up =   np.sqrt(np.square(temp_mass_df_unfolding["total_Up"]) + np.square(temp_mass_df_unfolding["stat_error"]))
                        temp_mass_total_down = np.sqrt(np.square(temp_mass_df_unfolding["total_Down"]) + np.square(temp_mass_df_unfolding["stat_error"]))

                        temp_pt_total_up =   np.sqrt(np.square(temp_pt_df_unfolding["total_Up"]) + np.square(temp_pt_df_unfolding["stat_error"]))
                        temp_pt_total_down = np.sqrt(np.square(temp_pt_df_unfolding["total_Down"]) + np.square(temp_pt_df_unfolding["stat_error"]))

                        mass_systematic=self.makeErrorNumpy(temp_mass_total_up, temp_mass_total_down)
                        pt_systematic=self.makeErrorNumpy(temp_pt_total_up, temp_pt_total_down)

                        current_facecolor = current_color
                        ax.errorbar(temp_mass_df_unfolding["mean"], temp_pt_df_unfolding["mean"], xerr=mass_systematic, yerr=pt_systematic, fmt=current_marker, color=current_color, mfc=current_facecolor, label=label_name, linewidth=1., ms = 4)

                else :
    
                    ax.plot(temp_mass_df["mean"], temp_pt_df["mean"], color=current_facecolor, linewidth=0.5)
                    ax.fill_between(temp_mass_df["mean"], temp_pt_df["mean"]-temp_pt_total_up, temp_pt_df["mean"]+temp_pt_total_up, facecolor=current_facecolor, alpha=0.2)

                print(label_name)
                #print(temp_mass_df["mean"].append(temp_pt_df["mean"], ignore_index = True))
                d = {'mean mass': temp_mass_df["mean"].values, "mean mass stat": temp_mass_df["stat_error"], "mean mass sys": temp_mass_df["total_Up"], 
                     'mean pt': temp_pt_df["mean"].values, "mean pt stat": temp_pt_df["stat_error"], "mean pt sys": temp_pt_df["total_Up"]}
                print(pd.DataFrame(data=d))

                #if do_linear_fit[index] :
                #    if isData :
                #        self.doLogLinearFit(ax, allMeanMassDF, allMeanPtDF, allErrorMassList, allErrorPtList, "black")

        if showCDF :
            
            # Muon
            cdf_muon_mass      = [47.72, 70.66, 90.99, 115.29, 243.33]
            cdf_muon_mass_stat = [0.05, 0.04, 0.01, 0.18, 1.63]
            cdf_muon_mass_sys  = [0.04, 0.07, 0.08, 0.14, 0.40]
            cdf_muon_mass_tot  = np.sqrt(np.square(cdf_muon_mass_stat)) + np.sqrt(np.square(cdf_muon_mass_sys))
            #cdf_muon_mass_tot  = self.makeErrorNumpy(cdf_muon_mass_tot, -1. * cdf_muon_mass_tot, False)

            cdf_electron_mass      = [47.83, 70.76, 90.98, 115.11, 245.46]
            cdf_electron_mass_stat = [0.05, 0.04, 0.01, 0.13, 1.29]
            cdf_electron_mass_sys  = [0.07, 0.04, 0.07, 0.14, 0.21]
            cdf_electron_mass_tot  = np.sqrt(np.square(cdf_electron_mass_stat)) + np.sqrt(np.square(cdf_electron_mass_sys))
            #cdf_electron_mass_tot  = self.makeErrorNumpy(cdf_electron_mass_tot, -1. * cdf_electron_mass_tot, False)

            # Electron
            cdf_muon_pt      = [9.12, 10.81, 11.84, 13.17, 16.18]
            cdf_muon_pt_stat = [0.09, 0.08, 0.03, 0.12, 0.61]
            cdf_muon_pt_sys  = [0.12, 0.14, 0.03, 0.12, 0.45]
            cdf_muon_pt_tot  = np.sqrt(np.square(cdf_muon_pt_stat)) + np.sqrt(np.square(cdf_muon_pt_sys))
            #cdf_muon_pt_tot  = self.makeErrorNumpy(cdf_muon_pt_tot, -1. * cdf_muon_pt_tot, False)

            cdf_electron_pt      = [9.10, 10.84, 11.79, 12.93, 16.41]
            cdf_electron_pt_stat = [0.13, 0.08, 0.02, 0.09, 0.56]
            cdf_electron_pt_sys  = [0.18, 0.10, 0.01, 0.09, 0.35]
            cdf_electron_pt_tot  = np.sqrt(np.square(cdf_electron_pt_stat)) + np.sqrt(np.square(cdf_electron_pt_sys))
            #cdf_electron_pt_tot  = self.makeErrorNumpy(cdf_electron_pt_tot, -1. * cdf_electron_pt_tot, False)

            if self.channel == "electron" :
                ax.errorbar(cdf_electron_mass, cdf_electron_pt, xerr=cdf_electron_mass_tot, yerr=cdf_electron_pt_tot, fmt="s", color="blue", mfc="none", label="CDF electron", linewidth=0.5)
                mass_systematic=self.makeErrorNumpy(cdf_electron_mass_tot, cdf_electron_mass_tot, fromDF=False) 
                pt_systematic=self.makeErrorNumpy(cdf_electron_pt_tot, cdf_electron_pt_tot, fromDF=False) 
                self.doLogLinearFit(ax, cdf_electron_mass, cdf_electron_pt, mass_systematic, pt_systematic, "blue", useDF=False, printPar=False) 
            else :
                ax.errorbar(cdf_muon_mass, cdf_muon_pt, xerr=cdf_muon_mass_tot, yerr=cdf_muon_pt_tot, fmt="s", color="blue", mfc="none", label="CDF muon", linewidth=0.5)
                mass_systematic=self.makeErrorNumpy(cdf_muon_mass_tot, cdf_muon_mass_tot, fromDF=False) 
                pt_systematic=self.makeErrorNumpy(cdf_muon_pt_tot, cdf_muon_pt_tot, fromDF=False) 
                self.doLogLinearFit(ax, cdf_muon_mass, cdf_muon_pt, mass_systematic, pt_systematic, "blue", useDF=False, printPar=False) 

        if showCMS : 

            # Muon
            cms_muon_mass      = [49.32, 73.62, 91.23, 117.33, 240.16, 427.27]
            cms_muon_mass_stat = [0.2, 0.01, 0.00, 0.03, 0.23, 1.75]
            cms_muon_mass_sys  = [0.05, 0.02, 0.08, 0.05, 0.49, 6.05]
            cms_muon_mass_tot  = np.sqrt(np.square(cms_muon_mass_stat)) + np.sqrt(np.square(cms_muon_mass_sys))
            #cms_muon_mass_tot  = self.makeErrorNumpy(cms_muon_mass_tot, -1. * cms_muon_mass_tot, False)

            cms_electron_mass      = [56.39, 73.69, 91.09, 116.89, 240.91, 428.18]
            cms_electron_mass_stat = [0.01, 0.01, 0.0, 0.04, 0.25, 1.77]
            cms_electron_mass_sys  = [0.03, 0.04, 0.06, 0.05, 0.22, 1.14]
            cms_electron_mass_tot  = np.sqrt(np.square(cms_electron_mass_stat)) + np.sqrt(np.square(cms_electron_mass_sys))
            #cms_electron_mass_tot  = self.makeErrorNumpy(cms_electron_mass_tot, -1. * cms_electron_mass_tot, False)

            # Electron
            cms_muon_pt      = [13.57, 16.62, 18.33, 20.04, 24.79, 26.80]
            cms_muon_pt_stat = [0.03, 0.03, 0.01, 0.03, 0.19, 0.34]
            cms_muon_pt_sys  = [0.11, 0.06, 0.04, 0.15, 0.39, 1.69]
            cms_muon_pt_tot  = np.sqrt(np.square(cms_muon_pt_stat)) + np.sqrt(np.square(cms_muon_pt_sys))
            #cms_muon_pt_tot  = self.makeErrorNumpy(cms_muon_pt_tot, -1. * cms_muon_pt_tot, False)

            cms_electron_pt      = [14.76, 16.85, 18.37, 20.39, 25.93, 27.66]
            cms_electron_pt_stat = [0.06, 0.04, 0.01, 0.03, 0.19, 0.41]
            cms_electron_pt_sys  = [0.14, 0.12, 0.05, 0.12, 0.42, 1.27]
            cms_electron_pt_tot  = np.sqrt(np.square(cms_electron_pt_stat)) + np.sqrt(np.square(cms_electron_pt_sys))
            #cms_electron_pt_tot  = self.makeErrorNumpy(cms_electron_pt_tot, -1. * cms_electron_pt_tot, False)

            cms_combined_pt      = [13.57, 14.76, 16.65, 18.35, 20.26, 25.31, 27.52]
            cms_combined_pt_stat = [0.03, 0.06, 0.03, 0.01, 0.02, 0.13, 0.35]
            cms_combined_pt_sys  = [0.11, 0.14, 0.07, 0.04, 0.11, 0.33, 1.27]
            cms_combined_pt_tot  = np.sqrt(np.square(cms_combined_pt_stat)) + np.sqrt(np.square(cms_combined_pt_sys))

            cms_combined_mass      = [49.32, 56.39,73.63, 91.10,117.18, 240.71, 428.09]
            cms_combined_mass_stat = [0.2, 0.01,0.01,0.00, 0.02, 0.19, 1.61]
            cms_combined_mass_sys  = [0.05, 0.03,0.02,0.06, 0.05, 0.22, 1.20]
            cms_combined_mass_tot  = np.sqrt(np.square(cms_combined_mass_stat)) + np.sqrt(np.square(cms_combined_mass_sys))

            ax.errorbar(cms_combined_mass, cms_combined_pt, xerr=cms_combined_mass_tot, yerr=cms_combined_pt_tot, fmt="s", color="black", mfc="none", label="CMS combined", linewidth=0.5)
            mass_systematic=self.makeErrorNumpy(cms_combined_mass_tot, cms_combined_mass_tot, fromDF=False) 
            pt_systematic=self.makeErrorNumpy(cms_combined_pt_tot, cms_combined_pt_tot, fromDF=False) 
            self.doLogLinearFit(ax, cms_combined_mass, cms_combined_pt, mass_systematic, pt_systematic, "black", useDF=False, printPar=True) 

            ax.errorbar(cms_electron_mass, cms_electron_pt, xerr=cms_electron_mass_tot, yerr=cms_electron_pt_tot, fmt="s", color="blue", mfc="none", label="CMS electron", linewidth=0.5)
            mass_systematic=self.makeErrorNumpy(cms_electron_mass_tot, cms_electron_mass_tot, fromDF=False) 
            pt_systematic=self.makeErrorNumpy(cms_electron_pt_tot, cms_electron_pt_tot, fromDF=False) 
            #self.doLogLinearFit(ax, cms_electron_mass, cms_electron_pt, mass_systematic, pt_systematic, "blue", useDF=False, printPar=False) 

            ax.errorbar(cms_muon_mass, cms_muon_pt, xerr=cms_muon_mass_tot, yerr=cms_muon_pt_tot, fmt="s", color="red", mfc="none", label="CMS muon", linewidth=0.5)
            mass_systematic=self.makeErrorNumpy(cms_muon_mass_tot, cms_muon_mass_tot, fromDF=False) 
            pt_systematic=self.makeErrorNumpy(cms_muon_pt_tot, cms_muon_pt_tot, fromDF=False) 
            #self.doLogLinearFit(ax, cms_muon_mass, cms_muon_pt, mass_systematic, pt_systematic, "blue", useDF=False, printPar=False) 

        if combined_result is not None :

            mass_err = np.sqrt(np.square(combined_result["mass"]["stat"]) + np.square(combined_result["mass"]["syst"]))
            pt_err = np.sqrt(np.square(combined_result["pt"]["stat"]) + np.square(combined_result["pt"]["syst"]))
            ax.errorbar(combined_result["mass"]["Nominal"], combined_result["pt"]["Nominal"], xerr=mass_err, yerr=pt_err, fmt="s", color="black", mfc="none", label="combined result", linewidth=0.5)

            mass_systematic=self.makeErrorNumpy(mass_err, mass_err, fromDF=False)
            pt_systematic=self.makeErrorNumpy(pt_err, pt_err, fromDF=False)

            self.doLogLinearFit(ax, combined_result["mass"]["Nominal"], combined_result["pt"]["Nominal"], mass_systematic,  pt_systematic, "black", useDF=False)

        ax.legend(loc='best', fontsize=15, fancybox=False, framealpha=0.0)

        outPdfName = self.outDirPath+self.plotPrefix+"_ISR_comparison.pdf"
        if outPdfPostfix is not None :
            outPdfName = self.outDirPath+self.plotPrefix+"_ISR_comparison.pdf"

        #plt.savefig(self.outDirPath+self.plotPrefix+"_"+outPdfPostfix + "_comparison.pdf", format="pdf", dpi=300)
        plt.savefig(self.outDirPath+self.plotPrefix+"_test_comparison.pdf", format="pdf", dpi=300)
        plt.close(fig)

    def make_error_boxes(self, ax, xdata, ydata, xerror, yerror,
                         showBox=True, showBar=False, facecolor='red', edgecolor='None', alpha=0.5, zorder=5, hatch_style=None, barFmt="None", barMfc='none'):

        # Loop over data points; create box from errors at each point
        try :
            if hatch_style is None :
                errorboxes = [Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum())
                              for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T)]
            else :
                errorboxes = [Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum(), fill=False, facecolor=None, 
                              edgecolor=None, hatch=hatch_style, linewidth=0.)
                              for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T)]
        except RuntimeWarning :
            for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T) :
                print("x, y, xe, ye: ", x, y, xe, ye)

        # Create patch collection with specified colour/alpha
        if showBox :
        
            if hatch_style is None :
                pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                                 edgecolor=edgecolor, zorder = zorder, linewidth=0.05)
            else :
                pc = PatchCollection(errorboxes, zorder = zorder, match_original=True, hatch=hatch_style, linewidth=0. )
            # Add collection to axes
            ax.add_collection(pc)

        if showBar :
            # Plot errorbars
            
            if barFmt != "None" :
                barFmt += facecolor[0]
            artists = ax.errorbar(xdata, ydata, xerr=xerror, yerr=yerror,
                                  fmt=barFmt, ecolor=facecolor, linewidth=0.1, mfc=barMfc)
        #return artists

    def makeErrorNumpy(self, Up, Down, fromDF=True) :

        if fromDF :
            Up=Up.values
        else :
            Up=Up
        Up[np.isinf(Up)] = 0
        Up=Up.reshape(1,len(Up))

        if fromDF :
            Down=Down.values
        else :
            Down=Down
        Down[np.isinf(Down)] = 0
        Down=Down.reshape(1,len(Down))

        UpDown=Up
        UpDown=np.append(UpDown, Down, axis = 0)

        return UpDown

    def getRawDict(self, dictName) :
        return self.rawHistsDict[dictName]

    def getDict(self, dictName="Data") :
        return self.dfs[dictName]
