import gc
import time
import os
import json
import pandas as pd
import numpy as np
import ROOT as rt
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, LogLocator, NullFormatter, LogFormatter, FixedLocator, FixedFormatter)
                               
from matplotlib.pyplot import cm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.ticker as mtick
from matplotlib.lines import Line2D

from scipy.odr import *
from scipy import stats

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

# root TH1 to DataFrame and plot using matplotlib
class ISRPlotter :
    
    def __init__(self, inputHistFilePath, jasonConfigFilePath, doSystematic=False, verbose=True, setQuantile=False, bin_names=None) :

        self.setQuantile = setQuantile

        # Set using self.binDef
        self.binDef={} # dictionary of TUnfoldBinning object
        self.massBins=[] # list of tuple, ex) [(40.,64.), (64., 81.), (81., 101.), (101., 200.), (200., 320.)]

        self.histTypes=["Data", "DataBkgMCSubtracted", "SigMC", "BkgMC", "MCTotal", "Hist"]

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

        self.lumi={"2016": 35.9, "2016a": 19.5, "2016b": 16.8, "2017": 41.5, "2018": 59.8}

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

            self.samples    = self.config['Samples']
            self.stackOrder = self.config['StackOrder']

            # Out directory
            self.outDirPath = "output/UltraLegacy/"+self.year+"/"+self.channel+"/"
            if "ChannelPostfix" in self.config :
                self.channelPostfix = self.config['ChannelPostfix']
                self.outDirPath     = "output/UltraLegacy/"+self.year+"/"+self.channel+"_"+self.channelPostfix+"/"
            if not os.path.exists(self.outDirPath) :
                os.makedirs(self.outDirPath)

        self.tunfoldBinNames=bin_names
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
            
            for var in ["[dimass-dipt]", "[dipt-dimass]"] :
                self.binDef[var]=self.inRootFile.Get(self.topDirName+"/[tunfold-bin]_"+var+"_"+self.tunfoldBinNames[var])

            # TODO use GetDistributionDimension()
            temp_tvecd=self.binDef["[dipt-dimass]"].GetDistributionBinning(1)
            temp_mass_bin_edges=temp_tvecd.GetMatrixArray()
            self.massBins.extend([ (temp_mass_bin_edges[index], temp_mass_bin_edges[index+1]) for index in range(temp_tvecd.GetNrows()-1)]) # list comprehension

        # get raw TH1 histograms
        # this code assume a histogram for a variable stored in a direcotry with the variable name as the directory name 
        for variable in self.variables : 
            varDir=variable
            
            if self.useTUnfoldBin : # FIXME use information of the second axis
                varDir=variable.split("__")[0] # 2D_dipt_dimass__0

            for combinedName in self.samples :
                # initialize
                hist_type = combinedName.split("_")[1]
                if variable not in self.rawHistsDict[hist_type] : # rawHistDict["Data"]["Pt__0"]
                    self.rawHistsDict[hist_type][variable] = dict()
                    temp_dict = self.rawHistsDict[hist_type][variable]
                    
                    if "BkgMC" in combinedName : self.bkgUsed = True
                else :
                    temp_dict = self.rawHistsDict[hist_type][variable]
               
                for mcFileName in self.samples[combinedName] :
                    temp_dict[mcFileName]=dict()

                    # Set nominal histograms
                    temp_dict[mcFileName]["Nominal"]=dict()
                    temp_dict[mcFileName]["Nominal"]["Nominal"]=dict()
                    self.set_rawHistsDict(varDir, mcFileName, "Nominal", "Nominal", variable, temp_dict)
                    
                    # "total" in temp_dict
                    if "total" in temp_dict :
                        self.set_rawHistsDict(varDir, mcFileName, "Nominal", "Nominal", variable, temp_dict, True)
                    else :
                        temp_dict["total"] = dict()
                        temp_dict["total"]["Nominal"] = dict()
                        temp_dict["total"]["Nominal"]["Nominal"] = dict()
                        self.set_rawHistsDict(varDir, mcFileName, "Nominal", "Nominal", variable, temp_dict, True)

                    # Set systemaitc histograms
                    for sysCategory in self.systematics.keys() :
                        for sysName, postfixs in self.systematics[sysCategory].items() :
                            temp_dict[mcFileName][sysName]=dict()

                            # sysName in temp_dict["total"]
                            if sysName not in temp_dict["total"] :
                                temp_dict["total"][sysName] = dict()
                            
                            for postfix in postfixs :
                                temp_dict[mcFileName][sysName][postfix]=dict()
                                self.set_rawHistsDict(varDir, mcFileName, sysName, postfix, variable, temp_dict)
                                
                                # postfix in temp_dict["total"]["sysName"]
                                if postfix not in temp_dict["total"][sysName] :
                                    temp_dict["total"][sysName][postfix] = dict()
                                    
                                self.set_rawHistsDict(varDir, mcFileName, sysName, postfix, variable, temp_dict, True)

        if self.bkgUsed :
            self.setBkgSubtractedDataHist()
        self.setMCTotalHists() # signal + background MC

        # convert histogram to DataFrame
        for histType in self.histTypes :
        
            if "DataBkgMCSubtracted" == histType or "BkgMC" == histType : 
                if self.bkgUsed == False : continue
            
            self.setDataFrameWithUnc(histType)

        # calculate uncertainty
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

    def set_rawHistsDict(self, varDir, inputHistFileName, sysName, sysPostfix, variable, targetDict, is_total = False) :
        if is_total :
            if "TH1" in targetDict["total"][sysName][sysPostfix] :
                targetDict["total"][sysName][sysPostfix]["TH1"].Add(targetDict[inputHistFileName][sysName][sysPostfix]["TH1"])
            else :
                targetDict["total"][sysName][sysPostfix]["TH1"] = targetDict[inputHistFileName][sysName][sysPostfix]["TH1"].Clone("total"+inputHistFileName+variable+sysName+sysPostfix)
        else :
            if sysName == "Nominal" and sysPostfix == "Nominal" :
                if self.plotPrefix == "detector" :
                    histToRead = self.topDirName+"/"+inputHistFileName+"/[tunfold-hist]_"+varDir+"_"+self.tunfoldBinNames[varDir] 
                else :
                    histToRead = self.topDirName+"/"+inputHistFileName+"/[tunfold-unfolded_fullphase_hist]_"+varDir+"_"+self.tunfoldBinNames[varDir] 
            else :
                histToRead = self.topDirName + "/"+varDir + "/" + self.histPrefix + inputHistFileName + '_' + sysName + sysPostfix
            # read histogram
            temp_TH1=self.inRootFile.Get(histToRead)
            
            # in case sysName+sysPostrix not exists, use nominal histogram
            if type(temp_TH1) != rt.TH1D :
                histToRead = self.topDirName + "/"+varDir + "/" + self.histPrefix + inputHistFileName
                temp_TH1=self.inRootFile.Get(histToRead)
            
            if self.useTUnfoldBin :
                temp_TH1 = self.binDef[varDir.split("__")[0]].ExtractHistogram(inputHistFileName + variable + sysName + sysPostfix, temp_TH1,
                                                                          0, True, self.steeringTUnfold[variable]) 
            targetDict[inputHistFileName][sysName][sysPostfix]["TH1"] = temp_TH1
            del temp_TH1

    def getOutBaseDir(self) :
        return self.outDirPath

    def loglinear_func(self, p, x):
        return 2.*p[0]*np.log(x)+p[1]
        
    def setBkgSubtractedDataHist(self) :

        temp_dict = self.rawHistsDict["DataBkgMCSubtracted"]
        for variable in self.variables :
            temp_dict[variable] = dict()
            temp_dict[variable]["total"] = dict()

            temp_dict[variable]["total"]["Nominal"] = dict()
            temp_dict[variable]["total"]["Nominal"]["Nominal"] = dict() 

            data_dict   = self.rawHistsDict["Data"][variable]["total"]["Nominal"]["Nominal"]
            mc_bkg_dict = self.rawHistsDict["BkgMC"][variable]["total"]["Nominal"]["Nominal"]

            temp_dict[variable]["total"]["Nominal"]["Nominal"]["TH1"] = data_dict["TH1"].Clone("data_bkg_subtracted")
            temp_dict[variable]["total"]["Nominal"]["Nominal"]["TH1"].Add(mc_bkg_dict["TH1"], -1)

            for sysCategory in self.systematics.keys() :
                for sysName, postfixs in self.systematics[sysCategory].items() :
                    temp_dict[variable]["total"][sysName]=dict()
                    for postfix in postfixs :
                        temp_dict[variable]["total"][sysName][postfix]=dict()

                        data_dict   = self.rawHistsDict["Data"][variable]["total"][sysName][postfix]
                        mc_bkg_dict = self.rawHistsDict["BkgMC"][variable]["total"][sysName][postfix]

                        temp_dict[variable]["total"][sysName][postfix]["TH1"] = data_dict["TH1"].Clone("data_bkg_subtracted")
                        temp_dict[variable]["total"][sysName][postfix]["TH1"].Add(mc_bkg_dict["TH1"], -1)

    def setMCTotalHists(self) :

        temp_dict = self.rawHistsDict["MCTotal"]
        for variable in self.variables :
            temp_dict[variable] = dict()
            temp_dict[variable]["total"] = dict()

            temp_dict[variable]["total"]["Nominal"] = dict()
            temp_dict[variable]["total"]["Nominal"]["Nominal"] = dict()

            # Lets combine MC histograms
            temp_dict[variable]["total"]["Nominal"]["Nominal"]["TH1"] = \
            self.rawHistsDict["SigMC"][variable]["total"]["Nominal"]["Nominal"]["TH1"].Clone("Clone_"+variable+"Nominal"+"Nominal")

            if self.bkgUsed :
                temp_dict[variable]["total"]["Nominal"]["Nominal"]["TH1"].Add(self.rawHistsDict["BkgMC"][variable]["total"]["Nominal"]["Nominal"]["TH1"])

            for sysCategory in self.systematics.keys() :
                for sysName, postfixs in self.systematics[sysCategory].items() :
                    temp_dict[variable]["total"][sysName]=dict()
                    for postfix in postfixs :
                        temp_dict[variable]["total"][sysName][postfix]=dict()
                        temp_dict[variable]["total"][sysName][postfix]["TH1"] = \
                        self.rawHistsDict["SigMC"][variable]["total"][sysName][postfix]["TH1"].Clone("Clone_"+variable+sysName+postfix)

                        if self.bkgUsed :
                            temp_dict[variable]["total"][sysName][postfix]["TH1"].Add(self.rawHistsDict["BkgMC"][variable]["total"][sysName][postfix]["TH1"])

    def setDataFrameWithUnc(self, dictName="Data") :

        in_dict=self.rawHistsDict[dictName]
        out_dict=self.dfs[dictName]

        for variable in in_dict.keys() :
            out_dict[variable]=dict()

            for sample in in_dict[variable].keys() : # Note loop over samples defined in the "dictionary!"
                
                out_dict[variable][sample]=dict()
                out_dict[variable][sample]["rawUnc"]    = self.convertTH1toDataFrame(in_dict[variable][sample]["Nominal"]["Nominal"]["TH1"])
                out_dict[variable][sample]["upDownUnc"] = self.convertTH1toDataFrame(in_dict[variable][sample]["Nominal"]["Nominal"]["TH1"])

                q = 0.5
                if self.useTUnfoldBin :
                    out_dict[variable][sample]["rawUnc_meanValue"]    = self.createMeanDataFrame(variable, in_dict[variable][sample]["Nominal"]["Nominal"]["TH1"])
                    out_dict[variable][sample]["upDownUnc_meanValue"] = self.createMeanDataFrame(variable, in_dict[variable][sample]["Nominal"]["Nominal"]["TH1"])

                    # quantile
                    if "[dipt-dimass]" in variable and self.setQuantile: 
                        #print("call createQuantileDataFrame")
                        out_dict[variable][sample]["rawUnc_qValue"]    = self.createQuantileDataFrame(variable, in_dict[variable][sample]["Nominal"]["Nominal"]["TH1"], q, True)
                        out_dict[variable][sample]["upDownUnc_qValue"] = self.createQuantileDataFrame(variable, in_dict[variable][sample]["Nominal"]["Nominal"]["TH1"], q, True)

                for sysCategory in self.systematics.keys() :
                    for sysName, postfixs in self.systematics[sysCategory].items() :

                        if sysName == "Nominal" : continue

                        # temp dictionary
                        temp_rawUnc_dict = {}
                        temp_rawUnc_meanValue_dict = {}
                        if self.setQuantile :
                            temp_rawUnc_qValue_dict = {}

                        for postfix in postfixs :

                            if "fsr" not in sysName :

                                temp_rawUnc_dict[sysName+"_"+postfix] = \
                                root_to_numpy(in_dict[variable][sample][sysName][postfix]["TH1"])[0][0]-root_to_numpy(in_dict[variable][sample]["Nominal"]["Nominal"]["TH1"])[0][0]


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

                                    if "[dipt-dimass]" in variable and self.setQuantile:
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

                        if "[dipt-dimass]" in variable and self.setQuantile:
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

                            if "[dipt-dimass]" in variable and self.setQuantile:
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
                                if "[dipt-dimass]" in variable and self.setQuantile:
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

        if "[dipt-dimass]" in variable and self.useTUnfoldBin:
            num_str      = variable.split("__")[1]
            nth_mass_bin = int(num_str)
            temp_mean    = TH1_hist.GetMean()
            temp_dict    = {"low mass cut": self.massBins[nth_mass_bin][0], "high mass cut": self.massBins[nth_mass_bin][1], "mean": temp_mean, "stat_error": TH1_hist.GetMeanError()}
            temp_df      = pd.DataFrame([temp_dict], columns=['low mass cut','high mass cut','mean', 'stat_error'])

            return temp_df

        else :
            if "[dimass-dipt]" in variable :
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
                if "[dipt-dimass]" not in variable :
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



    def convertTH1toDataFrame(self, temp_TH1, existing_df = None) :

        (temp_content, temp_binEdge), _ = root_to_numpy(temp_TH1)
        binWidth = (temp_binEdge[1:] - temp_binEdge[:-1])
        nBin = len(temp_binEdge)-1

        temp_stat_error = []
        for ibin in range(1, len(temp_content)+1) :
            temp_stat_error.append(temp_TH1.GetBinError(ibin))

        temp_stat_error_array = np.array(temp_stat_error)

        pd_series_binIndex    = pd.Series(range(1,nBin+1), range(1, nBin+1), name="bin_index")
        pd_series_binWidth    = pd.Series(binWidth, range(1, nBin+1), name="bin_width")
        pd_series_lowBinEdge  = pd.Series(temp_binEdge[0:-1], range(1, nBin+1), name="low_bin_edge")
        pd_series_highBinEdge = pd.Series(temp_binEdge[1:], range(1, nBin+1), name="high_bin_edge")

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

    def drawHistPlot(self, variable, divde_by_bin_width = False, setLogy=False, setLogx=False,
                     ratio_max=1.35, ratio_min=0.65,
                     figSize=(10,6), show_ratio=True, 
                     denominator="total_Data", setRatioLogy=False,
                     showNEvents=False, showChi2=False, outPdfPostfix=None, ratioName=None, top=0.9, bottom=0.1, xLabel="") :

        hep.style.use("CMS")
        n_columns = 1

        if show_ratio == True :
            num_rows = 2
            fig, axes = plt.subplots(num_rows, n_columns, sharex=False, figsize=figSize, gridspec_kw={'height_ratios':[1, 0.3]})

            top_axis = axes[0]
            bottom_axis = axes[1]
            top_axis.yaxis.get_major_ticks()[0].set_visible(False)
            bottom_axis.yaxis.get_minor_ticks()[0].set_visible(True)
        else :
            num_rows = 1
            fig, axes = plt.subplots(num_rows, n_columns, sharex=False, figsize=figSize)
            top_axis = axes[0] 

        plt.tight_layout()
        plt.subplots_adjust(left=0.15, right=0.95, bottom=bottom, top=top, hspace=0.05)

        denominator_hist_type = denominator.split("_")[1]
        denominator_hist_name = denominator.split("_")[0]
        denominator_label_name = denominator.split("_")[1] 
        if len(denominator.split("_")) > 2 :
            denominator_label_name = denominator.split("_")[2:]

        denominator_df=None
        denominator_hist=None

        nominator_hist_type="MCTotal"
        nominator_hist_name="total"
        nominator_hist=None

        minimum_content = 1
        df_filter = self.dfs[denominator_hist_type][variable][denominator_hist_name]["upDownUnc"]['content'] > minimum_content

        # set denominator
        denominator_df   = self.dfs[denominator_hist_type][variable][denominator_hist_name]["upDownUnc"].copy()
        denominator_hist = self.rawHistsDict[denominator_hist_type][variable][denominator_hist_name]["Nominal"]["Nominal"]["TH1"]

        # set nominators
        nominator_df      = self.dfs[nominator_hist_type][variable][nominator_hist_name]["upDownUnc"].copy()
        nominator_hist    = self.rawHistsDict[nominator_hist_type][variable][nominator_hist_name]["Nominal"]["Nominal"]["TH1"]
            
        if show_ratio : 
            ratio, ratio_error = divide(nominator_hist, denominator_hist, output_type='numpy')

        # bins
        bins = denominator_df.low_bin_edge.values
        last_edge = denominator_df.high_bin_edge.values[-1]
        bins = np.append(bins, last_edge)

        hep.cms.label("Preliminary", data=True, lumi=self.lumi[self.year], year=self.year, ax=top_axis, fontsize=20, loc=1)
        # draw
        if self.plotPrefix == "detector" :
            DY = self.dfs["SigMC"][variable]["total"]["upDownUnc"].copy()
            DY_tau = self.dfs["BkgMC"][variable]["DY_tau"]["upDownUnc"].copy()
            ttbar = self.dfs["BkgMC"][variable]["ttbar"]["upDownUnc"].copy()
            ww = self.dfs["BkgMC"][variable]["ww"]["upDownUnc"].copy()
            wz = self.dfs["BkgMC"][variable]["wz"]["upDownUnc"].copy()
            zz = self.dfs["BkgMC"][variable]["zz"]["upDownUnc"].copy()
            vv = ww + wz + zz
            top = self.dfs["BkgMC"][variable]["top"]["upDownUnc"].copy()
            antitop = self.dfs["BkgMC"][variable]["antitop"]["upDownUnc"].copy()
            stop = top + antitop

            hep.histplot([stop.content,ttbar.content,DY_tau.content, vv.content, DY.content], stack=True, bins=bins, ax=top_axis, binwnorm=1., histtype='fill', 
                         label=["top", "ttbar","tautau", "vv", "DY"], color=['blue', 'red', 'green', '#069AF3', "#C79FEF"])
            hep.histplot((nominator_df.content, bins), ax=top_axis, binwnorm=1., yerr=True, histtype='step', color='r')
        else :
            hep.histplot((nominator_df.content, bins), ax=top_axis, binwnorm=1., yerr=True, histtype='step', label="DY", color='r')
        hep.histplot((denominator_df.content, bins), ax=top_axis, binwnorm=1., yerr=denominator_df.stat_error.values, histtype='errorbar', label=denominator_label_name, color='black')
        top_axis.set_xlim(bins[0], bins[-1])
        top_axis.legend(fontsize=15, loc="upper right")

        hep.histplot((ratio[0], bins), ax=bottom_axis, yerr=ratio_error, histtype='step', label="ratio", color='red')
        bottom_axis.set_xlim(bins[0], bins[-1])
        bottom_axis.set_ylim(ratio_min, ratio_max)

        if setLogy : 
            top_axis.set_yscale("log")
            top_axis.set_ylim(top_axis.get_ylim()[0], top_axis.get_ylim()[1]*10)
        else :
            top_axis.set_ylim(top_axis.get_ylim()[0], top_axis.get_ylim()[1]*1.1)

        if setLogx :
            if top_axis.get_xlim()[0] == 0 : top_axis.set_xlim((1, top_axis.get_xlim()[1]))
            top_axis.set_xscale("log")
            if show_ratio : 
                if bottom_axis.get_xlim()[0] == 0 : bottom_axis.set_xlim((1, bottom_axis.get_xlim()[1]))
                bottom_axis.set_xscale("log")

        if show_ratio :
            top_axis.set_xticklabels([])
            bottom_axis.set_ylabel("MC/Data", fontsize=20) 
            bottom_axis.axhline(y=1, color='black', linestyle='--')
            bottom_axis.set_xlabel(xLabel, fontsize=20)

        top_axis.set_ylabel("Events/bin", fontsize=20)

        if "[dipt-dimass]" in variable :
            num_str = variable.split("__")[1]
            nth_bin = int(num_str)
            top_axis.text(0.03, .8, "{:.0f}".format(self.massBins[nth_bin][0])+"$<m^{\mathit{"+self.channel+"}}<$"+"{:.0f} GeV".format(self.massBins[nth_bin][1]),
                          fontsize='x-small', transform=top_axis.transAxes, horizontalalignment='left')

        # save plot as pdf
        outPdfName = self.outDirPath+self.plotPrefix+"_"+variable+"_"+self.channel+"_"+self.year+".pdf" 
        if outPdfPostfix is not None :
            outPdfName = self.outDirPath+self.plotPrefix+"_"+variable+"_"+outPdfPostfix+".pdf" 
        
        plt.savefig(outPdfName, format="pdf", dpi=300)
        plt.close(fig)

    def combinedPtDataFrame(self, name="Data") :

        in_df=self.dfs[name]

        combined_pt_df = None
        for nth_mass_bin in range(len(self.massBins)) :
            if nth_mass_bin == 0 :
                combined_pt_df=in_df["[dipt-dimass]__"+str(nth_mass_bin)]["total"]["upDownUnc_meanValue"].copy()
            else :
                combined_pt_df=pd.concat([combined_pt_df, in_df["[dipt-dimass]__"+str(nth_mass_bin)]["total"]["upDownUnc_meanValue"].copy()], ignore_index=True)
        return combined_pt_df

    def combinedQuantilePtDataFrame(self, name="Data") :

        in_df=self.dfs[name]

        combined_pt_df = None
        for nth_mass_bin in range(len(self.massBins)) :
            if nth_mass_bin == 0 :
                combined_pt_df=in_df["[dipt-dimass]__"+str(nth_mass_bin)]["total"]["upDownUnc_qValue"].copy()
            else :
                #combined_pt_df=combined_pt_df.append(in_df["[dipt-dimass]__"+str(nth_mass_bin)]["total"]["upDownUnc_qValue"].copy(), ignore_index=True)
                combined_pt_df=pd.concat([combined_pt_df, in_df["[dipt-dimass]__"+str(nth_mass_bin)]["total"]["upDownUnc_qValue"].copy()], ignore_index=True)
        return combined_pt_df


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

    def drawISRPlots(self, *objects_to_plot, names_in_objects, do_linear_fit=None, 
                    labels=None, markers=None, colors=None, facecolor_list = None, 
                    figSize=(8,8), ymin=13, ymax=30, xmin=30, xmax=4e2, outPdfPostfix=None, years = None, both_lepton = False, nominators = None,
                    show_quantile = False) :

        color=iter(cm.rainbow(np.linspace(0,1,len(names_in_objects))))
        isData=False

        hep.style.use("CMS")

        fig, ax = plt.subplots(figsize=figSize)
        hep.cms.label("Preliminary", data=True, lumi=self.lumi[self.year], year=self.year, ax=ax, fontsize=20) 

        plt.tight_layout()

        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.set_xscale("log")

        ax.grid(True, which='both', axis='x', color='black', linewidth=0.3, linestyle="--")
        ax.grid(True, axis='y', color='black', linewidth=0.3, linestyle="--")

        ax.tick_params(bottom=True, top=True, left=True, right=True, which='both', direction='in')
        ax.tick_params(length=10, which='major')
        ax.tick_params(length=5, which='minor')
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))

        minor_x_ticks = FixedLocator([50, 60, 70, 80, 90, 200, 300, 400, 500, 600])
        minor_x_ticklabels = ['50', '60', '70', '', '90', '200', '300', '400', '', '600']
        ax.xaxis.set_minor_locator(minor_x_ticks)
        ax.xaxis.set_minor_formatter(FixedFormatter(minor_x_ticklabels))

        major_x_ticks = FixedLocator([100])
        major_x_ticklabels = ['']
        ax.xaxis.set_major_locator(major_x_ticks)
        ax.xaxis.set_major_formatter(FixedFormatter(major_x_ticklabels))

        ax.set_ylabel("$p_{T}^{"+self.channel+"}$", fontsize=20) 
        ax.set_xlabel(r"$\langle m \rangle^{"+self.channel+"} (GeV)$", fontsize=20) 

        for index, name in enumerate(names_in_objects) :

            isData=False
            if name == "Data" :
                isData=True

            label_name=name
            # set label
            if labels is not None :
                label_name=labels[index]

            # set marker
            if markers is None :
                current_marker='o'
            else :
                current_marker=markers[index]

            # set color
            if colors is None :
                current_color="black"
            else :
                current_color=colors[index]

            # set face color
            if facecolor_list is not None :
                current_facecolor=facecolor_list[index]

            temp_mass_df = self.getDict(name)["[dimass-dipt]"]["total"]["upDownUnc_meanValue"]
            temp_pt_df   = self.combinedPtDataFrame(name)
            if show_quantile : temp_pt_df = self.combinedQuantilePtDataFrame(name)

            temp_mass_total_up   = np.sqrt(np.square(temp_mass_df["total_Up"]) + np.square(temp_mass_df["stat_error"]))
            temp_mass_total_down = np.sqrt(np.square(temp_mass_df["total_Down"]) + np.square(temp_mass_df["stat_error"]))

            temp_pt_total_up   = np.sqrt(np.square(temp_pt_df["total_Up"]) + np.square(temp_pt_df["stat_error"]))
            temp_pt_total_down = np.sqrt(np.square(temp_pt_df["total_Down"]) + np.square(temp_pt_df["stat_error"]))

            mass_systematic = self.makeErrorNumpy(temp_mass_total_up, temp_mass_total_down)
            pt_systematic   = self.makeErrorNumpy(temp_pt_total_up, temp_pt_total_down)

            key = "mean"
            if show_quantile : key = "q_0.5"
            ax.errorbar(temp_mass_df["mean"], temp_pt_df[key], xerr=mass_systematic, yerr=pt_systematic, 
                        fmt=current_marker, color=current_color, label=label_name, linewidth=0.5, ms = 4)

        ax.legend(loc='best', fontsize=15, fancybox=False, framealpha=0.0)

        outPdfName = self.outDirPath+self.plotPrefix+"_ISR_comparison.pdf"
        if outPdfPostfix is not None :
            outPdfName = self.outDirPath+self.plotPrefix+"_ISR_comparison.pdf"

        plt.savefig(self.outDirPath+self.plotPrefix+"_test_comparison.pdf", format="pdf", dpi=300)
        plt.close(fig)

    def makeErrorNumpy(self, Up, Down, fromDF=True) :

        if fromDF :
            Up=Up.values
            Down=Down.values
        else :
            Up=Up
            Down=Down

        Up[np.isinf(Up)] = 0
        Up = Up.reshape(1,len(Up))

        Down[np.isinf(Down)] = 0
        Down = Down.reshape(1,len(Down))

        UpDown = Up
        UpDown = np.append(UpDown, Down, axis = 0)

        return UpDown

    def getRawDict(self, dictName) :
        return self.rawHistsDict[dictName]

    def getDict(self, dictName="Data") :
        return self.dfs[dictName]
