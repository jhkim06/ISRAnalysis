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

from scipy.odr import *
from scipy import stats

# root TH1 to DataFrame and plot using matplotlib

class ISRPlotter :
    def __init__ (self, inputHistFilePath, jasonConfigFilePath, verbose=True) :
        # Set using self.binDef 
        self.binDef={} # dictionary of TUnfoldBinning object
        self.massBins=[] # list of tuple,  ex) [(40.,64.), (64., 81.), (81., 101.), (101., 200.), (200., 320.)]
        
        # Dictionary of raw histograms according DataFrames
        self.rawDataMeasured = {}
        self.rawDataBkgSubtracted = {}
        self.rawMCSignal = {}
        self.rawMCBackground = {}
        self.rawMCTotal = {}
        
        # Dictionary of DataFrames containing nominal and systematic histogram content
        self.Data = {}
        self.DataBkgSubtracted = {}
        self.MCSignal = {}
        self.MCBackground = {}
        self.MCTotal = {}
        
        self.bkgUsed = False
        
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
            self.systematics=self.config['Systematics']
            self.samples=self.config['Samples']
            self.stackOrder=self.config['StackOrder']
        
        if "PDF" in self.systematics["theory"] :
            n_max=int(((self.systematics["theory"]["PDF"][1]).split("_"))[1])
            pdf_prefix=((self.systematics["theory"]["PDF"][1]).split("_"))[0]
            
            self.systematics["theory"]["PDF"].pop()
            self.systematics["theory"]["PDF"]=[pdf_prefix+"{:0>3}".format(str(i)) for i in range(1,n_max+1)]
            
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
            for binName in self.tunfoldBinNames :
                varDirName = binName.split("_")[1]
                if len(self.variablePostfix) > 0 :
                    varDirName = varDirName + "_" + self.variablePostfix
                
                self.binDef[binName.split("_")[1]]=self.inRootFile.Get(self.topDirName+"/"+varDirName+"/"+binName)
                
            # Set massBins
            temp_tvecd=self.binDef["Pt"].GetDistributionBinning(1)
            temp_mass_bin_edges=temp_tvecd.GetMatrixArray()
            self.massBins.extend([ (temp_mass_bin_edges[index], temp_mass_bin_edges[index+1]) for index in range(temp_tvecd.GetNrows()-1)]) # list comprehension
        else :
            self.massBins.extend([(50, 64), (64, 81), (81, 101), (101, 200), (200, 320), (320,380), (380, 440), (440, 510), (510,600)])

        # Out directory
        self.outDirPath="output/"+self.year+"/"+self.channel+"/"
        if not os.path.exists(self.outDirPath) :
            os.makedirs(self.outDirPath)
        
        # Get raw TH1 histograms
        # Convert them into DataFrame
        for variable in self.variables :
            # dict[variable]
            self.rawDataMeasured[variable]=dict()
            self.rawMCSignal[variable]=dict()
            self.rawMCBackground[variable]=dict()
            self.rawMCTotal[variable]=dict()

            varDir=variable
            if self.useTUnfoldBin :
                varDir=variable.split("_")[0] # In case, TUnfoldBinning used
                if len(self.variablePostfix)>0 :
                    varDir=varDir+"_"+self.variablePostfix
            
            # dict[variable][sample]
            for combinedName in self.samples :
                if "Measured" in combinedName : temp_dict=self.rawDataMeasured[variable]
                if "Signal" in combinedName : temp_dict=self.rawMCSignal[variable]
                if "Background" in combinedName : 
                    temp_dict=self.rawMCBackground[variable]
                    self.bkgUsed = True
                
                temp_dict[combinedName]=dict() 
                first_sample=True
                for sampleName in self.samples[combinedName] :
                    temp_dict[sampleName]=dict()
                    
                    for sysCategory in self.systematics.keys() :  
                        for sysName, postfixs in self.systematics[sysCategory].items() :
                            temp_dict[sampleName][sysName]=dict()
                            if first_sample : temp_dict[combinedName][sysName]=dict()
                            for postfix in postfixs :
                                temp_dict[sampleName][sysName][postfix]=dict()
                                if first_sample : temp_dict[combinedName][sysName][postfix]=dict()
                                
                                # Get TH1 object!
                                if sysName == "Nominal" :
                                    temp_TH1=self.inRootFile.Get(self.topDirName+"/"+varDir+"/"+self.histPrefix+sampleName)
                                else :
                                    temp_TH1=self.inRootFile.Get(self.topDirName+"/"+varDir+"/"+self.histPrefix+sampleName+'_'+sysName+postfix)
                                    if type(temp_TH1) != rt.TH1D :
                                        temp_TH1=self.inRootFile.Get(self.topDirName+"/"+varDir+"/"+self.histPrefix+sampleName)
                                    
                                if self.useTUnfoldBin :
                                    temp_TH1=self.binDef[varDir.split("_")[0]].ExtractHistogram(sampleName+variable+sysName+postfix, temp_TH1, 0, True, self.steeringTUnfold[variable]) 
                                
                                temp_dict[sampleName][sysName][postfix]["TH1"]=temp_TH1
                                temp_dict[sampleName][sysName][postfix]["DataFrame"]=self.convertTH1toDataFrame(temp_TH1)
                                
                                if first_sample :
                                    temp_dict[combinedName][sysName][postfix]["TH1"]=temp_TH1.Clone("Clone_"+combinedName+sampleName)
                                    temp_dict[combinedName][sysName][postfix]["DataFrame"]=self.convertTH1toDataFrame(temp_TH1)
                                else :
                                    temp_dict[combinedName][sysName][postfix]["TH1"].Add(temp_TH1.Clone("Clone_"+combinedName+sampleName))
                                    temp_dict[combinedName][sysName][postfix]["DataFrame"].loc[:,"content":]= \
                                    temp_dict[combinedName][sysName][postfix]["DataFrame"].loc[:,"content":]+self.convertTH1toDataFrame(temp_TH1).loc[:,"content":]
                                    
                    
                    first_sample=False
                                
            for sysCategory in self.systematics.keys() :
                for sysName, postfixs in self.systematics[sysCategory].items() :
                    self.rawMCTotal[variable][sysName]=dict()
                    for postfix in postfixs :
                        self.rawMCTotal[variable][sysName][postfix]=dict()
                        
        self.combineHists(self.rawDataMeasured)
        self.combineHists(self.rawMCSignal)
        if self.bkgUsed : self.combineHists(self.rawMCBackground)
        self.setMCTotalHists()
        if self.bkgUsed : self.setBkgSubtractedDataHis()
        
        self.setDataFrames("Data")
        if self.bkgUsed : self.setDataFrames("DataBkgSub")
        self.setDataFrames("Signal")
        if self.bkgUsed : self.setDataFrames("Background")
        self.setDataFrames("TotalMC")
        
        self.calculateCombinedUnc("Data", "total")
        self.calculateCombinedUnc("Data", "theory")
        self.calculateCombinedUnc("Data", "measurement")
        
        self.calculateCombinedUnc("Data", "total", "refinedUnc_meanValue")
        self.calculateCombinedUnc("Data", "theory", "refinedUnc_meanValue")
        self.calculateCombinedUnc("Data", "measurement", "refinedUnc_meanValue")
        
        if len(self.rawMCBackground)>0 :
            self.calculateCombinedUnc("DataBkgSub", "total")
            self.calculateCombinedUnc("DataBkgSub", "theory")
            self.calculateCombinedUnc("DataBkgSub", "measurement")
        
            self.calculateCombinedUnc("Background", "total")
            self.calculateCombinedUnc("Background", "theory") 
            self.calculateCombinedUnc("Background", "measurement")
        
            self.calculateCombinedUnc("DataBkgSub", "total", "refinedUnc_meanValue")
            self.calculateCombinedUnc("DataBkgSub", "theory", "refinedUnc_meanValue")
            self.calculateCombinedUnc("DataBkgSub", "measurement", "refinedUnc_meanValue")
            
            self.calculateCombinedUnc("Background", "total", "refinedUnc_meanValue")
            self.calculateCombinedUnc("Background", "theory", "refinedUnc_meanValue")
            self.calculateCombinedUnc("Background", "measurement", "refinedUnc_meanValue")
    
        self.calculateCombinedUnc("Signal", "total")
        self.calculateCombinedUnc("Signal", "theory")
        self.calculateCombinedUnc("Signal", "measurement")
        
        self.calculateCombinedUnc("Signal", "total", "refinedUnc_meanValue")
        self.calculateCombinedUnc("Signal", "theory", "refinedUnc_meanValue")
        self.calculateCombinedUnc("Signal", "measurement", "refinedUnc_meanValue")
        
        self.calculateCombinedUnc("TotalMC", "total")
        self.calculateCombinedUnc("TotalMC", "theory")
        self.calculateCombinedUnc("TotalMC", "measurement")
        
        self.calculateCombinedUnc("TotalMC", "total", "refinedUnc_meanValue")
        self.calculateCombinedUnc("TotalMC", "theory", "refinedUnc_meanValue")
        self.calculateCombinedUnc("TotalMC", "measurement", "refinedUnc_meanValue")

    def loglinear_func(self, p, x):
        return 2.*p[0]*np.log(x)+p[1]
    
    def setBkgSubtractedDataHis(self) :
         
        for variable in self.variables :
            self.rawDataBkgSubtracted[variable]=dict()
            self.rawDataBkgSubtracted[variable]["total"]=dict()
            
            for sysCategory in self.systematics.keys() :
                for sysName, postfixs in self.systematics[sysCategory].items() :
                    self.rawDataBkgSubtracted[variable]["total"][sysName]=dict()
                    for postfix in postfixs :
                        self.rawDataBkgSubtracted[variable]["total"][sysName][postfix]=dict()
                        temp_dict=self.rawDataBkgSubtracted[variable]["total"][sysName][postfix]
                        
                        data_dict=self.rawDataMeasured[variable]["total"][sysName][postfix]
                        mc_bkg_dict=self.rawMCBackground[variable]["total"][sysName][postfix]
                        
                        temp_dict["TH1"]=data_dict["TH1"].Clone("data_bkg_subtracted")
                        temp_dict["TH1"].Add(mc_bkg_dict["TH1"], -1)
                        
                        temp_dict["DataFrame"]=data_dict["DataFrame"].copy()
                        temp_dict["DataFrame"].content=data_dict["DataFrame"].content-mc_bkg_dict["DataFrame"].content
                        
    def setMCTotalHists(self) : 

        temp_dict = self.rawMCTotal
        for variable in self.variables :
            # dict[variable]
            temp_dict[variable]=dict()
            temp_dict[variable]["total"]=dict()
            
            for sysCategory in self.systematics.keys() :
                for sysName, postfixs in self.systematics[sysCategory].items() :
                    temp_dict[variable]["total"][sysName]=dict()
                    for postfix in postfixs :
                        temp_dict[variable]["total"][sysName][postfix]=dict()
                        
                        # Lets combine MC histograms
                        temp_dict[variable]["total"][sysName][postfix]["TH1"] = \
                        self.rawMCSignal[variable]["total"][sysName][postfix]["TH1"].Clone("Clone_"+variable+sysName+postfix)
                        temp_dict[variable]["total"][sysName][postfix]["DataFrame"] = \
                        self.rawMCSignal[variable]["total"][sysName][postfix]["DataFrame"].copy()
                        
                        if self.bkgUsed :
                            temp_dict[variable]["total"][sysName][postfix]["TH1"].Add(self.rawMCBackground[variable]["total"][sysName][postfix]["TH1"])
                            temp_dict[variable]["total"][sysName][postfix]["DataFrame"].content= \
                            self.rawMCBackground[variable]["total"][sysName][postfix]["DataFrame"].content+temp_dict[variable]["total"][sysName][postfix]["DataFrame"].content 
                            
    def setDataFrames(self, dictName="Data") :
        
        if dictName == "Data" :
            in_dict=self.rawDataMeasured
            out_dict=self.Data
        if dictName == "DataBkgSub" :
            in_dict=self.rawDataBkgSubtracted
            out_dict=self.DataBkgSubtracted
        if dictName == "Signal" :
            in_dict=self.rawMCSignal
            out_dict=self.MCSignal
        if dictName == "Background" :
            in_dict=self.rawMCBackground
            out_dict=self.MCBackground
        if dictName == "TotalMC" :
            in_dict=self.rawMCTotal
            out_dict=self.MCTotal
        
        for variable in in_dict.keys() :
            out_dict[variable]=dict()
            
            for sample in in_dict[variable].keys() : # Note loop over samples defined in the "dictionary!"
                out_dict[variable][sample]=dict()
                out_dict[variable][sample]["rawUnc"]=in_dict[variable][sample]["Nominal"]["Nominal"]["DataFrame"].copy()
                out_dict[variable][sample]["refinedUnc"]=in_dict[variable][sample]["Nominal"]["Nominal"]["DataFrame"].copy()
                # TODO Create column for statistical uncertainty
                
                # DataFrame for mean values
                # low mass cut, high mass cut, mean value
                # For Mass, make a DataFrame for all the mass bins
                # For pT, make a DataFrame for a mass bin
                out_dict[variable][sample]["rawUnc_meanValue"]    = self.createMeanDataFrame(variable, in_dict[variable][sample]["Nominal"]["Nominal"]["TH1"])
                out_dict[variable][sample]["refinedUnc_meanValue"]= self.createMeanDataFrame(variable, in_dict[variable][sample]["Nominal"]["Nominal"]["TH1"])
                
                for sysCategory in self.systematics.keys() :  
                    for sysName, postfixs in self.systematics[sysCategory].items() :
                        if sysName == "Nominal" : continue
                        for postfix in postfixs :
                            if "fsr" not in sysName :
                                out_dict[variable][sample]["rawUnc"][sysName+"_"+postfix]= \
                                in_dict[variable][sample][sysName][postfix]["DataFrame"].content-in_dict[variable][sample]["Nominal"]["Nominal"]["DataFrame"].content
                        
                                temp_df=self.createMeanDataFrame(variable, in_dict[variable][sample][sysName][postfix]["TH1"])
                                out_dict[variable][sample]["rawUnc_meanValue"][sysName+"_"+postfix]= \
                                temp_df["mean"]-out_dict[variable][sample]["rawUnc_meanValue"]["mean"]
                            else :
                                current_index=postfixs.index(postfix)
                                the_other_postfix=None
                                if current_index == 0 : the_other_postfix=postfixs[1]
                                else : the_other_postfix=postfixs[0]
                                    
                                out_dict[variable][sample]["rawUnc"][sysName+"_"+postfix]= \
                                in_dict[variable][sample][sysName][postfix]["DataFrame"].content-in_dict[variable][sample][sysName][the_other_postfix]["DataFrame"].content
                        
                                temp_df1=self.createMeanDataFrame(variable, in_dict[variable][sample][sysName][postfix]["TH1"])
                                temp_df2=self.createMeanDataFrame(variable, in_dict[variable][sample][sysName][the_other_postfix]["TH1"])
                                out_dict[variable][sample]["rawUnc_meanValue"][sysName+"_"+postfix]= \
                                temp_df1["mean"]-temp_df2["mean"]
                                
                        if sysName != "PDF" :
                            out_dict[variable][sample]["refinedUnc"][sysName+'_Up']  =out_dict[variable][sample]["rawUnc"].filter(like=sysName).max(axis=1)
                            out_dict[variable][sample]["refinedUnc"][sysName+'_Down']=out_dict[variable][sample]["rawUnc"].filter(like=sysName).min(axis=1)
                            
                            out_dict[variable][sample]["refinedUnc_meanValue"][sysName+'_Up']  =\
                            out_dict[variable][sample]["rawUnc_meanValue"].filter(like=sysName).max(axis=1)
                            out_dict[variable][sample]["refinedUnc_meanValue"][sysName+'_Down']=\
                            out_dict[variable][sample]["rawUnc_meanValue"].filter(like=sysName).min(axis=1)
                        else :
                            out_dict[variable][sample]["refinedUnc"][sysName+'_Up']  =np.sqrt(out_dict[variable][sample]["rawUnc"].filter(like=sysName).var(axis=1))/2.
                            out_dict[variable][sample]["refinedUnc"][sysName+'_Down']= -1. * np.sqrt(out_dict[variable][sample]["rawUnc"].filter(like=sysName).var(axis=1))/2.
                            
                            out_dict[variable][sample]["refinedUnc_meanValue"][sysName+'_Up']  =\
                            np.sqrt(out_dict[variable][sample]["rawUnc_meanValue"].filter(like=sysName).var(axis=1))/ 2.
                            out_dict[variable][sample]["refinedUnc_meanValue"][sysName+'_Down']=\
                            -1. * np.sqrt(out_dict[variable][sample]["rawUnc_meanValue"].filter(like=sysName).var(axis=1))/ 2.
                            
    def createMeanDataFrame(self, variable, TH1_hist) :
        
        if "Pt" in variable :
            nth_mass_bin=int(variable.split("_")[1]) # 
            temp_mean=TH1_hist.GetMean()
            temp_dict={"low mass cut": self.massBins[nth_mass_bin][0], "high mass cut": self.massBins[nth_mass_bin][1], "mean": temp_mean}
            temp_df=pd.DataFrame([temp_dict], columns=['low mass cut','high mass cut','mean'])
            
            return temp_df 
        
        if "Mass" in variable :
            # Get mean mass for all mass bins
            #self.massBins
            # Create DataFrame
            row_list=[]
            for massBin in self.massBins :
                TH1_hist.GetXaxis().SetRangeUser(massBin[0], massBin[1])
                temp_mean=TH1_hist.GetMean()
                #print("low mass cut: {} high mass cut: {} mean : {}".format(massBin[0], massBin[1], temp_mean))
                temp_dict={"low mass cut": massBin[0], "high mass cut": massBin[1], "mean": temp_mean}
                row_list.append(temp_dict)
            temp_df = pd.DataFrame(row_list, columns=['low mass cut','high mass cut','mean'])
            
            return temp_df
        
    def calculateCombinedUnc(self, dictName="Data", sys_to_combine="total", column_name="refinedUnc") :
        
        if dictName == "Data" :
            in_dict=self.Data
        if dictName == "DataBkgSub" :
            in_dict=self.DataBkgSubtracted
        if dictName == "Signal" :
            in_dict=self.MCSignal
        if dictName == "Background" :
            in_dict=self.MCBackground
        if dictName == "TotalMC" :
            in_dict=self.MCTotal
        
        # total uncertainty (stat+sys), total theory, total measurement
        for variable in in_dict.keys() :
            for sample in in_dict[variable].keys() :
                
                combinedSys_Up = None
                combinedSys_Down = None
                first_sys=True
                for sysCategory in self.systematics.keys() :
                    if sys_to_combine != "total" :
                        if sysCategory != sys_to_combine : continue
                        
                    for sysName, postfixs in self.systematics[sysCategory].items() :
                        if sysName == "Nominal" : continue # Nominal don't have Up, Down...
                        if first_sys :
                            combinedSys_Up=np.square(in_dict[variable][sample][column_name][sysName+'_Up']).copy()
                            combinedSys_Down=np.square(in_dict[variable][sample][column_name][sysName+'_Down']).copy()
                            first_sys=False
                        else :
                            combinedSys_Up=combinedSys_Up+np.square(in_dict[variable][sample][column_name][sysName+'_Up'])
                            combinedSys_Down=combinedSys_Down+np.square(in_dict[variable][sample][column_name][sysName+'_Down'])
                                
                in_dict[variable][sample][column_name][sys_to_combine+"_Up"]=np.sqrt(combinedSys_Up)
                in_dict[variable][sample][column_name][sys_to_combine+"_Down"]=np.sqrt(combinedSys_Down)
            
    def combineHists(self, mc_dict) :
         
        for variable in mc_dict.keys() :
            mc_dict[variable]["total"]=dict()
            
            for sysCategory in self.systematics.keys() :
                for sysName, postfixs in self.systematics[sysCategory].items() :
                    mc_dict[variable]["total"][sysName]=dict()
                    for postfix in postfixs :
                        mc_dict[variable]["total"][sysName][postfix]=dict()
                        
                        # Lets combine MC histograms
                        first_mc=True
                        for sample in mc_dict[variable].keys() :
                            if sample in self.samples.keys(): continue
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
                                
                                
    def convertTH1toDataFrame(self, temp_TH1) :

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
            'bin_width': pd_series_binWidth, 
            'low_bin_edge': pd_series_lowBinEdge,
            'high_bin_edge': pd_series_highBinEdge,
            'content': temp_content,
            'stat_error': temp_stat_error_array,
        }

        pd_temp=pd.DataFrame(dict_temp, index=pd_series_binIndex)
        
        return pd_temp
   
    # Draw data and all MC and the ratio plot between data and total MC 
    def drawSubPlot(self, top_axis, bottom_axis, variable, divde_by_bin_width = False, setLogy=False, 
                    write_xaxis_title=False, write_yaxis_title=False, showMeanValue=False, setLogx = False) :
        # Get DataFrame
        dataName=""
        dataTH1 = None
        data_df={}
        totalMC=self.MCTotal[variable]["total"]["refinedUnc"].copy()
        totalMCTH1=self.rawMCTotal[variable]['total']['Nominal']['Nominal']['TH1'].Clone("total_mc")
        
        channelName = "e^{+}e^{-}"
        if self.channel == "muon" :
            channelName = "\mu^{+}\mu^{-}"
        
        for combinedName in self.samples :
            if "Data" in combinedName :
                dataName=combinedName
                data_df[combinedName]=self.Data[variable][combinedName]["refinedUnc"].copy()
                dataTH1=self.rawDataMeasured[variable][combinedName]['Nominal']['Nominal']['TH1'].Clone("data")
            if "Signal" in combinedName :
                data_df[combinedName]=self.MCSignal[variable][combinedName]["refinedUnc"].copy()
            if "Background" in combinedName :
                data_df[combinedName]=self.MCBackground[variable][combinedName]["refinedUnc"].copy()

        # Get basic histogram configuration
        x_bin_centers=data_df[dataName]['low_bin_edge'] + data_df[dataName]['bin_width']/2.
        bin_width=data_df[dataName]['bin_width']
        x_min=data_df[dataName]['low_bin_edge'].iloc[0]
        x_max=data_df[dataName]['high_bin_edge'].iloc[-1]
        nbins=len(x_bin_centers)

        if setLogx :
            if x_min == 0 :
                #x_min = x_bin_centers.values[0]
                x_min = 0.1
        
        binWidthnumpy = np.asarray([data_df[dataName]['bin_width']/2.])
        binWidthxerr = np.append(binWidthnumpy, binWidthnumpy, axis=0)
        
        if write_yaxis_title : top_axis.text(0., 1.1, "CMS Work in progress", fontsize=20, transform=top_axis.transAxes)
        if write_xaxis_title : top_axis.text(1., 1.1, "(13 TeV, " + self.year + ")", fontsize=20, transform=top_axis.transAxes, ha='right')
        if "Pt" in variable :
            nth_bin = int(variable.split("_")[1])
            top_axis.text(-0.05, 1.01, str(self.massBins[nth_bin][0])+"$ < M^{\mathit{"+channelName+"}} < $"+str(self.massBins[nth_bin][1]), 
                          fontsize=10, transform=top_axis.transAxes)
            
        if divde_by_bin_width: 
            if write_yaxis_title: top_axis.set_ylabel('Events/\n 1 GeV', fontsize=20, ha='right', y=1.0)
            totalMC["content"]=totalMC["content"]/data_df[combinedName]['bin_width']
            totalMC["stat_error"]=totalMC["stat_error"]/data_df[combinedName]['bin_width']
            totalMC["total_Up"]=totalMC["total_Up"]/data_df[combinedName]['bin_width']
            totalMC["total_Down"]=totalMC["total_Down"]/data_df[combinedName]['bin_width']

            for combinedName in self.samples :
                data_df[combinedName]["content"]=data_df[combinedName]["content"]/ data_df[combinedName]['bin_width']
                data_df[combinedName]["stat_error"]=data_df[combinedName]["stat_error"]/ data_df[combinedName]['bin_width']
                data_df[combinedName]["total_Up"]=data_df[combinedName]["total_Up"]/ data_df[combinedName]['bin_width']
                data_df[combinedName]["total_Down"]=data_df[combinedName]["total_Down"]/ data_df[combinedName]['bin_width']
        else :
            if write_yaxis_title: top_axis.set_ylabel('Events/\n Bin', fontsize=20, ha='right', y=1.0)
            
        color=iter(cm.rainbow(np.linspace(0,1,len(self.stackOrder))))
        
        #top_axis.scatter(x_bin_centers, data_df[dataName]["content"], marker='o', c='black', s=50, zorder=4, label='Data')
        #top_axis.errorbar(x_bin_centers, data_df[dataName]["content"], xerr=bin_width/2., yerr=data_df[dataName]["stat_error"], fmt='o', ecolor='black', zorder=4)
        top_axis.errorbar(x_bin_centers, data_df[dataName]["content"], xerr=bin_width/2., yerr=data_df[dataName]["stat_error"], fmt='ok', ecolor='black', zorder=4)

        if variable == "Mass" :
            print(data_df[dataName]["stat_error"])

        top_axis.set_xlim(x_min, x_max)
        data_abs_systematic=self.makeErrorNumpy(data_df[dataName].total_Up, data_df[dataName].total_Down)
        mc_abs_systematic=self.makeErrorNumpy(totalMC.total_Up, totalMC.total_Down)
        self.make_error_boxes(top_axis, x_bin_centers.values, data_df[dataName]["content"], binWidthxerr, data_abs_systematic, 
                              showBar=False, alpha=0.2, edgecolor='None', facecolor='black', zorder=2)
        self.make_error_boxes(top_axis, x_bin_centers.values, totalMC["content"], binWidthxerr, mc_abs_systematic, 
                              showBar=False, alpha=0.2, edgecolor='None', facecolor='red')
    
        for i, stack in enumerate(self.stackOrder) :
            if i==0 :
                if(len(self.stackOrder)==1) :
                    top_axis.errorbar(x_bin_centers, data_df[stack]['content'], xerr=bin_width/2., yerr=0, fmt='ro', ecolor='red')
                else :
                    top_axis.bar(x_bin_centers, data_df[stack]['content'], width = bin_width, color=next(color), label=stack)
                stacks=data_df[stack]['content']
            else :
                top_axis.bar(x_bin_centers, data_df[stack]['content'], width = bin_width, color=next(color), bottom=stacks, label='Drell-Yan')
                stacks=stacks+data_df[stack]['content']
                
        if setLogy :
            top_axis.set_yscale("log")
            top_axis.set_ylim(1e-1* data_df[dataName]["content"].min(), 10. * data_df[dataName]["content"].max())
            top_axis.yaxis.set_major_locator(LogLocator(10, numticks=14))
            top_axis.yaxis.set_minor_locator(LogLocator(10, subs=(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks=14))
        else :
            top_axis.set_ylim(0., 1.5 * data_df[dataName]["content"].max())

        if setLogx :
            top_axis.set_xscale("log")
            
        if "Mass" in variable :
            for bin_ in self.massBins[:-1] :
                top_axis.axvline(bin_[1], color='black', linewidth=0.5, zorder=3)
                
        if showMeanValue :
            temp_data_mean_df=self.Data[variable]["total"]["refinedUnc_meanValue"]
            if "Mass" in variable :
                for i in temp_data_mean_df.index :
                    #print("Mean: {}".format(data_df["Mass"]["total"]["refinedUnc_meanValue"]["mean"][i]))
                    top_axis.axvline(temp_data_mean_df["mean"][i], color='black', linewidth=1, linestyle=":")
            if "Pt" in variable :
                top_axis.axvline(temp_data_mean_df["mean"][0], color='black', linewidth=1, linestyle=":")
    
        bottom_axis.set_ylim(0.65,1.35)
            
        varName = "$p_{T}^{\mathit{"+channelName+"}}$"
        if variable == "Mass" :
            varName = "Mass$^{\mathit{"+channelName+"}}$"
            
        if write_xaxis_title : bottom_axis.set_xlabel(varName + " [GeV]", fontsize=20, ha='right', x=1.0)
        if write_yaxis_title : bottom_axis.set_ylabel("MC/\n Data", fontsize=20)

        dataTH1.Divide(dataTH1)

        temp_stat_error = []
        for ibin in range(1, dataTH1.GetNbinsX()+1) :
            temp_stat_error.append(dataTH1.GetBinError(ibin))
            
        temp_stat_error_array = np.array(temp_stat_error)
        if variable == "Mass" :
            print(temp_stat_error_array)
    
        ratio=totalMC.content/ data_df[dataName].content
        one_points=data_df[dataName].content/data_df[dataName].content
        #bottom_axis.scatter(x_bin_centers, ratio, facecolors='red', marker="_", edgecolor='red', s=40, zorder=2)
        bottom_axis.errorbar(x_bin_centers, ratio, xerr=bin_width/2., yerr=0, fmt='r,', ecolor='red', zorder=5)
        #bottom_axis.scatter(x_bin_centers, one_points, facecolors='black', edgecolor='black', s=40, zorder=3)
        bottom_axis.errorbar(x_bin_centers, one_points, xerr=0, yerr=temp_stat_error_array, fmt='r,', ecolor='black', zorder=6)
        bottom_axis.axhline(1., color='black', linewidth=1, zorder=1)
        if "Mass" in variable :
            for bin_ in self.massBins[:-1] :
                bottom_axis.axvline(bin_[1], color='black', linewidth=0.5, zorder=6)
        
        systematic=self.makeErrorNumpy(totalMC.total_Up/ data_df[dataName].content, totalMC.total_Down/ data_df[dataName].content)
        data_systematic=self.makeErrorNumpy(data_df[dataName].total_Up/data_df[dataName].content, data_df[dataName].total_Down/data_df[dataName].content)
        self.make_error_boxes(bottom_axis, x_bin_centers.values, ratio.values, binWidthxerr, systematic, 
                              showBar=False, alpha=0.2, edgecolor='None', facecolor='red', zorder=4)
        self.make_error_boxes(bottom_axis, x_bin_centers.values, one_points, binWidthxerr, data_systematic, 
                              showBar=False, alpha=0.2, edgecolor='None', facecolor='black', zorder=2)

    def drawSubPlotComparison(self, object_, nominator_name, denominator_name, top_axis, bottom_axis, variable, name, divde_by_bin_width = False, setLogy=False, 
                    write_xaxis_title=False, write_yaxis_title=False, showMeanValue=False, setLogx = False, marker_ = 'o', color_ = "black", ratio_max=1.35, ratio_min=0.65) :

        # Get DataFrame
        dataName=""
        denominator_df=None
        nominator_df=None

        channelName = "e^{+}e^{-}"
        if self.channel == "muon" :
            channelName = "\mu^{+}\mu^{-}"

        if "Data" in denominator_name :
            if object_ is None :
                nominator_df=self.Data[variable]["Data_Measured"]["refinedUnc"].copy()
                nominator_df.loc[:, "content":] = nominator_df.loc[:, "content":]/nominator_df["content"].sum()
            else :
                nominator_df = object_.Data[variable]["Data_Measured"]["refinedUnc"].copy()
                nominator_df.loc[:, "content":] = nominator_df.loc[:, "content":]/nominator_df["content"].sum()

                denominator_df = self.Data[variable]["Data_Measured"]["refinedUnc"].copy()
                denominator_df.loc[:, "content":] = denominator_df.loc[:, "content":]/denominator_df["content"].sum()
        else :
            if object_ is None :
                nominator_df=self.MCTotal[variable]["total"]["refinedUnc"].copy()
                nominator_df.loc[:, "content":] = nominator_df.loc[:, "content":]/nominator_df["content"].sum()

            else :
                nominator_df=object_.MCTotal[variable]["total"]["refinedUnc"].copy()
                nominator_df.loc[:, "content":] = nominator_df.loc[:, "content":]/nominator_df["content"].sum()

                denominator_df = self.MCTotal[variable]["total"]["refinedUnc"].copy()
                denominator_df.loc[:, "content":] = denominator_df.loc[:, "content":]/denominator_df["content"].sum()

        # Get basic histogram configuration
        x_bin_centers=nominator_df['low_bin_edge'] + nominator_df['bin_width']/2.
        bin_width=nominator_df['bin_width']
        x_min=nominator_df['low_bin_edge'].iloc[0]
        x_max=nominator_df['high_bin_edge'].iloc[-1]
        nbins=len(x_bin_centers)

        if setLogx :
            if x_min == 0 :
                #x_min = x_bin_centers.values[0]
                x_min = 0.1
        
        binWidthnumpy = np.asarray([nominator_df['bin_width']/2.])
        binWidthxerr = np.append(binWidthnumpy, binWidthnumpy, axis=0)
        
        if write_yaxis_title : top_axis.text(0., 1.1, "CMS Work in progress", fontsize=20, transform=top_axis.transAxes)
        if write_xaxis_title : top_axis.text(1., 1.1, "(13 TeV, " + self.year + ")", fontsize=20, transform=top_axis.transAxes, ha='right')
        if "Pt" in variable :
            nth_bin = int(variable.split("_")[1])
            top_axis.text(-0.05, 1.01, str(self.massBins[nth_bin][0])+"$ < M^{\mathit{"+channelName+"}} < $"+str(self.massBins[nth_bin][1]), 
                          fontsize=10, transform=top_axis.transAxes)
            
        if divde_by_bin_width: 
            if write_yaxis_title: top_axis.set_ylabel('A.U.', fontsize=20, ha='right', y=1.0)
            nominator_df["content"]=nominator_df["content"]/ nominator_df['bin_width']
            nominator_df["total_Up"]=nominator_df["total_Up"]/ nominator_df['bin_width']
            nominator_df["total_Down"]=nominator_df["total_Down"]/ nominator_df['bin_width']
            if object_ is not None :
                denominator_df["content"]=denominator_df["content"]/ denominator_df['bin_width']
                denominator_df["total_Up"]=denominator_df["total_Up"]/ denominator_df['bin_width']
                denominator_df["total_Down"]=denominator_df["total_Down"]/ denominator_df['bin_width']

        else :
            if write_yaxis_title: top_axis.set_ylabel('A.U.', fontsize=20, ha='right', y=1.0)
            
        top_axis.scatter(x_bin_centers, nominator_df["content"], marker=marker_, c=color_, s=50, zorder=4, label=name)
        top_axis.set_xlim(x_min, x_max)
        data_abs_systematic=self.makeErrorNumpy(nominator_df.total_Up, nominator_df.total_Down)
        self.make_error_boxes(top_axis, x_bin_centers.values, nominator_df["content"], binWidthxerr, data_abs_systematic, 
                              showBar=False, alpha=0.2, edgecolor='None', facecolor=color_, zorder=2)

        if write_yaxis_title :
            top_axis.legend(loc="lower left") 

        if setLogy :
            top_axis.set_yscale("log")
            top_axis.set_ylim(1e-1 * nominator_df["content"].min(), 10. * nominator_df["content"].max())
            top_axis.yaxis.set_major_locator(LogLocator(10, numticks=14))
            top_axis.yaxis.set_minor_locator(LogLocator(10, subs=(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks=14))
        else :
            top_axis.set_ylim(0., 1.5 * nominator_df["content"].max())

        if setLogx :
            top_axis.set_xscale("log")
            
        if "Mass" in variable :
            for bin_ in self.massBins[:-1] :
                top_axis.axvline(bin_[1], color='black', linewidth=0.5, zorder=3)
                
        if showMeanValue :
            if "Data" in denominator_name :
                if object_ is None :
                    temp_data_mean_df=self.Data[variable]["Data_Measured"]["refinedUnc_meanValue"]
                else :
                    temp_data_mean_df=object_.Data[variable]["Data_Measured"]["refinedUnc_meanValue"]
            else :
                if object_ is None :
                    temp_data_mean_df=self.MCTotal[variable]["total"]["refinedUnc_meanValue"]
                else :
                    temp_data_mean_df=object_.MCTotal[variable]["total"]["refinedUnc_meanValue"]

            if "Mass" in variable :
                for i in temp_data_mean_df.index :
                    #print("Mean: {}".format(nominator_df["Mass"]["total"]["refinedUnc_meanValue"]["mean"][i]))
                    top_axis.axvline(temp_data_mean_df["mean"][i], color=color_, linewidth=1, linestyle=":")
            if "Pt" in variable :
                top_axis.axvline(temp_data_mean_df["mean"][0], color=color_, linewidth=1, linestyle=":")
    
        bottom_axis.set_ylim(ratio_min, ratio_max)
            
        varName = "$p_{T}^{\mathit{"+channelName+"}}$"
        if variable == "Mass" :
            varName = "Mass$^{\mathit{"+channelName+"}}$"
            
        if write_xaxis_title : bottom_axis.set_xlabel(varName + " [GeV]", fontsize=20, ha='right', x=1.0)
        if write_yaxis_title and object_ is None : 
            bottom_axis.set_ylabel("/\n " + name, fontsize=20)
    
        if object_ is not None :
            ratio=nominator_df.content/ denominator_df.content
            #one_points=denominator_df.content/denominator_df.content
            bottom_axis.errorbar(x_bin_centers, ratio, xerr=bin_width/2., yerr=0, fmt=marker_, color = color_, ecolor=color_, zorder=5)
            #bottom_axis.scatter(x_bin_centers, one_points, facecolors='black', edgecolor='black', s=40, zorder=3)
            bottom_axis.axhline(1., color='black', linewidth=1, zorder=1)

            if "Mass" in variable :
                for bin_ in self.massBins[:-1] :
                    bottom_axis.axvline(bin_[1], color='black', linewidth=0.5, zorder=6)
            
            systematic=self.makeErrorNumpy(nominator_df.total_Up/ denominator_df.content, nominator_df.total_Down/ denominator_df.content)
            data_systematic=self.makeErrorNumpy(denominator_df.total_Up/denominator_df.content, denominator_df.total_Down/denominator_df.content)
            self.make_error_boxes(bottom_axis, x_bin_centers.values, ratio.values, binWidthxerr, systematic, 
                                  showBar=False, alpha=0.2, edgecolor='None', facecolor=color_, zorder=4)

            #self.make_error_boxes(bottom_axis, x_bin_centers.values, one_points, binWidthxerr, data_systematic, 
            #                      showBar=False, alpha=1.0, edgecolor='None', facecolor='white', zorder=2, hatch_style="/////")

    def drawHistPlotComparison(self, *object_list, nominator_name, denominator_name, variables, names, markers, colors, divde_by_bin_width = False, setLogy=False, showMeanValue=False, setLogx=False, ratio_max=1.35, ratio_min=0.65) :
        
        fig, axes = plt.subplots(2, len(variables), sharex=True, figsize=(10, 6), gridspec_kw={'hspace': 0, 'height_ratios':[1, 0.4]})
        fig.tight_layout()
        plt.subplots_adjust(left=0.12, right=0.97, bottom=0.15, top=0.85)
        
        write_xaxis_title = True
        write_yaxis_title = True

        for index, variable in enumerate(variables) :

            if len(variables) == 1:
                self.drawSubPlotComparison(None, nominator_name, denominator_name, axes[0], axes[1], variable, names[0], 
                                            divde_by_bin_width, setLogy, write_xaxis_title, write_yaxis_title,showMeanValue, setLogx, "o", "black", ratio_max, ratio_min)

            else :
                if index > 0 : 
                    write_yaxis_title=False
                if index < len(variables)-1 : 
                    write_xaxis_title=False
                else :
                    write_xaxis_title=True

                self.drawSubPlotComparison(None, nominator_name, denominator_name, axes[0][index], axes[1][index], variable, names[0], 
                                            divde_by_bin_width, setLogy, write_xaxis_title, write_yaxis_title,showMeanValue, setLogx, "o", "black", ratio_max, ratio_min) 

        write_xaxis_title = True
        write_yaxis_title = True

        for index, variable in enumerate(variables) :

            for object_index, object_ in enumerate(object_list) :

                if len(variables) == 1:
                    self.drawSubPlotComparison(object_, nominator_name, denominator_name, axes[0], axes[1], variable, names[object_index+1], 
                                                divde_by_bin_width, setLogy, write_xaxis_title, write_yaxis_title,showMeanValue, setLogx, markers[object_index], colors[object_index], ratio_max, ratio_min) 
                else :
                    if index > 0 : 
                        write_yaxis_title=False
                    if index < len(variables)-1 : 
                        write_xaxis_title=False
                    else :
                        write_xaxis_title=True

                    self.drawSubPlotComparison(object_, nominator_name, denominator_name, axes[0][index], axes[1][index], variable, names[object_index+1], 
                                                divde_by_bin_width, setLogy, write_xaxis_title, write_yaxis_title,showMeanValue, setLogx, markers[object_index], colors[object_index], ratio_max, ratio_min) 

        plt.savefig(self.outDirPath+self.plotPrefix+self.topDirName+"_Comparison"+variable+"_"+denominator_name+".pdf", format="pdf", dpi=300)
        plt.close(fig)
    
    def drawHistPlot(self, *variables, divde_by_bin_width = False, setLogy=False, showMeanValue=False, setLogx=False) :
        
        variable=variables[0]
        
        fig, axes = plt.subplots(2,len(variables), sharex=True, figsize=(10, 6), gridspec_kw={'hspace': 0, 'height_ratios':[1, 0.4]})
        fig.tight_layout()
        plt.subplots_adjust(left=0.12, right=0.97, bottom=0.15, top=0.85)
        
        write_xaxis_title = True
        write_yaxis_title = True

        for index, variable in enumerate(variables) :
            if len(variables) == 1:
                self.drawSubPlot(axes[0], axes[1], variable, divde_by_bin_width, setLogy, write_xaxis_title, write_yaxis_title,showMeanValue, setLogx)
            else :
                if index > 0 : 
                    write_yaxis_title=False
                if index < len(variables)-1 : 
                    write_xaxis_title=False
                else :
                    write_xaxis_title=True
                    
                self.drawSubPlot(axes[0][index], axes[1][index], variable, divde_by_bin_width, setLogy, write_xaxis_title, write_yaxis_title, showMeanValue, setLogx)
        
        plt.savefig(self.outDirPath+self.plotPrefix+self.topDirName+"_"+variable+".pdf", format="pdf", dpi=300)
        plt.close(fig)
        
    def combinedPtDataFrame(self, name="Data") :
        
        if name == "Data" :
            in_df=self.Data
        if name == "DataBkgSub" :
            in_df=self.DataBkgSubtracted
        if name == "Signal" :
            in_df=self.MCSignal
        if name == "Background" :
            in_df=self.MCBackground
        if name == "TotalMC" :
            in_df=self.MCTotal
        
        combined_pt_df = None
        for nth_mass_bin in range(len(self.massBins)) :
            if nth_mass_bin == 0 :
                combined_pt_df=in_df["Pt_"+str(nth_mass_bin)]["total"]["refinedUnc_meanValue"].copy()
            else :
                combined_pt_df=combined_pt_df.append(in_df["Pt_"+str(nth_mass_bin)]["total"]["refinedUnc_meanValue"].copy(), ignore_index=True)
        return combined_pt_df
    
    def doLogLinearFit(self, ax, mass_df, pt_df, mass_unc, pt_unc, line_color) :
        
        print("do log-linear fit!")
        x_err=[]
        y_err=[]
        for i in range(len(mass_unc.T)) :
            x_err.append(max(mass_unc.T[i]))
            y_err.append(max(pt_unc.T[i]))
            
        loglinear=Model(self.loglinear_func)
        data=RealData(mass_df["mean"], pt_df["mean"], sx=x_err, sy=y_err)
        odr=ODR(data, loglinear, beta0=[1.0, 0.0])
        out=odr.run()
        out.pprint()
        
        xn = np.linspace(min(mass_df["mean"]), max(mass_df["mean"]), 1000)
        yn = self.loglinear_func(out.beta, xn)
        ax.plot(xn, yn, color=line_color, linewidth=0.8)
        # prepare parameters for confidence interval curves
        nstd = 1. # to draw 1-sigma intervals
        popt_up = out.beta + nstd * out.sd_beta
        popt_dw = out.beta - nstd * out.sd_beta
        
        # calculate y values for 1 sigma
        fit_up = self.loglinear_func(popt_up, xn)
        fit_dw = self.loglinear_func(popt_dw, xn)
        
        # plot 1 sigma
        #ax.plot(xn, fit_up, '#0076D4', dashes=[9, 4.5], label='1 Sigma uncertainty', linewidth=0.8)
        #ax.plot(xn, fit_dw, '#0076D4', dashes=[9, 4.5], linewidth=0.8)
        
    def drawISRPlot(self, *list_to_plot, do_linear_fit=False) :
        print("draw isr plot, do_linear_fit {}".format(do_linear_fit))
        color=iter(cm.rainbow(np.linspace(0,1,len(list_to_plot))))
        isData=False
        
        fig, ax = plt.subplots(figsize=(10, 6))
        plt.subplots_adjust(left=0.12, right=0.97, bottom=0.15, top=0.9)
        ax.text(0., 1.05, "CMS Work in progress", fontsize=20, transform=ax.transAxes)
        ax.text(1., 1.05, "(13 TeV, " + self.year + ")", fontsize=20, transform=ax.transAxes, ha='right')

        ax.set_xlim(40,4e2) 
        ax.set_ylim(13,30) 
        ax.set_xscale("log")
        ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
        
        ax.set_xlabel("Mean $M^{ll}$ [GeV]", fontsize=20, ha='right', x=1.0)
        ax.set_ylabel("Mean $p_{T}^{ll}$ [GeV]", fontsize=20, ha='right', y=1.0)
        
        for index, name in enumerate(list_to_plot) :  
            if name == "Data" :
                isData=True
                temp_mass_df=self.Data["Mass"]["total"]["refinedUnc_meanValue"]
                temp_pt_df=self.combinedPtDataFrame("Data")
            if name == "DataBkgSub" :
                temp_mass_df=self.DataBkgSubtracted["Mass"]["total"]["refinedUnc_meanValue"]
                temp_pt_df=self.combinedPtDataFrame("DataBkgSub")
            if name == "Signal" :
                temp_mass_df=self.MCSignal["Mass"]["total"]["refinedUnc_meanValue"]
                temp_pt_df=self.combinedPtDataFrame("Signal")
            if name == "Background" :
                temp_mass_df=self.MCBackground["Mass"]["total"]["refinedUnc_meanValue"]
                temp_pt_df=self.combinedPtDataFrame("Background")
            if name == "TotalMC" :
                temp_mass_df=self.MCTotal["Mass"]["total"]["refinedUnc_meanValue"]
                temp_pt_df=self.combinedPtDataFrame("TotalMC")
        
            mass_systematic=self.makeErrorNumpy(temp_mass_df.total_Up, temp_mass_df.total_Down)
            pt_systematic=self.makeErrorNumpy(temp_pt_df.total_Up, temp_pt_df.total_Down)
            
            color_=next(color)
            
            if "Data" in name :
                color_ = 'black'

            ax.errorbar(temp_mass_df["mean"], temp_pt_df["mean"], xerr=mass_systematic, yerr=pt_systematic, fmt='o', color=color_, label=name)
            if do_linear_fit :
                if isData :
                    self.doLogLinearFit(ax, temp_mass_df, temp_pt_df, mass_systematic, pt_systematic, color_)
            isData=False
            
        ax.legend(loc='best', fontsize=15, fancybox=False, framealpha=0.0)
        
        plt.savefig(self.outDirPath+self.plotPrefix+"ISR.pdf", format="pdf", dpi=300)
        plt.close(fig)

    def drawISRPlots(self, *objects_to_plot, names_in_objects, do_linear_fit=None, labels=None, markers=None, colors=None, facecolor_list = None, ymin=13, ymax=30) :

        print("draw isr plot, do_linear_fit {}".format(do_linear_fit))
        color=iter(cm.rainbow(np.linspace(0,1,len(names_in_objects))))
        isData=False
        
        fig, ax = plt.subplots(figsize=(10, 6))
        plt.subplots_adjust(left=0.12, right=0.97, bottom=0.15, top=0.9)
        ax.text(0., 1.05, "CMS Work in progress", fontsize=20, transform=ax.transAxes)
        ax.text(1., 1.05, "(13 TeV, " + self.year + ")", fontsize=20, transform=ax.transAxes, ha='right')

        ax.set_xlim(30,4e2) 
        ax.set_ylim(ymin,ymax) 
        ax.set_xscale("log")
        ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
        
        ax.set_xlabel("Mean $M^{ll}$ [GeV]", fontsize=20, ha='right', x=1.0)
        ax.set_ylabel("Mean $p_{T}^{ll}$ [GeV]", fontsize=20, ha='right', y=1.0)
        
        for index, name in enumerate(names_in_objects) :  
            
            isData=False
            if name == "Data" :
                isData=True
                
            if index==0 :
                temp_mass_df=self.getDict(name)["Mass"]["total"]["refinedUnc_meanValue"]
                temp_pt_df=self.combinedPtDataFrame(name)
            else :
                temp_mass_df=objects_to_plot[index-1].getDict(name)["Mass"]["total"]["refinedUnc_meanValue"]
                temp_pt_df=objects_to_plot[index-1].combinedPtDataFrame(name)
            
            mass_systematic=self.makeErrorNumpy(temp_mass_df.total_Up, temp_mass_df.total_Down)
            pt_systematic=self.makeErrorNumpy(temp_pt_df.total_Up, temp_pt_df.total_Down)
            
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
                
            ax.errorbar(temp_mass_df["mean"], temp_pt_df["mean"], xerr=mass_systematic, yerr=pt_systematic, fmt=current_marker, color=current_color, mfc=current_facecolor, label=label_name)
            
            if do_linear_fit[index-1] :
                if isData :
                    self.doLogLinearFit(ax, temp_mass_df, temp_pt_df, mass_systematic, pt_systematic,color_)
            #self.doLogLinearFit(ax, temp_mass_df, temp_pt_df, mass_systematic, pt_systematic)
            isData=False
            
        ax.legend(loc='best', fontsize=15, fancybox=False, framealpha=0.0)
        
        plt.savefig(self.outDirPath+self.plotPrefix+"_ISR_comparison.pdf", format="pdf", dpi=300)
        plt.close(fig)
    
    def make_error_boxes(self, ax, xdata, ydata, xerror, yerror, 
                         showBox=True, showBar=False, facecolor='red', edgecolor='None', alpha=0.5, zorder=5, hatch_style=None):
    
        # Loop over data points; create box from errors at each point
        errorboxes = [Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum())
                      for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T)]
    
        # Create patch collection with specified colour/alpha
        if showBox :
            pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                                 edgecolor=edgecolor, zorder = zorder, linewidth=0.05, hatch=hatch_style)
            # Add collection to axes
            ax.add_collection(pc)

        if showBar :
            # Plot errorbars
            artists = ax.errorbar(xdata, ydata, xerr=xerror, yerr=yerror,
                                  fmt='None', ecolor='k')
    
        #return artists
    
    def makeErrorNumpy(self, Up, Down) :

        Up=Up.values
        Up=Up.reshape(1,len(Up))
        Down=Down.values
        Down=Down.reshape(1,len(Down))
        
        UpDown=Up
        UpDown=np.append(UpDown, Down, axis = 0)
        
        return UpDown
        
    def getRawDict(self, dictName) :
        if dictName == "Data" :
            return self.rawDataMeasured
        if dictName == "DataBkgSub" :
            return self.rawDataBkgSubtracted
        if dictName == "Signal" :
            return self.rawMCSignal
        if dictName == "Background" :
            return self.rawMCBackground
        if dictName == "TotalMC" :
            return self.rawMCTotal
        
    def getDict(self, dictName) :
        if dictName == "Data" :
            return self.Data
        if dictName == "DataBkgSub" :
            return self.DataBkgSubtracted
        if dictName == "Signal" :
            return self.MCSignal
        if dictName == "Background" :
            return self.MCBackground
        if dictName == "TotalMC" :
            return self.MCTotal
