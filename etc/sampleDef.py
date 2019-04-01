import copy
import os

class isrSample:
	def __init__(self, sName, path, isMC = False, isSig = False, isInc = False): # isInc is true if the signal sample contains tau events
        	self.path = []
                self.path.append(path)
        	self.name = sName
		self.isMC = isMC
 		self.isSig = isSig
 		self.isInc = isInc

    	def dump(self):
    	    print '**** name: %-*s ' % (100, self.name)
    	    print '  path    : ', self.path
    	    print '  is MC    : ', self.isMC
    	    print '  is Signal    : ', self.isSig

    	def add_sample(self, sample):
        	self.path.extend( sample.path )

    	def clone(self):
        	return copy.deepcopy(self)

myinputDir_2016Legacy = '/home/jhkim/data/Data/Legacy/2016/'
myinputDir_2017Legacy = '/home/jhkim/data/Data/Legacy/2017/'

ISR2016Legacy = {
    'DATA_electron'     : isrSample('DoubleEGamma',  myinputDir_2016Legacy + 'DoubleEG_All.root ', isMC = False, isSig = False),
    'DATA_muon'     : isrSample('DoubleMuon',  myinputDir_2016Legacy + 'DoubleMuon_All.root', isMC = False, isSig = False),
    'DY'       : isrSample('DYtoEE' ,       myinputDir_2016Legacy + 'DYJetsToLL_M-50.root', isMC = True, isSig = True, isInc = True),
    'DY10to50' : isrSample('DYtoEE10to50' , myinputDir_2016Legacy + 'DYJetsToLL_M-10to50.root', isMC = True, isSig = True, isInc = True),
    'TTbar'      : isrSample('TTbar' ,          myinputDir_2016Legacy + 'TT.root', isMC = True, isSig = False),
    'VV'      : isrSample('VV' ,          myinputDir_2016Legacy + 'VV.root', isMC = True, isSig = False),
    'Wjets'      : isrSample('Wjets' ,          myinputDir_2016Legacy + 'Wjets.root', isMC = True, isSig = False),
}

ISR2017Legacy = {
    'DATA_electron'     : isrSample('DoubleEGamma',  myinputDir_2017Legacy + 'DoubleEG_All.root ', isMC = False, isSig = False),
    'DATA_muon'     : isrSample('DoubleMuon',  myinputDir_2017Legacy + 'DoubleMuon_All.root', isMC = False, isSig = False),
    'DY'       : isrSample('DYtoEE' ,       myinputDir_2017Legacy + 'DYJetsToLL_M-50.root', isMC = True, isSig = True, isInc = True),
    'DY10to50' : isrSample('DYtoEE10to50' , myinputDir_2017Legacy + 'DYJetsToLL_M-10to50.root', isMC = True, isSig = True, isInc = True),
    'TTbar'      : isrSample('TTbar' ,          myinputDir_2017Legacy + 'TT.root', isMC = True, isSig = False),
    'VV'      : isrSample('VV' ,          myinputDir_2017Legacy + 'VV.root', isMC = True, isSig = False),
    'Wjets'      : isrSample('Wjets' ,          myinputDir_2017Legacy + 'Wjets.root', isMC = True, isSig = False),
}


# electron channel
samplesDef_electron2016Legacy = { 
    'data'   : ISR2016Legacy['DATA_electron'].clone(),
    'mcSig'  : ISR2016Legacy['DY'].clone(),
    'mcBkg1'  : ISR2016Legacy['TTbar'].clone(),
    'mcBkg2'  : ISR2016Legacy['VV'].clone(),
    'mcBkg3'  : ISR2016Legacy['Wjets'].clone(),
}

samplesDef_electron2016Legacy['mcSig'].add_sample(ISR2016Legacy['DY10to50'])

# muon channel
samplesDef_muon2016Legacy = {
    'data'   : ISR2016Legacy['DATA_muon'].clone(),
    'mcSig'  : ISR2016Legacy['DY'].clone(),
    'mcBkg1'  : ISR2016Legacy['TTbar'].clone(),
    'mcBkg2'  : ISR2016Legacy['VV'].clone(),
    'mcBkg3'  : ISR2016Legacy['Wjets'].clone(),
}

samplesDef_muon2016Legacy['mcSig'].add_sample(ISR2016Legacy['DY10to50'])


# electron channel
samplesDef_electron2017Legacy = { 
    'data'   : ISR2017Legacy['DATA_electron'].clone(),
    'mcSig'  : ISR2017Legacy['DY'].clone(),
    'mcBkg1'  : ISR2017Legacy['TTbar'].clone(),
    'mcBkg2'  : ISR2017Legacy['VV'].clone(),
    'mcBkg3'  : ISR2017Legacy['Wjets'].clone(),
}

samplesDef_electron2017Legacy['mcSig'].add_sample(ISR2017Legacy['DY10to50'])

# muon channel
samplesDef_muon2017Legacy = {
    'data'   : ISR2017Legacy['DATA_muon'].clone(),
    'mcSig'  : ISR2017Legacy['DY'].clone(),
    'mcBkg1'  : ISR2017Legacy['TTbar'].clone(),
    'mcBkg2'  : ISR2017Legacy['VV'].clone(),
    'mcBkg3'  : ISR2017Legacy['Wjets'].clone(),
}

samplesDef_muon2017Legacy['mcSig'].add_sample(ISR2017Legacy['DY10to50'])
