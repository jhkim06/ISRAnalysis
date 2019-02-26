import copy
import os

class tnpSample:
	def __init__(self, sName, path, isMC = False, isSig = False):
        	self.path = []
                self.path.append(path)
        	self.name = sName
		self.isMC = isMC
 		self.isSig = isSig

    	def dump(self):
    	    print '**** name: %-*s ' % (100, self.name)
    	    print '  path    : ', self.path
    	    print '  is MC    : ', self.isMC
    	    print '  is Signal    : ', self.isSig

    	def add_sample(self, sample):
        	self.path.extend( sample.path )

    	def clone(self):
        	return copy.deepcopy(self)

myinputDir = '/home/jhkim/data/Data/temp/'


ISR2016_electron = {
    'DATA' : tnpSample('DoubleEGamma', myinputDir + 'ISRee_unfolding_data_DoubleEG_cat_v8-0-7.root', isMC = False, isSig = False),
    'DY' : tnpSample('DYtoEE' , myinputDir + 'DYJets_cat_v8-0-7_10.root', isMC = True, isSig = True),
    'DY1' : tnpSample('DYtoEE1' , myinputDir + 'DYJets_cat_v8-0-7_11.root', isMC = True, isSig = True),
}

samplesDef_electron = {
    'data'   : ISR2016_electron['DATA'].clone(),
    'mcSig'  : ISR2016_electron['DY'].clone(),
}

#FIXME find better way to add files
samplesDef_electron['mcSig'].add_sample(ISR2016_electron['DY1'])
