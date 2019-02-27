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
    'DY' : tnpSample('DYtoEE' , myinputDir + 'DYJets_cat_v8-0-7_1.root', isMC = True, isSig = True),
    'DY10to50' : tnpSample('DYtoEE10to50' , myinputDir + 'ISRee_unfolding_DYJets_10to50_cat_v8-0-7.root', isMC = True, isSig = True),
    'BKG' : tnpSample('BKG' , myinputDir + 'ISRee_unfolding_background.root', isMC = True, isSig = False),
}

# for 50 DY ntuples
for i in range(2,51):
	ISR2016_electron['DY%s' % str(i)] = tnpSample('DYtoEE_%s' % str(i) , myinputDir + 'DYJets_cat_v8-0-7_%s.root' % str(i), isMC = True, isSig = True)

samplesDef_electron = {
    #'data'   : ISR2016_electron['DATA'].clone(),
    'mcSig'  : ISR2016_electron['DY'].clone(),
    #'mcBkg'  : ISR2016_electron['BKG'].clone(),
}

for i in range(2,51):
	samplesDef_electron['mcSig'].add_sample(ISR2016_electron['DY%s' % str(i)])

samplesDef_electron['mcSig'].add_sample(ISR2016_electron['DY10to50'])

