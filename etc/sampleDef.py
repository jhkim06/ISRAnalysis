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

myinputDir = '/home/jhkim/data/Data/temp/'
myinputDir_ = '/home/jhkim/data/Data/DY_MG/'

ISR2016_electron = {
    'DATA'     : isrSample('DoubleEGamma',  myinputDir + 'ISRee_unfolding_data_DoubleEG_cat_v8-0-7.root', isMC = False, isSig = False),
    'DY'       : isrSample('DYtoEE' ,       myinputDir + 'DYJets_cat_v8-0-7_1.root', isMC = True, isSig = True, isInc = True),
    'DY10to50' : isrSample('DYtoEE10to50' , myinputDir + 'ISRee_unfolding_DYJets_10to50_cat_v8-0-7.root', isMC = True, isSig = True, isInc = True),
    'DYMG'       : isrSample('DYMGtoEE' ,       myinputDir_ + 'DYJets_MG_cat_v8-0-7_1.root', isMC = True, isSig = True, isInc = True),
    'DYMG10to50' : isrSample('DYMGtoEE10to50' , myinputDir_ + 'ISRee_unfolding_DYJets_MG_10to50_cat_v8-0-7.root', isMC = True, isSig = True, isInc = True),
    'BKG'      : isrSample('BKG' ,          myinputDir + 'ISRee_unfolding_background.root', isMC = True, isSig = False),
    'TTbar'      : isrSample('TTbar' ,          myinputDir + 'ISRee_unfolding_SKTT_powheg_dilep_cat_v8-0-7.root', isMC = True, isSig = False),
    'VV'      : isrSample('VV' ,          myinputDir + 'ISRee_unfolding_VV.root', isMC = True, isSig = False),
    'Wjets'      : isrSample('Wjets' ,          myinputDir + 'ISRee_unfolding_SKWJets_dilep_cat_v8-0-7.root', isMC = True, isSig = False),
}

# for 50 DY ntuples
for i in range(2,51):
	ISR2016_electron['DY%s' % str(i)] = isrSample('DYtoEE_%s' % str(i) , myinputDir + 'DYJets_cat_v8-0-7_%s.root' % str(i), isMC = True, isSig = True, isInc = True)

for i in range(2,36):
        ISR2016_electron['DYMG%s' % str(i)] = isrSample('DYMGtoEE_%s' % str(i) , myinputDir_ + 'DYJets_MG_cat_v8-0-7_%s.root' % str(i), isMC = True, isSig = True, isInc = True)

samplesDef_electron = {
    'data'   : ISR2016_electron['DATA'].clone(),
    'mcSig'  : ISR2016_electron['DY'].clone(),
    #'mcSig'  : ISR2016_electron['DYMG'].clone(),
    #'mcBkg'  : ISR2016_electron['BKG'].clone(),
    'mcBkg1'  : ISR2016_electron['TTbar'].clone(),
    'mcBkg2'  : ISR2016_electron['VV'].clone(),
    'mcBkg3'  : ISR2016_electron['Wjets'].clone(),
}

for i in range(2,51):
	samplesDef_electron['mcSig'].add_sample(ISR2016_electron['DY%s' % str(i)])

samplesDef_electron['mcSig'].add_sample(ISR2016_electron['DY10to50'])

#for i in range(2,36):
#        samplesDef_electron['mcSig'].add_sample(ISR2016_electron['DYMG%s' % str(i)])
#
#samplesDef_electron['mcSig'].add_sample(ISR2016_electron['DYMG10to50'])

