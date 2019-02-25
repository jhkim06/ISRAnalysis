import copy
import os

class tnpSample:
	def __init__(self, sName, path):
        	self.path = []
        	self.name = sName

    	def add_sample(self, sample):
        	self.path.extend( sample.path )

    	def clone(self):
        	return copy.deepcopy(self)

myinputDir = '/home/jhkim/data/Data/temp/'

ISR2016_electron = {
    'DATA' : tnpSample('DoubleEGamma', myinputDir + 'ISRee_unfolding_data_DoubleEG_cat_v8-0-7.root'),
    'DY50plus' : tnpSample('DY50plus_aMC@NLO' , myinputDir + 'DYJets_cat_v8-0-7_10.root'),
}

samplesDef_electron = {
    'data'   : ISR2016_electron['DATA'].clone(),
    'mcSig'  : ISR2016_electron['DY50plus'].clone(),
}
