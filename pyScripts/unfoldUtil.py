import os
import sys
import ROOT as rt

# Main unfolding library
current_dir = os.getcwd()
rt.gSystem.Load(current_dir + "/lib/libisrunfold.so")

import gc
gc.collect()

