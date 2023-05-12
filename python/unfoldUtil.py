import os
import sys
import ROOT as rt

# Main unfolding library
current_dir = os.getcwd()
#os.environ['ISR_UNFOLD_WD']
rt.gSystem.Load(current_dir + "/lib/libisrunfold.so")
