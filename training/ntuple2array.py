import numpy as np
from numpy.lib.recfunctions import stack_arrays
from ROOT import *
from root_numpy import tree2array
import glob
import pandas as pd
import deepdish.io as io

qcd23 = TFile.Open('/home/sarakm0704/WORK/ntuples/eta24/IsoNtuple_QCD-20to30.root')
qcd23_tree = qcd23.Get('demo/tree')
qcd23_array = tree2array(qcd23_tree)
qcd23_df = pd.DataFrame(qcd23_array)
io.save('eta24_hdf/qcd-2030.h5',qcd23_df)

qcd35 = TFile.Open('/home/sarakm0704/WORK/ntuples/eta24/IsoNtuple_QCD-30to50.root')
qcd35_tree = qcd35.Get('demo/tree')
qcd35_array = tree2array(qcd35_tree)
qcd35_df = pd.DataFrame(qcd35_array)
io.save('eta24_hdf/qcd-3050.h5',qcd35_df)

qcd58 = TFile.Open('/home/sarakm0704/WORK/ntuples/eta24/IsoNtuple_QCD-50to80.root')
qcd58_tree = qcd58.Get('demo/tree')
qcd58_array = tree2array(qcd58_tree)
qcd58_df = pd.DataFrame(qcd58_array)
io.save('eta24_hdf/qcd-5080.h5',qcd58_df)

qcd812 = TFile.Open('/home/sarakm0704/WORK/ntuples/eta24/IsoNtuple_QCD-80to120.root')
qcd812_tree = qcd812.Get('demo/tree')
qcd812_array = tree2array(qcd812_tree)
qcd812_df = pd.DataFrame(qcd812_array)
io.save('eta24_hdf/qcd-80120.h5',qcd812_df)

dy = TFile.Open('/home/sarakm0704/WORK/ntuples/eta24/IsoNtuple_DY.root')
dy_tree = dy.Get('demo/tree')
dy_array = tree2array(dy_tree)
dy_df = pd.DataFrame(dy_array)
io.save('eta24_hdf/iso_dy.h5', dy_df)

"""
frames = [ttbb_df, ttbj_df, ttcc_df, ttLF_df, ttother_df]
result = pd.concat(frames, ignore_index=True)
io.save('ttbarJetCombinations.h5', result)
"""
