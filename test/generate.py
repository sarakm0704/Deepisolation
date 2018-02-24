import csv
import numpy as np

qcd_file='qcd_iso.txt'
qcd_file_obj=open(qcd_file,'r')
qcd_csv=csv.reader(qcd_file_obj,delimiter=' ')

qcd_iso = []
for row in qcd_csv:
  row=np.array(filter(None,row)).astype(float)
  qcd_iso.append(row) 

np.save('qcd_iso', qcd_iso)

