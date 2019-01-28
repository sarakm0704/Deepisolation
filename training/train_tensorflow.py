import numpy as np
import csv
from sklearn.utils import shuffle
import os
import pandas as pd
from ROOT import TFile, TTree
from root_numpy import array2tree, tree2array

os.environ["CUDA_VISIBLE_DEVICES"] = "0"

def weight_variable(shape):
  initial = tf.truncated_normal(shape, stddev=0.1)
  return tf.Variable(initial)

def bias_variable(shape):
  initial = tf.constant(0.1, shape=shape)
  return tf.Variable(initial)

input_variables=[]

p=pd.read_hdf("h5materials/iso_dy.h5").astype(np.float32)
n=pd.read_hdf("h5cut/qcd-20120_nocut.h5").astype(np.float32)

input_variables.extend(['pt', 'eta', 'phi', 'chIso0_005', 'chIso005_01', 'chIso01_015', 'chIso015_02', 'chIso02_025', 'chIso025_03', 'nhIso0_005', 'nhIso005_01', 'nhIso01_015', 'nhIso015_02', 'nhIso02_025', 'nhIso025_03', 'phIso0_005', 'phIso005_01', 'phIso01_015', 'phIso015_02', 'phIso02_025', 'phIso025_03', 'puiso0_005', 'puiso005_01', 'puiso01_015', 'puiso015_02', 'puiso02_025', 'puiso025_03', 'abIso', 'relIso', 'chN', 'nhN', 'phN', 'ch_dR', 'nh_dR', 'ph_dR', 'ch_dPhi', 'nh_dPhi', 'ph_dPhi', 'ch_dEta', 'nh_dEta', 'ph_dEta'])

p=p.filter(input_variables, axis=1).values
n=n.filter(input_variables, axis=1).values

train_p=p[:2800190] #pt 20to120
#train_p=p[:96023] #pt 20to30
#train_p=p[:504814] #pt 30to50
#train_p=p[:1032291] #pt 50to80
#train_p=p[:1167062] #pt 80to120
train_n=n

train_data=np.vstack((train_p,train_n))
train_out=np.array([1]*len(train_p)+[0]*len(train_n))

#for evaluation
#train_data=p
#train_out=np.array([1]*len(p))
#train_n=n
#train_out=np.array([0]*len(n))

numbertr=len(train_out)

#Shuffling
order=shuffle(range(numbertr),random_state=100)
train_out=train_out[order]
train_data=train_data[order,0::]

train_out = train_out.reshape( (numbertr, 1) )
trainnb=0.9 # Fraction used for training

#Splitting between training set and cross-validation set
valid_data=train_data[int(trainnb*numbertr):numbertr,0::]
valid_data_out=train_out[int(trainnb*numbertr):numbertr]

train_data_out=train_out[0:int(trainnb*numbertr)]
train_data=train_data[0:int(trainnb*numbertr),0::]

import tensorflow as tf
sess = tf.InteractiveSession()

x = tf.placeholder(tf.float32, shape=[None, 41]) #20
y_ = tf.placeholder(tf.float32, shape=[None, 1])

##### Model #####
a = 300
W1 = weight_variable( [41,a] )
b1 = bias_variable( [a] )
A1 = tf.nn.relu(tf.matmul(x, W1) + b1)
W2 = weight_variable( [a,a] )
b2 = bias_variable( [a] )
A2 = tf.nn.relu(tf.matmul(A1, W2) + b2)
W3 = weight_variable( [a,a] )
b3 = bias_variable( [a] )
A3 = tf.nn.relu(tf.matmul(A2, W3) + b3)
W4 = weight_variable( [a,a] )
b4 = bias_variable( [a] )
A4 = tf.nn.relu(tf.matmul(A3, W4) + b4)
W5 = weight_variable( [a,a] )
b5 = bias_variable( [a] )
A5 = tf.nn.relu(tf.matmul(A4, W5) + b5)
W6 = weight_variable( [a,a] )
b6 = bias_variable( [a] )
A6 = tf.nn.relu(tf.matmul(A5, W6) + b6)
W7 = weight_variable( [a,a] )
b7 = bias_variable( [a] )
A7 = tf.nn.relu(tf.matmul(A6, W7) + b7)
W8 = weight_variable( [a,a] )
b8 = bias_variable( [a] )
A8 = tf.nn.relu(tf.matmul(A7, W8) + b8)
W9 = weight_variable( [a,a] )
b9 = bias_variable( [a] )
A9 = tf.nn.relu(tf.matmul(A8, W9) + b9)
W10 = weight_variable( [a,a] )
b10 = bias_variable( [a] )
A10 = tf.nn.relu(tf.matmul(A9, W10) + b10)
W11 = weight_variable( [a,1] )
b11 = bias_variable( [1] )
y = tf.matmul(A10,W11) + b11
##################

cross_entropy = tf.reduce_mean(
    tf.nn.sigmoid_cross_entropy_with_logits(labels=y_, logits=y))

train_step = tf.train.AdamOptimizer(1e-4).minimize(cross_entropy)

ntrain = len(train_data)

#print ("ntrain = ", ntrain, " QCD = ", len(train_n) ," DY = ", len(train_p))

batch_size = 500
epoch_size = ntrain // batch_size
cur_id = 0

saver = tf.train.Saver()

model_output_name = "layer10"

loss_list = []

with tf.Session() as sess:
  sess.run(tf.global_variables_initializer())
  if os.path.exists('models/'+model_output_name+'/model_out.meta'):
      print ("Model file exists already!")
      saver.restore(sess, 'models/'+model_output_name+'/model_out')
  else:
    for i in range(20):
      for j in range(epoch_size):
        if j == epoch_size-1: cur_id = 0
        batch_data = train_data[cur_id:cur_id+batch_size]
        batch_data_out = train_data_out[cur_id:cur_id+batch_size]
        cur_id = cur_id+batch_size
        _, loss = sess.run([train_step, cross_entropy], feed_dict={x: batch_data, y_: batch_data_out})
        loss_list.append(loss)

    saver.save(sess,  'models/'+model_output_name+'/model_out')
    print ("Model saved!")

  prediction = tf.nn.sigmoid(y) 
  pred = prediction.eval( feed_dict={x: valid_data, y_: valid_data_out} )
  print (pred)  

  y = []
  signal_output = []
  signal_train_output = []
  background_output = []
  background_train_output = []

  #save output for later use
  with open('models/'+model_output_name+'/signal.csv', 'w') as f:
    writer = csv.writer(f, delimiter=" ")
    for i in range(len(valid_data)):
      x = valid_data_out[i] #valdiation label
      if x == 1: #signal output 
        signal_output.append( pred[i] ) 
        y = pred[i]
        writer.writerows( zip(y,x) )       

  with open('models/'+model_output_name+'/background.csv', 'w') as f:
    writer = csv.writer(f, delimiter=" ")
    for i in range(len(valid_data)):
      x = valid_data_out[i] #validation level
      if x == 0: #background output
        background_output.append( pred[i] ) 
        y = pred[i]
        writer.writerows( zip(y,x) )

"""
  outfile = TFile.Open(os.path.join('score_roots', 'score_dy.root'),'RECREATE')
  outtree = TTree("tree","tree")
  spectator = train_p.filter(['pt','eta'], axis=1)
  train_p.astype(np.float32)

  y = pred
  y.dtype = [('MLScore', np.float32)]
  array2tree(y, name='tree', tree=outtree)

  for colname, value in spectator.iteritems():
    spect = spectator[colname].values
    if colname == 'pt': branchname = 'muon_pt'
    elif colname == 'eta'    : branchname = 'muon_eta'
    else: branchname = colname

    spect.dtype = [(branchname, np.float32)]
    array2tree(spect, name='tree', tree=outtree)

  outtree.Fill()
  outfile.Write()
  outfile.Close()

  #convert to numpy array
  s_output = np.array(signal_output)
  b_output = np.array(background_output)
  threshold = 0.5 #deletey
  s = s_output > threshold
  b = b_output > threshold

  #signal count
  ns_sel = len(s_output[s]) # count only elements larger than threshold
  ns_total = len(signal_output) 
  
  #background count 
  nb_sel = len(b_output[b]) # count only elements larger than threshold
  nb_total = len(background_output)

  print ("signal : " , ns_sel ,  "/" , ns_total)
  print ("background : ", nb_sel ,  "/" , nb_total)
 
  #efficiency
  sig_eff = float(ns_sel)/float(ns_total) 
  bkg_eff = float(nb_sel)/float(nb_total)
 
  print ("signal eff = ", sig_eff, " background eff = ", bkg_eff)

"""
