import functools
import itertools
import os
import random
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import csv
import matplotlib.pyplot as plt
import numpy as np
from sklearn import metrics
from sklearn import preprocessing
from sklearn.model_selection import train_test_split, KFold, cross_val_score

from keras.layers import Dense, Activation, Flatten, Dropout, Reshape
from keras.layers import Conv1D,Conv2D, MaxPooling2D
from keras.models import Sequential,Model
from keras.utils.np_utils import to_categorical
from keras import optimizers
from keras.optimizers import Adam,SGD
from keras.layers.normalization import BatchNormalization
from keras.regularizers import l2
import copy
import sys

def predict_for_deepphos(train_file_name,db,sites,predictFrame = 'general',
                         hierarchy=None, kinase=None, prefix="deepphos"):
    '''

    :param train_file_name: input of your prdict file
                            it must be a .csv file and theinput format  is proteinName, postion,sites, shortseq
    :param sites: the sites predict: site = 'S','T' OR 'Y'
    :param predictFrame: 'general' or 'kinase'
    :param hierarchy: if predictFrame is kinse: you must input the hierarchy:
            group,family,subfamily,kinase to choose corresponding model
    :param kinase: kinase name
    :return:
     a file with the score
    '''


    win1 = 51
    win2 = 33
    win3 = 15
    from methods.dataprocess_predict import getMatrixInput
    [X_test1,y_test,ids,position] = getMatrixInput(train_file_name, sites,db, win1)
    [X_test2,_,_,_] = getMatrixInput(train_file_name, sites, db, win2)
    [X_test3,_,_,_]  = getMatrixInput(train_file_name, sites, db, win3)

#     print X_test1.shape
#     print len(position)

    from methods.model_n import model_net
    model = model_net(X_test1, X_test2, X_test3, y_test,nb_epoch = 0)

    #load model weight
    if predictFrame == 'general':
        outputfile = 'general_{:s}'.format(site)
        if site == 'ST':
            model_weight = './models/model_general_S,T.h5'
        if site == 'Y':
            model_weight = './models/model_general_Y.h5'


    if predictFrame == 'kinase':
        outputfile = 'kinase_{:s}_{:s}'.format(hierarchy, kinase)
        model_weight = './models/model_{:s}_{:s}.h5'.format(hierarchy, kinase)
#     print model_weight
    model.load_weights(model_weight)
    predictions_t = model.predict([X_test1, X_test2, X_test3])
    res = dict()
    res['protein'] = ids
    res['pos'] = position
    res['prob'] = predictions_t[:, 1]
    result = pd.DataFrame(res)
    
    dat = pd.read_table(train_file_name, header=0, sep="\t")
    dat['y_pred'] = result['prob']
    
    #results_ST = np.column_stack((ids, position,predictions_t[:, 1]))

    #result = pd.DataFrame(results_ST)
    
    dat.to_csv(outputfile + "_prediction_phosphorylation_" + str(prefix) + ".tsv", index=False, sep='\t')
    
if __name__ == '__main__':
    #train_file_name = 'test data.csv'
    #site = 'S','T'
    #predict_for_deepphos(train_file_name, site, predictFrame='kinase',
    #                     hierarchy='group', kinase='AGC')
    
    train_file_name = sys.argv[1]
    db = sys.argv[2]
    out_prefix = sys.argv[3]
    
    site = 'ST'
    dat = pd.read_table(train_file_name, header=0, sep="\t")
    
    sdata = dat.query('aa=="S"')
    tdata = dat.query('aa=="T"')
    
    stdata = pd.concat([sdata,tdata],axis=0)
    stdata.to_csv("st_input_data.tsv", index=False, sep='\t')
    predict_for_deepphos("st_input_data.tsv", db, site, predictFrame='general',prefix=out_prefix)
    site = 'Y'
    ydata = dat.query('aa=="Y"')
    ydata.to_csv("y_input_data.tsv", index=False, sep='\t')
    predict_for_deepphos("y_input_data.tsv", db, site, predictFrame='general',prefix=out_prefix)




