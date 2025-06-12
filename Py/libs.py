# =============================================================================
# MODULES AND ENVIRONMENT VARIABLES
# =============================================================================

# Python built-in
import os
import sys
import time
import logging
import json
import pickle
import re
import math
import importlib
import random
from datetime import datetime
from pprint import pformat
from os.path import exists

# Data wrangling
import pandas as pd
pd.options.mode.chained_assignment = None
pd.set_option('display.max_rows', 10)
import pyarrow.feather as ft
import numpy as np

# Predictions
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import pearsonr
import tensorflow as tf
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras import Model
from tensorflow.keras.layers import Dense, Dropout, Conv1D, add, Flatten, LSTM, AveragePooling1D, concatenate, MaxPool1D, BatchNormalization, Reshape, Activation
from tensorflow.keras.callbacks import EarlyStopping, Callback, ModelCheckpoint, ReduceLROnPlateau, TensorBoard
from tensorflow.keras.regularizers import l1, l2, L1L2
from tensorflow.keras.optimizers import Adam
from tensorflow.keras import backend as K
from tensorflow.keras import Input
import keras_tuner as kt

# Utility functions - implement missing functions as needed
def fit_model(final_model, params, train_x, val_x, train_y, val_y):
     
    # set variables
    tb_filepath, cp_filepath, b_size, epoch, vbs, sfl = [params['fit'][key] for key in ['tensorboard_fp', 'checkpoint_fp', 'batch_size', 'epochs', 'verbose', 'shuffle']]
    
    #set call backs
    tensorboard_cb = TensorBoard(tb_filepath)
    modelcheck_cb = ModelCheckpoint(filepath=cp_filepath,
                                    save_weights_only=True,
                                    monitor='val_loss',
                                    mode='min',
                                    save_best_only=True)
    model_cb = EarlyStopping(monitor='val_loss',
                                     min_delta=0.00001,
                                     patience=5,
                                     verbose=0,
                                     mode='min',
                                     baseline=None,
                                     restore_best_weights=True)
    final_model.fit(train_x, train_y, validation_data=(val_x, val_y),
                    batch_size = b_size,
                    epochs = epoch,
                    verbose = vbs,
                    shuffle = sfl,
                    callbacks=[modelcheck_cb, 
                               tensorboard_cb,
                               model_cb])
    
    final_model.load_weights(cp_filepath) # loads best weights
    return final_model

def predict_values (model, test_x, test_y, index, scaler):
    
    # perform predictions
    prediction = model.predict(test_x)
    
    # re-scale data
    obs = inverse_scale(scaler, test_y, verbose = False)
    pred = inverse_scale(scaler, prediction, verbose = False)
    out_data = pd.DataFrame([index, obs, pred], index=["index","obs","pred"]).T
    out_data["index"] = out_data["index"].astype('int')
    return out_data