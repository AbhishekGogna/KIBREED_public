from .libs import *

# Model 1 ---------------------------------------------------------------------------------------------------------------------
def tuner(hp):
    markers_no = 9797
    
    #### Define model
    model = Sequential()
    
    ### First layer
    model.add(Conv1D(filters = hp.Int('CNN_f_fl', min_value=64, max_value=512, step = 64, default = 512), 
                     kernel_size= hp.Int('CNN_ks_fl', min_value=3, max_value=36, step = 3), 
                     padding='valid', activation='relu', input_shape=(markers_no, 1), 
                     name="CNN_fl"))
    model.add(AveragePooling1D(pool_size =  hp.Int("CNN_ap_fl", min_value=2, max_value=32, default=16, step = 4), 
                               strides = 3, padding='same', 
                               name="CNN_ap_fl"))
    ### Variable layers
    for i in range(hp.Int("CNN_num_vl", min_value = 2, max_value = 4)):
        model.add(Conv1D(filters = hp.Int(f'CNN_f_vl_{i}', min_value=64, max_value=512, step = 32, default = 256), 
                         kernel_size = hp.Int(f'CNN_ks_vl_{i}', min_value=3, max_value=36, step = 3), 
                         padding='valid', activation='relu', 
                         name=f'CNN_Conv_{i}'))
        model.add(AveragePooling1D(pool_size = hp.Int(f'CNN_ap_vl_{i}', min_value=2, max_value=32, default=16, step = 4), 
                                   strides = 3, padding='same', 
                                   name=f'CNN_ap_{i}')) # pool size and strides
    
    ### Flattening layer
    model.add(Flatten(name="CNN_flatten"))
              
    ### Dense layers
    for i in range(hp.Int("CNN_num_dl", min_value = 1, max_value = 4)):
        model.add(Dense(units = hp.Int(f'CNN_unit_dl_{i}', min_value=32, max_value=256, step = 32, default = 128), 
                        activation = 'relu', 
                        name = f'CNN_dl_{i}'))
        model.add(Dropout(rate = hp.Float(f'CNN_drop_rate_dl_{i}', min_value=0.1, max_value=0.5, step = 0.01),
                          name=f'CNN_drop_dl_{i}'))
    ### Final layer
    model.add(Dense(1, activation='tanh',name="CNN_out"))
    
    ### hyperparameters
    l_rate = hp.Float("l_rate", min_value=1e-5, max_value=1e-2, sampling="log")
    beta_val_1 = hp.Float("beta_val_1", min_value=0, max_value=1)
    beta_val_2 = hp.Float("beta_val_2", min_value=0, max_value=1)
    
    ### Complie model
    model.compile(loss = 'mean_absolute_error', 
                  optimizer = Adam(learning_rate = l_rate, 
                                   beta_1 = beta_val_1, 
                                   beta_2 = beta_val_2),
                  metrics = ['mean_squared_error'])
    return(model)