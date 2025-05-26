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

#def tuner_expand(hp_out):
#    markers_no = 9796
#    
#    #### Define model
#    model = Sequential()
#    
#    ### First layer
#    model.add(Conv1D(filters = hp_out["CNN_f_fl"], kernel_size = hp_out["CNN_ks_fl"], 
#                     padding='valid', activation='relu', input_shape=(markers_no, 1), 
#                     name="CNN_fl"))
#    model.add(AveragePooling1D(pool_size =  hp_out["CNN_ap_fl"], strides = 3, padding='same', 
#                               name="CNN_ap_fl"))
#    
#    ### Variable layers
#    for i in range(hp_out["CNN_num_vl"]):
#        model.add(Conv1D(filters = hp_out[f'CNN_f_vl_{i}'], kernel_size = hp_out[f'CNN_ks_vl_{i}'], padding='valid', activation='relu', name=f'CNN_Conv_{i}'))
#        model.add(AveragePooling1D(pool_size = hp_out[f'CNN_ap_vl_{i}'], strides = 3, padding='same', name=f'CNN_ap_{i}'))
#    
#    ### Flattening layer
#    model.add(Flatten(name="CNN_flatten"))
#    
#    ### Dense layers
#    for i in range(hp_out["CNN_num_dl"]):
#        model.add(Dense(units = hp_out[f'CNN_unit_dl_{i}'], activation = 'relu', name = f'CNN_dl_{i}'))
#        model.add(Dropout(rate = hp_out[f'CNN_drop_rate_dl_{i}'], name=f'CNN_drop_dl_{i}'))
#    ### Final layer
#    model.add(Dense(1, activation='tanh',name="CNN_out"))
#    
#    ### hyperparameters
#    l_rate = hp_out["l_rate"]
#    beta_val_1 = hp_out["beta_val_1"]
#    beta_val_2 = hp_out["beta_val_2"]
#    
#    ### Complie model
#    model.compile(loss = 'mean_absolute_error', 
#                  optimizer = Adam(learning_rate = l_rate, 
#                                   beta_1 = beta_val_1, 
#                                   beta_2 = beta_val_2),
#                  metrics = ['mean_squared_error'])
#    return(model)

class tuner_2(kt.HyperModel):
    def __init__(self, marker_n):
        self.marker_n = marker_n # markers
        print(f'markers = {self.marker_n}')
    def build(self, hp):
        # generate network for g_a data
        X_list = CNN_net_flex(hp, n_in = self.marker_n, 
                          model_name = "g_a", base_layer = 7, max_layer = 3, tuning = self.tune)
        #net_g_d = CNN_g_data(hp, self.marker_n, "g_d")

        # split list for input and output
        input_list = [x[0] for x in X_list]
        output_list = [x[1] for x in X_list]
        
        # concat inputs
        model_concat = concatenate(output_list, name = "concat_in")
        
        model_compiled = dense_post_concat(hp, input_list, model_concat, tuning = not self.tune)
        
        return model_compiled