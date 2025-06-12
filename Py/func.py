from .libs import *

def scale_data(to_transform, pd_cols, pd_index):
    scaler = MinMaxScaler((0,1))
    data_scaled = scaler.fit_transform(to_transform)
    data_scaled_df = pd.DataFrame(data_scaled, columns = pd_cols, index = pd_index)
    return data_scaled_df, scaler

def create_train_val_data(index_train, index_test, index_val = None, prop = 0.1):
    if index_val is None:
        val_set = random.sample(index_train, int(len(index_train)*prop)) # cretes validation set from the remaining non_test set 
        train_set = list(set(index_train).difference(val_set))
    else:
        val_set = index_val
        train_set = index_train
    
    check = any(item in val_set for item in train_set)
    
    test_set = index_test
    
    if check:
        print("function failed since some elemets of val arer in the train set")
    else:
        return train_set, val_set, test_set

def CNN_net_fixed(hp, n_in, model_name, kernel_size = 2):
    # define model inputs
    model_input = Input(shape = (n_in, 1), name = f'cnn_{model_name}_in')
    # define first layer of CNN
    CNN = Conv1D(filters = 128,
                 kernel_size = kernel_size,
                 padding="valid", activation="relu",
                 name = f'CNN_{model_name}_fl')(model_input)
    # define variable layers of CNN
    CNN = Conv1D(filters = 64,
                 kernel_size = kernel_size,
                 padding="valid", activation="relu",
                 name = f'CNN_{model_name}_vl')(CNN)
    # define last layer of CNN
    CNN = Conv1D(filters = 32,
                 kernel_size = kernel_size,
                 padding="valid", activation="relu",
                 name = f'CNN_{model_name}_ll')(CNN)
    
    # flatten everything
    CNN_output = Flatten(name = f"flat_layer_{model_name}")(CNN)
    return model_input, CNN_output

def tuner_obj_int(hp, name, val):
        tuner_obj = hp.Int(name = name, min_value = math.ceil(val/2), max_value = math.ceil(val*2), step = math.ceil(val/4), default = val)
        return tuner_obj

def CNN_net_flex(hp, model_name, base_layer, max_layer, n_in = None, ec_layer = None, g_layer = None, tuning = False):
    # define inputs
    ## 1
    layer_units = [(2)**(x+1) if x == 0 else 2**(x) for x in range(base_layer, base_layer + max_layer, 1)][::-1]

    ##2
    kernel_size = layer_units[0]/64
    kernel_sizes = []
    for i in range(max_layer):
        if(len(kernel_sizes) == 0):
            kernel_sizes.append(int(kernel_size))
        elif(len(kernel_sizes) != 0) and (kernel_sizes[-1] > 2):
            kernel_sizes.append(math.ceil(kernel_sizes[-1] - 2))
        else:
            kernel_sizes.append(2)
    
    ##3
    pool_sizes = kernel_sizes
    
    ##4
    strides = [math.ceil(x/2) if x >= 2 else 2 for x in kernel_sizes]
    
    # modify input lists if tuning is to be done
    if tuning:
        layer_units = [tuner_obj_int(hp, f'l_u_d_{x}', x) for x in range(layer_units)]
        kernel_sizes = [tuner_obj_int(hp, f'k_s_d_{x}', x) for x in range(kernel_sizes)]
        pool_sizes = [tuner_obj_int(hp, f'p_s_d_{x}', x) for x in range(pool_sizes)]
        strides = [tuner_obj_int(hp, f's_d_{x}', x) for x in range(strides)]
    
    # define model inputs
    if n_in is not None:
        model_input = Input(shape = (n_in, 1), name = f'cnn_{model_name}_in')
    else:
        # concatenate layersbl
        concat_layer = concatenate([ec_layer, g_layer])
        model_input = Reshape((concat_layer.shape[1], 1))(concat_layer)
    
    # define first layer of CNN
    CNN = Conv1D(filters = layer_units[0],
                 kernel_size = kernel_sizes[0],
                 padding="valid", activation="relu",
                 name = f'CNN_{model_name}_fl')(model_input)
    CNN = AveragePooling1D(pool_size = pool_sizes[0],
                           strides = strides[0],
                           padding = "same",
                           name = f'CNN_{model_name}_fl_avg')(CNN)
    
    # define variable layers of CNN
    for i in range(len(layer_units[1:-1])):
        CNN = Conv1D(filters = layer_units[1:-1][i],
                     kernel_size = kernel_sizes[1:-1][i],
                     padding = "valid", activation = "relu",
                     name = f'CNN_{model_name}_vl_{i}')(CNN)
        CNN = AveragePooling1D(pool_size = pool_sizes[1:-1][i],
                               strides = strides[1:-1][i],
                               padding = "same",
                               name = f'CNN_{model_name}_vl_{i}_avg')(CNN)

    # define last layer of CNN
    CNN = Conv1D(filters = layer_units[-1],
                 kernel_size = kernel_sizes[-1],
                 padding="valid", activation="relu",
                 name = f'CNN_{model_name}_ll')(CNN)
    
    # flatten last layer
    CNN_output = Flatten(name=f'CNN_{model_name}_flat')(CNN)
    return model_input, CNN_output