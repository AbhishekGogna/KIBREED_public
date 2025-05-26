from .libs import *

# Common functions -----------------------------------------------------------
def get_random_string(length):
    # With combination of lower and upper case
    result_str = ''.join(random.choice(string.ascii_letters) for i in range(length))
    # print random string
    return result_str

def scale_data(to_transform, pd_cols, pd_index):
    scaler = MinMaxScaler((0,1))
    data_scaled = scaler.fit_transform(to_transform)
    data_scaled_df = pd.DataFrame(data_scaled, columns = pd_cols, index = pd_index)
    return data_scaled_df, scaler

def read_pkl(path):
    with open(path, "rb") as fp:   # Unpickling
        data = pickle.load(fp)
    return data

def write_pkl(data, path, verbose = False):
    with open(path, "wb") as fp:   # pickling
        pickle.dump(data, fp)
    if verbose:
        return print("Done")

def read_json(path):
    with open(path, encoding = "utf8") as json_file:
        data = json.load(json_file)
    return data
def write_json(data, path, verbose = False):
    with open(path, "w") as fp:   
        json.dump(data, fp)
    if verbose:
        return print("Done")

def print_function(func_name):
    lines = inspect.getsource(func_name)
    return print(lines)

def redo_argument_string(args):
    ret = args.split("\&")
    my_dict = {}
    for i in ret:
        key, value = i.split("=")
        my_dict[key] = value
    return my_dict

def scale_tanh(myyield, verbose): # from https://bitbucket.org/bucklerlab/maize_yield_pred/src/master/python_notebooks/5_Train_Test_replicate_val_train_tanh_final.ipynb
    #re-scale data between 0 and 1
    scaler = MinMaxScaler((0,1))
    myyield_tanh = scaler.fit_transform(myyield.reshape(-1, 1))
    myyield_tanh = myyield_tanh.flatten()
    if verbose:
        print(myyield.max(), myyield.min(), myyield.mean(), myyield.shape)
        print(scaler, myyield.shape, myyield_tanh.shape, myyield_tanh.max(), myyield_tanh.min(), myyield_tanh.mean())
    return scaler, myyield_tanh

def inverse_scale(scaler, myyield_tanh, verbose):
    #undo scalling
    myyield_tanh_inv = scaler.inverse_transform(myyield_tanh.reshape(-1, 1))
    myyield_tanh_inv = myyield_tanh_inv.flatten()
    if verbose:
        print(scaler, myyield_tanh.shape, myyield_tanh_inv.shape, myyield_tanh_inv.max(), myyield_tanh_inv.min(), myyield_tanh_inv.mean())
    return myyield_tanh_inv

def extract_dict_elements(dict_data, string, count = True):
    out = {key: value for key, value in dict_data.items() if string in key}
    out_len = len(out)
    if out_len == 0:
        print("Function failed, nothing to extract. Check manually")
    else:
        if count:
            return out, out_len
        else:
            return out

# Plots
def plot_fig(plot_data, save_at):
    ## plot parameters
    SMALL_SIZE = 16
    MEDIUM_SIZE = 20
    BIGGER_SIZE = 24
    # from - https://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    
    # plot 
    plt.figure()
    fig, ax = plt.subplots()
    ax.set_title(f'Correlations over {len(set(plot_data["run"].values))} runs')
    ax.set_ylim(0, 1)
    ax.set_xticklabels(['Hybrid', 'Non_hybrid'])
    ax.boxplot([plot_data[plot_data['type'] == 'Hybrid'].loc[:, 'Pred'].values, 
                 plot_data[plot_data['type'] == 'Non_hybrid'].loc[:, 'Pred'].values])
    plt.rcParams['figure.figsize'] = [20, 16]
    plt.savefig(save_at)
    return print("fig saved")

# Directory management
def set_dirs(base_dir_path, verbose = True, run_id = None):
    if run_id is None:
        run_id = time.strftime("run_%Y_%m_%d")
        base_folder = base_dir_path + '/' + run_id
    else:
        base_folder = base_dir_path + '/' + f'{str(run_id)}'
    cb_at = base_folder + '/callback_data'
    tb_cb = cb_at + '/tb_cb'
    mc_cb = cb_at + '/mc_cb/'
    pred_at = base_folder + '/pred'
    model_at = base_folder + '/model'
    tmp_at = base_folder + '/tmp_data'
    if(not os.path.isdir(base_folder)):
        os.system(f'mkdir -p  {base_folder} {pred_at} {model_at} {cb_at} {tb_cb} {mc_cb} {tmp_at}')
    if (verbose):
        print(f'base folder at {base_folder}, \ncallbacks at {cb_at}, \npredictions at {pred_at}, \nmodel at {model_at}, \ntmp at {tmp_at}')
    # output
    out = {}
    out['base_folder'] = base_folder
    out['tb_cb'] = tb_cb
    out['mc_cb'] = mc_cb
    out['pred_at'] = pred_at
    out['model_at'] = model_at
    out['tmp_at'] = tmp_at
    return out 
        
# Data wrangling
def create_train_val_data(index_train, index_test, index_val = None, prop = 0.1):
    if index_val is None:
        val_set = random.sample(index_train, int(len(index_train)*prop)) # cretes validation set from the remaining non_test set 
        train_set = list(set(index_train).difference(val_set))
    else:
        val_set = index_val
        train_set = index_train
    
    check = any(item in val_set for item in train_set)
    
    # adjust for index number in python
    train_set = [x - 1 for x in train_set]
    val_set = [x - 1 for x in val_set]
    test_set = [x - 1 for x in index_test]
    
    if check:
        print("function failed since some elemets of val arer in the train set")
    else:
        return train_set, val_set, test_set

def create_rn_frt_data(ecov_data, geno_data_a, idx, verbose = True):
    ec_data_sub = [x[idx] for x in ecov_data] # get this to work with train_data
    ec_data_raw_sub_reshaped = ec_data_sub[0]
    for i in ec_data_sub[1:]:
        ec_data_raw_sub_reshaped = np.concatenate([ec_data_raw_sub_reshaped, i], axis=1)
    ec_data_raw_sub_reshaped = ec_data_raw_sub_reshaped.astype('float32')
    
    g_data_reshaped = geno_data_a.reshape(geno_data_a.shape[0], geno_data_a.shape[1])
    
    out = np.concatenate([ec_data_raw_sub_reshaped, g_data_reshaped], axis = 1)
    
    return out

# Model fitting
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

def save_model(model, path, model_name = "kibreed_pred"):
    model_json = model.to_json()
    with open(path + '/' + model_name + ".json", "w") as json_file:
        json_file.write(model_json)
    # serialize weights to HDF5
    model.save_weights(path + '/' + model_name + ".h5")
    print("Saved model to disk")
    return

def stratified_sampling(data, test_prop = 0.2):
    test_size = test_prop * data.shape[0]
    fractions = np.rint(test_size * (data.groupby('Series').count().iloc[:, 0]/data.shape[0]))
    test_set_index = []
    for series in data['Series'].unique():
        target_index = data[data['Series'] == series].index.to_numpy()
        target_number = int(fractions[series])
        test_set_index.extend(np.random.choice(target_index, target_number, replace = False))
    return test_set_index

def get_sets(raw_p_data, to_test , val_pcent = 10,
             adj_year = True, remove_loc = True,
             remove_year = True, remove_geno = True,
             verbose = False, only_shapes = True, debug = False):
    # assign data
    p_data = raw_p_data
    words = to_test
    
    # subset required data
    base = r'^{}'
    expr = '(?=.*{})'  
    searchfor = base.format(''.join(expr.format(w) for w in words))
    
    # test set
    test_set = p_data.loc[p_data['Env'].str.contains(searchfor),:]
    test_idx = test_set.index
    
    # further filtering 
    Loc = test_set.Loc.unique().tolist()
    Year = test_set.Year.unique().tolist()
    Geno = test_set.Geno_new.unique().tolist()
    not_test = p_data[p_data.index.isin(test_idx) == False]
    
    if adj_year:
        if len(Year) > 1:
            yr_adj = np_min(Year)
        else:
            yr_adj = Year[0]
        yr_adj_idx = not_test[not_test.loc[:, 'Year'] > yr_adj].index # removes data points belonging to years after the test set
    
    if remove_loc:
        fil_loc = not_test[not_test.loc[:, 'Loc'].isin(Loc)].index
    else:
        fil_loc = [] # remove test locations from not_test data
        
    if remove_year:
        fil_year = not_test[not_test.loc[:, 'Year'].isin(Year)].index
    else:
        fil_year = []  # remove test years from not_test data
    
    if remove_geno:
        fil_geno = not_test[not_test.loc[:, 'Geno_new'].isin(Geno)].index
    else:
        fil_geno = [] # remove test geno from non test data
    
    fil = fil_loc.tolist() + fil_year.tolist() + fil_geno.tolist() + yr_adj_idx.to_list()
    fil = set(fil) #to get unique values.
    not_test_fil = not_test[not_test.index.isin(fil) == False]
    
    if debug:
        print(not_test[not_test.index.isin(fil) == True].loc[:, ['Env', 'BLUES_dt']].groupby('Env').count()) # for debugisng
    
    # validation set
    val_idx = random.sample(not_test_fil.index.to_list(),int(not_test_fil.shape[0]*(val_pcent/100))) # cretes validation set from the remaining non_test set 
    val_set = not_test_fil[not_test_fil.index.isin(val_idx) == True]
    
    #train set
    train_set = not_test_fil[not_test_fil.index.isin(val_idx) == False]
    
    if (train_set.shape[0] + val_set.shape[0] + test_set.shape[0] + len(fil)) == p_data.shape[0]:
        if verbose:
            print(f"Sets -> Train = {train_set.shape[0]}, Val = {val_set.shape[0]}, Test = {test_set.shape[0]}, and ignored = {len(fil)}")
    else:
        print('Something is wrong since the data taken in and excluded do not add up!. Check manually')
        print(f"Sets -> Train = {train_set.shape[0]}, Val = {val_set.shape[0]}, Test = {test_set.shape[0]}, and ignored = {len(fil)}")
        
    # create tensors
    #train_tensors = p_data[train_set.index].astype('float32'), raw_e_data[train_set.index].astype('float32'), raw_g_a_data[train_set.index].astype('float32'), raw_g_d_data[train_set.index].astype('float32')
    
    if only_shapes:
        return train_set.shape[0], val_set.shape[0], test_set.shape[0], len(fil)
    else:
        return train_set.index, val_set.index, test_set.index, fil

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

def dense_post_concat(hp, inputs, concatenated_model, base_layer = 6, max_layer = 5, dp_rate = 0.2, tuning = False):
    # define layer inputs
    layer_units = [(2)**(x+1) if x == 0 else 2**(x) for x in range(base_layer, base_layer + max_layer, 1)][::-1] # 6 layers
    dp_rate = 0.2 * np.ones(max_layer)
    
    if tuning:
        layer_units = [tuner_obj_int(hp, f'u_den_{x}', x) for x in layer_units]
        dp_rate = [hp.Float(f'dp_for_l_{layer_units[x]}', min_value=0.1, max_value=0.5, step = 0.01, default = dp_rate[x]) for x in range(len(dp_rate))]
    
    # define first layer
    model_concat = Dense(units = layer_units[0],
                         activation ='relu', 
                         name = f'concat_fl')(concatenated_model) # dense layer
    model_concat = Dropout(rate = dp_rate[0], 
                           name = f'concat_drop_fl')(model_concat) # dropout layer
    
    # define variable layers
    for i in range(len(layer_units[1:])):
        model_concat = Dense(units = layer_units[1:][i], 
                             activation ='relu', 
                             name = f'concat_vl_{i}')(model_concat) # dense layer
        model_concat = Dropout(rate = dp_rate[1:][i], 
                               name = f'concat_drop_vl_{i}')(model_concat) # dropout layer
    
    # define output layers
    model_concat_out = Dense(1, activation='relu',name="concat_out")(model_concat) # 1 unit since its a regression model
    
    # compile model
    compiled_model = Model(inputs=inputs, outputs = model_concat_out)
    compiled_model.compile(loss = 'mean_absolute_error',
                           optimizer = Adam(learning_rate = hp.Float("compile_l_rate",
                                                                     min_value=1e-5,
                                                                     max_value=1e-2,
                                                                     sampling="log"), #lr
                                            beta_1 = hp.Float("compile_beta_val_1",
                                                              min_value=0,
                                                              max_value=1), #beta1
                                            beta_2 = hp.Float("compile_beta_val_2",
                                                              min_value=0,
                                                              max_value=1)), # beta2
                           metrics = ['mean_squared_error'])
    # produce output
    return compiled_model
