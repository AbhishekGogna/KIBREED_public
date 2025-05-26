#!/proj/py_env/bin/python3
# sample command line input - pred_script.py 1 /proj/tmp/
# load functions -------------------------------------------------------------------------------------------
import sys
functions_at = '/proj/ext_dir/src'
sys.path.append(functions_at)
from Py.libs import *
from Py.func import *

# Define variables -------------------------------------------------------------------------------------------
## From command line
all_args = sys.argv[1:]
model_name = all_args[0]
model = importlib.import_module(f'Py.model_{model_name}')
cv_schema = all_args[1]
key = all_args[2]
model_inputs = read_json("/proj/model_args.json") # a josn formatted as string for model inputs if any
tune = model_inputs[f'{model_name}']["tune"]
#hparams_at = '/proj/ext_dir/results/M2/hp_tuning/best_params.pkl' #todo: add this if you want to add one hyparameter for all runs

## Constant variables
base_dir = '/proj' # where to save predictions
inputs_at = '/proj/ext_dir/results/Py/preprocessing'
tmp_at = '/proj/tmp_data'
os.system(f'export TMPDIR={tmp_at}')

## Run specifc
tb_cb = f'{base_dir}/callback_data/tb_cb'
mc_cb = f'{base_dir}/callback_data/mc_cb/model.ckpt'
tb_cb_tuning = f'{base_dir}/callback_data/tb_cb_tuning'
tuning_save_at = f'{base_dir}'
tune_dir = f'{model_name}/hp_tuning'
tuned_model_at = f'{base_dir}/model/model_tuned'
model_save_at = f'{base_dir}/model/model_fitted'
param_save_at = f'{base_dir}/model/best_params.json'
pred_save_at = f'{base_dir}/pred'
logs_at = f'{base_dir}/predictions.log'
path_to_pred_file = f'{pred_save_at}/output.csv'
cv_data = read_json(f'{inputs_at}/{cv_schema}')

# define logger -----------------------------------------------------------------------------------------------------
logging.basicConfig(filename=logs_at, level=logging.DEBUG, filemode='w')
logger = logging.getLogger(__name__)
logging.info(f'First attempt at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}')

# check for gpu -----------------------------------------------------------------------------------------------------
gpus = tf.config.list_physical_devices('GPU')
if gpus:
  # Restrict TensorFlow to only use the first GPU
  try:
    #tf.config.set_visible_devices(gpus[0], 'GPU')
    logical_gpus = tf.config.list_logical_devices('GPU')
    logging.info(f'{len(gpus)} Physical GPUs, {len(logical_gpus)} Logical GPU')
  except RuntimeError as e:
    # Visible devices must be set before GPUs have been initialized
    logging.info(e)

# conditional execution -------------------------------------------------------------------------------------------------
if not exists(path_to_pred_file):

    # read in data ------------------------------------------------------------------------------------------------------
    if "acr" in model_name:
        # set which data to load
        if "cv" in cv_schema:
            id_str = "cv"
        elif ("st" in cv_schema or "sce" in cv_schema):
            id_str = "st_sce"
        ## g_a data
        g_a_data = np.load(f'{inputs_at}/acr_g_a_{id_str}.npy')
        scaler_g_a = read_pkl(f'{inputs_at}/g_a.scl')
        
        ## p_data
        p_data = pd.read_csv(f'{inputs_at}/acr_p_{id_str}.csv')
        scaler_p = read_pkl(f'{inputs_at}/acr_p.scl')
        
        ## add further cols
        p_data["run_idx"] = re.sub(r"(r\_\d+)\_(f\_\d+)", r"\1", key)
        p_data["fold_idx"] = re.sub(r"(r\_\d+)\_(f\_\d+)", r"\2", key)
        p_data["cv_idx"] = None
        p_data["env"] = None
        
    elif "wtn" in model_name:
        ## g_a data
        g_a_data = np.load(f'{inputs_at}/wtn_g_a.npy')
        scaler_g_a = read_pkl(f'{inputs_at}/g_a.scl')

        ## g_d_data
        g_d_data = np.load(f'{inputs_at}/wtn_g_d.npy')
        scaler_d_a = read_pkl(f'{inputs_at}/g_d.scl')

        ## e_data
        #e_data = np.load(f'{inputs_at}/wtn_e.npy')
        #scaler_e = read_pkl(f'{inputs_at}/wtn_e.scl')
    
        ## ec_data
        ec_data = read_pkl(f'{inputs_at}/wtn_ec.pkl')
        scaler_e = read_pkl(f'{inputs_at}/wtn_ec.scl')
        
        ## g_s_data
        g_s_data = np.load(f'{inputs_at}/wtn_g_s.npy')
        scaler_g_s = read_pkl(f'{inputs_at}/wtn_g_s.scl')

        ## p_data
        p_data = pd.read_csv(f'{inputs_at}/wtn_p.csv')
        scaler_p = read_pkl(f'{inputs_at}/wtn_p.scl')
        
        ## add further cols
        p_data["run_idx"] = re.sub(r"(run\S+)\_(cv\d)", r"\1", key)
        p_data["fold_idx"] = None
        p_data["cv_idx"] = re.sub(r"(run\S+)\_(cv\d)", r"\2", key)
        p_data["idx_col"] = None
        p_data = p_data.rename(columns = {"Type" : "type", "BLUES_dt": "BLUEs_raw"})
        p_data['type'] = p_data['type'].apply(lambda x: 'Non_hybrid' if x != 'Hybrid' else x)
            
    ## further additions to p_data
    p_data = p_data.rename(columns = {"Geno_new" : "geno", "Env":"env"})
    p_data["model_idx"] = model_name
    p_data["run_type"] = re.sub(r"(\S+)\.json", r"\1", cv_schema)
    
    if "acr" in model_name:
        out_cols = ["series", "geno", "type", "idx_col", "run_type", "run_idx", "fold_idx", "model_idx", "cv_idx", "BLUEs_raw", "BLUEs_scaled"]
    else:
        out_cols = ["env", "geno", "type", "idx_col", "run_type", "run_idx", "fold_idx", "model_idx", "cv_idx", "BLUEs_raw", "BLUEs_scaled"]

    p_data = p_data.loc[:, out_cols]

    # create train, val and test sets -----------------------------------------------------------------------------------
    if "val" in cv_data[key].keys(): # the function substracts one from index value. 
        train_set, val_set, test_set = create_train_val_data(index_train = cv_data[key]["train"], index_test = cv_data[key]["test"], index_val = cv_data[key]["val"])
    else:
        train_set, val_set, test_set = create_train_val_data(index_train = cv_data[key]["train"], index_test = cv_data[key]["test"]) # makes a val set out of the training set
    
    ## get data as tensors 
    if "acr" in model_name:
        target_data = [p_data.loc[:, "BLUEs_scaled"].values.astype('float32'), \
                       g_a_data.astype('float32')]
    elif "wtn" in model_name:
        target_data = [p_data.loc[:, "BLUEs_scaled"].values.astype('float32'), \
                       ec_data.astype('float32'), \
                       g_a_data.astype('float32'), \
                       g_d_data.astype('float32'), \
                       g_s_data.astype('float32')]
        
    train_data = [x[train_set] for x in target_data]
    val_data = [x[val_set] for x in target_data]
    test_data = [x[test_set] for x in target_data]
    #if "wtn" in model_name:
    #    train_data.insert(1, [x[train_set].reshape(len(train_set), x.shape[1], 1).astype('float32') for x in ec_data])
    #    val_data.insert(1, [x[val_set].reshape(len(val_set), x.shape[1], 1).astype('float32') for x in ec_data])
    #    test_data.insert(1, [x[test_set].reshape(len(test_set), x.shape[1], 1).astype('float32') for x in ec_data]) # i dont have an elegant way to do this above so this chunk has to exist. 
    
    # customize input data and model tuner ------------------------------------------------------------------------------
    train_y = train_data[0]
    val_y = val_data[0]
    test_y = test_data[0]
    if "acr" in model_name:
        train_x = train_data[1]
        val_x = val_data[1]
        test_x = test_data[1]
        model_tuner = model.tuner
    elif "wtn" in model_name:
        if model_name == "wtn_CNN_LSTM": #modify this part to differentiate between M1 and M2 types
            train_x = [train_data[2], train_data[3]] # env_data and geno_data_a
            val_x = [val_data[2], val_data[3]]
            test_x = [test_data[2], test_data[3]]
            model_tuner = model.tuner(marker_n = train_x[1].shape[1], # number of markers
                                      env_n = train_x[0].shape[2], # number of env variables input
                                      days_n = train_x[0].shape[1]) # number of days in an env 
        elif model_name == "wtn_CNN_EC":
            train_x = [train_data[2], train_data[1]]
            val_x = [val_data[2], val_data[1]]
            test_x = [test_data[2], test_data[1]]
            model_tuner = model.tuner(marker_n = train_x[0].shape[1], # number of markers
                                      ec_n = train_x[1].shape[1],
                                      tune = tune) # model parameters
        elif model_name == "wtn_CNN_CGM_EC":
            train_x = [train_data[2], train_data[4]]
            val_x = [val_data[2], val_data[4]]
            test_x = [test_data[2], test_data[4]]
            model_tuner = model.tuner(marker_n = train_x[0].shape[1], # number of markers
                                      g_s_n = train_x[1].shape[1],
                                      tune = tune) # model parameters
        
    # Hparam Tuning -------------------------------------------------------------------------------------------------
    if not exists(tuned_model_at):
        start_time_tuning = time.time()
        stop_early = EarlyStopping(monitor='val_loss', patience=5, min_delta = 0.001)
        tb_cv_tuner = TensorBoard(tb_cb_tuning)
        tuner = kt.Hyperband(hypermodel=model_tuner,
                             objective=kt.Objective("val_mean_squared_error", direction="min"),
                             max_epochs=100,
                             factor=4,
                             hyperband_iterations=1,
                             overwrite = True,
                             directory=tuning_save_at,
                             project_name=tune_dir,
                             seed=30)
        tuner.search(train_x, train_y,
                     epochs=100,
                     validation_data=(val_x, val_y),
                     callbacks=[stop_early, tb_cv_tuner],
                     verbose=0)
        
        # save parameters
        for num_params in [3, 2, 1]:
            print(num_params)
            try:
                top3_params = tuner.get_best_hyperparameters(num_trials=num_params)
                if top3_params:
                    break  # If successful, exit the loop
            except tf.errors.NotFoundError as e:
                print("An error occurred:", e)
                if num_params == 1:
                    raise Exception("Error: Failed to retrieve best models with num_models=1. Script halted.")
        params = top3_params[0].values  # best hyperparameter values # can igonore warnings # https://stackoverflow.com/questions/58289342/tf2-0-translation-model-error-when-restoring-the-saved-model-unresolved-object
        write_json(params, param_save_at)
        
        # save model
        for num_models in [3, 2, 1]:
            print(num_models)
            try:
                top3_models = tuner.get_best_models(num_models=num_models)
                if top3_models:
                    break  # If successful, exit the loop
            except tf.errors.NotFoundError as e:
                print("An error occurred:", e)
                if num_models == 1:
                    raise Exception("Error: Failed to retrieve best models with num_models=1. Script halted.")
        best_model = top3_models[0]
        best_model.save(tuned_model_at)
        best_model = load_model(tuned_model_at) # loads a model with weights to run with test set directely
        #best_model = model_tuner(top3_params[0]) # alternative, but gives a weaker trend
 
        # clear space
        try:
            os.system(f'rm -rf {tuning_save_at}/{tune_dir}')
        except:
            logging.debug(f'Cannot delete {tuning_save_at}/{tune_dir}. Do it manually')
        
        # write Hparam tuning log
        end_time_tuning = time.time()
        elapsed_time_tuning = end_time_tuning - start_time_tuning
        logging.info(f'HP tuning took {elapsed_time_tuning} seconds, or {elapsed_time_tuning/60} minutes, or {elapsed_time_tuning/3600} hours')
        logging.info(f'parameters \n {pformat(params)}')
    else:
        logging.info(f'subsequent attempt at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}. Hparams were tuned earlier and were loaded for this run.')
        #os.system(f'rm -rf tune_dir')
        best_model = load_model(tuned_model_at) #todo: if log file exists then figure out a way to apeend lines rather than scratch all the previpous logs
    # Perform predictions -----------------------------------------------------------------------------------------------
    my_model = best_model
    if not exists(model_save_at):
        start_time_fit = time.time()
        fit_params = {'fit' : {'batch_size' : 32, # default is 32
                               'epochs' : 100,
                               'verbose' : 2,
                               'shuffle' : True,
                               'tensorboard_fp' : tb_cb,
                               'checkpoint_fp' : mc_cb}}
        my_model_fit = fit_model(final_model = my_model, params = fit_params, 
                                 train_x = train_x, 
                                 train_y = train_y, 
                                 val_x = val_x,
                                 val_y = val_y)
        my_model_fit.save(model_save_at)
        #save_model(model = my_model_fit, path = model_save_at, model_name = f'model_{key}')
        end_time_fit = time.time()
        elapsed_time_fit = end_time_fit - start_time_fit
        logging.info(f'Model fitting took {elapsed_time_fit} seconds, or {elapsed_time_fit/60} minutes, or {elapsed_time_fit/3600} hours')
    else:
        logging.info(f'subsequent attempt at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}. Model fit was done earlier and weights were loaded for this run.')
        my_model_fit = load_model(model_save_at) #todo: if log file exists then figure out a way to apeend lines rather than scratch all the previpous logs
    # Export data -------------------------------------------------------------------------------------------------------
    ## raw results
    pred_vals_test = predict_values(model = my_model_fit, 
                                    test_x = test_x, 
                                    test_y = test_y, 
                                    index = test_set, 
                                    scaler = scaler_p)
    pred_vals_test = pd.merge(pred_vals_test, 
                              p_data, 
                              how='left', 
                              left_on=['index'], 
                              right_index=True)
    pred_vals_test.to_csv(path_to_pred_file, index=False)

else:
    # define logger -----------------------------------------------------------------------------------------------------
    logging.basicConfig(filename=logs_at, level=logging.DEBUG, filemode='a')
    logger = logging.getLogger(__name__)
    logging.info(f'Followup attempt at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}. The whole script was not executed since pred/output.csv exists. You will need to delete it for script to work.')
