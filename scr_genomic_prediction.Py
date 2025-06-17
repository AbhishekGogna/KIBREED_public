# CNN Genomic Prediction Script - Standalone Version

# Import utilities and libraries
import sys
import os

# Add project directory to path and import utilities
if os.path.exists("/proj"):
    project_path = "/proj"
else:
    project_path = os.getcwd()

sys.path.append(project_path)
from Py.libs import *
from Py.func import *

# =============================================================================
# LOGGING SETUP
# =============================================================================

# Setup logging early
logs_at = os.path.join(project_path, "logs", "genomic_predictions_cnn.log")
os.makedirs(os.path.dirname(logs_at), exist_ok=True)

logging.basicConfig(filename=logs_at, level=logging.DEBUG, filemode='w')
logger = logging.getLogger(__name__)

logger.info("Starting CNN genomic prediction analysis")

# =============================================================================
# CONFIGURATION PARAMETERS
# =============================================================================

logger.info("Setting up configuration parameters")

# Use parameters from R script
PREDICTION_TYPE = "acr"  # Set to "acr" or "wtn" based on R script
RANDOM_SEED = 123
TEST_PROPORTION = 0.2

logger.info(f"Prediction type: {PREDICTION_TYPE}, Random seed: {RANDOM_SEED}")

# Model configuration
if PREDICTION_TYPE == "acr":
    model_name = "acr_CNN"
elif PREDICTION_TYPE == "wtn":
    model_name = "CNN_EV"
    
model = importlib.import_module(f'Py.model_{model_name}')

logger.info(f"Selected model: {model_name}")

# =============================================================================
# DIRECTORY SETUP
# =============================================================================

logger.info("Setting up directory structure")

# Base directories
base_dir = os.path.join(project_path, "results", "CNN", PREDICTION_TYPE)
results_dir = os.path.join(project_path, "results")
tmp_at = os.path.join(project_path, "tmp", "CNN", PREDICTION_TYPE)

# Model-specific directories
callback_dir = os.path.join(base_dir, "callback_data")
tb_cb = os.path.join(callback_dir, "tb_cb")
mc_cb = os.path.join(callback_dir, "mc_cb", "model.ckpt")
tb_cb_tuning = os.path.join(callback_dir, "tb_cb_tuning")
tuning_save_at = tmp_at # will be deleted since it takes a lot of space
tune_dir = f"{model_name}/hp_tuning"
tuned_model_at = os.path.join(base_dir, "model", "model_tuned")
model_save_at = os.path.join(base_dir, "model", "model_fitted")
param_save_at = os.path.join(base_dir, "model", "best_params.json")
pred_save_at = os.path.join(base_dir, "pred")
# Update logs path to use prediction-specific directory
logs_at = os.path.join(base_dir, "predictions.log")
path_to_pred_file = os.path.join(pred_save_at, "output.csv")

# Create necessary directories
dirs_to_create = ["model", "pred", "callback_data/tb_cb", "callback_data/mc_cb", "callback_data/tb_cb_tuning"]
for dir_name in dirs_to_create:
    dir_path = os.path.join(base_dir, dir_name)
    os.makedirs(dir_path, exist_ok=True)
os.makedirs(tmp_at, exist_ok=True)

logger.info(f"Created directories - Base: {base_dir}, Temp: {tmp_at}")

# =============================================================================
# GPU CONFIGURATION
# =============================================================================

logger.info("Checking GPU configuration")
gpus = tf.config.list_physical_devices('GPU')
if gpus:
    try:
        logical_gpus = tf.config.list_logical_devices('GPU')
        logger.info(f'{len(gpus)} Physical GPUs, {len(logical_gpus)} Logical GPU')
    except RuntimeError as e:
        logger.error(f"GPU configuration error: {e}")
else:
    logger.info("No GPU found - running on CPU")

# =============================================================================
# LOAD AND PREPROCESS DATA
# =============================================================================
try:
    logger.info("Loading data files...")
    
    if PREDICTION_TYPE == "acr":
        # Load phenotype data for ACR
        pheno_path = os.path.join(results_dir, "pheno_final_acr.feather")
        pheno_df = pd.read_feather(pheno_path)
        
        # Preprocess phenotype data
        p_data_mod, p_scl = scale_data(pheno_df.loc[:, "BLUEs_acr"].values.reshape(-1, 1), 
                                   ["BLUEs_acr"], pheno_df.index)
        p_data = pheno_df.merge(p_data_mod, how='left', left_index=True, right_index=True, 
                         sort=False, suffixes=('_raw', '_scaled'))
        
        # Load genomic data for ACR
        g_data_add_path = os.path.join(results_dir, "g_data_add_acr.feather")
        g_data_df = pd.read_feather(g_data_add_path).set_index('idx')
        
        # Preprocess genomic data
        g_a_data_scaled, g_a_scl = scale_data(g_data_df, g_data_df.columns, g_data_df.index)
        g_a_data_cv = np.stack([g_a_data_scaled[g_a_data_scaled.index == idx].iloc[0,:].values for idx in p_data['Geno_n']])
        
        # No environmental data for ACR
        ev_data_cv = None
        
    elif PREDICTION_TYPE == "wtn":
        # Load phenotype data for WTN
        pheno_path = os.path.join(results_dir, "pheno_final_wtn.feather")
        pheno_df = pd.read_feather(pheno_path)
        
        # Preprocess phenotype data
        p_data_mod, p_scl = scale_data(pheno_df.loc[:, "BLUEs_wtn"].values.reshape(-1, 1), 
                                   ["BLUEs_wtn"], pheno_df.index)
        p_data = pheno_df.merge(p_data_mod, how='left', left_index=True, right_index=True, 
                         sort=False, suffixes=('_raw', '_scaled'))

        # Load genomic data for WTN
        g_data_add_path = os.path.join(results_dir, "g_data_add_wtn.feather")
        g_data_df = pd.read_feather(g_data_add_path).set_index('idx')
        
        # Preprocess genomic data
        g_a_data_scaled, g_a_scl = scale_data(g_data_df, g_data_df.columns, g_data_df.index)
        g_a_data_cv_0 = np.stack([g_a_data_scaled[g_a_data_scaled.index == idx].iloc[0,:].values for idx in p_data['Geno_n']])
        g_a_data_cv = g_a_data_cv_0.reshape(g_a_data_cv_0.shape[0], g_a_data_cv_0.shape[1], 1)

        # Load environmental data for WTN (optional)
        ev_path = os.path.join(results_dir, "ev_data_wtn.feather")
        try:
            ev_df = pd.read_feather(ev_path).set_index('idx')
            ev_data = np.stack([ev_df[ev_df.index == idx].iloc[0,:].values for idx in p_data['Env_n']])
            ev_data_cv = ev_data.reshape(ev_data.shape[0], ev_data.shape[1], 1)
            logger.info(f"Environmental data loaded: {ev_df.shape[0]} x {ev_df.shape[1]}")
        except:
            ev_data_cv = None
            logger.info("Environmental data not found - continuing without it")
    
    else:
        raise ValueError(f"PREDICTION_TYPE must be 'acr' or 'wtn', got: {PREDICTION_TYPE}")
    
    logger.info(f"Data loaded successfully for {PREDICTION_TYPE.upper()}")
    logger.info(f"Genomic data: {g_a_data_scaled.shape[0]} genotypes x {g_a_data_scaled.shape[1]} markers")
    logger.info(f"Phenotype data: {len(pheno_df)} observations")
    
except Exception as e:
    logger.error(f"Failed to load data: {e}")
    raise

# =============================================================================
# CREATE DATA SPLITS
# =============================================================================
# Get indices for train/test splits from 'set' column
p_data_train = p_data[p_data['set'] == 'Train'].index.tolist()
p_data_test = p_data[p_data['set'] == 'Test'].index.tolist()
train_set, val_set, test_set = create_train_val_data(index_train=p_data_train, index_test=p_data_test)

## Prepare data arrays based on prediction type
if PREDICTION_TYPE == "acr":
    target_data = [p_data.loc[:, "BLUEs_acr_scaled"].values.astype('float32'), \
                   g_a_data_cv.astype('float32')]
elif PREDICTION_TYPE == "wtn":
    target_data = [p_data.loc[:, "BLUEs_wtn_scaled"].values.astype('float32'), \
                   ev_data_cv.astype('float32'), \
                   g_a_data_cv.astype('float32')]

# Split data arrays by train/val/test indices
train_data = [x[train_set] for x in target_data]
val_data = [x[val_set] for x in target_data]
test_data = [x[test_set] for x in target_data]

# Extract target variables (y)
train_y = train_data[0]
val_y = val_data[0]
test_y = test_data[0]

# Prepare input features (x) based on model type
if "acr" in model_name:
    train_x = train_data[1]
    val_x = val_data[1]
    test_x = test_data[1]
    model_tuner = model.tuner
elif PREDICTION_TYPE == "wtn":
    train_x = [train_data[2], train_data[1]]  # [genomic, environmental]
    val_x = [val_data[2], val_data[1]]
    test_x = [test_data[2], test_data[1]]
    model_tuner = model.tuner(marker_n=train_x[0].shape[1], # number of markers
                              ev_n=train_x[1].shape[1])

# =============================================================================
# HYPERPARAMETER TUNING
# =============================================================================

if not exists(tuned_model_at):
    logger.info("Starting hyperparameter tuning...")
    start_time_tuning = time.time()
    
    # Setup callbacks
    stop_early = EarlyStopping(monitor='val_loss', patience=5, min_delta=0.001)
    tb_cv_tuner = TensorBoard(tb_cb_tuning)
    
    # Initialize tuner
    tuner = kt.Hyperband(
        hypermodel=model_tuner,
        objective=kt.Objective("val_mean_squared_error", direction="min"),
        max_epochs=100,
        factor=4,
        hyperband_iterations=1,
        overwrite=True,
        directory=tuning_save_at,
        project_name=tune_dir,
        seed=RANDOM_SEED
    )
    
    # Run hyperparameter search
    tuner.search(
        train_x, train_y,
        epochs=100,
        validation_data=(val_x, val_y),
        callbacks=[stop_early, tb_cv_tuner],
        verbose=0
    )
    
    # Save best parameters
    for num_params in [3, 2, 1]:
        try:
            top3_params = tuner.get_best_hyperparameters(num_trials=num_params)
            if top3_params:
                break
        except tf.errors.NotFoundError as e:
            logger.error(f"Error retrieving parameters: {e}")
            if num_params == 1:
                raise Exception("Failed to retrieve best parameters")
    
    params = top3_params[0].values
    write_json(params, param_save_at)
    
    # Save best model
    for num_models in [3, 2, 1]:
        try:
            top3_models = tuner.get_best_models(num_models=num_models)
            if top3_models:
                break
        except tf.errors.NotFoundError as e:
            logger.error(f"Error retrieving models: {e}")
            if num_models == 1:
                raise Exception("Failed to retrieve best models")
    
    best_model = top3_models[0]
    best_model.save(tuned_model_at)
    best_model = load_model(tuned_model_at)
    
    # Cleanup
    try:
        tune_dir_path = os.path.join(tuning_save_at, tune_dir)
        if os.path.exists(tune_dir_path):
            os.system(f'rm -rf {tune_dir_path}')
    except:
        logger.warning(f'Cannot delete {tune_dir_path}')
    
    # Log tuning results
    end_time_tuning = time.time()
    elapsed_time_tuning = end_time_tuning - start_time_tuning
    logger.info(f'HP tuning took {elapsed_time_tuning/60:.2f} minutes')
    logger.info(f'Best parameters: \n{pformat(params)}')
else:
    logger.info("Loading previously tuned model...")
    best_model = load_model(tuned_model_at)

# =============================================================================
# MODEL TRAINING
# =============================================================================

if not exists(model_save_at):
    logger.info("Starting model training...")
    start_time_fit = time.time()
    
    fit_params = {
        'fit': {
            'batch_size': 32,
            'epochs': 100,
            'verbose': 2,
            'shuffle': True,
            'tensorboard_fp': tb_cb,
            'checkpoint_fp': mc_cb
        }
    }
    
    my_model_fit = fit_model(
        final_model=best_model,
        params=fit_params,
        train_x=train_x,
        train_y=train_y,
        val_x=val_x,
        val_y=val_y
    )
    
    my_model_fit.save(model_save_at)
    
    end_time_fit = time.time()
    elapsed_time_fit = end_time_fit - start_time_fit
    logger.info(f'Model fitting took {elapsed_time_fit/60:.2f} minutes')
else:
    logger.info("Loading previously fitted model...")
    my_model_fit = load_model(model_save_at)

# =============================================================================
# PREDICTIONS AND EXPORT
# =============================================================================

logger.info("Making predictions...")

# Generate predictions
pred_vals_test = predict_values(
    model=my_model_fit,
    test_x=test_x,
    test_y=test_y,
    index=test_set,
    scaler=p_scl
)

# Merge with phenotype data
pred_vals_test = pd.merge(
    pred_vals_test,
    p_data,
    how='left',
    left_on=['index'],
    right_index=True
)

# Save predictions
pred_vals_test.to_csv(path_to_pred_file, index=False)
logger.info(f"Predictions saved to {path_to_pred_file}")

# Clean up temporary files
temp_files = [f for f in os.listdir(tmp_at) if f.startswith(f"{PREDICTION_TYPE}_")]
if temp_files:
    for temp_file in temp_files:
        try:
            os.remove(os.path.join(tmp_at, temp_file))
        except:
            pass
    logger.info(f"Cleaned up {len(temp_files)} temporary files")

logger.info("Analysis completed successfully")

# Create a compact correlation summary
if PREDICTION_TYPE == "wtn":
   # Within environments - group by Env_n and Type
   within_env_corr = (pred_vals_test.groupby(['Env_n', 'Type'])
                      .apply(lambda x: x['obs'].corr(x['pred']).round(2))
                      .reset_index()
                      .rename(columns={0: 'correlation'}))
   
   # Mean across environments - group by Type and take mean of correlations
    correlations_by_type = (within_env_corr.groupby('Type')['correlation']
                            .mean()
                            .round(2)
                            .reset_index())

elif PREDICTION_TYPE == "acr":
    correlations_by_type = (pred_vals_test.groupby('Type')
                            .apply(lambda x: x['obs'].corr(x['pred']).round(2))
                            .reset_index()
                            .rename(columns={0: 'correlation'}))
                          
correlation_text = "; ".join([f"{row['Type']}: {row['correlation']}" 
                             for _, row in correlations_by_type.iterrows()])

logger.info(f"Analysis completed successfully - Correlations by type: {correlation_text}")
