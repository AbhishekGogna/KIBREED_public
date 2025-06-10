#!/usr/bin/env python3
# load functions -------------------------------------------------------------------------
import sys
all_args = sys.argv[1:]
source_code_at = all_args[0]
sys.path.append(source_code_at)
from Py.libs import *
from Py.func import *

# Define paths ----------------------------------------------------------------------------------------------------
input_paths = read_json(all_args[1])
save_at = all_args[2]
if not os.path.exists(save_at):
    os.makedirs(save_at, exist_ok = True)
task_name = all_args[3]
out_paths = input_paths # so that more keys-value pairs are appneded to an existing dictionary. should remove it and set out_paths to a dictionary if connection to the next step is needed
out_paths['base_folder'] = f'{save_at}'
out_paths['feature_imp_scores'] = f'{save_at}/feature_imp_scores.csv'
out_paths['feature_imp_model_fit_fig'] = f'{save_at}/feature_imp_model_fit_fig.png'
out_paths['feature_imp_model_fit.log'] = f'{save_at}/feature_imp_model_fit.log'

# define logger -----------------------------------------------------------------------------------------------------
logging.basicConfig(filename=out_paths['feature_imp_model_fit.log'], level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)
logging.info(f'First attempt at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}')

# Produce data ----------------------------------------------------------------------------------------------------
ec_data_raw = pd.read_feather(input_paths["ec_data"]).drop(['idx', 'env'], axis=1).rename(columns={"harvest_env": "idx"}).set_index('idx')
feature_imp_data = pd.read_csv(input_paths['clust_data'])

# format the data 
predictors = ec_data_raw.columns.values.tolist()

# Scaling is not needed per se
#feature_imp_data_to_scale = feature_imp_data.reset_index(drop = True)
#feature_imp_data_mod, feature_imp_data_scl = scale_data(feature_imp_data_to_scale.loc[:, predictors],
#                                                        predictors,
#                                                        feature_imp_data_to_scale.index)
#feature_imp_data_mod_2 = feature_imp_data_to_scale.loc[:, ['Env', 'Env_coded']].merge(feature_imp_data_mod, 
#                                                                                      how='left', left_index=True, right_index=True, 
#                                                                                      sort=False, suffixes=('_raw', '_scaled'))

## Drop duplicates
data_fil = feature_imp_data.loc[:, ["groups"] + predictors].drop_duplicates().reset_index(drop = True)

## Define data for classifier
y_0 = data_fil.loc[:, 'groups'].values.astype('int').tolist()
y = [x -1 for x in y_0] # adjust indices for python
X = data_fil.loc[:, predictors].astype('float32')

## change colnames
col_names = pd.DataFrame({"old" : X.columns,
                          "new" : ["var_" + str(i) for i in range(0, len(X.columns))]})
X.columns = col_names['new']

## create splits
train_x, test_x, train_y, test_y = train_test_split(X, y, test_size=0.2, random_state=42 )

## define eval set
evalset = [(train_x, train_y), (test_x, test_y)]

# Define classifier
importance_types = ['weight', 'gain', 'cover', 'total_gain', 'total_cover']

#‘weight’: the number of times a feature is used to split the data across all trees.
#‘gain’: the average gain across all splits the feature is used in.
#‘cover’: the average coverage across all splits the feature is used in.
#‘total_gain’: the total gain across all splits the feature is used in.
#‘total_cover’: the total coverage across all splits the feature is used in.

desired_type = importance_types[1]
model = XGBClassifier(n_jobs = 60, 
                      importance_type = desired_type,
                      learning_rate=0.1,
                      n_estimators=3000,
                      max_depth=30,
                      early_stopping_rounds = 100,
                      eval_metric='mlogloss',
                      seed=27,
                      verbosity = 0) # 3 for debug

# Training the model on the training data
model.fit(train_x, train_y, eval_set = evalset,  verbose = False) #https://machinelearningmastery.com/tune-xgboost-performance-with-learning-curves/

## Test_set accuracy
predictions_test = model.predict(test_x)

# Calculating accuracy
logging.info(f'Test set: Accurately predicted {accuracy_score(test_y, predictions_test, normalize=False)} out of {len(test_y)} instances')
logging.info(f'Test accuracy:, {round(accuracy_score(test_y, predictions_test), 2)}')

## Visualize model fit

# Retrieve performance metrics
results = model.evals_result()
all_vals = results['validation_0']['mlogloss'] + results['validation_1']['mlogloss']

# Plot learning curves
plt.plot(results['validation_0']['mlogloss'], label='train')
plt.plot(results['validation_1']['mlogloss'], label='test')
plt.vlines(model.best_iteration, 0, int(math.ceil(max(all_vals))), linestyles = "dashed", color = "black")

# Show the legend
plt.legend()

# Save the plot as an image (e.g., PNG)
plt.savefig(out_paths['feature_imp_model_fit_fig'])

# Extract feature importance scores
default_values = pd.DataFrame({f'{desired_type}' : model.feature_importances_},
                              index = train_x.columns.astype('str')).merge(col_names, left_index = True, right_on = 'new', how = "left")

#property feature_importances_: ndarray
#    Feature importances property, return depends on importance_type parameter. When model trained with
#    multi-class/multi-label/multi-target dataset, the feature importance is “averaged” over all targets. 
#    The “average” is defined based on the importance type. For instance, if the importance type is “total_gain”, then
#    the score is sum of loss change for each split from all trees.
#        Returns
#            • feature_importances_ (array of shape [n_features] except for multi-class)
#            • linear model, which returns an array with shape (n_features, n_classes)
#    feature_importances_ has the same values as get_score(), except that the values of the former are scaled so that they sum to 1.
#sum(default_values['gain'].values) # pretty close to 1

# Extract all
feature_df = None
for f_type in importance_types:
    importance_scores = model.get_booster().get_score(importance_type = f_type)
    score_row = pd.DataFrame(importance_scores, index = [f'{f_type}'])
    print(len(score_row.columns))
    if feature_df is None:
        feature_df = score_row
    else:
        feature_df = pd.concat([feature_df, score_row])
feature_df_transpose = feature_df.transpose() # Zero-importance features will not be included
scores_df = default_values.merge(feature_df_transpose, left_on = 'new', right_index = True, how = "left", suffixes = ['_scaled', '_raw']).drop(['new'], axis = 1)

# Produce output --------------------------------------------------------------------------------------------------
scores_df.to_csv(out_paths['feature_imp_scores'], index=False)
# Finish off the script -------------------------------------------------------------------------------------------
write_json(out_paths, f'{all_args[4]}')
print(f'{task_name} completed successfully. Paths written at {all_args[4]}')