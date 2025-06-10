#!/usr/bin/env python3
# load functions -------------------------------------------------------------------------
import sys
all_args = sys.argv[1:]
source_code_at = all_args[0]
sys.path.append(source_code_at)
from Py.libs import *
from Py.func import *

# Define variables -----------------------------------------------------------------------
## From command line
input_paths = read_json(all_args[1])
save_at = all_args[2]
if not os.path.exists(save_at):
    os.makedirs(save_at, exist_ok = True)
task_name = all_args[3]
print(save_at)

## Constant variables
elements_common = [
    'g_a.scl', 'g_d.scl', # genetic data scaling factors: source additive and dominance
    'preprocessing.log' # placeholder for storing logs from this file
]
elements_acr = [
    'acr_g_a_cv.npy', 'acr_g_d_cv.npy', 'acr_p_cv.csv', # genetic and phenotypic data for acr_cv
    'acr_g_a_st_sce.npy', 'acr_g_d_st_sce.npy', 'acr_p_st_sce.csv',  # genetic and phenotypic data for acr_st and acr_sce
    'acr_p.scl', # scaling factor for p_acr
    'acr_cv.json', 'acr_st.json', 'acr_sce.json' # respective cross validation schemes from R
]
elements_wtn = [
    'wtn_g_a.npy', 'wtn_g_d.npy', # genetic data: source additive and dominance
    'wtn_ec.pkl', 'wtn_ec.scl', # ec data and scaling factors
    'wtn_g_s.npy', 'wtn_g_s.scl', # g_s data and scaling factor
    'wtn_p.csv', 'wtn_p.scl', # phenodata and scaling factor
    'wtn_tra.json', 'wtn_LoO.json', 'wtn_cvL.json'  # cross validation scheme
]
elements_needed = elements_common + elements_acr + elements_wtn

out_paths = {x:f'{save_at}/{x}' for x in elements_needed}

# Define logger --------------------------------------------------------------------------
logging.basicConfig(filename=out_paths['preprocessing.log'], level=logging.DEBUG, filemode='w')
logger = logging.getLogger(__name__)
logging.info(f'Attempt at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}')

# Read data ------------------------------------------------------------------------------
p_data_acr = pd.read_feather(input_paths['p_data_BLUEs']).set_index('idx')
p_data_wtn = pd.read_feather(input_paths['p_data_raw']).set_index('idx')
g_data_a = pd.read_feather(input_paths['g_data_a']).set_index('idx')
g_data_d = pd.read_feather(input_paths['g_data_d']).set_index('idx')
e_data_raw = pd.read_feather(input_paths['e_data_raw']).set_index('idx')
ec_data_raw = pd.read_feather(input_paths['ec_data']).drop(['idx', 'env'], axis=1).rename(columns={"harvest_env": "idx"}).set_index('idx')
cgm_g_s_mat = pd.read_feather(input_paths['cgm_g_s_mat']).set_index('idx')

## cv scenarios
cv_acr_5f = read_json(input_paths['cv_acr_5f'])
#cv_acr_str = read_json(input_paths['cv_acr_str'])
cv_acr_sce = read_json(input_paths['cv_acr_sce'])
cv_wtn_tra = read_json(input_paths['cv_wtn_tra'])
cv_wtn_LoO = read_json(input_paths['cv_wtn_LoO'])
cv_wtn_cvL = read_json(input_paths['cv_wtn_cvL'])

#connect = pd.read_csv(input_paths['connect'], sep = " ")
#scenarios = pd.read_csv(input_paths['scenarios'])
#kws_old_new_cases = read_json(input_paths['kws_old_new'])
logging.info(f'Data read in')

# Scale data -----------------------------------------------------------------------------
## g data
g_a_data_scaled, g_a_scl = scale_data(g_data_a, g_data_a.columns, g_data_a.index)
g_d_data_scaled, g_d_scl = scale_data(g_data_d, g_data_d.columns, g_data_d.index)

## p data acr
p_data_mod, acr_p_scl = scale_data(p_data_acr.loc[:, "BLUEs"].values.reshape(-1, 1), 
                                   ["BLUEs"],
                                   p_data_acr.index)
acr_p = p_data_acr.merge(p_data_mod, how='left', left_index=True, right_index=True, 
                         sort=False, suffixes=('_raw', '_scaled'))

## p data wtn
p_data_wtn_scaled = p_data_wtn.dropna(subset = 'BLUES_dt') # drop rows where the column has missing values
wtn_p_scl = MinMaxScaler((0,1))
p_data_wtn_scaled['BLUEs_scaled'] = wtn_p_scl.fit_transform(p_data_wtn_scaled.loc[:, "BLUES_dt"].values.reshape(-1, 1))
wtn_p = p_data_wtn_scaled.reset_index(drop = True)

# g_s data
g_s_mat_scaled, g_s_mat_scl =  scale_data(cgm_g_s_mat, cgm_g_s_mat.columns, cgm_g_s_mat.index)

## ec_data
ec_data_scaled, wtn_ec_scl = scale_data(ec_data_raw, ec_data_raw.columns, ec_data_raw.index)

# Reshape data ---------------------------------------------------------------------------
## for acr_p_data: needs only genetic data
acr_p = acr_p.rename(columns = {"Geno_new" : "geno", "Type" : "type"})
acr_p["type"] = acr_p["type"].apply(lambda x: "Non_hybrid" if x != "Hybrid" else x)

### acr_cv
acr_p_cv = acr_p.drop(["Series", "idx_with_series"], axis = 1).drop_duplicates(subset = ["geno"]).rename(columns = {"unique_idx":"idx_col"}).reset_index(drop = True)
acr_p_cv["series"] = None

np_g_a_data_cv = np.stack([g_a_data_scaled[g_a_data_scaled.index == idx].iloc[0,:].values for idx in acr_p_cv['connect_geno_data']])
acr_g_a_cv = np_g_a_data_cv.reshape(np_g_a_data_cv.shape[0], np_g_a_data_cv.shape[1], 1)
np_g_d_data_cv = np.stack([g_d_data_scaled[g_d_data_scaled.index == idx].iloc[0,:].values for idx in acr_p_cv['connect_geno_data']])
acr_g_d_cv = np_g_d_data_cv.reshape(np_g_d_data_cv.shape[0], np_g_d_data_cv.shape[1], 1)

### acr_st_sce
acr_p_st_sce = acr_p.drop(["unique_idx"], axis = 1).drop_duplicates(subset = ["Series", "geno"]).rename(columns = {"idx_with_series":"idx_col", "Series":"series"}).reset_index(drop = True)

np_g_a_data_st_sce = np.stack([g_a_data_scaled[g_a_data_scaled.index == idx].iloc[0,:].values for idx in acr_p_st_sce['connect_geno_data']])
acr_g_a_st_sce = np_g_a_data_st_sce.reshape(np_g_a_data_st_sce.shape[0], np_g_a_data_st_sce.shape[1], 1)
np_g_d_data_st_sce = np.stack([g_d_data_scaled[g_d_data_scaled.index == idx].iloc[0,:].values for idx in acr_p_st_sce['connect_geno_data']])
acr_g_d_st_sce = np_g_d_data_st_sce.reshape(np_g_d_data_st_sce.shape[0], np_g_d_data_st_sce.shape[1], 1)

#acr_g_a_d = np.concatenate((acr_g_a, acr_g_d), axis=2)

## for wtn_p_data: genetic, g_s, and ec data
p_ec_data = np.stack([ec_data_scaled[ec_data_scaled.index == idx].iloc[0,:].values for idx in wtn_p['connect_climate']])
p_g_a_data = np.stack([g_a_data_scaled[g_a_data_scaled.index == idx].iloc[0,:].values for idx in wtn_p['connect_geno']])
p_g_d_data = np.stack([g_d_data_scaled[g_d_data_scaled.index == idx].iloc[0,:].values for idx in wtn_p['connect_geno']])
p_g_s_data = np.stack([g_s_mat_scaled[g_s_mat_scaled.index == idx].iloc[0,:].values for idx in wtn_p['connect_param']])

wtn_ec = p_ec_data.reshape(p_ec_data.shape[0], p_ec_data.shape[1], 1)
wtn_g_a = p_g_a_data.reshape(p_g_a_data.shape[0], p_g_a_data.shape[1], 1)
wtn_g_d = p_g_d_data.reshape(p_g_d_data.shape[0], p_g_d_data.shape[1], 1)
wtn_g_s = p_g_s_data.reshape(p_g_s_data.shape[0], p_g_s_data.shape[1], 1)

logging.info(f'Data preprocessing complete')

# Save data ------------------------------------------------------------------------------
## elements commom
write_pkl(g_a_scl, out_paths['g_a.scl'])
write_pkl(g_d_scl, out_paths['g_d.scl'])

## elements acr
np.save(out_paths['acr_g_a_cv.npy'], acr_g_a_cv)
acr_p_cv.to_csv(out_paths['acr_p_cv.csv'], index=False)
np.save(out_paths['acr_g_a_st_sce.npy'], acr_g_a_st_sce)
acr_p_st_sce.to_csv(out_paths['acr_p_st_sce.csv'], index=False)
write_pkl(acr_p_scl, out_paths['acr_p.scl'])
write_json(cv_acr_5f, out_paths["acr_cv.json"])
#write_json(cv_acr_str, out_paths["acr_st.json"])
write_json(cv_acr_sce, out_paths["acr_sce.json"])

#write_pkl(acr_g_d_scl, out_paths['acr_g_d.scl'])
#np.save(out_paths['acr_g_d.npy'], acr_g_d)
#np.save(out_paths['acr_g_a_d.npy'], acr_g_a_d)

## elements wtn
### genetic
np.save(out_paths['wtn_g_a.npy'], wtn_g_a)
np.save(out_paths['wtn_g_d.npy'], wtn_g_d) # their scalers will the same as are for acr data
### environment covariate
write_pkl(wtn_ec, out_paths['wtn_ec.pkl'])
write_pkl(wtn_ec_scl, out_paths['wtn_ec.scl'])
### geno_x_site
np.save(out_paths['wtn_g_s.npy'], wtn_g_s)
write_pkl(g_s_mat_scl, out_paths['wtn_g_s.scl'])
### phenotypic
wtn_p.to_csv(out_paths['wtn_p.csv'], index=False)
write_pkl(wtn_p_scl, out_paths['wtn_p.scl'])
### cv schema
write_json(cv_wtn_tra, out_paths["wtn_tra.json"])
write_json(cv_wtn_LoO, out_paths["wtn_LoO.json"])
write_json(cv_wtn_cvL, out_paths["wtn_cvL.json"])

### environment
#np.save(out_paths['wtn_e.npy'], wtn_e)
#write_pkl(wtn_e_scl, out_paths['wtn_e.scl'])

logging.info(f'Data saved')

# Finish off the script ------------------------------------------------------------------
write_json(out_paths, all_args[4])
logging.info(f'{task_name} completed successfully. Paths written at {all_args[4]}')
