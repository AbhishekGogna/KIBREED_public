#!/usr/bin/env python3
# load functions -------------------------------------------------------------------------
import sys
all_args = sys.argv[1:]
source_code_at = all_args[0]
sys.path.append(source_code_at)
from Py.libs import *
from Py.func import *

# Define paths ----------------------------------------------------------------------------------------------------------
input_paths = read_json(all_args[1])
save_at = all_args[2]
if not os.path.exists(save_at):
    os.makedirs(save_at, exist_ok = True)
task_name = all_args[3]
out_paths = input_paths # so that more keys-value pairs are appneded to an existing dictionary. should remove it and set out_paths to a dictionary if connection to the next step is needed
out_paths['base_folder'] = f'{save_at}'
out_paths['pred_dirs_paths'] = f'{save_at}/pred_dirs_paths.json'
out_paths['master_script_path'] = f'{save_at}/master_script_paths.txt'
parent_dir = "/qg-10/data/AGR-QG/Gogna"
project_name = "KIBREED_v2"

# Produce data ---------------------------------------------------------------------------------------------------------------------
chunk_1 = '''#!/usr/bin/env bash
# define variables
base_dir="/qg-10/data/AGR-QG/Gogna"
cont_dir="${base_dir}/computing_containers/containers/cc_jup_rst_v2.sif"
source_dir="%s" ######## input data
proj_dir="%s" ######## project_dir
slurm_id_at="%s" ######## file to save slurm id
instancename="%s" ######## instance name
python_env="${source_dir}/py_env"

# write slurm id for tracking
echo -e "${instancename}\t${SLURM_JOB_ID}\t$(hostname)" >> "${slurm_id_at}"

# run instance
singularity instance start --nv \\
        -W "/proj" \\
        -H "${proj_dir}:/proj" \\
        -B "%s:/proj/%s" \\
        -B "${python_env}:/proj/py_env" \\
        -B "${source_dir}:/proj/ext_dir" \\
        -B "/opt/Bio/cuda-toolkit/11.6/bin:/usr/local/cuda-11.6/bin" \\
        "${cont_dir}" "${instancename}"

# predictions
singularity exec --pwd "/proj" instance://"${instancename}" bash -c \\
        ""%s" "%s" "%s" "%s""

# stop instance
singularity instance stop "${instancename}"
exit'''

slurm_jobs = [("acr_CNN", "acr_cv.json"),
              #("acr_CNN", "acr_st.json"),
              ("acr_CNN", "acr_sce.json"),
              #("acr_CNN_BN", "acr_cv.json"),
              #("acr_CNN_BN", "acr_st.json"),
              #("acr_CNN_BN", "acr_sce.json"),
              ("wtn_CNN_EC", "wtn_tra.json"), #EC from climate data: run with cv1 to cv4
              ("wtn_CNN_CGM_EC", "wtn_tra.json"), #EC from CGM data: run with cv1 alone
              #("wtn_CNN_LSTM", "wtn_cv.json")
              ] # one line about what model to run and with what cv data
slurm_jobs_df = pd.DataFrame(slurm_jobs, columns=["model", "json_file"])
slurm_jobs_df.drop_duplicates(inplace=True)

# model inputs
model_inputs = {"wtn_CNN_EC": {"tune": False},
                "wtn_CNN_CGM_EC": {"tune": False},
                #"wtn_CNN_LSTM": {"tune": False},
                #"wtn_CNN_explicit": {"tune": False},
                "acr_CNN": {"tune": False},
                "acr_CNN_BN": {"tune": False}} # add more here

## delete master script path file if it exists
if os.path.exists(out_paths['master_script_path']):
    os.remove(out_paths['master_script_path'])

out_dirs = {}
for jobs in slurm_jobs_df.index:
    model = slurm_jobs_df.loc[jobs, "model"]
    json_file = slurm_jobs_df.loc[jobs, "json_file"]
    json_file_path = f'{input_paths[json_file]}'
    mem = str(20) if "acr" in model else str(100)
    #node_list = "" if "acr" in model else "slurm-gpu-01"
    node_list = "slurm-gpu-01"
    
    if os.path.isfile(json_file_path):
        data = read_json(json_file_path)
    else:
        sys.exit(f'No data found for type {jobs}')
    if model == "wtn_CNN_CGM_EC":
        to_remove = [x for x in list(data.keys()) if any(key in x for key in ["cv2", "cv3", "cv4"])]
        data = {key: value for key, value in data.items() if key not in to_remove}
    
    time_stamp = f'{get_random_string(10)}_{datetime.now().strftime("D_%d_%m_T_%H_%M")}'
    for run in range(len(data)):
        key = list(data.keys())[run]
        job_name = f'{model}@{json_file.replace(".json", "")}@{key}'
        # relative paths
        run_name = json_file.replace(".json", "")
        base_dir = f'{save_at}/{model}_{run_name}'
        sub_dirs = set_dirs(base_dir, run_id = key, verbose = False) # creates a series of folders in each directory
        out_dirs[f'{model}#{run_name}#{key}'] = sub_dirs
        master_script_path = f'{base_dir}/master_script.sh' # puts a script in the previous folder
        run_dir = f'{base_dir}/{key}'
        script_name = f'{run_dir}/run_script_{run}'
        script_path = f'{script_name}.sh' # puts a script in the previous folder    
        ext_hp_tuning = f'/qg-10/data/AGR-QG/temp/Gogna/KIBREED_v2/{model}_{run_name}/{key}'
        int_hp_tuning = ext_hp_tuning.replace("/qg-10/data/AGR-QG/temp/Gogna/KIBREED_v2", "/proj/tmp_data")
        
        if not os.path.exists(int_hp_tuning):
            os.makedirs(int_hp_tuning)
        #print(int_hp_tuning)
        
        # write model input file
        write_json(model_inputs, f'{run_dir}/model_args.json')
        
        # absolute paths for run_script
        project_dir_abs = f'{parent_dir}/{project_name}'
        
        # absolute path for slurm_id_file
        slurm_id_at = f'{project_dir_abs}/{base_dir[6:]}/{time_stamp}_srun_id.log'

        # absolute paths for run scripts and log as well as err file
        script_path_abs = f'{project_dir_abs}/{script_name[6:]}.sh'
        stdout = f'{project_dir_abs}/{script_name[6:]}.log'
        stderr =  f'{project_dir_abs}/{script_name[6:]}.err'        
        
        # parameters
        slurm_time = '1-0:0:0'
        proj_dir = f'{parent_dir}/{project_name}/{run_dir[6:]}'
        #hparams = f'{parent_dir}/{project_name}/{input_paths["best_param_rndm"][6:]}'
        pred_script = '/proj/ext_dir/run/Py/predictions.py'
        
        # write run script
        with open (script_path, 'w') as rsh:
            rsh.write(chunk_1 % (f'{project_dir_abs}',f'{proj_dir}', f'{slurm_id_at}' , f'ins_{key}', f'{ext_hp_tuning}', f'{model}', f'{pred_script}', f'{model}', f'{json_file}', f'{key}'))
        os.system(f'chmod +x {script_path}') 
        
        if run == 0:
            with open (master_script_path, 'w') as rsh:
                rsh.write('#!/usr/bin/env bash \n')
                rsh.write(f'if [ -f {slurm_id_at} ] ; then rm {slurm_id_at} ; fi \n')
                rsh.write(f'sbatch --auks=yes --job-name={job_name} --mem={mem}G -c 1 -p gpu -x "{node_list}" --gres=gpu:1 --time={slurm_time} --mail-type=FAIL --mail-user=gogna@ipk-gatersleben.de -o {stdout} -e {stderr} --wrap="{script_path_abs}" & \n')
        elif (run >0) and (run < (len(data)-1)):
            with open (master_script_path, 'a+') as rsh:
                rsh.write(f'sbatch --auks=yes --job-name={job_name} --mem={mem}G -c 1 -p gpu -x "{node_list}" --gres=gpu:1 --time={slurm_time} --mail-type=FAIL --mail-user=gogna@ipk-gatersleben.de -o {stdout} -e {stderr} --wrap="{script_path_abs}" & \n')
        elif run == (len(data)-1):
            with open (master_script_path, 'a+') as rsh:
                rsh.write(f'sbatch --auks=yes --job-name={job_name} --mem={mem}G -c 1 -p gpu -x "{node_list}" --gres=gpu:1 --time={slurm_time} --mail-type=FAIL --mail-user=gogna@ipk-gatersleben.de -o {stdout} -e {stderr} --wrap="{script_path_abs}" & \nwait \necho "All scripts Done!"')
        os.system(f'chmod +x {master_script_path}')
        
    ## write master script paths at save_at
    with open (out_paths['master_script_path'], 'a') as rsh:
        rsh.write(f'{project_dir_abs}/{master_script_path[6:]}\t{len(data)}\n')

# Produce output --------------------------------------------------------------------------------------------------
write_json(out_dirs, out_paths['pred_dirs_paths'])
# Finish off the script --------------------------------------------------------------------------------------------------
write_json(out_paths, f'{all_args[4]}')
print(f'{task_name} completed successfully. Paths written at {all_args[4]}')

