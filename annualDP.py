
import os
import subprocess
import data_partition as dp
from time import process_time
#ls_output=subprocess.Popen(["sleep", "30"])

number_of_subdomains = 4200  # 4200 for 700 Summit nodes
i_timesteps = -1   # -1 means the full dataset
input_path= '/home/7xw/data/GWSP3_DayMet/2014NA/'
#input_path = os.getcwd()
print(input_path)
files = os.listdir(input_path) 

files.sort() 
output_path = input_path + 'subdomain_batch/'

file_no =0

files_xls = [f for f in files if (f[-2:] == 'nc')] 
print("total " + str(len(files_xls)) + " files need to be processed")
print('number of subdomains: '+ str(number_of_subdomains))


for f in files_xls: 
    var_name = f[20:-11]
    period = f[-10:-3]
    print('processing '+ var_name +' in the file ' + f )
    dp.launch_job(input_path, f, output_path, var_name, period, number_of_subdomains, i_timesteps)

