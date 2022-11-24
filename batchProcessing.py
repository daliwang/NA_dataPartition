import os, sys, shlex, subprocess
#import dpFunction

#  make sure the path ends with "/"
input_path = "/gpfs/alpine/cli144/proj-shared/NA_dataset/"
output_path = "/gpfs/alpine/cli144/scratch/wangd/subdomains/"
test = 1
fullrun = 0 
AOI_mask = "/gpfs/alpine/cli144/proj-shared/wangd/Python4data/AOI_mask/AOI_mask_AK.nc"

def get_dirs(input_path):
    print(input_path)
    files = os.listdir(input_path)

    dirs_no =0
    files.sort()

    dirs = [f for f in files if (os.path.isdir(os.path.join(input_path,f)))]

    print("total " + str(len(dirs)) + " direcories need to be processed")
    return dirs


def main():
    number_of_subdomains = 5 
    if fullrun:
        i_timesteps = -1
    else:
        i_timesteps = 1

    dirs = get_dirs(input_path)
    print(dirs)
    if test:
        dirs = ["2000"]
    for d in dirs:
        print('processing folder ' + d + ' in ' + input_path)
        #dp.launch_job(input_path, f, output_path, var_name, period, number_of_subdomains, i_timesteps)
        command = "python3 DataPartition.py " + input_path  + d + " " + output_path + " " +str(number_of_subdomains) + " " + str(i_timesteps) + " " + AOI_mask
        print(command)
        args= shlex.split(command)
        p_output = subprocess.Popen(args)
        print('result in '+ output_path  + '(' + str(number_of_subdomains) + ') with timesteps (' + str(i_timesteps) + ')')
        ls_output=subprocess.Popen(["sleep", "60"])

if __name__ == '__main__':
    main()
