import os, sys
import dpFunction as dp


# input_path = /gpfs/alpine/cli144/proj-shared/NA_dataset
# output_path = /gpfs/alpine/cli144/scratch/wangd/subdomains


def get_files(input_path):
    print(input_path)
    files = os.listdir(input_path) 

    files.sort() 

    file_no =0

    files_nc = [f for f in files if (f[-2:] == 'nc')] 
    print("total " + str(len(files_nc)) + " files need to be processed")
    return files_nc

    
def main():
    args = sys.argv[1:]
    input_path = args[0]
    output_path = args[1]
    number_of_subdomains = int(args[2])
    i_timesteps = int(args[3])
    AOI_mask = args[4]
    
    files_nc = get_files(input_path)

    for f in files_nc: 
        var_name = f[20:-11]
        period = f[-10:-3]
        print('processing '+ var_name + '(' + period + ') in the file ' + f + 'in the folder ' + input_path)
        dp.launch_job(input_path, f, output_path, var_name, period, number_of_subdomains, i_timesteps, AOI_mask)
        print('result in '+ output_path  + '(' + str(number_of_subdomains) + ') with timesteps (' + str(i_timesteps) + ')')

if __name__ == '__main__':
    main()

