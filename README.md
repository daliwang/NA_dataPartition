# NA_dataPartition


batchProccessing.py    batch processing all the data with pipeline
DataPartion.py             scripts to launch actual python jobs
dpFunction.py    python function read/partition/save netcdf files
 
Usage:   
 
batch processing 

       python3 batchProcessing.py 
       
command line:

      python3 DataPartition.py <input_path> <output_path> <number_of_subdomains> <i_timesteps>
 
i_timestep = -1     all the timesteps in a file
i_timestep >= 1    process few timesteps of data in a file
