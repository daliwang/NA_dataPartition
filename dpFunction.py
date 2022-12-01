# code for whole NA domain using Kao's dataset
# 4th version with array_split/memory profiling
import os 
import netCDF4 as nc
import numpy as np
from itertools import cycle
#from pprint import pprint
from time import process_time
from memory_profiler import profile

def data_read(file_name, var_name, times):
    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    r_nc_fid = nc.Dataset(file_name, 'r', format='NETCDF4')

    total_cols = r_nc_fid.dimensions['x'].size
    total_rows = r_nc_fid.dimensions['y'].size
    total_timesteps = r_nc_fid.dimensions['time'].size

    data = r_nc_fid[var_name][0:times, :, :] # read (timestep, y, x) format
    return total_rows, total_cols, data

def gridId_data(total_rows, total_cols, times, data, AOI_mask):
    total_gridcells = total_rows * total_cols
    grid_ids = np.linspace(0, total_gridcells-1, total_gridcells, dtype=int)

    # create a mask for land grid_ids (1)
    if (AOI_mask =="NA") :
        mask = data[0]    # FSDS is in (time, Y, X) format
        mask = np.where(~np.isnan(mask), 1, np.nan)
    else:
        r_nc_fid = nc.Dataset(AOI_mask, 'r', format='NETCDF4')
        mask = r_nc_fid['Mask'][:, :] # read mask(y, x) format
        mask = np.where((mask==30), 1, np.nan)

    # create an flattened list of land gridID and reduce the size of gridIDs array
    grid_ids = grid_ids.reshape(total_rows,total_cols)
    grid_ids = np.multiply(mask,grid_ids)
    grid_ids = grid_ids[~np.isnan(grid_ids)]

    # use the size of land gridcells to resize the FSDS matrix

    landcells = len(grid_ids)    
    #data1 = np.empty([times, landcells],dtype = float)
    data = np.multiply(mask, data)
    data = data[~np.isnan(data)]
    data = np.reshape(data,(times,landcells))
    return grid_ids, data

def data_partition_RR(number_of_subdomains, grid_ids, data):
    debug = 0
    # cyclic (round-robin) partition
    domains = [[] for _ in range(number_of_subdomains)]
    for element, domain in zip(grid_ids, cycle(domains)):
        domain.append(element)

    grid_id_domains = domains.copy()
    
    landcells = len(grid_ids) 
    # partition the data over landcells
    # landcell_idx is alse the column_idx of FSDS
    landcell_idx = np.linspace(0, landcells-1, landcells, dtype=int)

    domains = [[] for _ in range(number_of_subdomains)]
    for element, domain in zip(landcell_idx, cycle(domains)):
        domain.append(element)
    
    # save the boundaries of each subdomain (for array_split)
    size_of_subdomains = [ len(domain) for domain in domains]

    # partitioned landcells_idx in subdomains 
    arranged_grid_idx = np.concatenate(domains).ravel()
    if debug:
        print('the frist 5 grid_idx are: '+ str(arranged_grid_idx[0:5]))

    # find the original index of landcells for column swap
    np.sort(arranged_grid_idx)
    grid_swap_idx = (np.argsort(arranged_grid_idx))

    if debug:
        print('the frist 5 grid_swap_idx are: '+ str(grid_swap_idx[0:5]))

    # create swap index and arrange data
    idx = np.empty_like(grid_swap_idx)
    idx[grid_swap_idx] = np.arange(len(grid_swap_idx))
    data = data[:,idx]

    # split the FSDS into subdomains using the boundary index
    subdomain_idx = size_of_subdomains
    for i in range(0,len(subdomain_idx)-1):
        subdomain_idx[i+1] +=subdomain_idx[i]
    
    #FSDS_list = np.hsplit(FSDS,subdomain_idx[:-1])
    data = np.hsplit(data,subdomain_idx[:-1])
    
    return grid_id_domains, data

def data_save(number_of_subdomains, grid_id_domains, subdomain_path, i_timesteps, var_name, period, data):
    for i in range(number_of_subdomains):
        # convert local grid_id_lists into an array
        grid_id_arr = np.array(grid_id_domains[i])

        #data_arr = np.array(FSDS_list[i])
        data_arr = np.array(data[i])
        file_name = subdomain_path + var_name + '_' + period +'_subdomain'+ str(i) + '.nc'

        # Open a new NetCDF file to write the data to. For format, you can choose from
        # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
        w_nc_fid = nc.Dataset(file_name, 'w', format='NETCDF4')
        w_nc_fid.title = 'The ELM domain files on individudal process: '+str(i)

        # create the gridIDs variable
        x_dim = w_nc_fid.createDimension('x_dim', grid_id_arr.size)
        time_dim = w_nc_fid.createDimension('time_dim', i_timesteps)
        w_nc_var = w_nc_fid.createVariable('gridIDs', np.int32, ('x_dim',))
        w_nc_var.long_name = 'gridIds in the subdomain'    
        w_nc_fid.variables['gridIDs'][:] = grid_id_arr.reshape(grid_id_arr.size)

        w_nc_var = w_nc_fid.createVariable(var_name, np.float32, ('time_dim', 'x_dim'))
        w_nc_var.long_name = 'FSDS in the subdomain'    
        w_nc_fid.variables[var_name][:] =data_arr.reshape(i_timesteps,grid_id_arr.size)
        w_nc_fid.close()  # close the new file

def launch_job(input_path, file_name, output_path, var_name, period, \
               number_of_subdomains, i_timesteps, AOI_mask):
    file_name = input_path + "/" + file_name  
    start = process_time()
    [total_rows, total_cols, data] = data_read(file_name, var_name, i_timesteps)
    end = process_time()
    print("Reading " + var_name + " takes  {}".format(end-start))
    
    start = process_time()  
    [grid_ids, data] = gridId_data(total_rows, total_cols, i_timesteps, data, AOI_mask)
    end = process_time()
    print("Creating dense " + var_name + " takes  {}".format(end-start))

    start = process_time()
    [grid_id_domains, data] = data_partition_RR(number_of_subdomains, grid_ids, data)
    end = process_time()
    print("Partitioning data/GridID takes  {}".format(end-start))
    
    start = process_time() 
    data_save(number_of_subdomains, grid_id_domains, output_path, i_timesteps, var_name, period, data)
    end = process_time()
    print("Saving data/GridID takes  {}".format(end-start))
          



