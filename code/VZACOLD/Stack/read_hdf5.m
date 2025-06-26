function data = read_hdf5(file_id,data_name,dims,startpxl)
mem_space_id = H5S.create_simple(2, dims, []);
dset_id = H5D.open(file_id, data_name);
file_space_id = H5D.get_space(dset_id);
H5S.select_hyperslab(file_space_id, 'H5S_SELECT_SET', startpxl, [], [], dims);
data = (H5D.read(dset_id, 'H5ML_DEFAULT', mem_space_id, file_space_id, 'H5P_DEFAULT'))';

H5D.close(dset_id)
end