function initialradiance(dir_BRDF, dir_ref, tile_name, ...
    n_rst_data, n_rst_save, date_start, row1, row2, col1, col2, varargin)

warning('off', 'all')

%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'num_vza', 4); 

% request user's input
parse(p, varargin{:});
num_vza = p.Results.num_vza;

yr_s = floor(date_start/1000);
doy_s = mod(date_start, 1000);
date_start = datenum(yr_s, 1, doy_s);

% Initilize the gradual change record
InitialRadiance = zeros(row2-row1+1, col2-col1+1, num_vza);
for irow_ids =  row1:row2 
    %% Read and merge mapped changes.
    filename_rec = sprintf('%s_changemetric%04d.mat', tile_name, irow_ids);
    filename_rec_null = sprintf('%s_changemetric%04d_null.mat', tile_name, irow_ids);
    % For the individual VZA stratified COLD methods
    if isfile(fullfile(dir_BRDF, n_rst_data, filename_rec))
        load(fullfile(dir_BRDF, n_rst_data, filename_rec), 'rec_cg');
    elseif isfile(fullfile(dir_BRDF, n_rst_data, filename_rec_null))
        continue
    end

    %% Calculate the overall values
    check_col = 0;
    for i = 1:length(rec_cg)
        icol_ids = rec_cg(i).pos;
        if icol_ids == check_col
            continue
        end

        num_vza = size(rec_cg(i).coefs, 2);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
        if num_vza == 1
            InitialRadiance(irow_ids, icol_ids, end) = 0.1*(rec_cg(i).coefs(1, :)+rec_cg(i).coefs(2, :)*date_start);
        else
            InitialRadiance(irow_ids, icol_ids, :) = 0.1*(rec_cg(i).coefs(1, :)+rec_cg(i).coefs(2, :)*date_start);
        end
        check_col = icol_ids;
    end
end 
InitialRadiance(isnan(InitialRadiance)) = 0;

% Save the record
[~, R] = readgeoraster(fullfile(dir_ref, sprintf('%s.tif', tile_name)));
for n_vza = 1:num_vza
    fname_save = sprintf('%s_initialRadiance_VZAint%d.tif', tile_name, n_vza);
    geotiffwrite(fullfile(dir_BRDF, n_rst_save, fname_save), InitialRadiance, R)
end
end