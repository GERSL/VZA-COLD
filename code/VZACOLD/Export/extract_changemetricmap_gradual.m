function extract_changemetricmap_gradual(var, fname, yr_s, row1, row2, col1, col2, ...
    dir_rowchg, dir_chgmap, dir_tiffref, tile_name, yr, num_vza, replace_code)

%% Get the directories
% For the change metric records
n_rst_rowchg = sprintf('ChangeMetricMap_%s', tile_name);
% For gradual change
n_rst_dailygrdchg = sprintf('Daily_GradualChangeMap_%s', tile_name);
if ~isfolder(fullfile(dir_chgmap, n_rst_dailygrdchg))
    mkdir(fullfile(dir_chgmap, n_rst_dailygrdchg))
end

%% Generate themaps
% Initialize the daily maps
nrows = row2-row1+1;
ncols = col2-col1+1;
n_yr = yr-yr_s+2;
sdate = datenum(yr, 1, 1):datenum(yr, 12, 31);
num_dates = length(sdate);

% For gradual change
Daily_grdmag = zeros(nrows, ncols, num_dates, 'double'); % gradual change magnitude
for irow_ids = row1:row2
    %% Read in and summerize changemetric
    filename_data = sprintf('%s_changemetric%04d.mat', tile_name, irow_ids);
    pathdata = fullfile(dir_rowchg, n_rst_rowchg, filename_data);
    if ~isfile(pathdata)
        continue
    end

    if contains(var, 'net')
        load(pathdata, 'GradualMag_net_row');
        GradualMag_row = GradualMag_net_row;
    elseif contains(var, 'pos')
        load(pathdata, 'GradualMag_pos_row');
        GradualMag_row = GradualMag_pos_row;
    elseif contains(var, 'neg')
        load(pathdata, 'GradualMag_neg_row');
        GradualMag_row = GradualMag_neg_row;
    end
    Daily_grdmag(irow_ids, :, :) = GradualMag_row(:, :, n_yr, 4, 1:num_dates);
end
Daily_grdmag(isnan(Daily_grdmag)) = 0;

%% Save and extract the maps
% Get the georeference
[~, R] = readgeoraster(fullfile(dir_tiffref, sprintf('%s.tif', tile_name)));
% For the annual records
for i_date = 1:num_dates
    file_datenum = sdate(i_date);
    file_doy = day(datetime(file_datenum, 'ConvertFrom', 'datenum'), 'dayofyear');
    fname_save = sprintf(fname, tile_name, yr, file_doy, num_vza);
    filename_save_tif = sprintf('%s.tif', fname_save);
    % Skip if existed
    if isfile(fullfile(dir_chgmap, n_rst_dailygrdchg, filename_save_tif))
        continue
    end

    %% For gradual
    % Save the daily net/positive/negative magnitude map
    tmp = Daily_grdmag(:, :, i_date);
    if sum(tmp(:) ~= 0) == 0
        filename_save_tif = sprintf('%s_null.mat', fname_save);
    end
    clear tmp
    path_save = fullfile(dir_chgmap, n_rst_dailygrdchg, filename_save_tif);
    SaveGeotiff(replace_code, path_save, Daily_grdmag(:, :, i_date), R);
end
clear Daily_grdmag
end

function SaveGeotiff(replace_code, path_save, data, R)
if isfile(path_save)
    if replace_code == 0
        return
    else
        delete(path_save)
    end  
end

if contains(path_save, 'null')
    null_data = [];
    save(path_save, 'null_data')
else
    geotiffwrite(path_save, data, R)
end
end