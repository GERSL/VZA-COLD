function error_stack = find_stackerror(tile_name, year, doy, dir_data_raw, ...
    dir_data_brdf, dir_stack, row1, row2, nlines_read, varargin)
% This script will check for the storage of target data within the input
% period and tile specified, and create a .txt list of the missing/broken
% dates
%
% Inputs:
% n_prdct: the product number, 1=Raw Black Marble (VNP46A1);
%     2=BRDF-corrected Black_Marble (VNP46A2);
%     3=L2_SurfaceReflectance (VNP09GA).
% varargin('date_start'): start date in 'yyyyddd' format, default 2012019.
% varargin('date_end'): end date in 'yyyyddd' format, default 2021240.
% varargin('tile_list'): a string array list of tile name, with the format
%     'h*_v*'.
%
% Output:
% 'list_*.txt': the .txt record of error dates listed in 'yyyyddd' format,
%     which is named by the input start and end time
%
% Example: find_laads_download_errordate(2, 'h10v04')


%% Check if data available for the selected date & tile
v_prod = 2;
name_date = [num2str(year), sprintf('%03d', doy)];
% Read the .csv record for the raw data
name_csv = ['allData5000', 'VNP46A1', name_date, '.csv'];
if ~isfile(fullfile(dir_data_raw, name_csv))
    error_stack = [];
    return
end
rec_csv = readtable(fullfile(dir_data_raw, name_csv));
header_name = rec_csv.Properties.VariableNames;
if ~isequal(header_name{2}, 'name')
    % Change column name
    rec_csv.Properties.VariableNames{2} = 'name';
    rec_csv.Properties.VariableNames{3} = 'last_modified';
    rec_csv.Properties.VariableNames{4} = 'size';
    rec_csv.Properties.VariableNames{5} = 'checkfile';
    % Write to CSV file
    writetable(rec_csv, fullfile(dir_data_raw, name_csv))
end

tilelist_raw = regexpi(rec_csv.name, 'h(\w*)v(\w*)', 'match');
tilelist_raw = [tilelist_raw{:}];
tilelist_raw = string(tilelist_raw);
checkfile = rec_csv.checkfile;
% Get the successfully downloaded list
tilelist_raw = tilelist_raw(checkfile == 2);

% Read the .csv record for the Lunar-BRDF-corrected data
name_csv = ['allData5000', 'VNP46A2', name_date, '.csv'];
if ~isfile(fullfile(dir_data_brdf, name_csv))
    %     error_stack = [];
    %     return
    v_prod = 1;
    tilelist_brdf = "";
else
    rec_csv = readtable(fullfile(dir_data_brdf, name_csv));
    header_name = rec_csv.Properties.VariableNames;
    if ~isequal(header_name{2}, 'name')
        % Change column name
        rec_csv.Properties.VariableNames{2} = 'name';
        rec_csv.Properties.VariableNames{3} = 'last_modified';
        rec_csv.Properties.VariableNames{4} = 'size';
        rec_csv.Properties.VariableNames{5} = 'checkfile';
        % Write to CSV file
        writetable(rec_csv, fullfile(dir_data_brdf, name_csv))
    end
    tilelist_brdf = regexpi(rec_csv.name, 'h(\w*)v(\w*)', 'match');
    tilelist_brdf = [tilelist_brdf{:}];
    tilelist_brdf = string(tilelist_brdf);
    checkfile = rec_csv.checkfile;
    % Get the successfully downloaded list
    tilelist_brdf = tilelist_brdf(checkfile == 2);
end

% If no valid data available, skip
if ~ismember(tile_name, tilelist_raw) % ||
    error_stack = [];
    return
elseif ~ismember(tile_name, tilelist_brdf)
    v_prod = 1;
end

%% Check for the stack data for the tile list of target date
% The saving row array
irows = row1+nlines_read.*(0:ceil((row2-row1)/nlines_read)-1);
error_stack = 0;

% Checking for the last irows
n_row1 = irows(end);
n_row2 = n_row1+nlines_read-1;
n_rst = fullfile(dir_stack, sprintf('%04d%04d', n_row1, n_row2));
switch v_prod
    case 2
        filename_row = sprintf('map_tsdata_buffer_%04d%04d_%s', ...
            n_row1, n_row2, name_date);
    case 1
        filename_row = sprintf('map_tsdata_buffer_%04d%04d_%s_raw', ...
            n_row1, n_row2, name_date);
end

% Continue if exist
if isfile(fullfile(n_rst, [filename_row, '.mat']))
    % Delete the file if previously identified as error/no data
    var_error = {'errordata'};
    var_nodata = {'nodata'};
    var_nulldata = {'nulldata'};
    variableInfo = who('-file', fullfile(n_rst, [filename_row, '.mat']));
    if isequal(variableInfo, var_error)
        fprintf('Errordata file found for %s %s\n', tile_name, name_date);
        delete(fullfile(n_rst, [filename_row, '.mat']))
    elseif isequal(variableInfo, var_nodata)
        fprintf('Nodata file found for %s lines %04d%04d\n', tile_name, n_row1, n_row2);
        delete(fullfile(n_rst, [filename_row, '.mat']))
    elseif isequal(variableInfo, var_nulldata)
        movefile(fullfile(n_rst, [filename_row, '.mat']), fullfile(n_rst, [filename_row, '_null.mat']))
    else
        return
    end
end

filename_row_null = [filename_row, '_null.mat'];
if isfile(fullfile(n_rst, filename_row_null))
    return
end

error_stack = v_prod;
    
    
end