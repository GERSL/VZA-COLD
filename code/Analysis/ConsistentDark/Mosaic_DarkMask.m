function Mosaic_DarkMask(varargin)

%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17);
addParameter(p, 'num_row', 2400);
addParameter(p, 'num_col', 2400);

% request user's input
parse(p, varargin{:});
H_list = p.Results.H_list;
V_list = p.Results.V_list;
num_row = p.Results.num_row;
num_col = p.Results.num_col;

% Get the georeference information for the global map
path_tiffref = '/shared/zhulab/Tian/VIIRS_NTL/Tiff_ref/global_mosaic.tif'; % Arcpy version
[~, R] = readgeoraster(path_tiffref);
dir_l = '/shared/zhulab/Tian/Prod_09082023/ChangeMetricMap_l_20132023/';
dir_l = fullfile(dir_l, '/Analysis/DarkPixel/DarkMask');
% Initialize the maps
nrows = num_row*length(V_list);
ncols = num_col*length(H_list);
GlobalMap = zeros(nrows, ncols, 'uint8');
for i_H = H_list
    for i_V = V_list(3:end-3)
        tile_name = sprintf('h%02dv%02d', i_H, i_V);
        row_list = i_V*num_row+1:(i_V+1)*num_row;
        col_list = i_H*num_col+1:(i_H+1)*num_col;

        fname_data = sprintf('DarkMask_%s.tif', tile_name);
        [tmp, ~] = readgeoraster(fullfile(dir_l, fname_data));
        GlobalMap(row_list, col_list) = tmp;
    end
end
fname_save = 'DarkMask_Global.tif';
path_save = fullfile(dir_l, fname_save);
% Save the global map
geotiffwrite(path_save, GlobalMap, R);
end