function Mosaic_initialradiancemap(varargin)
%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17);
addParameter(p, 'num_row', 2400);
addParameter(p, 'num_col', 2400);
addParameter(p, 'n_VZA', 4);

% request user's input
parse(p, varargin{:});
H_list = p.Results.H_list;
V_list = p.Results.V_list;
num_row = p.Results.num_row;
num_col = p.Results.num_col;
n_VZA = p.Results.n_VZA;

dir_l = '/shared/zhulab/Tian/Prod_09082023/';
dir_l = fullfile(dir_l, 'ChangeMetricMap_l_20132023', 'BlackMarble_BRDFcorrected_new');
dir_l_save = fullfile(dir_l, 'Mosaic_IniRad');
if ~isfolder(dir_l_save)
    mkdir(dir_l_save)
end

% Get the map list
n_rst = "InitialRadianceMap_%s";
fname = "%s_initialRadiance_VZAint%d.tif";
nrows = num_row*length(V_list);
ncols = num_col*length(H_list);

%% Get the total task number, task input parameters and allocated the tasks to the jobs
% Initialize the maps
GlobalMap = zeros(nrows, ncols, 'double');
for i_H = H_list
    for i_V = V_list(3:end-3) % Only include +70 to -60 degree latitude
        tile_name = sprintf('h%02dv%02d', i_H, i_V);
        dir_data = fullfile(dir_l, tile_name, sprintf(n_rst, tile_name));
        fname_data = sprintf(fname, tile_name, 4);
        path_data = fullfile(dir_data, fname_data);
        if ~isfile(path_data)
            continue
        end
        row_list = i_V*num_row+1:(i_V+1)*num_row;
        col_list = i_H*num_col+1:(i_H+1)*num_col;
        [tmp, ~] = readgeoraster(path_data);

        GlobalMap(row_list, col_list) = tmp(:, :, n_VZA);
    end
end

R = generate_globaltif_Ref('save_code', 0);
% Save the global map
fname_save = sprintf(fname, 'Global', n_VZA);
path_save = fullfile(dir_l_save, fname_save);
geotiffwrite(path_save, GlobalMap, R, TiffType="bigtiff");
end