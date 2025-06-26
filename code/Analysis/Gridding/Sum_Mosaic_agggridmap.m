function Sum_Mosaic_agggridmap(varargin)

%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17); % +70 to -60degree latitude
addParameter(p, 'num_row', 2400);
addParameter(p, 'num_col', 2400);
addParameter(p, 'yr_list', 2014:2022);
addParameter(p, 'grid_size', 120); %0.05/0.05-degree

% request user's input
parse(p, varargin{:});
H_list = p.Results.H_list;
V_list = p.Results.V_list;
num_row = p.Results.num_row;
num_col = p.Results.num_col;
yr_list = p.Results.yr_list;
grid_size = p.Results.grid_size;

dir_l = '/shared/zhulab/Tian/Prod_09082023/ChangeMetricMap_l_20132023/Analysis/Grid_AreaAdj';
dir_l_data = fullfile(dir_l, sprintf('%dx%d', grid_size, grid_size), 'Mosaic');
var_list = [...
    "ttlarea_pos", ...
    "ttlarea_neg", ...
    "ttlmag_pos", ...
    "ttlmag_neg", ...
    "ttlint_pos", ...
    "ttlint_neg", ...
    ];

% Get the georeference information for the global map
R = generate_globaltif_Ref('save_code', 0);
num_row_agggrid = num_row/grid_size;
num_col_agggrid = num_col/grid_size;
R.RasterSize = [num_row_agggrid*length(V_list), num_col_agggrid*length(H_list)];

nrows = num_row_agggrid*length(V_list);
ncols = num_col_agggrid*length(H_list);
for var = var_list
    fname_save = sprintf('GlobalAggGrid_%s_%d.tif', var, yr_list(1)*10000+yr_list(end));
    path_save = fullfile(dir_l_data, fname_save);

    if contains(var, 'int')
        % Initialize the maps
        SumGlobalMap = zeros(nrows, ncols, length(yr_list), 'double');
        AreaGlobalMap = zeros(nrows, ncols, length(yr_list), 'double');
        for i_yr = 1:length(yr_list)
            yr = yr_list(i_yr);

            fname_data = sprintf('GlobalAggGrid_%s_%d_mean.tif', var, yr);
            [tmp, ~] = readgeoraster(fullfile(dir_l_data, fname_data));
            SumGlobalMap(:, :, i_yr) = tmp;
            
            switch var
                case "ttlint_pos"
                    fname_data = sprintf('GlobalAggGrid_ttlarea_pos_%d_mean.tif', yr);
                case "ttlint_neg"
                    fname_data = sprintf('GlobalAggGrid_ttlarea_neg_%d_mean.tif', yr);
            end
            [tmp, ~] = readgeoraster(fullfile(dir_l_data, fname_data));
            AreaGlobalMap(:, :, i_yr) = tmp;

        end

        AreaGlobalMap = AreaGlobalMap./sum(AreaGlobalMap, 3);
        AreaGlobalMap(isnan(AreaGlobalMap)) = 0;
        SumGlobalMap = SumGlobalMap.*AreaGlobalMap;
        SumGlobalMap = sum(SumGlobalMap, 3);
    else
        % Initialize the maps
        SumGlobalMap = zeros(nrows, ncols, 'double');
        for yr = yr_list
            fname_data = sprintf('GlobalAggGrid_%s_%d.tif', var, yr);
            [tmp, ~] = readgeoraster(fullfile(dir_l_data, fname_data));
            SumGlobalMap = SumGlobalMap+tmp;
        end
    end

    % Save the global map
    geotiffwrite(path_save, SumGlobalMap, R);
end
end