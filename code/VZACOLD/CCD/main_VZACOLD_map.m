function rec_cg = main_VZACOLD_map(dir_BRDF_new, tile_name, irow, ...
    n_irows, nlines_read, date_start, date_end, varargin)
% Matlab code for Continous Change Detection (for standalone version)
% Inputs:
% conse:  number of consecutive observations

%% get parameters from inputs
p = inputParser;

% default values.
addParameter(p, 'col1', 1); % Dimension of the map
addParameter(p, 'col2', 2400); 
addParameter(p, 'save_code', 1); 
addParameter(p, 'vzaints', [0 20; 20 40; 40 60]);
addParameter(p, 'update_model', 0.02); 
addParameter(p, 'Tmax_cg', 1); % maximum number of coefficients
addParameter(p, 'T_const', 1); % Threshold for cloud, shadow, and snow outliers.
addParameter(p, 'T_cg', 0.75); % change threshold
addParameter(p, 'conse', 14); 
addParameter(p, 'num_toler', 1); % noise-tolerant number %floor(conse*0.2); % assume 20% comparing observaitons are noises
addParameter(p, 'num_c', 4); 
addParameter(p, 'NTL_lmt', 1); % Threshold for consistent dark pixels (that will not be runned)
addParameter(p, 'run_dir', 1); % 1=ascending; 2=descending

% request user's input
parse(p, varargin{:});
col1 = p.Results.col1;
col2 = p.Results.col2;
save_code = p.Results.save_code;
vzaints = p.Results.vzaints;
update_model = p.Results.update_model;
Tmax_cg = p.Results.Tmax_cg; % Treshold of noise % Tmax_cg = 1-1e-5;
T_cg = p.Results.T_cg;
T_const = p.Results.T_const;
conse = p.Results.conse;
num_toler = p.Results.num_toler;
num_c = p.Results.num_c;
NTL_lmt = p.Results.NTL_lmt;
run_dir = p.Results.run_dir;

% folder name of all CCDC results
cmd_l = '/home/til19015/GlobalNTLAnalyze/VZACOLD_master/Stable/Map';
addpath(cmd_l);  
n_rst = sprintf('TSFitMap_%s', tile_name);

% make TSFitMap folder for storing coefficients
if ~isfolder(fullfile(dir_BRDF_new, n_rst))
    fprintf('Generating directory: %s\n', fullfile(dir_BRDF_new, n_rst));
    mkdir(fullfile(dir_BRDF_new, n_rst));
end
fprintf('!! Processing row%d of the %s: CI = %f, conse = %d\n', n_irows, ...
    tile_name, T_cg, conse);

rec_cg = TrendSeasonalFit_map_VZACOLD( ...
    dir_BRDF_new, n_rst, nlines_read, irow, n_irows, col1, ...
    col2, tile_name, T_cg, Tmax_cg, conse, num_c, vzaints, ...
    date_start, date_end, num_toler, update_model, NTL_lmt);

end % end of function 