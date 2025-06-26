function loop_extractchangemetricmap_gradual(job, njobs, varargin)

%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'replace_code', 0); 
addParameter(p, 'row1', 1);
addParameter(p, 'row2', 2400);
addParameter(p, 'col1', 1);
addParameter(p, 'col2', 2400);
addParameter(p, 'yr_s', 2014);
addParameter(p, 'yr_e', 2022);
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17);
addParameter(p, 'num_vza', 4); 
addParameter(p, 'check_task', 0);  

% request user's input
parse(p, varargin{:});
replace_code = p.Results.replace_code;
row1 = p.Results.row1;
row2 = p.Results.row2;
col1 = p.Results.col1;
col2 = p.Results.col2;
yr_s = p.Results.yr_s;
yr_e = p.Results.yr_e;
H_list = p.Results.H_list;
V_list = p.Results.V_list;
num_vza = p.Results.num_vza;
check_task = p.Results.check_task;

%% Get the total task number, task input parameters and allocated the tasks to the jobs 
dir_l = '/shared/zhulab/Tian/Prod_09082023/ChangeMetricMap_l_20132023';
dir_l_rowchg = fullfile(dir_l, 'BlackMarble_BRDFcorrected_daily');
dir_l_chgmap = fullfile(dir_l, 'Result_ChangeMap');
dir_tiffref = '/shared/zhulab/Tian/VIIRS_NTL/Tiff_ref/';

% Get the tile list
tile_list = [];
for i_H = H_list
    for i_V = V_list(3:end-3)
        tile_name = sprintf('h%02dv%02d', i_H, i_V); 
        tile_list = [tile_list, tile_name];
    end
end

var_list = [...
    "GradualMag_net_row", ...
    "GradualMag_pos_row", ...
    "GradualMag_neg_row", ...
    ];
fname_list = [...
    "%s_GradualChangeMagnitude_net_%04d%03d_VZAint%d", ...
    "%s_GradualChangeMagnitude_pos_%04d%03d_VZAint%d", ...
    "%s_GradualChangeMagnitude_neg_%04d%03d_VZAint%d", ...
    ];
n_rst = 'Daily_GradualChangeMap_%s';

filename_task = 'tmp_task_extractchangemetricmap.mat';
path_task = fullfile(dir_l, filename_task);
if check_task == 0 && isfile(path_task)
    load(path_task, 'tasks');
else    
    % Set the input parameters for each of the task
    tasks = [];
    itask = 0;
    for i_tile = 1:length(tile_list)
        tile_name = tile_list(i_tile);
        dir_rowchg = fullfile(dir_l_rowchg, tile_name);
        dir_chgmap = fullfile(dir_l_chgmap, tile_name);

        for i_var = 1:length(var_list)
            fname = fname_list(i_var);
            var = var_list(i_var);
            for yr = yr_s:yr_e
                fname_save = sprintf(fname, tile_name, yr, day(datetime(yr, 12, 31), 'dayofyear'), num_vza);
                path_save = fullfile(dir_chgmap, sprintf(n_rst, ...
                    tile_name), sprintf('%s.tif', fname_save));
    
                % Skip if all files for the target tile has been processed
                if isfile(path_save) && replace_code == 0
                    continue
                end
    
                itask = itask+1;
                tasks(itask).tile = tile_name;
                tasks(itask).dir_chgmap = dir_chgmap;
                tasks(itask).dir_rowchg = dir_rowchg;
                tasks(itask).var = var;
                tasks(itask).fname = fname;
                tasks(itask).yr = yr;
            end
        end
    end
    save(path_task, 'tasks')
end

% The total task number
total_task = length(tasks);
rng(1);
task_ids = 1:total_task;

% Allocate the tasks to each jobs based on the total job number
ntasks_job = ones(njobs, 1);
num_add = mod(total_task, njobs);
ntasks_job = ntasks_job.*floor(total_task/njobs);
ntasks_job(1:num_add) = ntasks_job(1:num_add)+1;
task_s = sum(ntasks_job(1:job-1))+1;
task_e = sum(ntasks_job(1:job));
task_list = task_ids(task_s:task_e);

% The start and end saving line array irows for the map
for itask = task_list
    % Get the input parameters for the target task
    tile_name = tasks(itask).tile;
    dir_rowchg = tasks(itask).dir_rowchg;
    dir_chgmap = tasks(itask).dir_chgmap;
    var = tasks(itask).var;
    fname = tasks(itask).fname;
    yr = tasks(itask).yr;
    fprintf('Extracting for the tile %s, year %d\n', tile_name, yr)

    extract_changemetricmap_gradual(var, fname, yr_s, row1, row2, col1, col2, ...
        dir_rowchg, dir_chgmap, dir_tiffref, tile_name, yr, num_vza, replace_code);
end

end