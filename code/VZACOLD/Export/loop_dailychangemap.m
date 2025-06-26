function loop_dailychangemap(job, njobs, varargin)

%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'replace_code', 1); 
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17);
addParameter(p, 'row1', 1);
addParameter(p, 'row2', 2400);
addParameter(p, 'col1', 1);
addParameter(p, 'col2', 2400);
addParameter(p, 'yr_s', 2013);
addParameter(p, 'yr_e', 2023);

% request user's input
parse(p, varargin{:});
replace_code = p.Results.replace_code;
H_list = p.Results.H_list;
V_list = p.Results.V_list;
row1 = p.Results.row1;
row2 = p.Results.row2;
col1 = p.Results.col1;
col2 = p.Results.col2;
yr_s = p.Results.yr_s;
yr_e = p.Results.yr_e;

filepath_m = mfilename('fullpath');
[cmd_l, ~] = fileparts(filepath_m);
addpath(cmd_l);

% Get the tile list
tile_list = [];
for i_H = H_list
    for i_V = V_list(3:end-3)
        tile_name = sprintf('h%02dv%02d', i_H, i_V); 
        tile_list = [tile_list, tile_name];
    end
end

%% Get the total task number, task input parameters and allocated the tasks to the jobs 
dir_l = '/shared/zhulab/Tian/Prod_09082023/';
dir_l = fullfile(dir_l, 'ChangeMap_l_MQbuf2_gap_75DQ');
dir_tiffref = '/shared/zhulab/Tian/VIIRS_NTL/Tiff_ref/';

% Set the input parameters for each of the task
tasks = [];
itask = 0;
for i_tile = 1:length(tile_list)   
    tile_name = tile_list(i_tile);
    dir_rowchg = fullfile(dir_l, 'BlackMarble_BRDFcorrected', tile_name);
    dir_changemap = fullfile(dir_l, 'Result_ChangeMap', tile_name);
    n_rst_accmap = sprintf('AccumulateChangeMap_%s', tile_name);
    filename_save_color = sprintf('%s_AccuChangeMap.clr', tile_name, yr_e);      
    if replace_code == 0
        dir_save = fullfile(dir_changemap, n_rst_accmap, filename_save_color);
        if isfile(dir_save)
            continue
        end
    end

    itask = itask+1;
    tasks(itask).tile = tile_name;
    tasks(itask).dir_changemap = dir_changemap;
    tasks(itask).dir_rowchg = dir_rowchg;
end

% The total task number
total_task = length(tasks);
rng(1);
% task_ids = randperm(total_task);
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
    dir_changemap = tasks(itask).dir_changemap;

    dailychangemap(yr_s, yr_e, row1, row2, col1, col2, dir_rowchg, ...
        dir_changemap, dir_tiffref, tile_name);
end

end