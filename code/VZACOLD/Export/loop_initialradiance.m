function loop_initialradiance(job, njobs, varargin)
% Updated on 10/10/2024:
% Extract the overal NTL radiance at the beginning of the study period.

%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'date_start', 2014001); % start year of the obs. time-series
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17);
addParameter(p, 'row1', 1);
addParameter(p, 'row2', 2400);
addParameter(p, 'col1', 1);
addParameter(p, 'col2', 2400)
addParameter(p, 'check_task', 0);
addParameter(p, 'random_code', 0);

% request user's input
parse(p, varargin{:});
date_start = p.Results.date_start;
H_list = p.Results.H_list;
V_list = p.Results.V_list;
row1 = p.Results.row1;
row2 = p.Results.row2;
col1 = p.Results.col1;
col2 = p.Results.col2;
check_task = p.Results.check_task;
random_code = p.Results.random_code;

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
dir_l = '/shared/zhulab/Tian/Prod_09082023/ChangeMetricMap_l_20132023/BlackMarble_BRDFcorrected_new';
dir_ref = '/shared/zhulab/Tian/Analysis/BlackMarbleTiles/Tiff_ref/';
filename_task = 'tmp_task_changemetric.mat';
dir_task = fullfile(dir_l, filename_task);
if check_task == 0 && isfile(dir_task)
    load(dir_task, 'tasks');
else
    % Set the input parameters for each of the task
    tasks = [];
    itask = 0;
    for i_tile = 1:length(tile_list)
        tile_name = tile_list(i_tile);
        n_rst_save = sprintf('InitialRadianceMap_%s', tile_name);
        fname_save = sprintf('%s_initialRadiance.tif', tile_name);
        dir_BRDF = fullfile(dir_l, tile_name);
        if ~isfolder(fullfile(dir_BRDF, n_rst_save))
            mkdir(fullfile(dir_BRDF, n_rst_save));
        end

        if isfile(fullfile(dir_BRDF, n_rst_save, fname_save))
            continue
        end

        itask = itask+1;
        tasks(itask).tile = tile_name;
        tasks(itask).n_rst_save = n_rst_save;
        tasks(itask).n_rst_data = sprintf('ChangeMetricRec_%s', tile_name);
        tasks(itask).dir_BRDF = dir_BRDF;
    end

    save(fullfile(dir_l, filename_task), 'tasks');
end

% The total task number
total_task = length(tasks);
incmplt_tasks = 1:total_task;
if random_code == 1
    rng(1);
    task_ids = randperm(total_task);
    task_ids = incmplt_tasks(task_ids);
else
    task_ids = 1:total_task;
end

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
    n_rst_save = tasks(itask).n_rst_save;
    n_rst_data = tasks(itask).n_rst_data;
    dir_BRDF = tasks(itask).dir_BRDF; 
    
    % Run COLD for the time-series in rows
    initialradiance(dir_BRDF, dir_ref, tile_name, ...
        n_rst_data, n_rst_save, date_start, row1, row2, col1, col2)
end

end