function loop_mosaicagggrid(job, njobs, varargin)

%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17);
addParameter(p, 'num_row', 2400);
addParameter(p, 'num_col', 2400);
addParameter(p, 'year_list', 2014:2022);
addParameter(p, 'grid_size', 120); %0.5/0.05-degree
addParameter(p, 'replace_code', 0); 
addParameter(p, 'check_task', 0); 

% request user's input
parse(p, varargin{:});
H_list = p.Results.H_list;
V_list = p.Results.V_list;
num_row = p.Results.num_row;
num_col = p.Results.num_col;
year_list = p.Results.year_list;
grid_size = p.Results.grid_size;
replace_code = p.Results.replace_code;
check_task = p.Results.check_task;

tile_list = [];
for i_H = H_list
    for i_V = V_list
        tile_list = [tile_list; sprintf("h%02dv%02d", i_H, i_V)];
    end
end

% Get the directories
dir_l = fullfile('/shared/zhulab/Tian/Prod_09082023/ChangeMetricMap_l_20132023/Analysis/Grid_AreaAdj', sprintf('%dx%d', grid_size, grid_size));
dir_l_save = fullfile(dir_l, 'Mosaic');
if ~isfolder(dir_l_save)
    mkdir(dir_l_save)
end

var_list = [...
    "ttlarea_pos", ...
    "ttlarea_neg", ...
    "ttlrad_pos", ...
    "ttlrad_neg", ...
    "ttlint_pos", ...
    "ttlint_neg", ...
    ];

%% Get the total task number, task input parameters and allocated the tasks to the jobs
filename_task = 'tmp_task_mosaicagggrid.mat';
path_task = fullfile(dir_l, filename_task);
if check_task == 0 && isfile(path_task)
    load(path_task, 'tasks');
else
    % Set the input parameters for each of the task
    tasks = [];
    itask = 0;
    for year = year_list
        for var = var_list
            filename_save = sprintf('GlobalAggGrid_%s_%d.tif', var, year);
            path_save = fullfile(dir_l_save, filename_save);
            if isfile(path_save) && replace_code == 0
                continue
            end

            itask = itask + 1;
            tasks(itask).var = var;
            tasks(itask).year = year;
            tasks(itask).path_save = path_save;
        end
    end
    save(path_task, 'tasks')
end
% The total incompleted task number
total_task = length(tasks);
task_ids = 1:total_task;

% Allocate the tasks to each jobs based on the total job number
ntasks_job = ones(njobs, 1);
num_add = mod(total_task, njobs);
ntasks_job = ntasks_job.*floor(total_task/njobs);
ntasks_job(1:num_add) = ntasks_job(1:num_add)+1;
task_s = sum(ntasks_job(1:job-1))+1;
task_e = sum(ntasks_job(1:job));
task_list = task_ids(task_s:task_e);

% Get the georeference information for the global map
R = generate_globaltif_Ref('save_code', 0);
num_row_agggrid = num_row/grid_size;
num_col_agggrid = num_col/grid_size;
R.RasterSize = [num_row_agggrid*length(V_list), num_col_agggrid*length(H_list)];
% The start and end saving line array irows for the map
for itask = task_list
    % Get the input parameters for the target task
    var = tasks(itask).var;
    year = tasks(itask).year;
    path_save = tasks(itask).path_save;
    
    Mosaic_agggridmap(H_list, V_list, num_row_agggrid, num_col_agggrid, ...
        dir_l, var, year, path_save, R)
end
end