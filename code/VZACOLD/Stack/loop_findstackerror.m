function loop_findstackerror(job, njobs, varargin)

%% get parameters from inputs
p = inputParser;

% default values.
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17);

% request user's input
parse(p, varargin{:});

H_list = p.Results.H_list;
V_list = p.Results.V_list;

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

% Read in the default tile list
% Request user's input
parse(p, varargin{:});
date_start = p.Results.date_start;
date_end = p.Results.date_end;
tile_list = p.Results.tile_list;
replace_code = p.Results.replace_code;

%% Get the total task number, task input parameters and allocated the tasks to the jobs 
% The directories of the source and desitination
dir_l_stack = '/shared/zhulab/Tian/Stacked_new';
% Get the start and end year/doy
yr_s = floor(date_start/1000);
doy_s = mod(date_start, 1000);
yr_e = floor(date_end/1000);
doy_e = mod(date_end, 1000);

% Set the input parameters for each of the task
tasks = [];
i_task = 0;
for i_tile = 1:length(tile_list)    
    % Get the check stack result
    tile_name = tile_list(i_tile);
    dir_stack = fullfile(dir_l_stack, tile_name);

    for yr = yr_s:yr_e
        % Get the check stack list directories
        if yr == yr_s
            doy_start = doy_s;
        else
            doy_start = 1;
        end
        
        if yr == yr_e
            doy_end = doy_e;
        elseif mod(yr, 4) == 0
            doy_end = 366;
        else
            doy_end = 365;
        end        
        filename_checkstack = sprintf('checkstack_%d%03d_%d%03d.mat', yr, doy_start, yr, doy_end);
        
        % Skip if existing
        if replace_code == 0 && isfile(fullfile(dir_stack, filename_checkstack))
            continue
        end
        
        i_task = i_task+1;
        tasks(i_task).tile = tile_name;
        tasks(i_task).yr = yr;   
        tasks(i_task).doy_start = doy_start;   
        tasks(i_task).doy_end = doy_end;   
    end
end

%% Allocate the tasks to each jobs based on the total job number
% Allocate the tasks randomly
total_task = length(tasks);% The total task number
rng(1);
tasks = tasks(randperm(total_task));

% Get the range of tasks for the job
ntasks_job = ones(njobs, 1);
num_add = mod(total_task, njobs);
ntasks_job = ntasks_job.*floor(total_task/njobs);
ntasks_job(1:num_add) = ntasks_job(1:num_add)+1;

% The start/end task ids for the target job
task_s = sum(ntasks_job(1:job-1))+1;
task_e = sum(ntasks_job(1:job));

%% Loop running all the tasks for the target job
for itask = task_e:-1:task_s %task_s:task_e
    % Get the input parameters for the target task
    tile_name = tasks(itask).tile;
    yr = tasks(itask).yr;
    doy_start = tasks(itask).doy_start;
    doy_end = tasks(itask).doy_end;
    
    % Check the stack list for the target tile & year
    find_stackerror_tile(dir_l_stack, tile_name, yr, doy_start, doy_end);
end
end
