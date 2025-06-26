function loop_DarkMask(job, njobs, varargin)

%% get parameters from inputs           
p = inputParser;
% default values.
addParameter(p, 'check_task', 0); 
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17);
addParameter(p, 'num_row', 2400);
addParameter(p, 'num_col', 2400);
addParameter(p, 'nlines_read', 60); 
addParameter(p, 'irows_s', 1);
addParameter(p, 'irows_e', 40);
addParameter(p, 'Dark_lmt', 1);
addParameter(p, 'replace_code', 0); 

% request user's input
parse(p, varargin{:});

check_task = p.Results.check_task;
H_list = p.Results.H_list;
V_list = p.Results.V_list;
nlines_read = p.Results.nlines_read;
irows_s = p.Results.irows_s;
irows_e = p.Results.irows_e;
num_row = p.Results.num_row;
num_col = p.Results.num_col;
Dark_lmt = p.Results.Dark_lmt;
replace_code = p.Results.replace_code;

filepath_m = mfilename('fullpath');
[cmd_l, ~] = fileparts(filepath_m);
addpath(cmd_l);

%% Get the total task number, task input parameters and allocated the tasks to the jobs
% Get the tile list
tile_list = [];
for i_H = H_list
    for i_V = V_list(3:end-3)
        tile_name = sprintf('h%02dv%02d', i_H, i_V); 
        tile_list = [tile_list, tile_name];
    end
end

dir_l = '/shared/zhulab/Tian/Prod_09082023/ChangeMetricMap_l_20132023/';
dir_l_data = fullfile(dir_l, '/Analysis/DarkPixel/');
dir_l_save = fullfile(dir_l_data, 'DarkMask');
if ~isfolder(dir_l_save)
    mkdir(dir_l_save)
end
dir_l_tif = '/shared/zhulab/Tian/VIIRS_NTL/Tiff_ref';
dir_l_task = '/home/til19015/GlobalNTLAnalyze/Analyze/DarkPixel';
if ~isfolder(dir_l_task)
    mkdir(dir_l_task)
end

filename_task = 'tmp_task_darkmask.mat';
dir_task = fullfile(dir_l_task, filename_task);
if check_task == 0 && isfile(dir_task)
    load(dir_task, 'tasks');
else
    % Set the input parameters for each of the task
    tasks = [];
    itask = 0;
    for i_tile = 1:length(tile_list)
        tile_name = tile_list(i_tile);
        fname_save = sprintf('DarkMask_%s.tif', tile_name);
        dir_data_row = fullfile(dir_l_data, tile_name);

        if replace_code == 0
            dir_save = fullfile(dir_l_save, fname_save);
            if isfile(dir_save)
                continue
            end
        end

        itask = itask+1;
        tasks(itask).tile = tile_name;
        tasks(itask).dir_data_row = dir_data_row;
    end

    save(fullfile(dir_l_task, filename_task), 'tasks');
end 
 
% The total incompleted task number
total_task = length(tasks);
rng(1);
task_ids = randperm(total_task);

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
    dir_data_row = tasks(itask).dir_data_row;

    DarkMask(nlines_read, irows_s, irows_e, dir_l_tif, num_row, num_col, ...
        tile_name, dir_l_save, dir_data_row, Dark_lmt);
end

end