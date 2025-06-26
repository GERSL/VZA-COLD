function loop_grid(job, njobs, varargin)
% Pixel area in km2 applied

%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17);
addParameter(p, 'year_start', 2014);
addParameter(p, 'year_end', 2022);
addParameter(p, 'grid_size', 120); %0.05/0.5-degree
addParameter(p, 'replace_code', 0); 
addParameter(p, 'check_task', 0); 
addParameter(p, 'n_vza', 4); 

% request user's input
parse(p, varargin{:});
H_list = p.Results.H_list;
V_list = p.Results.V_list;
year_start = p.Results.year_start;
year_end = p.Results.year_end;
grid_size = p.Results.grid_size;
replace_code = p.Results.replace_code;
check_task = p.Results.check_task;
n_vza = p.Results.n_vza;

tile_list = [];
for i_H = H_list
    for i_V = V_list
        tile_list = [tile_list; sprintf("h%02dv%02d", i_H, i_V)];
    end
end

dir_l = '/shared/zhulab/Tian/Prod_09082023/ChangeMetricMap_l_20132023';
dir_l_save = fullfile(dir_l, '/Analysis/Grid', sprintf('%dx%d', grid_size, grid_size));
dir_data = fullfile(dir_l, 'Result_ChangeMap_new');
dir_PixelArea = '/shared/zhulab/Tian/Analysis/PixelArea/'; % pixel area 
dir_wm = '/shared/zhulab/Tian/Analysis/BlackMarbleTiles/WaterMask/Geo/';
fname_wm = '%s_WaterMask_2014_geo.tif';

var_list = [...
    "ttlarea", ...
    "ttlarea_pos", ...
    "ttlarea_neg", ...
    "ttlrad", ...
    "ttlrad_pos", ...
    "ttlrad_neg", ...
    "ttlint_pos", ...
    "ttlint_neg", ...
    "posfreq", ...
    "negfreq", ...
    ];

fname_data_list = {...
    "%s_TotalChangeMagnitude_%d_VZAint%d.tif", ...
    "%s_TotalChangeMagnitude_pos_%d_VZAint%d.tif", ...
    "%s_TotalChangeMagnitude_neg_%d_VZAint%d.tif", ...
    "%s_TotalChangeMagnitude_%d_VZAint%d.tif", ...
    "%s_TotalChangeMagnitude_pos_%d_VZAint%d.tif", ...
    "%s_TotalChangeMagnitude_neg_%d_VZAint%d.tif", ...
    "%s_TotalChangeMagnitude_pos_%d_VZAint%d.tif", ...
    "%s_TotalChangeMagnitude_neg_%d_VZAint%d.tif", ...
    "%s_TotalChangeMagnitude_pos_%d_VZAint%d.tif", ...
    "%s_TotalChangeMagnitude_neg_%d_VZAint%d.tif", ...
    };

n_rst_list = {...
    "Annual_TotalChangeMap_%s", ...
    "Annual_TotalChangeMap_%s", ...
    "Annual_TotalChangeMap_%s", ...
    "Annual_TotalChangeMap_%s", ...
    "Annual_TotalChangeMap_%s", ...
    "Annual_TotalChangeMap_%s", ...
    "Annual_TotalChangeMap_%s", ...
    "Annual_TotalChangeMap_%s", ...
    "Annual_TotalChangeMap_%s", ...
    "Annual_TotalChangeMap_%s", ...
    };

%% Get the total task number, task input parameters and allocated the tasks to the jobs
filename_task = 'tmp_task_agggrid.mat';
path_task = fullfile(dir_l_save, filename_task);
if check_task == 0 && isfile(path_task)
    load(path_task, 'tasks');
else
    % Set the input parameters for each of the task
    tasks = [];
    itask = 0;
    for i_tile = 1:length(tile_list)
        tile_name = tile_list(i_tile, :);
        dir_save = fullfile(dir_l_save, tile_name);
        if ~isfolder(dir_save)
            mkdir(dir_save)
        end

        for i_var = 1:length(var_list)
            for year = year_start:year_end
                filename_save = sprintf('AggGrid_%s_%s_%d.tif', var_list(i_var), tile_name, year);
                path_save = fullfile(dir_save, filename_save);
                if isfile(path_save) && replace_code == 0
                    continue
                end

                itask = itask + 1;
                tasks(itask).tile_name = tile_name;
                tasks(itask).var = var_list(i_var); 
                tasks(itask).year = year;
                tasks(itask).n_rst = n_rst_list{i_var};
                tasks(itask).fname_data = fname_data_list{i_var};
                tasks(itask).path_save = path_save;
            end
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

% The start and end saving line array irows for the map
for itask = task_list
    % Get the input parameters for the target task
    tile_name = tasks(itask).tile_name; 
    var = tasks(itask).var;
    year = tasks(itask).year;
    n_rst = tasks(itask).n_rst;
    fname_data = tasks(itask).fname_data;
    path_save = tasks(itask).path_save;
    
    AggregateGrid(tile_name, var, year, dir_data, n_rst, fname_data, ...
        dir_PixelArea, dir_wm, fname_wm, path_save, grid_size, n_vza);
end
end