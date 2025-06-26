function loop_imshow_posneg_heatmap(job, njobs, varargin)

%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17);
addParameter(p, 'year_list', 2014:2022);
addParameter(p, 'num_row', 2400);
addParameter(p, 'num_col', 2400);
addParameter(p, 'grid_size', 120); %0.05/0.5-degree
addParameter(p, 'replace_code', 1); 
addParameter(p, 'check_task', 0);

% request user's input
parse(p, varargin{:});
H_list = p.Results.H_list;
V_list = p.Results.V_list;
year_list = p.Results.year_list;
num_row = p.Results.num_row;
num_col = p.Results.num_col;
grid_size = p.Results.grid_size;
replace_code = p.Results.replace_code;
check_task = p.Results.check_task;

dir_l = '/shared/zhulab/Tian/Prod_09082023/ChangeMetricMap_l_20132023/Analysis/Grid_AreaAdj';
dir_l_data = fullfile(dir_l, sprintf('%dx%d', grid_size, grid_size), 'Mosaic');
dir_l_save = fullfile(dir_l_data, 'Imshow');
if ~isfolder(dir_l_save)
    mkdir(dir_l_save)
end
var_list = [...
    "ttlarea", ...
    "ttlrad", ...
    "ttlint", ...
    ];

%% Get the total task number, task input parameters and allocated the tasks to the jobs
filename_task = 'tmp_task_imshowheatmap.mat';
path_task = fullfile(dir_l, filename_task);
if check_task == 0 && isfile(path_task)
    load(path_task, 'tasks');
else
    % Set the input parameters for each of the task
    tasks = [];
    itask = 0;
    for i_var = 1:length(var_list)
        var = var_list(i_var);

        year = year_list(1)*10000+year_list(end);
        filename_savetif = sprintf('HeatMap_PosNeg%s_%d.tif', var, year);
        path_savetif = fullfile(dir_l_save, filename_savetif);
        if ~isfile(path_savetif) || replace_code == 1
            itask = itask + 1;
            tasks(itask).var = var;
            tasks(itask).year = year;
            tasks(itask).path_savetif = path_savetif;
        end
    end

    save(path_task, 'tasks')
end


% The total incompleted task number
total_task = length(tasks);
task_ids = 1:total_task;

% Get the georeference
% [~, Img_R] = readgeoraster(path_tiffref);
Img_R = generate_globaltif_Ref('save_code', 0);
num_row_agggrid = num_row/grid_size;
num_col_agggrid = num_col/grid_size;
Img_R.RasterSize = [num_row_agggrid*length(V_list), num_col_agggrid*length(H_list)];

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
    var = tasks(itask).var;
    year = tasks(itask).year;
    path_savetif = tasks(itask).path_savetif;

    imshow_posneg_heatmap(var, dir_l_data, year, path_savetif, Img_R)
end
end