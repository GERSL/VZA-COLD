function loop_mosaicchangemetricmap(job, njobs, varargin)

%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17);
addParameter(p, 'num_row', 2400);
addParameter(p, 'num_col', 2400);
addParameter(p, 'yr_s', 2014);
addParameter(p, 'yr_e', 2022);
addParameter(p, 'replace_code', 0); 
addParameter(p, 'check_task', 0); 
addParameter(p, 'n_VZA', 4); 

% request user's input
parse(p, varargin{:});
H_list = p.Results.H_list;
V_list = p.Results.V_list;
num_row = p.Results.num_row;
num_col = p.Results.num_col;
yr_s = p.Results.yr_s;
yr_e = p.Results.yr_e;
replace_code = p.Results.replace_code;
check_task = p.Results.check_task;
n_VZA = p.Results.n_VZA;

sdate = datenum(yr_s, 1, 1):datenum(yr_e, 12, 31);
dir_l = '/shared/zhulab/Tian/Prod_09082023/';
dir_l = fullfile(dir_l, 'ChangeMetricMap_l_20132023', 'Result_ChangeMap_new');
dir_l_save = fullfile(dir_l, 'Mosaic_DailyGradual');
if ~isfolder(dir_l_save)
    mkdir(dir_l_save)
end

% Get the map list
n_rst_list = [...
    "Daily_GradualChangeMap_%s", ... % Daily gradual change map
    "Daily_GradualChangeMap_%s", ...
    "Daily_GradualChangeMap_%s", ...
    ];

fname_list = [...
    "_GradualChangeMagnitude_net_%04d%03d_VZAint%d.tif", ... 
    "_GradualChangeMagnitude_pos_%04d%03d_VZAint%d.tif", ... 
    "_GradualChangeMagnitude_neg_%04d%03d_VZAint%d.tif", ... 
    ];

%% Get the total task number, task input parameters and allocated the tasks to the jobs
filename_task = 'tmp_task_mosaicchangemetricmap.mat';
path_task = fullfile(dir_l_save, filename_task);
if check_task == 0 && isfile(path_task)
    load(path_task, 'tasks');
else
    % Initialize the task
    tasks = [];
    itask = 0;
    for i_map = length(n_rst_list)
        n_rst = n_rst_list(i_map);

        for i_date = sdate
            yr_date = year(datetime(i_date, 'ConvertFrom', 'datenum'));
            doy_date = day(datetime(i_date, 'ConvertFrom', 'datenum'), 'dayofyear');
            
            fname = sprintf(fname_list(i_map), yr_date, doy_date, n_VZA);
            fname_save = sprintf('%s%s', 'Global', fname);
            if isfile(fullfile(dir_l_save, fname_save)) && replace_code == 0
                continue
            end

            itask = itask+1;
            tasks(itask).n_rst = n_rst;
            tasks(itask).fname = fname;
        end
    end
    save(path_task, 'tasks')
end
% The total task number
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

% Get the georeference information for the global map
R = generate_globaltif_Ref('save_code', 0);
% The start and end saving line array irows for the map
for itask = task_list   
    % Get the input parameters for the target task
    n_rst = tasks(itask).n_rst;
    fname = tasks(itask).fname;
    fprintf('Mosaicing the global daily %s\n', fname) 

    Mosaic_changemetricmap(H_list, V_list, num_row, num_col, dir_l, ...
        dir_l_save, n_rst, fname, R);
end

end