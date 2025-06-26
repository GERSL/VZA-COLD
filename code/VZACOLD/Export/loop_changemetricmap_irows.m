function loop_changemetricmap_irows(job, njobs, varargin)

%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'date_start', 2013001); % start year of the obs. time-series
addParameter(p, 'date_end', 2023365); 
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17);
addParameter(p, 'irows_s', 1);
addParameter(p, 'irows_e', 40);
addParameter(p, 'nlines_read', 60); 
addParameter(p, 'col1', 1);
addParameter(p, 'col2', 2400)
addParameter(p, 'sgn', 0.05)
addParameter(p, 'rerun', 0);
addParameter(p, 'check_task', 0);
addParameter(p, 'random_code', 0); 
addParameter(p, 'T_cg', 0.75); % change threshold

% request user's input
parse(p, varargin{:});
date_start = p.Results.date_start;
date_end = p.Results.date_end;
H_list = p.Results.H_list;
V_list = p.Results.V_list;
irows_s = p.Results.irows_s;
irows_e = p.Results.irows_e;
col1 = p.Results.col1;
col2 = p.Results.col2;
nlines_read = p.Results.nlines_read;
sgn = p.Results.sgn;
rerun = p.Results.rerun;
check_task = p.Results.check_task;
random_code = p.Results.random_code;
T_cg = p.Results.T_cg;

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
% adjust threshold based on normal distribution
T_cg = norminv(1-(1-T_cg)/2);
dir_l = '/shared/zhulab/Tian/Prod_09082023/';
dir_chgmap = fullfile(dir_l, 'ChangeMetricMap_l_20132023');
dir_chgrec = fullfile(dir_l, 'COLDResult_l_20132023');
if ~isfolder(dir_chgmap)
    mkdir(dir_chgmap)           
end
n_rst_cold = 'BlackMarble_BRDFcorrected';
n_rst_save = 'BlackMarble_BRDFcorrected_new_daily';
filename_task = 'tmp_task_changemetric.mat';
dir_task = fullfile(dir_chgmap, filename_task);

if check_task == 0 && isfile(dir_task)
    load(dir_task, 'tasks');
else
    % Set the input parameters for each of the task
    tasks = [];
    itask = 0;
    for i_tile = 1:length(tile_list)
        tile_name = tile_list(i_tile);

        dir_scn_BRDF = fullfile(n_rst_save, tile_name);
        dir_BRDF_new = fullfile(dir_chgmap, dir_scn_BRDF);
        n_rst_rec = sprintf('ChangeMetricRec_%s', tile_name);
        n_rst_map = sprintf('ChangeMetricMap_%s', tile_name);

        for i_irows = irows_s:irows_e 
            check_row = 0;
            for irow_ids = 1:nlines_read
                n_row = irow_ids+(i_irows-1)*nlines_read;

                filename_save = sprintf('%s_changemetric%04d.mat', tile_name, n_row);
                filename_save_null = sprintf('%s_changemetric%04d_null.mat', tile_name, n_row);
                if ~isfile(fullfile(dir_BRDF_new, n_rst_map, filename_save)) && ...
                        ~isfile(fullfile(dir_BRDF_new, n_rst_map, filename_save_null))
                    check_row = 1;
                    break
                end

                if ~isfile(fullfile(dir_BRDF_new, n_rst_rec, filename_save)) && ...
                        ~isfile(fullfile(dir_BRDF_new, n_rst_rec, filename_save_null))
                    check_row = 1;
                    break
                end
            end
            % Skip if all lines have been run
            if check_row == 0 && rerun == 0
                continue
            end

            itask = itask+1;
            tasks(itask).tile = tile_name;
            tasks(itask).n_rst_rec = n_rst_rec;
            tasks(itask).n_rst_map = n_rst_map;
            tasks(itask).n_irows = i_irows;
        end
    end
    
    save(fullfile(dir_chgmap, filename_task), 'tasks');
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
    n_rst_rec = tasks(itask).n_rst_rec;
    n_rst_map = tasks(itask).n_rst_map;
    n_irows = tasks(itask).n_irows;
    
    % Run COLD for the time-series in rows
    changemetricmap_Row(dir_chgmap, dir_chgrec, n_rst_cold, n_rst_save, ...
        tile_name, n_rst_rec, n_rst_map, date_start, date_end, n_irows, ...
        nlines_read, col1, col2, sgn, T_cg)
end
end