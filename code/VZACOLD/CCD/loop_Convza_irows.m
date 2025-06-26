function loop_Convza_irows(job, njobs, varargin)

%% get parameters from inputs           
p = inputParser;
% default values.
addParameter(p, 'date_start', 2013001); 
addParameter(p, 'date_end', 2023365); 
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17);
addParameter(p, 'nlines_read', 60); 
addParameter(p, 'irows_s', 1);
addParameter(p, 'irows_e', 40);
addParameter(p, 'check_task', 0); 
addParameter(p, 'run_dir', 1); 

filepath_m = mfilename('fullpath');
[cmd_l, ~] = fileparts(filepath_m);
addpath(cmd_l);

% request user's input
parse(p, varargin{:});
date_start = p.Results.date_start;
date_end = p.Results.date_end;
H_list = p.Results.H_list;
V_list = p.Results.V_list;
nlines_read = p.Results.nlines_read;
irows_s = p.Results.irows_s;
irows_e = p.Results.irows_e;
check_task = p.Results.check_task;
run_dir = p.Results.run_dir;

% Get the tile list
tile_list = [];
for i_H = H_list
    for i_V = V_list(3:end-3)
        tile_name = sprintf('h%02dv%02d', i_H, i_V); 
        tile_list = [tile_list, tile_name];
    end
end

%% Get the total task number, task input parameters and allocated the tasks to the jobs
dir_l_new = '/shared/zhulab/Tian/Prod_09082023/COLDResult_l_20132023';
dir_l_task = '/home/til19015/GlobalNTLAnalyze/VZACOLD_master/Stable/Map';
dir_stack = '/shared/zhulab/Tian/Stacked_new/';
filename_task = 'tmp_task_cold.mat';
dir_task = fullfile(dir_l_task, filename_task);
if check_task == 0 && isfile(dir_task)
    load(dir_task, 'tasks');
else
    % Set the input parameters for each of the task
    tasks = [];
    itask = 0;
    for i_tile = 1:length(tile_list)
        tile_name = tile_list(i_tile);
        dir_BRDF_new = fullfile(dir_l_new, 'BlackMarble_BRDFcorrected', tile_name);
        n_rst = sprintf('TSFitMap_%s', tile_name);
        for i_irows = irows_s:irows_e
            % Skip if all lines have been run
            check_row = 0;
            for irow_ids = 1:nlines_read
                n_row = irow_ids+(i_irows-1)*nlines_read;

                filename_save = sprintf('%s_record_change%04d.mat', tile_name, n_row);
                filename_save_null = sprintf('%s_record_change%04d_null.mat', tile_name, n_row);
                if ~isfile(fullfile(dir_BRDF_new, n_rst, filename_save)) && ...
                        ~isfile(fullfile(dir_BRDF_new, n_rst, filename_save_null))
                    check_row = 1;
                    break
                end
            end
            if check_row == 0
                continue
            end

            n_row1 = nlines_read*(i_irows-1)+1;
            n_row2 = nlines_read*i_irows;
            folderpath_stack_rows = fullfile(dir_stack, tile_name, sprintf('%04d%04d', n_row1, n_row2));
            dir_datanum = fullfile(folderpath_stack_rows, 'DataNum');
            filename_obsnum = sprintf('DataNum.%04d%04d.*.txt', n_row1, n_row2);
            imf_datanum = dir(fullfile(dir_datanum, filename_obsnum));
            if isempty(imf_datanum)
                num_data = 0;
            else
                imf_datanum = {imf_datanum.name};
                imf_datanum = string(imf_datanum);
                num_data = split(imf_datanum, '.');
                num_data = num_data(:, :, 4);
                num_data = str2double(num_data);
                num_data = sum(num_data);
            end

            if num_data == 0
                continue
            end

            itask = itask + 1;
            tasks(itask).tile = tile_name;
            tasks(itask).n_irows = i_irows;
            tasks(itask).dir_BRDF_new = dir_BRDF_new;
            tasks(itask).num_data = num_data;
        end
    end
    save(fullfile(dir_l_task, filename_task), 'tasks');
end 
 
% The total incompleted task number
total_task = length(tasks);
% Sort the tasks by data number
[~, index_sort] = sort([tasks.num_data]);
tasks = tasks(index_sort);
% The start and end saving line array irows for the map
for itask = job:njobs:total_task
    % Get the input parameters for the target task
    tile_name = tasks(itask).tile;
    n_irows = tasks(itask).n_irows;
    dir_BRDF_new = tasks(itask).dir_BRDF_new;
    
    % Run COLD for the time-series in rows
    main_VZACOLD_map(dir_BRDF_new, tile_name, 0, n_irows, nlines_read, ...
        date_start, date_end, 'run_dir', run_dir)
end

end