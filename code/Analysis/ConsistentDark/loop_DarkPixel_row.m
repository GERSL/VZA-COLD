function loop_DarkPixel_row(job, njobs, varargin)

%% get parameters from inputs           
p = inputParser;
% default values.
addParameter(p, 'date_start', 2013001); 
addParameter(p, 'date_end', 2023365); 
addParameter(p, 'nlines_read', 60); 
addParameter(p, 'irows_s', 1);
addParameter(p, 'irows_e', 40);
addParameter(p, 'check_task', 0); 
addParameter(p, 'col1', 1); % Dimension of the map
addParameter(p, 'col2', 2400); 
addParameter(p, 'num_c', 4); 
addParameter(p, 'n_times', 12); % initial value for 16-day Landsat imagery
addParameter(p, 'conse', 14); 
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17);

% request user's input
parse(p, varargin{:});
date_start = p.Results.date_start;
date_end = p.Results.date_end;
nlines_read = p.Results.nlines_read;
irows_s = p.Results.irows_s;
irows_e = p.Results.irows_e;
check_task = p.Results.check_task;
col1 = p.Results.col1;
col2 = p.Results.col2;
num_c = p.Results.num_c;
n_times = p.Results.n_times;
conse = p.Results.conse;
H_list = p.Results.H_list;
V_list = p.Results.V_list;

filepath_m = mfilename('fullpath');
[cmd_l, ~] = fileparts(filepath_m);
addpath(cmd_l);

tile_list = [];
for i_H = H_list
    for i_V = V_list(3:end-3)
        tile_list = [tile_list; sprintf("h%02dv%02d", i_H, i_V)];
    end
end

%% Get the total task number, task input parameters and allocated the tasks to the jobs
dir_l_new = '/shared/zhulab/Tian/Prod_09082023/ChangeMetricMap_l_20132023/Analysis/DarkPixel';
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

        folderpath_stack = sprintf('/shared/zhulab/Tian/Stacked_new/%s/', tile_name);
        dir_DarkPixel = fullfile(dir_l_new, tile_name);
        if ~isfolder(dir_DarkPixel)
            mkdir(dir_DarkPixel)
        end

        for i_irows = irows_s:irows_e
            % Skip if all lines have been run
            check_row = 0;
            for irow_ids = 1:nlines_read
                n_row = irow_ids+(i_irows-1)*nlines_read;

                filename_save = sprintf('%s_DarkPixel_%04d.mat', tile_name, n_row);
                if ~isfile(fullfile(dir_DarkPixel, filename_save))
                    check_row = 1;
                    break
                end
            end
            if check_row == 0
                continue
            end

            n_row1 = nlines_read*(i_irows-1)+1;
            n_row2 = nlines_read*i_irows;
            folderpath_stack_rows = fullfile(folderpath_stack, sprintf('%04d%04d', n_row1, n_row2));
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
            tasks(itask).dir_DarkPixel = dir_DarkPixel;
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
    dir_DarkPixel = tasks(itask).dir_DarkPixel;
    
    % Run COLD for the time-series in rows
    DarkPixel_row(num_c, n_times, date_start, date_end, dir_DarkPixel, ...
        tile_name, n_irows, nlines_read, col1, col2, conse)
end

end