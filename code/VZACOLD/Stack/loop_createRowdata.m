function loop_createRowdata(job, njobs, tile_name, ...
    folderpath_NTLRaw, folderpath_NTLBRDF, folderpath_stack, varargin)
%% get parameters from inputs
p = inputParser;

% default values.
addParameter(p, 'date_start', 2013001);
addParameter(p, 'date_end', 2023365);
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17);
addParameter(p, 'row1', 1);% position limits
addParameter(p, 'row2', 2400);
addParameter(p, 'nlines_read', 60); % number of lines per stack files
addParameter(p, 'replace_code', 0); 

% request user's input
parse(p, varargin{:});

date_start = p.Results.date_start;
date_end = p.Results.date_end;
H_list = p.Results.H_list;
V_list = p.Results.V_list;
row1 = p.Results.row1;
row2 = p.Results.row2;
nlines_read = p.Results.nlines_read;
replace_code = p.Results.replace_code;

filepath_m = mfilename('fullpath');
[cmd_l, ~] = fileparts(filepath_m);
addpath(cmd_l);

%% Get the data directories
folderpath_NTLRaw = sprintf('/shared/cn449/VIIRS_NTL/BlackMarble_Raw/%s/', tile_name);
folderpath_NTLBRDF = sprintf('/shared/cn449/VIIRS_NTL/BlackMarble_BRDFcorrected/%s/', tile_name);
folderpath_stack = sprintf('/shared/zhulab/Tian/Stacked_new/%s/', tile_name);

%% Get the yeardoy list for the Black Marble data
% Dates of available raw Black Marble data
imf_raw = dir(fullfile(folderpath_NTLRaw, 'VNP46A1*.h5'));
yeardoy_raw = regexpi({imf_raw.name}, 'A20(\w*)', 'match');
yeardoy_raw = [yeardoy_raw{:}];
yeardoy_raw = vertcat(yeardoy_raw{:});
yeardoy_raw = string(yeardoy_raw(:, 2:8));

% Dates of available lunar-BRDF-corrected Black Marble data
imf_brdf = dir(fullfile(folderpath_NTLBRDF, 'VNP46A2*.h5'));
yeardoy_brdf = regexpi({imf_brdf.name}, 'A20(\w*)', 'match');
yeardoy_brdf = [yeardoy_brdf{:}];
yeardoy_brdf = vertcat(yeardoy_brdf{:});
yeardoy_brdf = string(yeardoy_brdf(:, 2:8));
% Only process the ones with corresponding raw data available
index_vld = ismember(yeardoy_brdf, yeardoy_raw) ...
    & double(yeardoy_brdf) >= date_start ...
    & double(yeardoy_brdf) <= date_end;
yeardoy_brdf = yeardoy_brdf(index_vld);
imf_brdf = imf_brdf(index_vld);

%% Get the yeardoy list for the stacking record
% Dates of the existing stack non-null and null data
folderpath_stack_rows = fullfile(folderpath_stack, sprintf('%04d%04d', ...
    row2-nlines_read+1, row2));
filename_stack = sprintf('map_tsdata_buffer_%04d%04d_*.mat', ...
    row2-nlines_read+1, row2);
imf_stack = dir(fullfile(folderpath_stack_rows, filename_stack));
yeardoy_stack = char({imf_stack.name}');
if ~isempty(yeardoy_stack)
    yeardoy_stack = string(yeardoy_stack(:, 28:34));
end

% Dates of the existing data number record
dir_datanum = fullfile(folderpath_stack_rows, 'DataNum');
if ~isfolder(dir_datanum)
    mkdir(dir_datanum)
end
filename_datanum = sprintf('DataNum.%04d%04d.*.txt', row2-nlines_read+1, row2);
imf_datanum = dir(fullfile(dir_datanum, filename_datanum));
yeardoy_datanum = char({imf_datanum.name}');
if ~isempty(yeardoy_datanum)
    yeardoy_datanum = string(yeardoy_datanum(:, 18:24));
end

%% Get all tasks
% If replacement is needed, check if the stacking & data number file exist
if replace_code == 0
    index_stack = ~ismember(yeardoy_brdf, yeardoy_stack);
    index_sum = ismember(yeardoy_brdf, yeardoy_stack) ...
        & ~ismember(yeardoy_brdf, yeardoy_datanum);
    
    yeardoy_brdf = [yeardoy_brdf(index_stack); yeardoy_brdf(index_sum)];
    imf_brdf = [imf_brdf(index_stack); imf_brdf(index_sum)];
    process_list = [repmat("stack", [sum(index_stack), 1]); repmat("summary", [sum(index_sum), 1])];
else 
    yeardoy_brdf = [yeardoy_brdf; yeardoy_brdf];
    imf_brdf = [imf_brdf; imf_brdf];
    process_list = [repmat("stack", [length(yeardoy_brdf), 1]); repmat("summary", [length(yeardoy_brdf), 1])];
end


% Assign the task list
i_task = 0;
tasks = [];
for i = 1:length(yeardoy_brdf)
    yeardoy = yeardoy_brdf(i);
    index_raw = yeardoy_raw == yeardoy;

    % Append the tasks
    i_task = i_task+1;
    tasks(i_task).yeardoy = yeardoy;
    tasks(i_task).filename_brdf = fullfile(folderpath_NTLBRDF, imf_brdf(i).name);
    tasks(i_task).filename_raw = fullfile(folderpath_NTLRaw, imf_raw(index_raw).name);
    tasks(i_task).process = process_list(i);
end

%% Allocate the tasks to each jobs based on the total job number
% Allocate the tasks randomly
total_task = length(tasks);% The total task number
if total_task > 0
    fprintf('%s: %d tasks left\n', tile_name, total_task)    
else
    fprintf('%s: done!!\n', tile_name)
end

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
for itask = task_s:task_e
    % Get the input parameters for the target task
    yeardoy = tasks(itask).yeardoy;
    filename_brdf = tasks(itask).filename_brdf;
    filename_raw = tasks(itask).filename_raw;
    process = tasks(itask).process;

    switch process
        case 'stack'
            % Read and stack the data into rows
            createRowdata(yeardoy, filename_brdf, filename_raw, ...
                folderpath_stack, nlines_read);

        case 'summary'
            % Summarize the valid observation number for each stacked rows
            summaryRowdata(yeardoy, folderpath_stack, nlines_read)
    end
end
end