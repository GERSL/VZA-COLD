function tasks_tile = find_stackerror_tile(dir_l_stack, tile_name, year, doy_start, doy_end, varargin)

%% Get input and defalut values for the directories of the source and destination. 
p = inputParser;

% Default inputs.
addParameter(p, 'nlines_read', 60); 
addParameter(p, 'row1', 1);% position limits
addParameter(p, 'row2', 2400);

% Read in the default tile list
% Request user's input
parse(p, varargin{:});
nlines_read = p.Results.nlines_read;
row1 = p.Results.row1;
row2 = p.Results.row2;

%% Get the default directories
cmd_l = '/home/til19015/GlobalNTLAnalyze/StackData';
cd(cmd_l)
addpath(cmd_l)

% The top directories of the source and desitination
dir_l_data = '/shared/cn449/VIIRS_NTL/';

% The storage directory for the download list
dir_data_raw = fullfile(dir_l_data, 'BlackMarble_Raw', 'Download_list', 'Checked');
dir_data_brdf = fullfile(dir_l_data, 'BlackMarble_BRDFcorrected', 'Download_list', 'Checked');
dir_stack = fullfile(dir_l_stack, tile_name);

% Create the stack folder for tcheck_taskshe records if not existing
if ~isfolder(dir_stack)
    mkdir(dir_stack)
end

%% Get the total task number, task input parameters and allocated the tasks to the jobs 
% Set the input parameters for each of the task
tasks = [];
i_task = 0;
for doy = doy_start:doy_end
    i_task = i_task+1;
    tasks(i_task).year = year;
    tasks(i_task).doy = doy;
end

%% Check the stack data for all the tasks
% Initialized the check error number
tasks_tile = [];
i_tasktile = 0;
% Loop running all the tasks
for i_task = 1:length(tasks)
    % Get the input parameters for the target task
    year = tasks(i_task).year;
    doy = tasks(i_task).doy;
    
    % Check for the stack data of the target date
    error_stack = find_stackerror(tile_name, year, doy, dir_data_raw, ...
        dir_data_brdf, dir_stack, row1, row2, nlines_read);
    
    % If no valid data, skip
    if isempty(error_stack)
        continue
    end
    
    i_tasktile = i_tasktile+1;
    tasks_tile(i_tasktile).tile = tile_name;
    tasks_tile(i_tasktile).t = sprintf('%d%03d', year, doy);
    tasks_tile(i_tasktile).checkstack = error_stack;
end

%% Save the check stack result 
filename_checkstack = sprintf('checkstack_%d%03d_%d%03d', year, doy_start, year, doy_end);
save(fullfile(dir_stack, filename_checkstack), 'tasks_tile')

if isempty(tasks_tile)
    fprintf('All tasks fully stacked for year %d:D\n', year)
else
    check_tasks = sum([tasks_tile.checkstack] ~= 0);
    if check_tasks == 0
        fprintf('All tasks fully stacked for year %d:D\n', year)
    else
        fprintf('%d tasks to be stack...\n', check_tasks)
    end
end
end
