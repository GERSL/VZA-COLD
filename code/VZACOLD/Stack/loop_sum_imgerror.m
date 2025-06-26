function loop_sum_imgerror(job, njobs, varargin)
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

% Get the directory
dir_stack = '/shared/zhulab/Tian/DataQuality';

tasks = [];
i_task = 0;
for tile_name = tile_list'
    filename_data = sprintf('DataQualityRec_%s.mat', tile_name);
    if ~isfile(fullfile(dir_stack, filename_data))
        continue
    end

    filename_save = sprintf('%s_sum.csv', tile_name);
    dir_save = fullfile(dir_stack, filename_save);
    if isfile(dir_save)
        % continue
    end

    i_task = i_task+1;
    tasks(i_task).tile_name = tile_name;
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
for itask = task_s:task_e
    % Get the input parameters for the target task
    tile_name = tasks(itask).tile_name;
    sum_imgerror(dir_stack, tile_name)
end
end

function sum_imgerror(dir_stack, tile_name)

filename_data = sprintf('DataQualityRec_%s.mat', tile_name);
load(fullfile(dir_stack, filename_data), 'T')

num_range = [T.clr_rangepix];
[~, index_num_range] = sort(num_range);

total_clr = [T.clear_pix]';
num_clr = [T.clr_rangepix_raw];
num_clr = reshape(num_clr, [10, length(T)]);
num_clr_4050 = (num_clr(5, :))';
num_clr_4060 = (num_clr(5, :)+num_clr(6, :))';
num_clr_3060 = (num_clr(4, :)+num_clr(5, :)+num_clr(6, :))';
[~, index_num_clr_4050] = sort(num_clr_4050);
[~, index_num_clr_4060] = sort(num_clr_4060);
[~, index_num_clr_3060] = sort(num_clr_3060);
prc_clr_4050 = num_clr_4050./total_clr;
prc_clr_4050(isnan(prc_clr_4050) | isinf(prc_clr_4050)) = 0;
[~, index_prc_clr_4050] = sort(prc_clr_4050);
prc_clr_4060 = num_clr_4060./total_clr;
prc_clr_4060(isnan(prc_clr_4060) | isinf(prc_clr_4060)) = 0;
[~, index_prc_clr_4060] = sort(prc_clr_4060);
prc_clr_3060 = num_clr_3060./total_clr;
prc_clr_3060(isnan(prc_clr_3060) | isinf(prc_clr_3060)) = 0;
[~, index_prc_clr_3060] = sort(prc_clr_3060);

total_all = [T.total_pix]';
num_all = [T.total_rangepix_raw];
num_all = reshape(num_all, [10, length(T)]);
num_all_4050 = (num_all(5, :))';
num_all_4060 = (num_all(5, :)+num_all(6, :))';
num_all_3060 = (num_all(4, :)+num_all(5, :)+num_all(6, :))';
[~, index_num_all_4050] = sort(num_all_4050);
[~, index_num_all_4060] = sort(num_all_4060);
[~, index_num_all_3060] = sort(num_all_3060);
prc_all_4050 = num_all_4050./total_all;
prc_all_4050(isnan(prc_all_4050) | isinf(prc_all_4050)) = 0;
[~, index_prc_all_4050] = sort(prc_all_4050);
prc_all_4060 = num_all_4060./total_all;
prc_all_4060(isnan(prc_all_4060) | isinf(prc_all_4060)) = 0;
[~, index_prc_all_4060] = sort(prc_all_4060);
prc_all_3060 = num_all_3060./total_all;
prc_all_3060(isnan(prc_all_3060) | isinf(prc_all_3060)) = 0;
[~, index_prc_all_3060] = sort(prc_all_3060);

for i = 1:length(T)
    if ~isfield(T, 'dip_clr') || isempty(T(i).dip_clr)
        T(i).dip_clr = [99, 99, 99, 99];
    end

    if ~isfield(T, 'dip_all') || isempty(T(i).dip_all)
        T(i).dip_all = [99, 99, 99, 99];
    end
end
dip_all = reshape([T.dip_all], [4, length(T)]);
p_value_all = dip_all(1, :);
[~, index_p_value_all] = sort(p_value_all);
dip_clr = reshape([T.dip_clr], [4, length(T)]);
p_value_clr = dip_clr(1, :);
[~, index_p_value_clr] = sort(p_value_clr);

save(fullfile(dir_stack, sprintf('summary_DataQuality_%s.mat', tile_name)), ...
    'T', 'index_num_range', ...
    'num_clr_4050', 'num_clr_4060', 'num_clr_3060', ...
    'index_num_clr_4050', 'index_num_clr_4060', 'index_num_clr_3060', ...
    'prc_clr_4050', 'prc_clr_4060', 'prc_clr_3060', ...
    'index_prc_clr_4050', 'index_prc_clr_4060', 'index_prc_clr_3060', ...
    'num_all_4050', 'num_all_4060', 'num_all_3060', ...
    'index_num_all_4050', 'index_num_all_4060', 'index_num_all_3060', ...
    'prc_all_4050', 'prc_all_4060', 'prc_all_3060', ...
    'index_prc_all_4050', 'index_prc_all_4060', 'index_prc_all_3060', ...
    'p_value_all', 'p_value_clr', 'index_p_value_all', 'index_p_value_clr')

Name = [T.name]';
ImgMean = [T.mean]';
AvgMean = mean(ImgMean, 'omitnan')*ones(length(ImgMean), 1);
StdMean = std(ImgMean, 'omitnan')*ones(length(ImgMean), 1);
DeltaMean = ImgMean-AvgMean;

ImgMedian = [T.median]';
AvgMedian = mean(ImgMedian, 'omitnan')*ones(length(ImgMedian), 1);
StdMedian = std(ImgMedian, 'omitnan')*ones(length(ImgMedian), 1);
DeltaMedian = ImgMedian-AvgMedian;

rec_tile = table(Name, ImgMean, AvgMean, StdMean, DeltaMean, ImgMedian, ...
    AvgMedian, StdMedian, DeltaMedian, total_clr, num_clr_3060, ...
    prc_clr_3060, total_all, num_all_3060, prc_all_3060);
writetable(rec_tile, fullfile(dir_stack, sprintf('%s_sum.csv', tile_name)))
end