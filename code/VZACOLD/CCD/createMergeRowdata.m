function [sdate, line_tbrdf_mrg, line_tvzagrp_mrg, DeltaMedian, StdMedian, DeltaMean, StdMean] = ...
    createMergeRowdata(tile_name, n_irows, sdate, varargin)

%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'nlines_read', 60); % default 10 lines
addParameter(p, 'T_cg_error', norminv(1-(1-0.75)/2)); % 

% request user's input
parse(p, varargin{:});
nlines_read = p.Results.nlines_read;
T_cg_error = p.Results.T_cg_error;

% Set the upper directories for torage
dir_l_stack = '/shared/zhulab/Tian/Stacked_new/';
dir_stack = fullfile(dir_l_stack, tile_name);
dir_qa = '/shared/zhulab/Tian/DataQuality';

% Get the default position limits
[row1, row2, col1, col2] = deal(1, 2400, 1, 2400);
ncols = col2-col1+1;
% The start and end saving line array irows for the map
irows = row1+nlines_read.*(0:ceil((row2-row1)/nlines_read)-1);
% For the target row array
n_row1 = irows(n_irows);
n_row2 = n_row1+nlines_read-1;
n_rst_stack = fullfile(dir_stack, sprintf('%04d%04d', n_row1, n_row2));

% Get the data quality issue list
filename_qa = sprintf('%s_sum.csv', tile_name);
rec_tile = readtable(fullfile(dir_qa, filename_qa));
StdMedian = rec_tile.StdMedian(1);
DeltaMedian = [rec_tile.DeltaMedian];
StdMean = rec_tile.StdMean(1);
DeltaMean = [rec_tile.DeltaMean];

yeardoy_data = [rec_tile.Name];
yeardoy_data = split(yeardoy_data, '.');
yeardoy_data = yeardoy_data(:, 2);
[yeardoy_data, i_uni, ~] = unique(yeardoy_data, 'last');
DeltaMedian = DeltaMedian(i_uni);
DeltaMean = DeltaMean(i_uni);
index_error = DeltaMean >= T_cg_error*StdMean;

yeardoy_data = char(yeardoy_data);
yeardoy_data = yeardoy_data(:, 2:8);
yeardoy_error = yeardoy_data(index_error, :);

% Get the actual available good quality data date list
name_date = datetime(sdate, 'ConvertFrom', 'datenum', 'Format', 'yyyyDDD');
name_date = string(name_date);
num_t = length(sdate);
for n_t = 1:num_t
    filename_row = sprintf('map_tsdata_buffer_%04d%04d_%s.mat', ...
        n_row1, n_row2, name_date(n_t));
    
    % Skip if no stack data
    if ~isfile(fullfile(n_rst_stack, filename_row))
        sdate(n_t) = nan;
    end
    
    % Skip if error data
    if ismember(name_date(n_t), yeardoy_error)
        sdate(n_t) = nan;
    end
end
name_date(isnan(sdate)) = [];
sdate(isnan(sdate)) = [];
DeltaMedian = DeltaMedian(ismember(string(yeardoy_data), name_date));
DeltaMean = DeltaMean(ismember(string(yeardoy_data), name_date));

% Preparing the default values and casts for the stack variables
num_t = length(sdate);
line_tbrdf_mrg = zeros(nlines_read, ncols, num_t, 'uint16'); % BRDF-corrected DNB
line_tvzagrp_mrg = zeros(nlines_read, ncols, num_t, 'int16'); % VZAs
% For each of the date
for n_t = 1:num_t
    filename_row = sprintf('map_tsdata_buffer_%04d%04d_%s.mat', ...
        n_row1, n_row2, name_date(n_t));
   
    % Check for stack data verions
    variableInfo = who('-file', fullfile(n_rst_stack, filename_row));
    var_null = {'nulldata'};
    % Skip if null data
    if isequal(variableInfo, var_null) 
        continue
    % Load the stack data: summarized stack data
    elseif ismember('line_tvzagrp', variableInfo) 
        load(fullfile(n_rst_stack, filename_row), 'line_tbrdf', 'line_tvzagrp')

    end
    line_tbrdf_mrg(:, :, n_t) = line_tbrdf;
    line_tvzagrp_mrg(:, :, n_t) = line_tvzagrp;
    
    clear line_tbrdf line_tvza line_tclr line_tvzagrp
end
end