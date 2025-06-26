function summaryRowdata(yeardoy, folderpath_stack, nlines_read, varargin)
% This script aims to summary the stacked data's observaton number for each
% row
% created on May, 02, 2024

%% get parameters from inputs
p = inputParser;

% default values.
addParameter(p, 'row1', 1);% position limits
addParameter(p, 'row2', 2400);

% request user's input
parse(p, varargin{:});
row1 = p.Results.row1;
row2 = p.Results.row2;

%% Add seaching paths
path_mfile = fileparts(mfilename('fullpath'));
addpath(fileparts(path_mfile)); % add parent folder, including sets.m

%% Loop stack for the row array
% Get the row array
irows = row1+nlines_read.*(0:ceil((row2-row1)/nlines_read)-1);
for i_rows = 1:length(irows)
    % Read the stack data
    n_row1 = irows(i_rows);
    n_row2 = n_row1+nlines_read-1;
    folderpath_stack_rows = fullfile(folderpath_stack, sprintf('%04d%04d', n_row1, n_row2));
    filename_stack = sprintf('map_tsdata_buffer_%04d%04d_%s.mat', ...
        n_row1, n_row2, yeardoy);
    filename_stack_null = sprintf('map_tsdata_buffer_%04d%04d_%s_null.mat', ...
        n_row1, n_row2, yeardoy);

    if isfile(fullfile(folderpath_stack_rows, filename_stack))
        % Summarize the data number file
        load(fullfile(folderpath_stack_rows, filename_stack), 'line_tbrdf');
        num_data = line_tbrdf < 65535;
        num_data = sum(num_data(:));
    elseif isfile(fullfile(folderpath_stack_rows, filename_stack_null))
        num_data = 0;
    else
        continue
    end

    % Save the available data number information
    dir_datanum = fullfile(folderpath_stack_rows, 'DataNum');
    if ~isfolder(dir_datanum)
        mkdir(dir_datanum)
    end
    filename_obsnum = sprintf('DataNum.%04d%04d.%s.%09d.txt', n_row1, n_row2, yeardoy, num_data);
    fid = fopen(fullfile(dir_datanum, filename_obsnum), 'w');
    fclose(fid);
end
clear line_tbrdf line_tvza line_tclr line_tsn DNBblk_BRDF ...
    SFblk MQblk QFDNBblk CMblk VZAblk VZA_grp_blk CM_vld CM_nght CM_cld_c...
    CM_cld_p CM_crs CM_snw CM_clr Phy_blk Snw_blk CM_cld_blk Clr_blk
end

