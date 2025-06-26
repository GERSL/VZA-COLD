function createRowdata(yeardoy, filename_brdf, filename_raw, ...
    folderpath_stack, nlines_read, varargin)
% This script aims to stack the Black Marble NTL data for the target tile &
% date
% created on Aug, 22, 2023

%% get parameters from inputs
p = inputParser;

% default values.
addParameter(p, 'row1', 1);% position limits
addParameter(p, 'row2', 2400);
addParameter(p, 'col1', 1);
addParameter(p, 'col2', 2400);
addParameter(p, 'win_buf', 2); % dialog size for the cloud buffer window
addParameter(p, 'vza_lmt', [0, 20; 20, 40; 40, 60]);

% request user's input
parse(p, varargin{:});
row1 = p.Results.row1;
row2 = p.Results.row2;
col1 = p.Results.col1;
col2 = p.Results.col2;
win_buf = p.Results.win_buf;
vza_lmt = p.Results.vza_lmt;

%% Add seaching paths
path_mfile = fileparts(mfilename('fullpath'));
addpath(fileparts(path_mfile)); % add parent folder, including sets.m

%% Get the default variables for the .h5 data
% group directory for the hdf5 file datasets
fld = '/HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/';
nlines = row2-row1+1;
% The HDF5 library uses C-style ordering for multidimensional
% arrays, while MATLAB uses FORTRAN-style ordering
% dimention of extracted data
dims = [nlines, col2-col1+1];
% starting pixel location
startpxl = [row1-1, col1-1];

%% Read in VNP46A2 observations
file_id = H5F.open(fullfile(filename_brdf), 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

% BRDF Corrected DNB Radiance
data_name = [fld, 'DNB_BRDF-Corrected_NTL'];
DNBblk_BRDF = read_hdf5(file_id, data_name, dims, startpxl);
% Snow Flag
data_name = [fld, 'Snow_Flag'];
SFblk = read_hdf5(file_id, data_name, dims, startpxl);
% Mandatory quality Flag
data_name = [fld, 'Mandatory_Quality_Flag'];
MQblk = read_hdf5(file_id, data_name, dims, startpxl); % 0=permanent NTL; 1=ephermeral NTL; 2=poor quality retrieval

H5F.close(file_id)

%% Read in VNP46A1 observations
file_id = H5F.open(filename_raw, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

% DNB quality flag
data_name = [fld, 'QF_DNB'];
QFDNBblk = read_hdf5(file_id, data_name, dims, startpxl); % 1=subsitute calculation; >0=poor quality
% Cloud Mask Status
data_name = [fld, 'QF_Cloud_Mask'];
CMblk = read_hdf5(file_id, data_name, dims, startpxl);
% Sensor Zenith
data_name = [fld, 'Sensor_Zenith'];
VZAblk = read_hdf5(file_id, data_name, dims, startpxl);

H5F.close(file_id)

%% Extract the useful records
% Initialize the records
CM_nght = zeros(size(CMblk), 'logical');
CM_cld_c = zeros(size(CMblk), 'logical');
CM_cld_p = zeros(size(CMblk), 'logical');
CM_crs = zeros(size(CMblk), 'logical');
CM_snw = zeros(size(CMblk), 'logical');

% Get the valid VCM records
CM_vld = CMblk < 65535;
% Day/Night (0=Night, 1=Day) - Find Night time observation
CM_nght(CM_vld & bitget(CMblk, 1, 'uint16') == 0) = 1;
% Cloud Detection Results & Confidence Indicator (00=Condifent Clear,...
% 01=Probably Clear, 10=Probably Cloudy, 11=Confident Cloudy)...
% - 'Confident Cloudy' observations
CM_cld_c(CM_vld & bitget(CMblk, 8, 'uint16') == 1 & bitget(CMblk, 7, 'uint16') == 1) = 1;
% - 'Probably Cloudy' observations
CM_cld_p(CM_vld & bitget(CMblk, 8, 'uint16') == 1 & bitget(CMblk, 7, 'uint16') == 0) = 1;
% 'Snow/ice' pixels
CM_snw(CM_vld & bitget(CMblk, 11,'uint16') == 1) = 1;
% - Cirrus cloud % Cirrus Detection (IR) (BTM15-BTM16) (1=Cloud, 0=No
% Cloud)
CM_crs(CM_vld & bitget(CMblk, 10, 'uint16') == 1) = 1;

%% Summarize the useful masks
% Valid pixel mask
Phy_blk = DNBblk_BRDF < 30000 ... % valid Lunar-BRDF-corrected DNB values (excluding saturation obervations >3000 nW*cm-2*sr-1)
    & MQblk < 2 ... % good DNB quality
    & QFDNBblk == 0 ... % good Lunar-BRDF-correction retrieval quality
    & VZAblk >= vza_lmt(1, 1)*100 & VZAblk < vza_lmt(end, 2)*100; % valid VZA range
% Snow mask (1=snow, 0=no snow)
Snw_blk = CM_snw | SFblk; % Snow pixels labled in both Raw and Lunar-BRDF-corrected Black Marble products

% Clear pixel mask
% Cloud/no cloud mask (1=cloud, 0=no cloud)
CM_cld_blk = CM_nght & CM_vld & (CM_cld_c | CM_cld_p | CM_crs); % Including confidence cloudy, probably cloudy, snow/ice, and cirrus cloud pixels.
% Cloud & snow buffered mask (Buffering the edges of cloud & snow)
buf_win = 1+2*win_buf; % Get the buffer window
buf_blk = CM_cld_blk | Snw_blk == 1 | MQblk == 2;
buf_blk = imdilate(buf_blk, true(buf_win));
buf_blk = ~buf_blk;
% Mask for the valid obs. for COLD
Clr_blk = Phy_blk & buf_blk;
DNBblk_BRDF(~Clr_blk) = 65535;

VZA_grp_blk = zeros(size(VZAblk));
for i_grp = 1:size(vza_lmt, 1)
    VZA_grp_blk(VZAblk >= vza_lmt(i_grp, 1)*100 & VZAblk < vza_lmt(i_grp, 2)*100) = i_grp;
end

% Loop stack for the row array
% Get the row array
irows = row1+nlines_read.*(0:ceil((row2-row1)/nlines_read)-1);
for i_rows = 1:length(irows)
    %% Get the saveing directories
    n_row1 = irows(i_rows);
    n_row2 = n_row1+nlines_read-1;
    folderpath_stack_rows = fullfile(folderpath_stack, sprintf('%04d%04d', n_row1, n_row2));
    % make the storage folder for the rows
    if ~isfolder(folderpath_stack_rows)
        mkdir(folderpath_stack_rows);
    end
    filename_stack = sprintf('map_tsdata_buffer_%04d%04d_%s.mat', ...
        n_row1, n_row2, yeardoy);
    % Skip if stack data already exist
    if isfile(fullfile(folderpath_stack_rows, filename_stack))
        % Summarize the data number file if not existing
        load(fullfile(folderpath_stack_rows, filename_stack), 'line_tbrdf');
        num_data = line_tbrdf < 65535;
        num_data = sum(num_data(:));

        % Save the available data number information
        dir_datanum = fullfile(folderpath_stack_rows, 'DataNum');
        if ~isfolder(dir_datanum)
            mkdir(dir_datanum)
        end
        filename_obsnum = sprintf('DataNum.%04d%04d.%s.%09d.txt', n_row1, n_row2, yeardoy, num_data);
        fid = fopen(fullfile(dir_datanum, filename_obsnum), 'w');
        fclose(fid);

        continue
    end

    %% Extract records for the target rows
    blk_numrow = n_row1:n_row2;
    line_tbrdf = uint16(DNBblk_BRDF(blk_numrow, :)); % Lunar-BRDF-corrected DNB records (uint16)
    line_tvzagrp = uint8(VZA_grp_blk(blk_numrow, :)); % Viewing zenith angle group record (uint8)

    %% Save the stack file
    % If no clear data at all, save an empty 'nulldata' file
    num_data = line_tbrdf < 65535;
    num_data = sum(num_data(:));

    % Save the available data number information
    dir_datanum = fullfile(folderpath_stack_rows, 'DataNum');
    if ~isfolder(dir_datanum)
        mkdir(dir_datanum)
    end
    filename_obsnum = sprintf('DataNum.%04d%04d.%s.%09d.txt', n_row1, n_row2, yeardoy, num_data);
    fid = fopen(fullfile(dir_datanum, filename_obsnum), 'w');
    fclose(fid);

    if num_data == 0
        nulldata = [];
        filename_stack = sprintf('map_tsdata_buffer_%04d%04d_%s_null.mat', ...
            n_row1, n_row2, yeardoy);
        save(fullfile(folderpath_stack_rows, filename_stack), 'nulldata');
        
        continue
    end
    
    % For normal stack data
    save(fullfile(folderpath_stack_rows, [filename_stack, '.temp']), ...
        'line_tbrdf', 'line_tvzagrp', '-v7', '-nocompression');
    movefile(fullfile(folderpath_stack_rows, [filename_stack, '.temp']), ...
        fullfile(folderpath_stack_rows, filename_stack));
end
clear line_tbrdf line_tvza line_tclr line_tsn DNBblk_BRDF ...
    SFblk MQblk QFDNBblk CMblk VZAblk VZA_grp_blk CM_vld CM_nght CM_cld_c...
    CM_cld_p CM_crs CM_snw CM_clr Phy_blk Snw_blk CM_cld_blk Clr_blk
end

