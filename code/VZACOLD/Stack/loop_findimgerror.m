function loop_findimgerror(job, njobs, varargin)
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

% The top directories of the source and desitination
dir_save = '/shared/zhulab/Tian/DataQuality';
if ~isfolder(dir_save)
    mkdir(dir_save)
end

% Set the input parameters for each of the task
tasks = [];
i_task = 0;
for tile_name = tile_list
    filename_save = sprintf('DataQualityRec_%s.mat', tile_name);
    if isfile(fullfile(dir_save, filename_save))
        continue
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
    % Get dates of available lunar-BRDF-corrected Black Marble data
    folderpath_NTLBRDF = sprintf('/shared/cn449/VIIRS_NTL/BlackMarble_BRDFcorrected/%s/', tile_name);
    imf_brdf = dir(fullfile(folderpath_NTLBRDF, 'VNP46A2*.h5'));
    fname_list = {imf_brdf.name}';
    if isempty(fname_list)
        continue
    end
    filename_save = sprintf('DataQualityRec_%s.mat', tile_name);
    dir_save = fullfile(dir_save, filename_save);
    
    fprintf('Processing %s\n', dir_save)

    % Check the stack list for the target tile & year
    find_imgerror_doy(dir_save, fname_list);
end

function find_imgerror_doy(dir_save, fname_list)
    if isfile(dir_save)
        return
    end
    range_raw = [0 10; 10 20; 20 30; 30 40; 40 50; 50 60; 60 70; 70 80; 80 90; 90 100];
    range = [20 70];
    T = [];
    i_T = 0;
    fld = '/HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/';
    dims = [2400, 2400];% dimention of extracted data
    startpxl = [0, 0];% starting pixel location
    for i = 1:length(fname_list)
        fname = fname_list(i);
        tilename = regexpi(fname, 'h(\w*)v(\w*)', 'match');
        tilename = [tilename{:}];
        tilename = string(tilename);
        dir_data_tile = fullfile('/shared/cn449/VIIRS_NTL/BlackMarble_BRDFcorrected/', tilename);
        path_data = fullfile(dir_data_tile, fname);

        % Skip if data is broken
        try
            file_id = H5F.open(path_data, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
            % BRDF Corrected DNB Radiance
            data_name = [fld, 'DNB_BRDF-Corrected_NTL'];
            DNBblk_BRDF = read_hdf5(file_id, data_name, dims, startpxl);
            % Snow Flag
            data_name = [fld, 'Snow_Flag'];
            SFblk = read_hdf5(file_id, data_name, dims, startpxl);
            % Mandatory quality Flag
            data_name = [fld, 'Mandatory_Quality_Flag'];
            MQblk = read_hdf5(file_id, data_name, dims, startpxl);
            % Cloud Mask Status
            data_name = [fld, 'QF_Cloud_Mask'];
            CMblk = read_hdf5(file_id, data_name, dims, startpxl);

            H5F.close(file_id)
        catch
            fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
            fprintf('error data %s\n', path_data)
            fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
            continue
        end
        fprintf('Checking %s\n', tilename)

        % Initialize the records
        CM_nght = zeros(size(CMblk), 'logical');
        CM_cld_c = zeros(size(CMblk), 'logical');
        CM_cld_p = zeros(size(CMblk), 'logical');
        CM_crs = zeros(size(CMblk), 'logical');
        %     CM_sdw = zeros(size(CMblk), 'logical');
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

        total_pix = DNBblk_BRDF > 0 & DNBblk_BRDF < 30000;
        cld_pix = CM_nght & CM_vld & (CM_cld_c | CM_cld_p | CM_crs);
        cldsnw_pix = cld_pix | SFblk == 1 | CM_snw;
        clear_pix = total_pix & MQblk < 2 & ~cldsnw_pix;
        total_rangepix = total_pix & DNBblk_BRDF >= range(1) & ...
            DNBblk_BRDF <= range(2);
        clr_rangepix = clear_pix & DNBblk_BRDF >= range(1) & ...
            DNBblk_BRDF <= range(2);
        
        i_T = i_T+1;
        T(i_T).name = fname;
        T(i_T).dir_data = dir_data_tile;
        T(i_T).total_pix = sum(total_pix(:));
        T(i_T).total_rangepix = sum(total_rangepix(:));
        T(i_T).total_rangepix_raw = [];
        for i_range = 1:size(range_raw, 1)
            tmp = total_pix & DNBblk_BRDF >= range_raw(i_range, 1) & ...
                DNBblk_BRDF < range_raw(i_range, 2);
            T(i_T).total_rangepix_raw = [T(i_T).total_rangepix_raw sum(tmp(:))];
        end

        T(i_T).clear_pix = sum(clear_pix(:));
        T(i_T).clr_rangepix = sum(clr_rangepix(:));
        T(i_T).clr_rangepix_raw = [];
        for i_range = 1:size(range_raw, 1)
            tmp = clear_pix & DNBblk_BRDF >= range_raw(i_range, 1) & ...
                DNBblk_BRDF < range_raw(i_range, 2);
            T(i_T).clr_rangepix_raw = [T(i_T).clr_rangepix_raw sum(tmp(:))];
        end

        DNB = DNBblk_BRDF(DNBblk_BRDF <= 100);
        DNB = double(DNB(:));
        % Skip if too short to perform the unimodal/bimodal test
        if length(DNB) >= 15
            [p_value, dip, xl, xu] = dipTest(DNB, 1);
            T(i_T).dip_all = [p_value, dip, xl, xu];
        end
        
        DNB_clr = DNBblk_BRDF(clear_pix & DNBblk_BRDF <= 100);
        DNB_clr = double(DNB_clr(:));
        % Skip if too short to perform the unimodal/bimodal test
        if length(DNB_clr) >= 15
            [p_value, dip, xl, xu] = dipTest(DNB_clr, 1);
            T(i_T).dip_clr = [p_value, dip, xl, xu];
        end
        
        DNB = DNBblk_BRDF(~cldsnw_pix & DNBblk_BRDF <= 100 & DNBblk_BRDF > 0);
        DNB = double(DNB(:));
        T(i_T).mean = mean(DNB);

        DNB = DNBblk_BRDF(~cldsnw_pix & DNBblk_BRDF < 30000 & DNBblk_BRDF > 0);
        DNB = double(DNB(:));
        T(i_T).median = median(DNB);
    end
    
    save(dir_save, 'T')
    fprintf('done !\n')
end

end
