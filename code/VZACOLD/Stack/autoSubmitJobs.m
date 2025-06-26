function autoSubmitJobs(clusterName)
% This is to submit the job for processing NTL tiles
% DATE: Aug. 22, 2023

%% Testing setups
clusterName = 'general'; % general priority
tiles = readTileList;% VIIRS tiles
for i =  1:length(tiles)
    % Obtain the tile name of data
    tile = tiles{i};

    % Submit a dependent job, which starts stack jobs and after all the jobs are finished well, continuedly start the COLD process.
    folderpath_NTLRaw = sprintf(globalsets.FolderpathNTLRaw, tile); % locate the direct folder of the raw Black Marble images
    folderpath_NTLBRDF = sprintf(globalsets.FolderpathNTLBRDF, tile); % locate the direct folder of the Lunar-BRDF-corrected Black Marble images
    folderpath_stack  = sprintf(globalsets.FolderpathStack, tile); % locate to <StackData>

    fprintf('\r>>> Start to submit jobs to stack\r');
    folderpath_code = globalsets.FolderpathStackCode;
    submitJobStackUCONN(clusterName, tile, folderpath_code, ...
        folderpath_NTLRaw, folderpath_NTLBRDF, folderpath_stack, ...
        globalsets.NumberCoresStack);
    fprintf('>>> Finish submitting jobs to stack with %0.2f mins\r', toc/60);
end % end of tile
end

function jobid = submitJobStackUCONN(clusterName, tile, folderpath_code, ...
    folderpath_NTLRaw, folderpath_NTLBRDF, folderpath_stack, totalcores)
%% submitJobStackTTU is to submit job of stacking the Landsat images on TTU HPCC system
%
% Inputs:
% % tile: the name of the path/row
% folderpath_landsat: the directory of the landsat dataset
% filepath_reference: the file path of the reference image
% totalcores: total number of cores that will be used
%
% Return:
% jobid: The ID of the job submitted by this function
    
    % Template .sh files
    FileJobs = 'ArrayJobs2Stack.sh';
    filepath_job_template = fullfile(folderpath_code, 'job', FileJobs); % the template of the jobs
    fid_in = fopen(filepath_job_template,'a+');
    shfile=fscanf(fid_in,'%c',inf);

    % Target .sh files
    folderpath_job_record = fullfile(folderpath_code, 'job', 'log'); % the outputs will be stored in this directory
    if ~isfolder(folderpath_job_record)
        mkdir(folderpath_job_record);
    end
    FileJobsReady = sprintf('ArrayJobs2Stack_%s.sh', tile);
    path_job_ready = fullfile(folderpath_job_record, FileJobsReady);
    
    % create the job file and submit it.
    shfile_ready = strrep(shfile, 'JobName', sprintf('S_%s', tile));
    shfile_ready = strrep(shfile_ready, 'TileName', tile);
    shfile_ready = strrep(shfile_ready, 'ClusterName', clusterName);
    shfile_ready = strrep(shfile_ready, 'TotalCores', num2str(totalcores));
    shfile_ready = strrep(shfile_ready, 'FolderpathCode', folderpath_code);
    shfile_ready = strrep(shfile_ready, 'FolderpathNTLRaw', folderpath_NTLRaw);
    shfile_ready = strrep(shfile_ready, 'FolderpathNTLBRDF', folderpath_NTLBRDF);
    shfile_ready = strrep(shfile_ready, 'FolderpathStack', folderpath_stack);
    shfile_ready = strrep(shfile_ready, '%', '%%'); % print %

    %open file identifier
    fid = fopen(path_job_ready,'w');
    fprintf(fid, shfile_ready);
    fclose(fid);
    record_folder = pwd; % recording the current folder
    cd(folderpath_job_record);
    [~, jobid] = system(sprintf('sbatch %s', FileJobsReady));
    cd(record_folder); % go back to the code package's folder
    
    jobid= split(jobid, ' ');
    jobid = str2double(jobid{end});
    
    fprintf('Stacking jobs for %s with %d cores (ID: %d) has been submited successfully\r', tile, totalcores, jobid);
end
