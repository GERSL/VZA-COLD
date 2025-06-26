%average daily DNB, read by block, DNB multipel 10 to change it's range to
%1-10,000 to fit the ramda value
function rec_cg = TrendSeasonalFit_map_VZACOLD(...
    dir_BRDF_new, n_rst, nlines_read, irow, n_irows, col1, ...
    col2, tile_name, T_cg, Tmax_cg, conse, num_c, vzaints, ...
    date_start, date_end, num_toler, update_model, NTL_lmt)

warning ('off', 'all');

% number of VZA groups
num_vza = size(vzaints, 1);
num_mod = num_vza+1;
% number of clear observation / number of coefficients
n_times = 12; 
% minimum year for model intialization
mini_yrs = 1;
lmt_yr = 4;
% threshold (degree) of mean included angle
nsign = 45;
% Tmasking of noise
Tmax_cg = norminv(1-(1-Tmax_cg)/2);
% adjust threshold based on normal distribution
T_cg = norminv(1-(1-T_cg)/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Get ready for Xs & Ys

%% Read the n*row time series data
% The year and doy for the start and end time
yr_s = floor(date_start/1000);
doy_s = mod(date_start, 1000);
yr_e = floor(date_end/1000);
doy_e = mod(date_end, 1000);
date_start = datenum(yr_s, 1, doy_s);
date_end = datenum(yr_e, 1, doy_e);
sdate = (date_start: date_end)';

tic
[sdate, line_tbrdf_mrg, line_tvzagrp_mrg] = createMergeRowdata(tile_name, n_irows, sdate);
fprintf('!! Finished reading the %f th 10 rows for tile %s in %0.2f mins\n', ...
    n_irows, tile_name, toc/60);

% Loop running the line(s)
for irow_ids = 1:nlines_read
    n_row = irow_ids+(n_irows-1)*nlines_read;
    % check if row already exist
    filename_save = sprintf('%s_record_change%04d.mat', tile_name, n_row);
    filename_save_null = sprintf('%s_record_change%04d_null.mat', tile_name, n_row);
    if (isfile(fullfile(dir_BRDF_new, n_rst, filename_save)) || ...
            isfile(fullfile(dir_BRDF_new, n_rst, filename_save_null))) && ...
            irow == 0
        continue
    end

    if isempty(sdate)
        rec_cg = [];
        save(fullfile(dir_BRDF_new, n_rst, [filename_save_null, '.temp']), 'rec_cg');
        movefile(fullfile(dir_BRDF_new, n_rst, [filename_save_null, '.temp']), fullfile(dir_BRDF_new, n_rst, filename_save_null));
        continue
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Get ready for Xs & Ys
    % initialize NUM of Functional Curves
    num_fc = 0;
    % initialize the struct data of RECording of ChanGe (rec_cg)
    rec_cg = struct('t_start', [], 't_end', [], 't_break', [], 'coefs', [], ...
        'rmse', [], 'pos', [], 'change_prob', [], 'num_obs', [], 'category', ...
        [], 'magnitude', [], 'adj_rmse', [], 'rec_rmse', [], 'rec_normdif', [], ...
        'clrx', [], 'rec_magdif', [], 'clrx_sdate', [], 'num_yrs', []);
    % Loop COLD for the columns
    for icol_ids = col1:col2
        % Get the data 
        y_all = line_tbrdf_mrg(irow_ids, icol_ids, :);
        y_all = y_all(:);
        vzagrp_all = line_tvzagrp_mrg(irow_ids, icol_ids, :);
        vzagrp_all = vzagrp_all(:);

        % Obtaining the valid obs. ids for the modeling
        idgood = y_all < 65535; % clear pixel IDs
        %% Case 1: All time series data do not have enough clear observations for change detection
        if sum(idrange_all) < n_times*num_c % not enough for all data
            continue;
        end

        %% Case 2: Stratified zenith angle data to detect changes
        % Xs & Ys for computation after removing noises such as cloud
        clrx = sdate(idgood);
        clry = double(y_all(idgood));
        clrvzaints = double(vzagrp_all(idgood)); % 1, 2, 3 means different angle range's data
        % Skip if not enough clear pixels, all data time series
        if length(clrx) < n_times*num_c
            continue;
        end  
        
        %% Skip if is consistent dark pixel
        check_dark = 1;
        for i_y = 1:length(clry)-conse+1
            if median(clry(i_y:i_y+conse-1)) > NTL_lmt*10
                check_dark = 0;
                break
            end
        end

        if check_dark == 1
            continue
        end

        %% Find out the annual period
        clrx_sdate = clrx; % keep a record of the actual dates
        clrx_doy = day(datetime(clrx_sdate, 'ConvertFrom', 'datenum'), 'dayofyear');
        ann_doy = unique(clrx_doy);
        num_yrs = length(ann_doy);
        if num_yrs > 365
            num_yrs = 365.25;
        end
        
        doy_lack = clrx(1):clrx(end);
        index_date = ismember(doy_lack, clrx);
        doy_lack = day(datetime(doy_lack, 'ConvertFrom', 'datenum'), 'dayofyear');
        index_lack = ~ismember(doy_lack, ann_doy);
        doy_lack = double(index_lack)*-1;
        doy_lack = cumsum(doy_lack);
        doy_lack = doy_lack(index_date);
        clrx = clrx+doy_lack';
        
        % !!! Start here !!!
        % !!!
        % !!!
        % !!!
        
        % Initial the model parameters
        qa = 0;
        fit_cft = nan(num_c, num_mod);
        rmse = nan(1, num_mod);
        fit_num = nan(1, num_mod);
        pre_end = nan(1, num_mod);
        % one variable, all data list, which will benefit code organizing
        clr_data = [clrx, clry, clrvzaints, clrx_sdate];
        rec_normdif  = zeros(size(clr_data, 1), num_mod);
        rec_magdif = rec_normdif;
        rec_rmse  = zeros(size(clr_data, 1), num_mod, 2); % [mod_rmse, tmp_rmse]

        % adjust RMSE for each stratifed data
        adj_rmse = zeros(num_mod, 1); % +1 means all data
        for nvza = 1:num_vza
            clry_vza = clr_data(clr_data(:, 3) == nvza, 2);
            if length(clry_vza) < 2
                continue
            end
            adj_rmse(nvza) = median(abs(clry_vza(2:end)-clry_vza(1:end-1)), 1);
        end

        % adjust RMSE for all data for Tmask
        adj_rmse_all = median(abs(clr_data(2:end, 2)-clr_data(1:end-1, 2)), 1);
        adj_rmse(end) = adj_rmse_all;
        clear clry_vza;

        % the first observation for TSFit for each zenith angle data
        i_start = ones(num_mod, 1); % the all data is 1
        for nvza = 1:num_vza
            try
                i_start(nvza) = find(clr_data(:, 3) == nvza, 1);
            catch
            end
        end

        % record the start of the model initialization (0=>initial;1=>done)
        BL_train = zeros(num_mod, 1);
        i_count = zeros(num_mod, 1);
        % identified and move on for the next curve
        num_fc = num_fc+1; % NUM of Fitted Curves (num_fc)
        % record the num_fc at the beginning of each pixel
        rec_fc = num_fc;
        % initialize i_dense
        i_dense = ones(num_mod, 1);
        % start with the miminum requirement of clear obs for all data
        i = n_times*num_c; % first i  goes further

        % record the difference of conse observations, all records same format
        % dimension 1: v_dif_mag % record the magnitude of change
        % dimension 2: v_dif     % value of difference for conse obs
        % dimension 3: vec_mag   % chagne vector magnitude
        % dimension 4: index of observation in orginal datalist
        rec_dif = nan(4, conse, num_mod);

        % while loop - process till the last clear observation - conse
        while i <= size(clr_data, 1)-conse
            % Get the ID of vza interval of the i's observation
            n_id = [num_mod clr_data(i, 3)]; % include all data & the zenith data's demension in "clr_data"
            changeconfirmed = 0; % label change confirmed, and there is no change at default

            % "all data" and "zenith angle data" will be checked at the same time
            for n = n_id % n indicates the ID of stratified zenith data, and all data has higher priority
                % if change has been confirmed by all data, do not repeat zenith angle data
                if changeconfirmed > 0
                    continue;
                end

                % set input dataset all the same between "all data" and
                % "zenith angle data", which can share same process code below
                if n > num_vza
                    clr_data_n = clr_data;
                    clr_data_n(:, 3) = n; % all in the 4th stratifed data
                else
                    clr_data_n = clr_data; % do not change the orginal zenith angle ID
                end

                % remaining observations: if remaining observations less than conse, next process
                if sum(clr_data_n(i+1:end, 3) == n) < conse
                    continue;
                end

                % basic requrirements: 1) enough observations; 2) enough time
                if ~basicRequire(clr_data_n, i_start, i, n, num_c, mini_yrs, n_times, lmt_yr)
                    continue; % if it fails, we will move to next observation
                end

                % initializing model
                if BL_train(n) == 0
                    [clrx_n, clry_n, clrx_sdate_n] = findStratifiedData(clr_data_n, n, i_start(n), i, lmt_yr);  % select the stratified data x, accordingly
                    % check max time difference
                    if max(clrx_sdate_n(2:end)-clrx_sdate_n(1:end-1)) > num_yrs % max time difference
                        i_start(n)  = moveStartIndex(clr_data_n, i_start, n); % start from next clear obs ONLY from the stratified data list
                        i_dense(n) = i_start(n); % i that is dense enough
                        continue;
                    end

                    % initial model fit
                    [fit_cft(:, n), rmse(n), rec_v_dif, ~] = ...
                        autoRobustFit_lmtyr_dynT(clrx_n, clry_n, num_c, num_yrs);
                    fit_num(n) = length(clrx_n);
                    pre_end(n) = i;

                    % record the ith RMSEs and normalized differences for check
                    rec_rmse(i+1, n, 1) = rmse(n); 
                    [isstable, rec_normdif(i+1, n)] = isStableModel(rec_v_dif, fit_cft(:, n), adj_rmse(n), rmse(n), clrx_n, T_cg);
                    rec_magdif(i+1, n) = rec_v_dif(end);

                    % stable model requirement
                    if isstable
                        BL_train(n) = 1;
                    else
                        fit_cft(:, n) = nan(num_c, 1); % withdraw model
                        rmse(n) = nan; % withdraw model's rmse
                        i_start(n)  = moveStartIndex(clr_data_n, i_start, n); % start from next clear obs
                        continue;
                    end

                    % backward change detection if model ready
                    i_count(n) = 0;

                    % find the previous break point
                    if num_fc == rec_fc
                        i_break = 1; % first curve
                    else
                        i_break = find(clr_data_n(:, 4) >= rec_cg(num_fc-1).t_break, 1); % after the first curve
                    end

                    % backward change detection
                    changeconfirmed_back = 0;
                    if i_start(n) > i_break % if observations are available between i_break and i_start(n)
                        % loop observations from i_start(n) -1 back to i_break to confirm its change

                        % find out the remaining conse observations
                        ids_check = find(clr_data_n(:, 3) == n);
                        ids_check = ids_check(i_break <= ids_check & ids_check <= i_start(n)-1); % back to real ID in all datalist, this code can select the real IDs that we need to check
                        if ~isempty(ids_check) % sometimes no data
                            ids_check = ids_check(end:-1:1); % backfarward e.g., 8, 7, 6, 5, 4, 3, 2, 1
                            for itmp = 1:length(ids_check) % loop the observation
                                i_ini = ids_check(itmp); % real ID in all datalist
                                ids_conse = ids_check(itmp:length(ids_check)); % observations that will be compared
                                ids_conse = ids_conse(1:min(conse, length(ids_conse))); % not exceed conse

                                % compare observation and predication
                                rec_dif_back = diffConseObs(n, clr_data_n, fit_cft, rmse, adj_rmse, ids_conse, num_mod, conse, num_yrs);
                                rec_normdif(i_ini-min(conse, length(ids_conse))+1:i_ini, n) = rec_dif_back(3, 1:min(conse, length(ids_conse)), n);

                                % confirm change
                                rec_rmse(i_ini-min(conse, length(ids_conse))+1:i_ini, n, 1)  = rmse(n); 
                                overall = 0.1*(fit_cft(1, n)+fit_cft(2, n)*i_start(n));
                                [changeconfirmed_back, v_dif_mag, rec_normdif(i_ini, n)] = ...
                                    confirmChange(n, rec_dif_back, num_toler, T_cg, Tmax_cg, NTL_lmt, overall);
                                rec_magdif(i_ini-min(conse, length(ids_conse))+1:i_ini, n) = rec_dif_back(1, 1:min(conse, length(ids_conse)), n);

                                if changeconfirmed_back > 0
                                    break; % stop detecting change backward
                                elseif changeconfirmed_back == -1
                                    clr_data(i_ini, :) = []; % remove noise
                                    rec_normdif(i_ini, :) = [];
                                    rec_rmse(i_ini, :, :) = [];
                                    rec_magdif(i_ini, :) = [];

                                    rec_dif(4, :, :) = rec_dif(4, :, :)-1; % remove index in the record of difference in the orginal datalist
                                    i = i-1; % stay & check again after noise removal
                                end
                                % update new_i_start if i_ini is not a confirmed break
                                i_start(n) = i_ini; % may other stratified data varied
                                %% end of confirming change
                            end
                        end
                    end

                    % fit the fisrt curve if enough observations remaining
                    % only fit first curve if 1) have more than
                    % conse obs 2) previous obs is more than a year
                    if changeconfirmed_back > 0 && num_fc == rec_fc ...
                            && clr_data_n(i_start(n)-1, 4)-clr_data_n(i_dense(n), 4)+1 >= 365.25 % previous obs is more than a year
                        % move the i_dense to within 4 years data
                        if clr_data_n(i_dense(n), 4) < clr_data_n(i_start(n), 4)-lmt_yr*365.25
                            tmp_idense = find(clr_data_n(:, 4) >= clr_data_n(i_start(n), 4)-lmt_yr*365.25);
                            tmp_idense = tmp_idense(1);
                        else
                            tmp_idense = i_dense(n);
                        end

                        % find out the remaining conse observations
                        ids_check = find(clr_data_n(:, 3) == n);

                        % 1) Skip if have less than conse obs
                        % if sum(i_dense(n) <= ids_check & ids_check <= i_start(n)) >= conse
                        if sum(tmp_idense <= ids_check & ids_check <= i_start(n)-1) < conse
                            continue
                        end

                        fit_cft_disturb = nan(num_c, num_mod);
                        rmse_disturb = nan(1, num_mod);
                        qa = 10;
                        [fit_cft_disturb(:, end), rmse_disturb(end)] = ...
                            autoRobustFit_lmtyr_dynT(clr_data_n(tmp_idense:i_start(n)-1, 1), ...
                            clr_data_n(tmp_idense:i_start(n)-1, 2), num_c, num_yrs);
                        
                        % check the data at each stratified
                        for nvza =1:num_vza
                            [clrx_n, clry_n] = findStratifiedData(clr_data, nvza, i_dense(n), i_start(n)-1, lmt_yr);
                            if length(clrx_n) < conse % if less than conse obs, no model fitting
                                continue;
                            end
                            [fit_cft_disturb(:, nvza), rmse_disturb(nvza)] = ...
                                autoRobustFit2_lmtyr_dynT(clrx_n, clry_n, num_c, num_yrs);
                        end

                        % if model available, just record it or them
                        if sum(~isnan(fit_cft_disturb(:))) > 0
                            % record time of curve end
                            rec_cg(num_fc).t_end = clr_data_n(i_start(n)-1, 4);
                            % record postion of the pixel
                            rec_cg(num_fc).pos = icol_ids;
                            % record fitted coefficients
                            rec_cg(num_fc).coefs = fit_cft_disturb;
                            % record rmse of the pixel
                            rec_cg(num_fc).rmse = rmse_disturb;
                            % record break time
                            rec_cg(num_fc).t_break = clr_data(i_start(n), 4);
                            % record change probability
                            rec_cg(num_fc).change_prob = changeconfirmed_back;
                            % record time of curve start
                            rec_cg(num_fc).t_start = clr_data(i_dense(n), 4);
                            % record fit category
                            rec_cg(num_fc).category = qa+num_c;
                            % record number of observations
                            rec_cg(num_fc).num_obs = i_start(n)-i_dense;
                            % record change magnitude
                            rec_cg(num_fc).magnitude = -median(v_dif_mag, 'omitnan');
                            % record the RMSEs and normalized differences
                            rec_cg(num_fc).adj_rmse = adj_rmse;
                            rec_cg(num_fc).rec_rmse = rec_rmse;
                            rec_cg(num_fc).rec_normdif = rec_normdif;
                            rec_cg(num_fc).clrx = clr_data(:, 1);
                            rec_cg(num_fc).clrx_sdate = clr_data(:, 4);
                            rec_cg(num_fc).rec_magdif = rec_magdif;
                            rec_cg(num_fc).num_yrs = num_yrs;

                            % identified and move on for the next functional curve
                            num_fc = num_fc+1;
                            i_start(:) = i_start(n); % unique i_start because of new change is confirmed
                            i_dense(:) = i; % update to record i_dense
                        end
                        clear fit_cft_disturb rmse_disturb;
                    end
                end % end of inialization

                % continous change detection
                [clrx_n, clry_n] = findStratifiedData(clr_data_n, n, i_start(n), i, lmt_yr);
                update_num_c = num_c;
                % initial model fit when there are not many obs
                if  i_count(n) == 0 || length(clrx_n) >= num_c*n_times && length(clrx_n) > num_c
                    % update i_count at each interation
                    i_count(n) = sum(clr_data_n(i_start(n):i, 3) == n);
                    % update i_count at each interation
                    [fit_cft(:, n), rmse(n), ~] = ...
                        autoRobustFit_lmtyr_dynT(clrx_n, clry_n, num_c, num_yrs); % update the time series model
                    
                    fit_num(n) = length(clrx_n);
                    pre_end(n) = i;

                    tmpcg_rmse = 0;
                    rec_rmse(i, n, 2)  = tmpcg_rmse;

                    %% start of updating record
                    % record time of curve start
                    rec_cg(num_fc).t_start = clr_data(min(i_start), 4);
                    % record time of curve end
                    rec_cg(num_fc).t_end = clr_data(i, 4);
                    % record break time
                    rec_cg(num_fc).t_break = 0; % no break at the moment
                    % record postion of the pixel
                    rec_cg(num_fc).pos = icol_ids;
                    % record fitted coefficients
                    rec_cg(num_fc).coefs = fit_cft;
                    % record rmse of the pixel
                    rec_cg(num_fc).rmse = rmse;
                    % record change probability
                    rec_cg(num_fc).change_prob = 0;
                    % record number of observations
                    rec_cg(num_fc).num_obs = i-min(i_start)+1;
                    % record fit category
                    rec_cg(num_fc).category = qa+0; % no change at default
                    % record change magnitude
                    rec_cg(num_fc).magnitude = 0;
                    % record the RMSEs and normalized differences
                    rec_cg(num_fc).adj_rmse = adj_rmse;
                    rec_cg(num_fc).rec_rmse = rec_rmse;
                    rec_cg(num_fc).rec_normdif = rec_normdif;
                    rec_cg(num_fc).clrx = clr_data(:, 1);
                    rec_cg(num_fc).clrx_sdate = clr_data(:, 4);
                    rec_cg(num_fc).rec_magdif = rec_magdif;
                    rec_cg(num_fc).num_yrs = num_yrs;

                    %% end of updating record
                    % find out the remaining conse observations
                    ids_conse = find(clr_data_n(:, 3) == n);
                    ids_conse = ids_conse(find(ids_conse > i, conse)); % the ID after i at the current dataset

                    if length(ids_conse) == conse % enough conse data
                        % compare observation and predication
                        rec_dif_n = diffConseObs(n, clr_data_n, fit_cft, rmse, adj_rmse, ids_conse, num_mod, conse, num_yrs);
                        rec_dif(:, :, n) = rec_dif_n(:, :, n); % update to record
                    end
                else
                    % Skip for n% of obs., close to lower value
                    num_skip = fix(fit_num(n)*update_model); 
                    if sum(clr_data_n(pre_end(n):i, 3) == n)-1 >= num_skip % Update model if the skip condition is matched
                        % update i_count at each interation
                        i_count(n) = sum(clr_data_n(i_start(n):i, 3) == n);
                        % update the model fitting
                        [fit_cft(:, n), rmse(n)] = ...
                            autoRobustFit_lmtyr_dynT(clrx_n, clry_n, num_c, num_yrs);
                        
                        fit_num(n) = length(clrx_n);
                        pre_end(n) = i;

                        % record fitted coefficients
                        rec_cg(num_fc).coefs = fit_cft;
                        % record rmse of the pixel
                        rec_cg(num_fc).rmse = rmse;
                        % record number of observations
                        rec_cg(num_fc).num_obs =  i-min(i_start)+1;
                        % record fit category
                        rec_cg(num_fc).category = qa+0;
                    end

                    % record time of curve end
                    rec_cg(num_fc).t_end = clr_data_n(i, 4);
                    % record the RMSEs and normalized differences
                    rec_cg(num_fc).adj_rmse = adj_rmse;
                    rec_cg(num_fc).rec_rmse = rec_rmse;
                    rec_cg(num_fc).rec_normdif = rec_normdif;
                    rec_cg(num_fc).clrx = clr_data(:, 1);
                    rec_cg(num_fc).clrx_sdate = clr_data(:, 4);
                    rec_cg(num_fc).rec_magdif = rec_magdif;
                    rec_cg(num_fc).num_yrs = num_yrs;

                    % move the ith col to i-1th col at zenith angle data
                    i_end_n = rec_dif(4, end, n); % record the index of observation
                    % find out next conse observation in the zenith angle data
                    rec_dif(:, 1:end-1, n) = rec_dif(:, 2:end, n);
                    rec_dif(:, end, n) = nan;
                    if ~isnan(i_end_n) % no comparation at last
                        i_end_n = find(clr_data_n(:, 1) > clr_data_n(i_end_n, 1) & clr_data_n(:, 3) == n, 1); % move the last observation in this zenith angle datalist to the next one
                    end
                    if ~isempty(i_end_n) && ~isnan(i_end_n) % to the end
                        % all data new comparision
                        % absolute difference
                        rec_dif(1, end, n) = clr_data_n(i_end_n, 2)-...
                            autoTSPred_OLS_dynT(clr_data_n(i_end_n, 1), ...
                            fit_cft(:, clr_data_n(i_end_n, 3)), num_yrs); % i_conse may be different from n

                        % minimum rmse
                        mini_rmse = max([adj_rmse(clr_data_n(i_end_n, 3)), rmse(clr_data_n(i_end_n, 3)), tmpcg_rmse]); %  must be the same as i + conse
                        % z-scores
                        rec_dif(2, end, n) = rec_dif(1, end, n)/mini_rmse;
                        rec_dif(3, end, n) = norm(rec_dif(2, end, n)); 
                        rec_dif(4, end, n) = i_end_n;
                    end
                end

                % confirm change
                rec_rmse(i+1:i+conse, n, 1)  = rmse(n); 
                rec_normdif(i+1:i+conse, n) = rec_dif(3, :, n);
                overall = 0.1*(fit_cft(1, n)+fit_cft(2, n)*clr_data_n(i, 1));
                [changeconfirmed, v_dif_mag, rec_normdif(i+1, n)] = ...
                    confirmChange(n, rec_dif, num_toler, T_cg, Tmax_cg, NTL_lmt, overall);
                rec_magdif(i+1:i+conse, n) = rec_dif(1, :, n);

                if changeconfirmed > 0
                    % record break time
                    rec_cg(num_fc).t_break = clr_data_n(i+1, 4);
                    % record change probability
                    rec_cg(num_fc).change_prob = 10+changeconfirmed;
                    %1X: change prob,
                    % X1: confirmed by continous model,
                    % X2: confirmed by the 1st zenith model
                    % X3: confirmed by the 2st zenith model
                    % X4: confirmed by the 3st zenith model
                    % record change magnitude
                    rec_cg(num_fc).magnitude = median(v_dif_mag, 'omitnan');
                    % record the RMSEs and normalized differences
                    rec_cg(num_fc).adj_rmse = adj_rmse;
                    rec_cg(num_fc).rec_rmse = rec_rmse;
                    rec_cg(num_fc).rec_normdif = rec_normdif;
                    rec_cg(num_fc).clrx = clr_data(:, 1);
                    rec_cg(num_fc).clrx_sdate = clr_data(:, 4);
                    rec_cg(num_fc).rec_magdif = rec_magdif;
                    rec_cg(num_fc).num_yrs = num_yrs;

                    % identified and move on for the next functional curve
                    num_fc = num_fc+1;
                    % start from i+1 for the next functional curve
                    i_start(:) = i+1; % all i_start new

                    % start training again
                    BL_train(1: num_mod) = 0;
                    i_count(1: num_mod) = 0;
                    % initialize the variables
                    rec_dif = nan(4, conse, num_mod); % record:
                    fit_cft = nan(num_c, num_mod);
                    rmse = nan(1, num_mod);
                    i_dense(:) = i; % update to record i_dense
                elseif changeconfirmed == -1
                    % remove noise
                    clr_data(rec_dif_back(4, 1, n), :) = [];
                    rec_normdif(rec_dif_back(4, 1, n), :) = [];
                    rec_rmse(rec_dif_back(4, 1, n), :, :) = [];
                    rec_magdif(rec_dif_back(4, 1, n), :, :) = [];

                    % remove index for diff variable
                    rec_dif(4, :, n) = rec_dif(4, :, n)-1;% index remove to index -1
                    i = i-1; % stay & check again after noise removal
                end
            end
            % move forward to the i+1th clear observation
            i = i+1;
        end % end of while iterative

        % Two ways for processing the end of the time series
        n_id = [num_mod 1:num_vza];
        for i_n = n_id % priority all data, small angle, midel angle, and large angle
            if BL_train(i_n) == 1 && i_n == num_mod %  prob based on all data
                % 1) if no break find at the end of the time series
                % define probability of change based on conse
                rec_dif_n_last = rec_dif(:, :, i_n);
                rec_dif_n_last(:, isnan(rec_dif_n_last(4, :))) = [];
                if isempty(rec_dif_n_last)
                    continue;
                end

                id_last = 0;
                for i_conse = size(rec_dif_n_last, 2):-1:1
                    % sign of change vector
                    max_angl = mean(angl(rec_dif_n_last(2, i_conse:size(rec_dif_n_last, 2))));
                    if rec_dif_n_last(3, i_conse) <= T_cg || max_angl >= nsign
                        % the last stable id
                        id_last = i_conse;
                        break;
                    end
                end

                % update change probability
                rec_cg(num_fc).change_prob = (size(rec_dif_n_last, 2)-id_last)/size(rec_dif_n_last, 2);
                % update end time of the curve
                rec_cg(num_fc).t_end = clr_data(end-size(rec_dif_n_last, 2)+id_last, 4);
                % record the RMSEs and normalized differences
                rec_cg(num_fc).adj_rmse = adj_rmse;
                rec_cg(num_fc).rec_rmse = rec_rmse;
                rec_cg(num_fc).rec_normdif = rec_normdif;
                rec_cg(num_fc).clrx = clr_data(:, 1);
                rec_cg(num_fc).clrx_sdate = clr_data(:, 4);
                rec_cg(num_fc).rec_magdif = rec_magdif;
                rec_cg(num_fc).num_yrs = num_yrs;
                if size(rec_dif, 2) > id_last % > 1
                    % update time of the probable change
                    rec_cg(num_fc).t_break = clr_data(end-size(rec_dif_n_last, 2)+id_last+1, 4);
                    % update magnitude of change
                    rec_cg(num_fc).magnitude = median(rec_dif_n_last(1, id_last+1:size(rec_dif_n_last, 2))); % based on all data
                end
            elseif BL_train(i_n) == 0
                % 2) if break find close to the end of the time series
                % Use [conse, num_c*n_times+conse) to fit curve
                % end of fit qa = 20
                qa = 20;
                if num_fc == rec_fc
                    % first curve
                    i_start = 1;
                else
                    i_start = find(clr_data(:, 4) >= rec_cg(num_fc-1).t_break);
                    i_start = i_start(1);
                end

                if i_n > num_vza
                    clr_data_n = clr_data;
                    clr_data_n(:, 3) = i_n; % all in the 4th stratifed data
                else
                    clr_data_n = clr_data; % do not change the orginal zenith angle ID
                end
                [clrx_n, clry_n] = findStratifiedData(clr_data_n, i_n, i_start, size(clr_data, 1), lmt_yr);
                clear clr_data_n;

                if length(clrx_n) >= conse
                    [fit_cft(:, i_n), rmse(i_n)] = ...
                        autoRobustFit_lmtyr_dynT(clrx_n, clry_n, num_c, num_yrs);
                    fit_num(i_n) = length(clrx_n);
                    pre_end(i_n) = size(clr_data, 1);

                    % record time of curve start
                    rec_cg(num_fc).t_start = clrx_sdate(i_start);
                    % record time of curve end
                    rec_cg(num_fc).t_end = clrx_sdate(end);
                    % record break time
                    rec_cg(num_fc).t_break = 0;
                    % record postion of the pixel
                    rec_cg(num_fc).pos = icol_ids;
                    % record fitted coefficients
                    rec_cg(num_fc).coefs = fit_cft;
                    % record rmse of the pixel
                    rec_cg(num_fc).rmse = rmse;
                    % record change probability
                    rec_cg(num_fc).change_prob = 0;
                    % record number of observations
                    rec_cg(num_fc).num_obs = length(clrx(i_start:end));
                    % record fit category
                    rec_cg(num_fc).category = qa+num_c;
                    % record change magnitude
                    rec_cg(num_fc).magnitude = 0;
                    % record the RMSEs and normalized differences
                    rec_cg(num_fc).adj_rmse = adj_rmse;
                    rec_cg(num_fc).rec_rmse = rec_rmse;
                    rec_cg(num_fc).rec_normdif = rec_normdif;
                    rec_cg(num_fc).clrx = clr_data(:, 1);
                    rec_cg(num_fc).clrx_sdate = clr_data(:, 4);
                    rec_cg(num_fc).rec_magdif = rec_magdif;
                    rec_cg(num_fc).num_yrs = num_yrs;
                end
            end
        end
    end%end of for icol_ids loop

    save(fullfile(dir_BRDF_new, n_rst, [filename_save, '.temp']), 'rec_cg');
    movefile(fullfile(dir_BRDF_new, n_rst, [filename_save, '.temp']), ...
        fullfile(dir_BRDF_new, n_rst, filename_save));
    clear clr_vza

    fprintf('Finishing row %d for tile %s in %0.2f mins\n', n_row, tile_name);
end%end of irow_ids loop
end % end of function

% function to caculate included angle between ajacent pair of change vectors
function y = angl(v_dif)
v_dif(all(isnan(v_dif), 2), :) = []; % remove nan rows
row = size(v_dif, 1);
y = zeros(row-1, 1);
if row > 1
    for i = 1:row-1
        a = v_dif(i, :);
        b = v_dif(i+1, :);
        b(isnan(a)) = [];
        a(isnan(a)) = [];

        % y measures the opposite of cos(angle)
        y(i) = acos(a*b'/(norm(a)*norm(b)));
    end
else
    y = 0;
end

% convert angle from radiance to degree
y = rad2deg(y);
end

% function to find the next i_start for the target VZA interval
function i_start_next = moveStartIndex(clr_data, i_start, n)
i_start_next = find(clr_data(:, 1) > clr_data(i_start(n), 1) & clr_data(:, 3) == n, 1);
end

function [clrx_n, clry_n, clrx_sdate_n]= findStratifiedData(clr_data, n, i_start, i_end, lmt_yr)
% This function is to find out the stratified data from the orginal data
% list, according to indexes of start and end
%
% clr_data: the 3-dimensional clear dataset
% n: the index of the stratified level, such as 1, 2, 3
% i_start: the starting index in all data list
% i_end: the ending index in all data list

% update selecting the stratified data x and y, accordingly
clrx_n = clr_data(clr_data(:, 3) == n, 1);
clrx_sdate_n = clr_data(clr_data(:, 3) == n, 4);
clry_n = clr_data(clr_data(:, 3) == n, 2);

if ~exist('i_start', 'var') && ~exist('i_end', 'var')
    return; % all data at stratified n
else
    ids_n = find(clr_data(i_start, 1) <= clrx_n & clrx_n <= clr_data(i_end, 1)); % IDs of the current period within the stratified data
    clrx_n = clrx_n(ids_n); % pick up the final data
    clrx_sdate_n = clrx_sdate_n(ids_n); % pick up the final data
    clry_n = clry_n(ids_n); % pick up the final data
end
index_lmtyr = max(clrx_sdate_n)-clrx_sdate_n <= lmt_yr*365.25;
clrx_n = clrx_n(index_lmtyr);
clrx_sdate_n = clrx_sdate_n(index_lmtyr);
clry_n = clry_n(index_lmtyr);
end


function isenough = basicRequire(clr_data, istart, i, n, num_c, mini_yrs, n_times, lmt_yr)

% examine enough observations
if n == length(istart)
    clrx_n = clr_data(istart(n):i, 1);
    clrx_sdate_n = clr_data(istart(n):i, 4);
    index_lmtyr = clrx_sdate_n(end)-clrx_sdate_n <= lmt_yr*365.25;

    clrx_n = clrx_n(index_lmtyr);
    
    i_span = length(clrx_n);
    time_span = (clrx_n(end)-clrx_n(1))./365.25;
else
    clrx_n = findStratifiedData(clr_data, n, istart(n), i, lmt_yr);
    if isempty(clrx_n) % no enough data at stratried dataset
        isenough = 0;
        return;
    end
    i_span = length(clrx_n);
    time_span = (clrx_n(end)-clrx_n(1))./365.25;
end
% 1) not enough observations;   % 2) not enough time
isenough = i_span >= n_times*num_c && time_span >= mini_yrs;
end

function [isstable, v_dif] = isStableModel(rec_v_dif, fit_cft, adj_rmse, rmse, clrx, T_cg)
% mini rmse
mini_rmse = max(adj_rmse, rmse);
% compare the first clear obs
v_start = rec_v_dif(1)/mini_rmse;
% compare the last clear observation
v_end = rec_v_dif(end)/mini_rmse;
% anormalized slope values
v_slope = fit_cft(2)*(clrx(end)-clrx(1))/mini_rmse;
% differece in model intialization
v_dif = abs(v_slope)+max(abs(v_start), abs(v_end));
v_dif = norm(v_dif); %^2;
isstable = v_dif <= T_cg; % true is stable model, false is not
end

%% function of comparing observations and predications according to conse
function rec_dif_back = diffConseObs(n, clr_data_n, fit_cft, rmse, adj_rmse, ids_conse, num_mod, conse, num_yrs)
% initialize the difference variables, that are same as normal change detection below
rec_dif_back = nan(4, conse, num_mod); % record:
% dimension 1: v_dif_mag % record the magnitude of change
% dimension 2: v_dif     % value of difference for conse obs
% dimension 3: vec_mag   % chagne vector magnitude
% dimension 4: index of observation in orginal datalist
for i = 1:length(ids_conse)
    id_now = ids_conse(i);
    v_dif_mag_tmp = clr_data_n(id_now, 2)-...
        autoTSPred_OLS_dynT(clr_data_n(id_now, 1), fit_cft(:, n), num_yrs);
    v_dif_tmp= v_dif_mag_tmp/max(adj_rmse(n), rmse(n)); % z-scores
    vec_mag_tmp = norm(v_dif_tmp); %^2;
    % record zenith angle comparisions
    rec_dif_back(1, i, n) = v_dif_mag_tmp;
    rec_dif_back(2, i, n) = v_dif_tmp;
    rec_dif_back(3, i, n) = vec_mag_tmp;
    rec_dif_back(4, i, n) = id_now;
end
end

%% function of change confirmation
function [changeconfirmed_back, v_dif_mag, v_dif] = ...
    confirmChange(n, rec_dif_back, num_toler, T_cg, Tmax_cg, NTL_lmt, overall)
% confirm change
changeconfirmed_back = 0;
v_dif_mag = [];
vec_mag_sort = rec_dif_back(3, :, n); % zenith angle data first
vec_mag_sort = sort(vec_mag_sort);
v_dif = rec_dif_back(3, 1, n);
if v_dif > T_cg && vec_mag_sort(min(1 + num_toler, length(vec_mag_sort))) > T_cg && ...
        (abs(median(rec_dif_back(1, :, n))*0.1) > NTL_lmt || overall > NTL_lmt)
        changeconfirmed_back = n; % change identified by zenith angle
        v_dif_mag = rec_dif_back(1, :, n); % record change magnitude
elseif v_dif > Tmax_cg
    changeconfirmed_back = -1; % false change identified by all data
end
end
