function changemetricmap_Row(dir_chgmap, dir_chgrec, n_rst_cold, n_rst_save, ...
    tile_name, n_rst_rec, n_rst_map, date_start, date_end, n_irows, ...
    nlines_read, col1, col2, sgn, T_cg, varargin)
% Count abrupt and transition changes

warning('off', 'all')

%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'num_c', 4); 
addParameter(p, 'NTL_lmt', 1); 
addParameter(p, 'Grad_lmt', 0.1); 
addParameter(p, 'num_vza', 4); 
addParameter(p, 'n_times', 12); % number of clear observation / number of coefficients
addParameter(p, 'num_DOY', 366); % Maximum number of DOY

% request user's input
parse(p, varargin{:});
num_c = p.Results.num_c;
NTL_lmt = p.Results.NTL_lmt;
Grad_lmt = p.Results.Grad_lmt;
num_vza = p.Results.num_vza;
n_times = p.Results.n_times;
num_DOY = p.Results.num_DOY;

% Read the mannual interpreted table
dir_BRDF_old = fullfile(dir_chgrec, n_rst_cold, tile_name);
dir_BRDF_new = fullfile(dir_chgmap, n_rst_save, tile_name);
n_rst_vza = sprintf('TSFitMap_%s', tile_name); 
if ~isfolder(fullfile(dir_BRDF_new, n_rst_rec))
    mkdir(fullfile(dir_BRDF_new, n_rst_rec));
end
if ~isfolder(fullfile(dir_BRDF_new, n_rst_map))
    mkdir(fullfile(dir_BRDF_new, n_rst_map));
end

yr_s = floor(date_start/1000);
doy_s = mod(date_start, 1000);
yr_e = floor(date_end/1000);
doy_e = mod(date_end, 1000);
date_start = datenum(yr_s, 1, doy_s);
date_end = datenum(yr_e, 1, doy_e);
sdate = (date_start:date_end)';
yr_list = yr_s:yr_e;
num_yr = length(yr_list);
[sdate, line_tbrdf_mrg, line_tvzagrp_mrg] = createMergeRowdata(tile_name, n_irows, sdate);
num_col = length(col1:col2);

% Loop for the rows
for irow_ids = 1:nlines_read
    n_row = irow_ids+(n_irows-1)*nlines_read;
    filename_save = sprintf('%s_changemetric%04d.mat', tile_name, n_row);
    filename_save_null = sprintf('%s_changemetric%04d_null.mat', tile_name, n_row);
    % Skip if exist
    if (isfile(fullfile(dir_BRDF_new, n_rst_map, filename_save)) || ...
            isfile(fullfile(dir_BRDF_new, n_rst_map, filename_save_null))) && ...
            (isfile(fullfile(dir_BRDF_new, n_rst_rec, filename_save)) || ...
            isfile(fullfile(dir_BRDF_new, n_rst_rec, filename_save_null)))
        continue 
    end

    %% Read and merge mapped changes.
    filename_rec = sprintf('%s_record_change%04d.mat', tile_name, n_row);
    filename_rec_null = sprintf('%s_record_change%04d_null.mat', tile_name, n_row);

    % For the individual VZA stratified COLD methods
    if isfile(fullfile(dir_BRDF_old, n_rst_vza, filename_rec))
        load(fullfile(dir_BRDF_old, n_rst_vza, filename_rec), 'rec_cg');
    elseif isfile(fullfile(dir_BRDF_old, n_rst_vza, filename_rec_null))
        null_chg = [];

        save(fullfile(dir_BRDF_new, n_rst_rec, [filename_save_null, '.temp']), 'null_chg');
        movefile(fullfile(dir_BRDF_new, n_rst_rec, [filename_save_null, '.temp']), ...
            fullfile(dir_BRDF_new, n_rst_rec, filename_save_null));

        save(fullfile(dir_BRDF_new, n_rst_map, [filename_save_null, '.temp']), 'null_chg');
        movefile(fullfile(dir_BRDF_new, n_rst_map, [filename_save_null, '.temp']), ...
            fullfile(dir_BRDF_new, n_rst_map, filename_save_null));

        fprintf(sprintf('%s done\n', filename_save))
        continue
    else
        null_chg = [];

        save(fullfile(dir_BRDF_new, n_rst_rec, [filename_save_null, '.temp']), 'null_chg');
        movefile(fullfile(dir_BRDF_new, n_rst_rec, [filename_save_null, '.temp']), ...
            fullfile(dir_BRDF_new, n_rst_rec, filename_save_null));

        save(fullfile(dir_BRDF_new, n_rst_map, [filename_save_null, '.temp']), 'null_chg');
        movefile(fullfile(dir_BRDF_new, n_rst_map, [filename_save_null, '.temp']), ...
            fullfile(dir_BRDF_new, n_rst_map, filename_save_null));

        fprintf(sprintf('%s done\n', filename_save))
        continue
    end

    % continue if there is no model available
    l_pos = length(rec_cg);
    if l_pos == 0
        null_chg = [];
        save(fullfile(dir_BRDF_new, n_rst_rec, [filename_save_null, '.temp']), 'null_chg');
        movefile(fullfile(dir_BRDF_new, n_rst_rec, [filename_save_null, '.temp']), ...
            fullfile(dir_BRDF_new, n_rst_rec, filename_save_null));
        save(fullfile(dir_BRDF_new, n_rst_map, [filename_save_null, '.temp']), 'null_chg');
        movefile(fullfile(dir_BRDF_new, n_rst_map, [filename_save_null, '.temp']), ...
            fullfile(dir_BRDF_new, n_rst_map, filename_save_null));
        fprintf(sprintf('%s done\n', filename_save))
        continue
    end

    % initialize change row indices for the detected change
    [rec_cg(:).QA] = deal(0); % time segment QA
    [rec_cg(:).overall_start] = deal(nan(1, num_vza)); % start time overall
    [rec_cg(:).overall_end] = deal(nan(1, num_vza)); % end time overall
    [rec_cg(:).pValues] = deal(zeros(num_c, num_vza)); % significance of the slope term
    [rec_cg(:).check_notdark] = deal(false); %check dark
    [rec_cg(:).break_yr] = deal(0); % Break year
    [rec_cg(:).t_drn] = deal(0); % duration time
    [rec_cg(:).chgmag_drn] = deal(nan(1, num_vza)); % during change's magnitude
    [rec_cg(:).chgmag_aft] = deal(nan(1, num_vza)); % after change's magintude
    [rec_cg(:).chg_rate] = deal(nan(1, num_vza)); % 100*after change's magintude/end time overall magnitude
    
    %% Refit the harmonic model
    i = 1;
    while i <= length(rec_cg)
        icol_ids = rec_cg(i).pos;
        num_vza = size(rec_cg(i).coefs, 2);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
        % Get the clear observations 
        [clrx, clrx_sdate, clry, clr_vza, num_yrs, ~] = ...
            Find_ClrObs(line_tbrdf_mrg, line_tvzagrp_mrg, irow_ids, icol_ids, sdate);

        %% For the VZA-COLD identified time segment 
        % Get the valid ramge 
        [clrx_mod, clrx_sdate_mod, clry_mod, clr_vza_mod] = ...
            Find_RngObs(clrx, clrx_sdate, clry, clr_vza, rec_cg(i).t_start, rec_cg(i).t_end);
        % Refit the harmonic models
        [rec_cg(i).coefs, rec_cg(i).pValues] = RefitModel(rec_cg(i), ...
            num_vza, clrx_mod, clry_mod, clr_vza_mod, n_times, num_c, num_yrs);

        if length(clrx_mod) > n_times*num_c 
            for n_vza = 1:num_vza
                if n_vza == num_vza
                    index_vza = ones(size(clr_vza_mod), 'logical');
                else
                    index_vza = clr_vza_mod == n_vza;
                end
                if sum(index_vza) > n_times*num_c
                    rec_cg(i).QA = 1;
                    [rec_cg(i).coefs(:, n_vza), rec_cg(i).rmse(n_vza), ~, ~, rec_cg(i).pValues(:, n_vza)] = ...
                        autoRobustFit_lmtyr_dynT(clrx_mod(index_vza), clry_mod(index_vza), num_c, num_yrs);
                else
                    rec_cg(i).coefs(:, n_vza) = nan;
                    rec_cg(i).rmse(n_vza) = nan;
                    rec_cg(i).pValues(:, n_vza) = nan;
                end
            end

            % Calculate the overall values
            rec_cg(i).overall_start = 0.1*(rec_cg(i).coefs(1, :)+rec_cg(i).coefs(2, :)*rec_cg(i).t_start);
            rec_cg(i).overall_end = 0.1*(rec_cg(i).coefs(1, :)+rec_cg(i).coefs(2, :)*rec_cg(i).t_end);
        else
            rec_cg(i).QA = 2;
            rec_cg(i).coefs = nan(num_c, num_vza);
            rec_cg(i).rmse = nan(1, num_vza);
            rec_cg(i).pValues = nan(num_c, num_vza);
            rec_cg(i).overall_start = nan(1, num_vza);
            rec_cg(i).overall_end = nan(1, num_vza);
            rec_cg(i).overall_start(end) = median(0.1*clry_mod);
            rec_cg(i).overall_end(end) = median(0.1*clry_mod);
        end
        rec_cg(i).clrx = clrx_mod;
        rec_cg(i).clrx_sdate = clrx_sdate_mod;
        rec_cg(i).num_yrs = num_yrs;
        rec_cg(i).num_obs = length(clrx_mod);

        % move to the next record
        i = i+1;

        %% if is no end model fitted after the last break
        if rec_cg(i-1).t_break ~= 0 && (i > length(rec_cg) || icol_ids ~= rec_cg(i).pos)
            % Add a mandatory fit model after the break
            rec_cg(i+1:end+1) = rec_cg(i:end);
            % Initialize all fields for the mandatory end model
            fn_reccg = fieldnames(rec_cg);
            for i_fn = 1:length(fn_reccg)
                rec_cg(i).(fn_reccg{i_fn}) = 0;
            end
            rec_cg(i).coefs = nan(num_c, num_vza);
            rec_cg(i).rmse = nan(num_vza, 1);
            rec_cg(i).pValues = nan(num_c, num_vza);

            % get the valid data
            [clrx_mod, clrx_sdate_mod, clry_mod, ~] = ...
                Find_RngObs(clrx, clrx_sdate, clry, clr_vza, rec_cg(i-1).t_end+1, sdate(end));  

            % if not enough data to fit the robust model
            if length(clrx_mod) <= n_times*num_c 
                rec_cg(i).QA = 5;
                rec_cg(i).t_start = clrx_sdate_mod(1);
                rec_cg(i).overall_start = nan(1, num_vza);
                rec_cg(i).overall_end = nan(1, num_vza);
                % Use the previous model to calculate the differences
                % between obs. and predictions
                rec_cg(i).overall_start(end) = median(0.1*clry_mod);
                rec_cg(i).overall_end(end) = median(0.1*clry_mod);

            % Conduct the initial stable test
            else
                for i_aftchg = 1:length(clrx_mod)-n_times*num_c+1
                    % adjust RMSE for each stratifed data
                    adj_rmse = median(abs(clry_mod(i_aftchg+1:end)-clry_mod(i_aftchg:end-1)));
                    [fit_cft, rmse, rec_v_dif, ~] = autoRobustFit_lmtyr_dynT(clrx_mod(i_aftchg:end), clry_mod(i_aftchg:end), num_c, num_yrs);
                    % record the ith RMSEs and normalized differences for check
                    [isstable, ~] = isStableModel(rec_v_dif, fit_cft, adj_rmse, rmse, clrx_mod(i_aftchg:end), T_cg);
                    if isstable
                        rec_cg(i).QA = 3;
                        rec_cg(i).t_start = clrx_sdate_mod(i_aftchg);
                        [rec_cg(i).coefs(:, end), rec_cg(i).rmse(end), ~, ~, rec_cg(i).pValues(:, end)] = ...
                            autoRobustFit_lmtyr_dynT(clrx_mod(i_aftchg:end), clry_mod(i_aftchg:end), num_c, num_yrs);
                        break
                    elseif i_aftchg == length(clrx_mod)-n_times*num_c+1 
                        rec_cg(i).QA = 4;
                        rec_cg(i).t_start = clrx_sdate_mod(1);
                        [rec_cg(i).coefs(:, end), rec_cg(i).rmse(end), ~, ~, rec_cg(i).pValues(:, end)] = ...
                            autoRobustFit_lmtyr_dynT(clrx_mod, clry_mod, num_c, num_yrs);
                    else
                        continue
                    end
                end
                rec_cg(i).overall_start = 0.1*(rec_cg(i).coefs(1, :)+rec_cg(i).coefs(2, :)*rec_cg(i).t_start);
                rec_cg(i).overall_end = 0.1*(rec_cg(i).coefs(1, :)+rec_cg(i).coefs(2, :)*rec_cg(i).t_end);
            end
            rec_cg(i).clrx = clrx_mod;
            rec_cg(i).clrx_sdate = clrx_sdate_mod;
            rec_cg(i).num_yrs = num_yrs;
            rec_cg(i).t_end = clrx_sdate_mod(end);
            rec_cg(i).pos = icol_ids;
            rec_cg(i).num_obs = length(clrx_mod);
            
            % Move to the next record 
            i = i+1; 
        end
    end

    %% Calculate the change indices
    num_vza = 4;
    % postions
    pos = [rec_cg.pos];
    % change probability
    change_prob = [rec_cg.change_prob];
    % change identified VZA group
    change_vza = mod(change_prob, 10);
    % break time
    t_break = [rec_cg.t_break];
    % change vector magnitude
    mag = 0.1*[rec_cg.magnitude];
    % start and end time
    t_start = [rec_cg.t_start];
    t_end = [rec_cg.t_end];
    QA = [rec_cg.QA];

    % Initilize the gradual change record
    Gradual_pValue_row = 99*ones(1, num_col, num_yr); 
    GradualSlope_all_row = zeros(1, num_col, num_yr); 
    GradualSlope_sgn_row = zeros(1, num_col, num_yr); 
    GradualMag_net_row = zeros(1, num_col, num_yr, num_vza, num_DOY);
    GradualMag_pos_row = zeros(1, num_col, num_yr, num_vza, num_DOY);
    GradualMag_neg_row = zeros(1, num_col, num_yr, num_vza, num_DOY);
    for i = 1:length(rec_cg)
        icol_ids = pos(i);
        % Get the clear observations
        [clrx, clrx_sdate, clry, clr_vza, num_yrs, ann_doy] = ...
            Find_ClrObs(line_tbrdf_mrg, line_tvzagrp_mrg, irow_ids, icol_ids, sdate);

        % For abrupt change if applied
        if change_prob(i) >= 1
            % Get the valid range for during change observations
            [clrx_mod, ~, clry_mod, clr_vza_mod] = ...
                Find_RngObs(clrx, clrx_sdate, clry, clr_vza, t_end(i)+1, t_start(i+1)-1);
            chgmag_drn = nan(num_vza, 1);
            
            for n_vza = 1:num_vza
                index_vza = clr_vza_mod == n_vza;
                if sum(index_vza) > 0
                    chgmag_drn(n_vza) = median(clry_mod(index_vza)-...
                        autoTSPred_OLS_dynT(clrx_mod(index_vza), rec_cg(i).coefs(:, n_vza), num_yrs)); 
                end
            end
            rec_cg(i).chgmag_aft = rec_cg(i+1).overall_start-rec_cg(i).overall_end;
            rec_cg(i).chg_rate = rec_cg(i).chgmag_aft./rec_cg(i).overall_end;
            rec_cg(i).t_drn = t_start(i+1)-t_end(i)+1; % change duration
            rec_cg(i).chg_vzaint = change_vza(i); % change identified VZA interval
            rec_cg(i).break_yr = year(datetime(t_break(i), 'ConvertFrom', 'datenum')); % change year
            rec_cg(i).chgmag_drn = chgmag_drn*0.1; 
            % Remove the consistant dark pixels, when chenge, before
            % change overall, and after change overall magnitudes all 
            % less than threshold (1=pass, not dark pixel; 0=fail, dark pixel)
            rec_cg(i).check_notdark = rec_cg(i).overall_end(change_vza(i)) >= NTL_lmt ...
                || abs(mag(i)) >= NTL_lmt || rec_cg(i+1).overall_start(change_vza(i)) >= NTL_lmt;
        end

        if QA(i) ~= 1
            continue
        end
        % For gradual change
        [Gradual_pValue_row, GradualSlope_all_row, GradualSlope_sgn_row, ...
            GradualMag_net_row, GradualMag_pos_row, GradualMag_neg_row] = ...
            Cal_GrdChg(yr_s, t_start(i), t_end(i), clrx_sdate, clrx, ...
            rec_cg(i).coefs(2, :), rec_cg(i).pValues(2, :), icol_ids, sgn, ...
            Grad_lmt, Gradual_pValue_row, GradualSlope_all_row, ...
            GradualSlope_sgn_row, GradualMag_net_row, ...
            GradualMag_pos_row, GradualMag_neg_row, ann_doy, num_vza);
    end
    % Save the record
    save(fullfile(dir_BRDF_new, n_rst_rec, [filename_save, '.temp']), 'rec_cg');
    movefile(fullfile(dir_BRDF_new, n_rst_rec, [filename_save, '.temp']), ...
        fullfile(dir_BRDF_new, n_rst_rec, filename_save));

    %% For the abrupt changes summary
    % Initialize the abrupt change record
    map_change = struct('pos', [], 'break_point', [], 't_drn', [], ...
        'chg_vzaint', [], 'chgmag_drn', [], 'chgmag_aft', [], 'chg_rate', ...
        [], 'overall_bfr', [], 'overall_aft', []);
    % initialize change map indices for the detected change
    col_list = num2cell(col1:col2);
    [map_change(col1:col2).pos] = col_list{:}; % column postition
    % Collect and merge mapped breaks of the stratified models
    % Get the IDs for the confirmed changes
    check_notdark = [rec_cg.check_notdark];
    index_chg = find(check_notdark == 1);
    for i = index_chg
        icol_ids = pos(i);
        map_change(icol_ids).break_point = [map_change(icol_ids).break_point, t_break(i)];
        map_change(icol_ids).t_drn = [map_change(icol_ids).t_drn, rec_cg(i).t_drn];
        map_change(icol_ids).chg_vzaint = [map_change(icol_ids).chg_vzaint, rec_cg(i).chg_vzaint];
        map_change(icol_ids).chgmag_drn = [map_change(icol_ids).chgmag_drn; rec_cg(i).chgmag_drn];
        map_change(icol_ids).chgmag_aft = [map_change(icol_ids).chgmag_aft; rec_cg(i).chgmag_aft];
        map_change(icol_ids).chg_rate = [map_change(icol_ids).chg_rate; rec_cg(i).chg_rate];
        map_change(icol_ids).overall_bfr = [map_change(icol_ids).overall_bfr; rec_cg(i).overall_end];
        map_change(icol_ids).overall_aft = [map_change(icol_ids).overall_aft; rec_cg(i+1).overall_start];
    end

    %% Save the merged VZA stratified change records
    save(fullfile(dir_BRDF_new, n_rst_map, [filename_save, '.temp']), ...
        'Gradual_pValue_row', 'GradualSlope_all_row', 'GradualSlope_sgn_row', ...
        'GradualMag_net_row', 'GradualMag_pos_row', 'GradualMag_neg_row', ...
        'map_change');
    movefile(fullfile(dir_BRDF_new, n_rst_map, [filename_save, '.temp']), ...
        fullfile(dir_BRDF_new, n_rst_map, filename_save));
    fprintf(sprintf('%s done\n', filename_save))
end
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
v_dif = norm(v_dif); 
isstable = v_dif <= T_cg; % true is stable model, false is not
end

function [Gradual_pValue_row, GradualSlope_all_row, GradualSlope_sgn_row, ...
    GradualMag_net_row, GradualMag_pos_row, GradualMag_neg_row] = ...
    Cal_GrdChg(yr_s, t_start, t_end, clrx_sdate, clrx, coefs_slp, pValues, ...
    icol_ids, sgn, Grad_lmt, Gradual_pValue_row, GradualSlope_all_row, ...
    GradualSlope_sgn_row, GradualMag_net_row, GradualMag_pos_row, ...
    GradualMag_neg_row, ann_doy, num_vza)
%% For the grdual changes
tmp_yr_s = year(datetime(t_start, 'ConvertFrom', 'datenum'));
tmp_yr_e = year(datetime(t_end, 'ConvertFrom', 'datenum'));
for yr = tmp_yr_s:tmp_yr_e
    n_yr = yr-yr_s+1;
    tmp_start = max(t_start, datenum(yr, 1, 1));
    tmp_end = min(t_end, datenum(yr, 12, 31));
    index_ann = find(clrx_sdate >= tmp_start & clrx_sdate <= tmp_end);
    if isempty(index_ann)
        continue
    end
    grad_chgmag = (0.1*coefs_slp*(clrx(index_ann(end))-clrx(index_ann(1))+1))';

    % Record the lowest pValues per year
    if pValues(end) < Gradual_pValue_row(1, icol_ids, n_yr)
        Gradual_pValue_row(1, icol_ids, n_yr) = pValues(end);
    end

    % Record the largest slope per year - based on absolute value
    if abs(coefs_slp(end)) > abs(GradualSlope_all_row(1, icol_ids, n_yr))
        GradualSlope_all_row(1, icol_ids, n_yr) = coefs_slp(end);
    end

    % Skip if slope is not significant (Null hypothesis: coefficient is 0)
    if pValues(end) >= sgn
        continue
    end

    % Set to zero if too dark
    if max(abs(grad_chgmag)) < Grad_lmt
        continue
    end

    vza_pos = find(coefs_slp > 0);
    vza_neg = find(coefs_slp < 0);
    index_DOY = day(datetime(clrx_sdate(index_ann), 'ConvertFrom', 'datenum'), 'dayofyear');
    for i_DOY = index_DOY(1):index_DOY(end)
        % Skip if constant lack of data for this DOY (with corresponding
        % model period shrink)
        if ~ismember(i_DOY, ann_doy)
            continue
        end
        % Net magnitude
        for n_vza = 1:num_vza
            GradualMag_net_row(1, icol_ids, n_yr, n_vza, i_DOY) = ...
                GradualMag_net_row(1, icol_ids, n_yr, n_vza, i_DOY)+0.1*coefs_slp(n_vza);
        end

        % Positive magnitude
        if ~isempty(vza_pos)
            for n_vza = vza_pos
                GradualMag_pos_row(1, icol_ids, n_yr, n_vza, i_DOY) = ...
                    GradualMag_pos_row(1, icol_ids, n_yr, n_vza, i_DOY)+0.1*coefs_slp(n_vza);
            end
        end

        % Negative magnitude
        if ~isempty(vza_neg)
            for n_vza = vza_neg
                GradualMag_neg_row(1, icol_ids, n_yr, n_vza, i_DOY) = ...
                    GradualMag_neg_row(1, icol_ids, n_yr, n_vza, i_DOY)+(0.1*coefs_slp(n_vza));
            end
        end
    end

    %% Record the gradual change slopes/magnitudes
    % Recotd the largest slope per year - based on absolute value
    if abs(coefs_slp(end)) > abs(GradualSlope_sgn_row(1, icol_ids, n_yr))
        GradualSlope_sgn_row(1, icol_ids, n_yr) = coefs_slp(end);
    end
    
end
end


function [rec_cg_coefs, rec_cg_pValues] = RefitModel(rec_cg, num_vza, clrx, clry, clr_vza, n_times, num_c, num_yrs)
for n_vza = 1:num_vza
    if n_vza == num_vza
        index_vza = ones(size(clr_vza), 'logical');
    else
        index_vza = clr_vza == n_vza;
    end

    if sum(index_vza) > n_times*num_c 
        [rec_cg.coefs(:, n_vza), ~, ~, ~, rec_cg.pValues(:, n_vza)] = ...
            autoRobustFit2_lmtyr_dynT(clrx(index_vza), clry(index_vza), num_c, num_yrs);
    end
end
rec_cg_coefs = rec_cg.coefs; 
rec_cg_pValues = rec_cg.pValues;
end


function [clrx, clrx_sdate, clry, clr_vza, num_yrs, ann_doy] = ...
    Find_ClrObs(line_tbrdf_mrg, line_tvzagrp_mrg, irow_ids, icol_ids, sdate)

clry = line_tbrdf_mrg(irow_ids, icol_ids, :);
clry = double(clry(:));
clrvza = line_tvzagrp_mrg(irow_ids, icol_ids, :);
clr_vza = double(clrvza(:));
% get the valid data
id_range = clry < 65535;
clry = clry(id_range);
clr_vza = clr_vza(id_range);
clrx_sdate = sdate(id_range);
% Find out the annual period
clrx_doy = day(datetime(clrx_sdate, 'ConvertFrom', 'datenum'), 'dayofyear');
ann_doy = unique(clrx_doy);
num_yrs = min(length(ann_doy), 365.25); % number of available obs. per year
doy_lack = clrx_sdate(1):clrx_sdate(end);
index_date = ismember(doy_lack, clrx_sdate);
doy_lack = day(datetime(doy_lack, 'ConvertFrom', 'datenum'), 'dayofyear');
index_lack = ~ismember(doy_lack, ann_doy);
doy_lack = double(index_lack)*-1;
doy_lack = cumsum(doy_lack);
doy_lack = doy_lack(index_date);
clrx = clrx_sdate+doy_lack';
end  

function [clrx_mod, clrx_sdate_mod, clry_mod, clr_vza_mod] = ...
    Find_RngObs(clrx, clrx_sdate, clry, clr_vza, t_start, t_end)
irange_date = find(clrx_sdate >= t_start & clrx_sdate <= t_end);
clrx_mod = clrx(irange_date);
clrx_sdate_mod = clrx_sdate(irange_date);
clry_mod = clry(irange_date);
clr_vza_mod = clr_vza(irange_date);
end
