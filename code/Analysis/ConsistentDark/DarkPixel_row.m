function DarkPixel_row(num_c, n_times, date_start, date_end, dir_DarkPixel, ...
    tile_name, n_irows, nlines_read, col1, col2, conse)

%% Read the n*row time series data
% The year and doy for the start and end time
yr_s = floor(date_start/1000);
doy_s = mod(date_start, 1000);
yr_e = floor(date_end/1000);
doy_e = mod(date_end, 1000);
date_start = datenum(yr_s, 1, doy_s);
date_end = datenum(yr_e, 1, doy_e);
sdate = (date_start: date_end)';
[sdate, line_tbrdf_mrg, ~] = createMergeRowdata(tile_name, n_irows, sdate);
irow_list = 1:nlines_read;

% Loop running the line(s)
for irow_ids = irow_list
    n_row = irow_ids+(n_irows-1)*nlines_read;
    % check if row already exist
    filename_save = sprintf('%s_DarkPixel_%04d.mat', tile_name, n_row);
    if isfile(fullfile(dir_DarkPixel, filename_save))
        continue
    end

    % Initialize the record
    RollMedian = zeros(col2-col1+1, 1);
    for icol_ids = col1:col2
        % Get the data
        y_all = line_tbrdf_mrg(irow_ids, icol_ids, :);
        y_all = y_all(:);
        % Obtaining the valid obs. ids for the modeling
        idgood = y_all < 65535; % clear pixel IDs
        %% Case 1: All time series data do not have enough clear observations for change detection
        if sum(idgood) < n_times*num_c % not enough for all data
            continue
        end

        % Xs & Ys for computation after removing noises such as cloud
        clrx = sdate(idgood);
        clry = double(y_all(idgood));
        % Skip if not enough clear pixels, all data time series
        if length(clrx) < n_times*num_c
            continue
        end

        % Fix rolling conse.
        for i_y = 1:length(clry)-conse+1
            RollMedian(icol_ids) = max(RollMedian(icol_ids), median(clry(i_y:i_y+conse-1)));
        end
    end
    save(fullfile(dir_DarkPixel, filename_save), 'RollMedian');
end
end
