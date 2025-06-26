function dailychangemap(yr_s, yr_e, row1, row2, col1, col2, ...
    dir_rowchg, dir_changemap, dir_tiffref, tile_name)

% Set the label names of the legend for different stratification methods
% Read the mannual interpreted table
% Add path for the coding folder
filepath_m = mfilename('fullpath');
[cmd_l, ~] = fileparts(filepath_m);

n_rst_rowchg = sprintf('ChangeMap_%s', tile_name);
n_rst_annmap = sprintf('AnnualChangeMap_%s', tile_name);
n_rst_annVZAmap = sprintf('AnnualChangeVZAintMap_%s', tile_name);
n_rst_annmagmap = sprintf('AnnualChangeMagnitudeMap_%s', tile_name);
n_rst_annChgRateDismap = sprintf('AnnualChangePctDisMap_%s', tile_name);
n_rst_annChgRateStbmap = sprintf('AnnualChangePctStbMap_%s', tile_name);
n_rst_accmap = sprintf('AccumulateChangeMap_%s', tile_name);
dir_ColorMap_Accu = fullfile(cmd_l, sprintf('ColorMap_Accu_%d%d.clr', yr_s, yr_e));
dir_ColorMap_VZAint = fullfile(cmd_l, 'ColorMap_VZAint.clr');

if ~isfolder(fullfile(dir_changemap, n_rst_annmap))
    mkdir(fullfile(dir_changemap, n_rst_annmap))
end

if ~isfolder(fullfile(dir_changemap, n_rst_annVZAmap))
    mkdir(fullfile(dir_changemap, n_rst_annVZAmap))
end

if ~isfolder(fullfile(dir_changemap, n_rst_annmagmap))
    mkdir(fullfile(dir_changemap, n_rst_annmagmap))
end

if ~isfolder(fullfile(dir_changemap, n_rst_annChgRateDismap))
    mkdir(fullfile(dir_changemap, n_rst_annChgRateDismap))
end

if ~isfolder(fullfile(dir_changemap, n_rst_annChgRateStbmap))
    mkdir(fullfile(dir_changemap, n_rst_annChgRateStbmap))
end

if ~isfolder(fullfile(dir_changemap, n_rst_accmap))
    mkdir(fullfile(dir_changemap, n_rst_accmap))
end

nrows = row2-row1+1;
ncols = col2-col1+1;
nyrs = yr_e-yr_s+1;
iyrs = yr_s:yr_e;
AnnMap = zeros(nrows, ncols, nyrs, 'int16');
AnnMap_vzaint = zeros(nrows, ncols, nyrs, 'int16');
AnnMap_mag = zeros(nrows, ncols, nyrs, 2, 'double');
AnnMap_chgscore = zeros(nrows, ncols, nyrs, 1, 'double');
AnnMap_chgduration = zeros(nrows, ncols, nyrs, 1, 'double');
AnnMap_slope_bfr = zeros(nrows, ncols, nyrs, 1, 'double');
AnnMap_slope_aft = zeros(nrows, ncols, nyrs, 1, 'double');
for irow_ids = row1:row2
    filename_save_row = sprintf('%s_changemap%04d.mat', tile_name, irow_ids);
    filename_save_row_null = sprintf('%s_changemap%04d_null.mat', tile_name, irow_ids);
    if ~isfile(fullfile(dir_rowchg, n_rst_rowchg, filename_save_row)) && ...
            ~isfile(fullfile(dir_rowchg, n_rst_rowchg, filename_save_row_null))
        continue
    end

    if ~isfile(fullfile(dir_rowchg, n_rst_rowchg, filename_save_row)) && ...
            isfile(fullfile(dir_rowchg, n_rst_rowchg, filename_save_row_null))
        continue
    end

    load(fullfile(dir_rowchg, n_rst_rowchg, filename_save_row), ...
        'DailyMap_row', 'chgvzaint_row', 'chgmag_row', '',...
        'chgscore_row', 'chgduration_row', 'chgslope_bfr_row', ...
        'chgslope_aft_row');
    tmp_map = DailyMap_row;
    AnnMap(irow_ids, :, :) = tmp_map(:, :, 1:nyrs);

    tmp_map = chgvzaint_row;
    AnnMap_vzaint(irow_ids, :, :) = tmp_map(:, :, 1:nyrs);

    tmp_map = chgmag_row;
    tmp_map(tmp_map < 0) = 0;
    AnnMap_mag(irow_ids, :, :, 1) = tmp_map(:, :, 1:nyrs);
    tmp_map = chgmag_row;
    tmp_map(tmp_map > 0) = 0;
    AnnMap_mag(irow_ids, :, :, 2) = tmp_map(:, :, 1:nyrs);

    tmp_map = chgmag_row;
    tmp_map(tmp_map < 0) = 0;
    AnnMap_mag(irow_ids, :, :, 1) = tmp_map(:, :, 1:nyrs);
    tmp_map = chgmag_row;
    tmp_map(tmp_map > 0) = 0;
    AnnMap_mag(irow_ids, :, :, 2) = tmp_map(:, :, 1:nyrs);

    tmp_map = chgscore_row;
    AnnMap_chgscore(irow_ids, :, :) = tmp_map(:, :, 1:nyrs);

    tmp_map = chgduration_row;
    AnnMap_chgduration(irow_ids, :, :) = tmp_map(:, :, 1:nyrs);

    tmp_map = chgslope_bfr_row;
    AnnMap_slope_bfr(irow_ids, :, :) = tmp_map(:, :, 1:nyrs);
    tmp_map = chgslope_aft_row;
    AnnMap_slope_aft(irow_ids, :, :) = tmp_map(:, :, 1:nyrs);
end
AnnMap(isnan(AnnMap)) = 0;
AnnMap_vzaint(isnan(AnnMap_vzaint)) = 0;
AnnMap_mag(isnan(AnnMap_mag)) = 0;
AccuMap = zeros(nrows, ncols, 'int16');
for n_yr = 1:length(iyrs)
    tmp_map = squeeze(AnnMap(:, :, n_yr));
    AccuMap(tmp_map > 0) = tmp_map(tmp_map > 0);
end
AccuMap(isnan(AccuMap)) = 0;
[~, R] = readgeoraster(fullfile(dir_tiffref, sprintf('%s.tif', tile_name)));

for show_year = yr_s:yr_e
    show_nyear = show_year-yr_s+1;
    AccuMap(AnnMap(:, :, show_nyear) > 0) = show_year;

    % Save the annual change map
    filename_save_map = sprintf('%s_DailyChangeMap_%d', tile_name, show_year);
    filename_save_tif = sprintf('%sAnnMap.tif', filename_save_map);
    geotiffwrite(fullfile(dir_changemap, n_rst_annmap, filename_save_tif), AnnMap(:, :, show_nyear), R);

    % Save the annual change VZA interval map
    filename_save_map = sprintf('%s_DailyChangeVZAintMap_%d', tile_name, show_year);
    filename_save_tif = sprintf('%sAnnMap.tif', filename_save_map);
    filename_save_color = sprintf('%sAnnMap.clr', filename_save_map);
    geotiffwrite(fullfile(dir_changemap, n_rst_annVZAmap, filename_save_tif), AnnMap_vzaint(:, :, show_nyear), R);
    copyfile(dir_ColorMap_VZAint, fullfile(dir_changemap, n_rst_annVZAmap, filename_save_color))

    % Save the annual change magnitudde map
    filename_save_map = sprintf('%s_DailyChangeMagnitudeMap_%d', tile_name, show_year);
    filename_save_tif = sprintf('%sAnnMap_pos.tif', filename_save_map);
    geotiffwrite(fullfile(dir_changemap, n_rst_annmagmap, filename_save_tif), AnnMap_mag(:, :, show_nyear, 1), R);

    % Save the annual change magnitudde map
    filename_save_map = sprintf('%s_DailyChangeMagnitudeMap_%d', tile_name, show_year);
    filename_save_tif = sprintf('%sAnnMap_neg.tif', filename_save_map);
    geotiffwrite(fullfile(dir_changemap, n_rst_annmagmap, filename_save_tif), AnnMap_mag(:, :, show_nyear, 2), R);
end

% Save the accumulative annual change map
filename_save_map = sprintf('%s_AccuChangeMap_%d', tile_name, show_year);
filename_save_tif = sprintf('%sAccuMap.tif', filename_save_map);
filename_save_color = sprintf('%sAccuMap.clr', filename_save_map);
geotiffwrite(fullfile(dir_changemap, n_rst_accmap, filename_save_tif), AccuMap, R);
copyfile(dir_ColorMap_Accu, fullfile(dir_changemap, n_rst_accmap, filename_save_color))
end