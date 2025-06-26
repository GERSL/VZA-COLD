function DarkRec(nlines_read, irows_s, irows_e, num_row, ...
    num_col, tile_name, dir_save, dir_data_row, dir_l_chg)

DarkRec_RollMedian = zeros(num_row, num_col);
for i_irows = irows_s:irows_e
    for irow_ids = 1:nlines_read
        n_row = irow_ids+(i_irows-1)*nlines_read;
        filename_data = sprintf('%s_DarkPixel_%04d.mat', tile_name, n_row);
        if ~isfile(fullfile(dir_data_row, filename_data))
            continue
        end
        load(fullfile(dir_data_row, filename_data), 'RollMedian');
        DarkRec_RollMedian(n_row, :) = RollMedian;
    end
end 

% Only include the mapped change pixels
dir_chg = fullfile(dir_l_chg, tile_name, sprintf('Accumulate_ChangeMap_%s', tile_name));
fname_chg = sprintf('%s_LatestAbruptChangeYear_20142022.tif', tile_name);
[LocChg, ~] = readgeoraster(fullfile(dir_chg, fname_chg));
id_chg = LocChg > 0;
RollMedian = DarkRec_RollMedian(id_chg);

save(dir_save, 'RollMedian')
end