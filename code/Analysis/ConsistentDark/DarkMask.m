function DarkMask(nlines_read, irows_s, irows_e, dir_l_tif, num_row, ...
    num_col, tile_name, dir_l_save, dir_data_row, Dark_lmt)

DarkMask_RollMedian = zeros(num_row, num_col, 'logical');
for i_irows = irows_s:irows_e
    for irow_ids = 1:nlines_read
        n_row = irow_ids+(i_irows-1)*nlines_read;
        filename_data = sprintf('%s_DarkPixel_%04d.mat', tile_name, n_row);
        if ~isfile(fullfile(dir_data_row, filename_data))
            continue
        end
        load(fullfile(dir_data_row, filename_data), 'RollMedian');
        DarkMask_RollMedian(n_row, :) = RollMedian >= Dark_lmt*10;
    end
end

[~, R] = readgeoraster(fullfile(dir_l_tif, sprintf('%s.tif', tile_name)));
fname_save = sprintf('DarkMask_RollMedian_%s.tif', tile_name);
geotiffwrite(fullfile(dir_l_save, fname_save), DarkMask_RollMedian, R)
end