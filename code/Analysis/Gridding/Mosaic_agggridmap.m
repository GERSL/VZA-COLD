function Mosaic_agggridmap(H_list, V_list, num_row_agggrid, num_col_agggrid, ...
    dir_l, var, year, path_save, R)

% Initialize the maps
nrows = num_row_agggrid*length(V_list);
ncols = num_col_agggrid*length(H_list);
GlobalMap = zeros(nrows, ncols, 'double');
for i_H = H_list
    for i_V = V_list
        tile_name = sprintf('h%02dv%02d', i_H, i_V);
        row_list = i_V*num_row_agggrid+1:(i_V+1)*num_row_agggrid;
        col_list = i_H*num_col_agggrid+1:(i_H+1)*num_col_agggrid;

        fname_data = sprintf('AggGrid_%s_%s_%d.tif', var, tile_name, year);
        [tmp, ~] = readgeoraster(fullfile(dir_l, tile_name, fname_data));
        [~, ~, num_vza] = size(tmp);
        if num_vza > 1
            tmp = tmp(:, :, end);
        end
        GlobalMap(row_list, col_list) = tmp;
    end
end

% Save the global map
geotiffwrite(path_save, GlobalMap, R);
end