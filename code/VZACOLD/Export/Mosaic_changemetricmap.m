function Mosaic_changemetricmap(H_list, V_list, num_row, num_col, ...
    dir_l, dir_l_save, n_rst, fname, R)

% Get the directories
fname_save = sprintf('%s%s', 'Global', fname);
path_save = fullfile(dir_l_save, fname_save);
% Initialize the maps
nrows = num_row*length(V_list);
ncols = num_col*length(H_list);
datatype = 'double';
GlobalMap = zeros(nrows, ncols, datatype);

for i_H = H_list
    for i_V = V_list(3:end-3) % Only include +70 to -60 degree latitude
        tile_name = sprintf('h%02dv%02d', i_H, i_V);
        dir_data = fullfile(dir_l, tile_name, sprintf(n_rst, tile_name));
        fname_data = sprintf('%s%s', tile_name, fname);
        path_data = fullfile(dir_data, fname_data);
        if ~isfile(path_data)
            continue
        end
        row_list = i_V*num_row+1:(i_V+1)*num_row;
        col_list = i_H*num_col+1:(i_H+1)*num_col;
        [tmp, ~] = readgeoraster(path_data);

        GlobalMap(row_list, col_list) = tmp;
    end
end

% Save the global map
geotiffwrite(path_save, GlobalMap, R);
end