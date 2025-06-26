function AggregateGrid(tile_name, var, year, dir_data, n_rst, fname_data, ...
    dir_PixelArea, dir_wm, fname_wm, path_save, grid_size, n_vza)
% This function aims to aggregate the NTL change result into equal
% latitude/longitude degree grids.

% Read in the pixel area (km2)
fname_PixelArea = sprintf('%s.tif', tile_name);
[PixelArea, ~] = readgeoraster(fullfile(dir_PixelArea, fname_PixelArea));
PixelArea = PixelArea./1000000;

% Read in the water mask
[WaterMask, ~] = readgeoraster(fullfile(dir_wm, sprintf(fname_wm, tile_name)));
loc_mask = ismember(WaterMask, 1:5);

% Read in the data
[ChgRec, R] = readgeoraster(fullfile(dir_data, tile_name, sprintf(n_rst, tile_name), sprintf(fname_data, tile_name, year, n_vza)));
ChgRec(~loc_mask) = 0; % exclude the ocean pixels
if contains(var, 'int')
    ChgRec = ChgRec.*(PixelArea/mean(PixelArea(:)));
elseif contains(var, 'area')
    ChgRec = double(ChgRec ~= 0).*PixelArea;
else
    ChgRec = ChgRec.*PixelArea;
end

%% this is aggregate to grid with sum
scale = 1/grid_size;
if contains(var, 'int')
    ChgRec(ChgRec == 0) = NaN;
    ChgRec(ChgRec == Inf) = NaN;
    [ChgRec_agg, R_grid] = agr2grid(ChgRec, R, scale, 'mean');
else
    [ChgRec_agg, R_grid] = agr2grid(ChgRec, R, scale, 'sum');
end

%% saving image using geosave
geotiffwrite(path_save, ChgRec_agg, R_grid);
end

function [grid, R_grid] = agr2grid(image, R, scale, method)
    % Use the georesize to modify the R
    [~, R_grid] = georesize(image, R, scale); % 0.05*2400 pixels = 120 pixels, each grid has 20 pixels by 20 pixels
    
    % Use blockproc to aggregate the image
    switch method
        case 'sum'
            fun_bock_sum = @(block) sum(block.data(:)); % this is sum, and similarily
            % we can also define other functions
        case 'mean'
            fun_bock_sum = @(block) mean(block.data(:), 'omitmissing'); 
    end
    grid = blockproc(image, [int32(1/scale), int32(1/scale)], fun_bock_sum); % set up the 20 pixels by 20 pixels as a grid
end