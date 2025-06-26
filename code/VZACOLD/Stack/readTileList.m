function tile_list = readTileList(varargin)
%% get parameters from inputs
p = inputParser;

% default values.
addParameter(p, 'H_list', 0:35);
addParameter(p, 'V_list', 0:17);

% request user's input
parse(p, varargin{:});

% Get the tile list
tile_list = [];
for i_H = H_list
    for i_V = V_list(3:end-3)
        tile_name = sprintf('h%02dv%02d', i_H, i_V); 
        tile_list = [tile_list, tile_name];
    end
end
end