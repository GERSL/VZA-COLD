function R_new = generate_globaltif_Ref(varargin)

%% get parameters from inputs
p = inputParser;
% default values.
addParameter(p, 'H_list', 0:35); 
addParameter(p, 'V_list', 0:17); 
addParameter(p, 'num_row', 2400); 
addParameter(p, 'num_col', 2400); 
addParameter(p, 'save_code', 1); 

% request user's input
parse(p, varargin{:});
H_list = p.Results.H_list;
V_list = p.Results.V_list;
num_row = p.Results.num_row;
num_col = p.Results.num_col;
save_code = p.Results.save_code;

scn_Nlat_1 = 90-10*V_list(end)-10;
scn_Nlat_2 = 90-10*V_list(1);
scn_Wlon_1 = 10*H_list(1)-180;
scn_Wlon_2 = 10*H_list(end)-180+10;

path_data = '/shared/zhulab/Tian/Analysis/BlackMarbleTiles/Tiff_ref/h15v02.tif';
[~, R] = geotiffread(path_data);
R_new = R;
R_new.LatitudeLimits = [scn_Nlat_1 scn_Nlat_2];
R_new.LongitudeLimits = [scn_Wlon_1 scn_Wlon_2];
R_new.RasterSize = [num_row*length(V_list) num_col*length(H_list)];

% Save the global geotiff file if needed
if save_code == 1
    Img = zeros(R_new.RasterSize);
    fname_save = 'Global.tif';
    geotiffwrite(fname_save, Img, R_new)
end
end