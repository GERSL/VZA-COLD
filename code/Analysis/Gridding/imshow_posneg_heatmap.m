function imshow_posneg_heatmap(var, dir_l_data, year, path_savetif, ...
    Img_R)

% Get the directories
[ChgRec_pos, ~] = readgeoraster(fullfile(dir_l_data,  sprintf("GlobalAggGrid_%s_pos_%d.tif", var, year)));
[ChgRec_neg, ~] = readgeoraster(fullfile(dir_l_data, sprintf("GlobalAggGrid_%s_neg_%d.tif", var, year)));
if contains(var, 'int')
    ChgRec_pos(isnan(ChgRec_pos)) = 0;
    ChgRec_neg(isnan(ChgRec_neg)) = 0;
end
RGB_ChgRec = zeros(size(ChgRec_pos, 1), size(ChgRec_pos, 2), 3, 'uint8');

% Get the color map setting
val1 = 1;
val2 = 0.1;
R=[val2 1
   0 val1];
G=[val2 1
   0 val2];
B=[val1 1
   0 val2];
R = interp2(R,8);
G = interp2(G,8);
B = interp2(B,8);
color_map = uint8(255*cat(3,R,G,B));

if contains(var, 'int')
    lmt_up = 1;
elseif contains(var, 'area')
    lmt_up = 5000; 
elseif contains(var, 'rad')
    lmt_up = 2000; 
end
lmt_low = 0;

% Assign the range limit to pos/neg maps
loc_pos_color = round(size(color_map, 2)*(ChgRec_pos-lmt_low)/(lmt_up-lmt_low));
loc_pos_color(loc_pos_color < 1) = 1;
loc_pos_color(loc_pos_color > size(color_map, 2)) = size(color_map, 2);
if contains(var, 'area')
    loc_neg_color = size(color_map, 1)-round(size(color_map, 1)*...
        (ChgRec_neg-lmt_low)/(lmt_up-lmt_low))+1;
else
    loc_neg_color = size(color_map, 1)-round(size(color_map, 1)*...
        (lmt_low-ChgRec_neg)/(lmt_low+lmt_up))+1;
end
loc_neg_color(loc_neg_color < 1) = 1;
loc_neg_color(loc_neg_color > size(color_map, 1)) = size(color_map, 1);

% Assign the heat map color
for row = 1:size(RGB_ChgRec, 1)
    for col = 1:size(RGB_ChgRec, 2)
        RGB_ChgRec(row, col, :) = color_map(loc_neg_color(row, col), loc_pos_color(row, col), :);
    end
end

% Save geotiff image
geotiffwrite(path_savetif, RGB_ChgRec, Img_R);
end