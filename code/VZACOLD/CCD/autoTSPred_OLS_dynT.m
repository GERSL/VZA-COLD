function outfity = autoTSPred_OLS_dynT(outfitx, fit_cft, num_yrs)
% Auto Trends and Seasonal Predict
% INPUTS:
% outfitx - Julian day [1; 2; 3];
% fit_cft - fitted coefficients;
% OUTPUTS:
% outfity - predicted reflectances [0.1; 0.2; 0.3];
% General model TSModel:
% f(x) =  a0 + b0*x + a1*cos(x*w) + b1*sin(x*w) 

w=2*pi/num_yrs; % anual cycle 
outfity = [ones(size(outfitx)), outfitx, ...% overall ref + trending
        cos(w*outfitx), sin(w*outfitx), ...% add seasonality
        ]*fit_cft; % add trimodal seasonality
end