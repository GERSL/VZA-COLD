function [fit_cft, rmse, v_dif, weight_robustfit, p_robustfit] = ...
    autoRobustFit_lmtyr_dynT(x, y, df, num_yrs)
% INPUTS:
% x - Julian day [1; 2; 3];
% y - predicted reflectances [0.1; 0.2; 0.3];
% w - cofficient related with T=365d;
% yr - number of years ceil([end(day)-start(day)]/365);
% outfitx - fitted x values;

% OUTPUTS:
% outfity - fitted y values;
% fit_cft - fitted coefficients;
% General fitting model: f(x) =  a0 + a1*cos(x*w/N) + b1*sin(x*w/N) + ...
% a2*cos(x*w/1) + b2*sin(x*w/1)

n = length(x); % number of clear pixels
w = 2*pi/num_yrs; % anual cycle

% % build X
X = zeros(n, df);
X(:, 1) = ones(n, 1);
X(:, 2) = x;
X(:, 3) = cos(w*x);
X(:, 4) = sin(w*x);

% Robust fitting
[fit_cft, stats] = robustfit(X, y, [], [], 'off');
weight_robustfit = stats.w;
p_robustfit = stats.p;
ids_stats_w = stats.w;
ids_stats_w = ids_stats_w > 0;
yhat = autoTSPred_OLS_dynT(x, fit_cft, num_yrs);
v_dif = y-yhat;
rmse = RMSE(y(ids_stats_w), yhat(ids_stats_w), df);
end

function r = RMSE(data, estimate, df)
    % Function to calculate root mean square error from a data vector or matrix 
    % and the corresponding estimates.
    % Usage: r=rmse(data,estimate)
    % Note: data and estimates have to be of same size
    % Example: r=rmse(randn(100,100),randn(100,100));
    % delete records with NaNs in both datasets first
    I = ~isnan(data) & ~isnan(estimate); 
    data = data(I); 
    estimate = estimate(I);
    r = sqrt(sum((data(:)-estimate(:)).^2)/(numel(data))-df);
end

