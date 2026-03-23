function auc = area_under_curve(data, xspan, cfg)
%AREA_UNDER_CURVE This function calculates area under each curve in the
%"data" array, over a given interval, using Matlab trapz function.
%
%   auc = area_under_curve(data, xspan, cfg)
%
% Input:
%  data:    N x F, N = number of samples, F = number of features.
%           For example, autocorrelation functions of N regions stacked in
%           an N-by-L matrix, where L is the number of time lags.
%  xspan:   1 x F, coordinates associated with each column in "data"
%           For example, the time lags for autocorrelation functions.
%  cfg:     struct of configurations
%   cfg.xrange: interval to integrate the curves; required
% Output:
%  auc:     N x 1, area under curve

if ~isfield(cfg, 'xrange'), error('Please specify an interval!'); else, xrange = cfg.xrange; end

idx = (xspan > xrange(1) & xspan < xrange(2));

data2int = data(:, idx);
xspan2int = xspan(idx);

auc = trapz(xspan2int, data2int, 2);

end