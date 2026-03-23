function [exponent, fit_results] = fit_power_law(data, xspan, cfg)
%FIT_POWER_LAW This function fits a power law of the form y(x) = b * x^(-a)
%to each row of the input "data" array using Matlab polyfit function and
%returns the exponent a.
%
%   exponent = fit_power_law(data, xspan, cfg)
%   [exponent, fit_results] = fit_power_law(data, xspan, cfg)
%
% Input:
%  data:    N x F, N = number of samples, F = number of features.
%           For example, power spectrums of N regions stacked in an N-by-F
%           matrix, where F is the number of frequency points.
%  xspan:   1 x F, coordinates associated with each column in "data"
%           For example, frequency points associated with spectrums.
%  cfg:     struct of configurations
%   cfg.xrange: interval over which the data are assumed to demonstrate a
%               power law relation; required
% Output:
%  exponent:    exponent in the power law fitted to each row, N x 1
%  fit_results: other fitting results, a struct including
%   fit_results.xspan2fit:  coordinates selected for fitting, 1 x F'
%   fit_results.data_pred:  predicted data values, N x F'
%   fit_results.rsquared:   R^2 in linear regression, N x 1

if ~isfield(cfg, 'xrange'), error('Please specify an interval!'); else, xrange = cfg.xrange; end
fprintf("Interval to fit the power law: [%.4f, %.4f].\n", xrange(1), xrange(2));

% Get data to fit and the associated coordinates
idx = (xspan > xrange(1) & xspan < xrange(2));
data2fit = data(:, idx);
xspan2fit = xspan(idx);

% Pre-allocate arrays to store results
[N, F] = size(data2fit);
exponent = zeros(N,1);
fit_results = [];
fit_results.xspan2fit = xspan2fit;
fit_results.data_pred = zeros(N,F);
fit_results.rsquared = zeros(N,1);

% Loop over samples
for i = 1:N

    % Fit log y = log b + (-a) * log x, equivalent to y = b * x^(-a)
    [p, S] = polyfit(log(xspan2fit), log(data2fit(i,:)), 1);
    
    % Obtain exponent a
    % (Note: p(1) = -a, p(2) = log b)
    exponent(i) = -p(1);
    
    % Calculate predicted data values
    data_pred = polyval(p, log(xspan2fit), S);
    data_pred = exp(data_pred);
    
    % Export other results
    fit_results.data_pred(i,:) = data_pred;
    fit_results.rsquared(i) = S.rsquared;
end

end