function [relax_time, fit_results] = fit_exponetial_decay(data, xspan, cfg)
%FIT_EXPONENTIAL_DECAY This function fits an exponential function of the
%form y(t) = A * exp(-t/t0) to each row of the input "data" array using
%Matlab polyfit function and returns the relaxation time t0.
%
%   relax_time = fit_exponetial_decay(data, xspan, cfg)
%   [relax_time, fit_results] = fit_exponetial_decay(data, xspan, cfg)
%
% Input:
%  data:    N x F, N = number of samples, F = number of features.
%           For example, autocorrelation functions of N regions stacked in
%           an N-by-L matrix, where L is the number of time lags.
%  xspan:   1 x F, coordinates associated with each column in "data"
%           For example, the time lags for autocorrelation functions.
%  cfg:     struct of configurations
%   cfg.xrange: interval over which the data are assumed to demonstrate an
%               exponential decay; required
% Output:
%  relax_time:  relaxation time in the decaying exponential function fitted
%               to each row, N x 1
%  fit_results: other fitting results, a struct including
%   fit_results.xspan2fit:  coordinates selected for fitting, 1 x F'
%   fit_results.data_pred:  predicted data values, N x F'
%   fit_results.rsquared:   R^2 in linear regression

if ~isfield(cfg, 'xrange'), error('Please specify an interval!'); else, xrange = cfg.xrange; end
% Get data to fit and the associated coordinates
idx = (xspan > xrange(1) & xspan < xrange(2));
data2fit = data(:, idx);
xspan2fit = xspan(idx);

% Pre-allocate arrays to store results
[N, F] = size(data2fit);
relax_time = zeros(N,1);
fit_results = [];
fit_results.xspan2fit = xspan2fit;
fit_results.data_pred = zeros(N,F);
fit_results.rsquared = zeros(N,1);
fit_results.delta = zeros(N,1);

% Loop over samples
for i = 1:N

    % Fit log y = log A + (-1/t0) * t, equivalent to y = A * exp(-t/t0)
    [p, S] = polyfit(xspan2fit, log(data2fit(i,:)), 1);
    
    % Obtain exponent t0
    % (Note: p(1) = -1/t0, p(2) = log A)
    relax_time(i) = -1/p(1);
    
    % Calculate predicted data values
    data_pred = polyval(p, log(xspan2fit), S);
    data_pred = exp(data_pred);
    
    % Export other results
    fit_results.data_pred(i,:) = data_pred;
    fit_results.rsquared(i) = S.rsquared;
end

end