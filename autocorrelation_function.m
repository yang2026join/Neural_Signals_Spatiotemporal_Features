function [acf, lag] = autocorrelation_function(signal, time, fs, cfg)
%AUTOCORRELATION_FUNCTION This function calculates the autocorrelation
%funtions of multiple signals.
%
% Autocorrelation function is the correlation of a signal with a delayed
% copy of itself, measuring the self-similarity of the signal over
% different time lags. This function calculates autocorrelation function
% for each signal in the input "signal".
%
%   [acf, lag] = autocorrelation_function(signal, time, fs, cfg)
%
% Input:
%  signal:  N x T, N = number of signals, T = number of time points
%  time:    1 x T, time points associated with the signals, in second
%  fs:      sampling frequency of the signals, in Hz
%  cfg:     struct of configurations, must contain the following fields
%   cfg.max_lag:    largest lag to evaluate, in second; required
%   cfg.method:     method for correlation coefficient, must be one of the
%                   "Pearson", "Spearman", "Kendall"; default: Pearson
% Output:
%  acf:     N x M autocorrelation data, M = number of lags
%  lag:     lags in second
%
%                                                    YANG Hao, 2025-26 Fall
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
if ~isfield(cfg, 'max_lag'), error('Please specify an upper bound of lag to evaluate.'); else, max_lag = cfg.max_lag; end
if ~isfield(cfg, 'method'), method = 'Pearson'; else, method = cfg.method;  end

% Prepare time lags
lag_idx = 0:floor(max_lag*fs);
lag = lag_idx / fs;

% Prepare indices of sample time points
time_idx_x = 1:(length(time)-length(lag));
time_idx_xlag = time_idx_x + lag_idx';

% Loop over signals
acf = zeros(size(signal,1),length(lag));
for i = 1:size(signal,1)
    sig_temp = signal(i,:);
    x = sig_temp(time_idx_x);    % x[t], row vector
    x_lag = sig_temp(time_idx_xlag);    % x[t + k], for k = 0 to K, stacked in a matrix
    acf(i,:) = corr(x', x_lag', "Type", method, "Rows", "pairwise");
end

end