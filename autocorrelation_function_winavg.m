function [acf, lag] = autocorrelation_function_winavg(signal, time, fs, cfg)
%AUTOCORRELATION_FUNCTION_WINAVG This function calculates the
%autocorrelation funtions of multiple signals by averaging ACF across
%multiple windows. This function depends on function
%AUTOCORRELATION_FUNCTION.
%
% Autocorrelation function is the correlation of a signal with a delayed
% copy of itself, measuring the self-similarity of the signal over
% different time lags. This function separates the given signals into
% multiple slices based on the given time window and step, calculate
% autocorrelation function for each slice of signals, and then average the
% ACF across time windows.
%
%   [acf, lag] = autocorrelation_function_winavg(signal, time, fs, cfg)
%
% Input:
%  signal:  N x T, N = number of signals, T = number of time points
%  time:    1 x T, time points associated with the signals, in second
%  fs:      sampling frequency of the signals, in Hz
%  cfg:     struct of configurations, must contain the following fields
%   cfg.max_lag:    largest lag to evaluate, in second; required
%   cfg.method:     method for correlation coefficient, must be one of the
%                   'Pearson', 'Spearman', 'Kendall'; default: Pearson
%   cfg.win:        length of time window, in second; required
%   cfg.step:       length of moving step, in second; required
% Output:
%  acf:     N x M autocorrelation data, M = number of lags
%  lag:     lags in second
%
%                                                    YANG Hao, 2025-26 Fall
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
if ~isfield(cfg, 'win'), error('Please specify the window length!'); else, win = cfg.win; end
if ~isfield(cfg, 'step'), error('Please specify the step length!'); else, step = cfg.step; end
if ~isfield(cfg, 'max_lag'), error('Please specify the maximal lag to evaluate.'); else, max_lag = cfg.max_lag; end
if ~isfield(cfg, 'method'), method = 'Pearson'; else, method = cfg.method; end

% Prepare time lags
lag_idx = 0:floor(max_lag*fs);
lag = lag_idx / fs;

% Prepare time indices of each window
n_time = length(time);
win_tpnts = round(win*fs);
step_tpnts = round(step*fs);
n_wins = floor((n_time-win_tpnts)/step_tpnts) + 1;
time_idx = (1:win_tpnts) + step_tpnts*(0:(n_wins-1))';
fprintf("Window length turns out to be %.4f s.\n", win_tpnts/fs)
fprintf("Step length turns out to be %.4f s.\n", step_tpnts/fs)
fprintf("Maximal lag to evaluate ACF turns out to be %.4f s.\n", lag(end))

% Prepare for looping over time windows
cfgIn = [];    % configurations passed to the function for a single time window
cfgIn.max_lag = max_lag;
cfgIn.method = method;
acf = zeros(size(signal,1),length(lag));    % pre-allocate array to store results
n_skips = 0;

% Loop over time windows
for i = 1:n_wins
    sig_temp = signal(:,time_idx(i,:));
    time_temp = time(time_idx(i,:));
    if any(diff(time_temp) > 1.01/fs)    % check discontinuous time points
        n_skips = n_skips + 1;    % record and do not calculate
    else
        [acf_temp, ~] = autocorrelation_function(sig_temp, time_temp, fs, cfgIn);
        acf = acf + acf_temp;
    end
end
acf = acf / (n_wins - n_skips);

end