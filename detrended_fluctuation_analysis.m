function [Fval, win_lens] = detrended_fluctuation_analysis(signal, time, fs, cfg)
%DETRENDED_FLUCTUATION_ANALYSIS This function performs detrended
%fluctuation analysis for each signal in the "signal".
%
% Detrended fluctuation analysis (DFA) is a technique to estimate
% long-range temporal correlation in time series. DFA includes fitting
% a polynomial to each time slice of a signal and therefore accounts for
% non-stationarity in the signal ("detrending"). This function performs DFA
% for each signal in the input "signal" and returns the fluctuation
% function "Fval" of time window length. The steps are described below.
% Given N x T signals (N = number of signals, T = number of time points),
%   1. Calculate the profile of each signal, i.e., the cumulative sum over
%      time points.
%   2. For each window length W, segment the profile time series into
%      floor(T/W) slices.
%      - For each window, fit a polynomial of order n to the profile and
%        subtract the profile time series from the fitted polynomial,
%        resulting in a detrended profile for each time window.
%      - Calculate the standard deviations of detrended profile in each
%        window and average the standard deviations across windows.
%      - The averaged standard deviation is the value of fluctuation
%        function at this window length.
%   3. Loop over window lengths, resulting in the values of fluctuation
%      function at each window length ("Fval" and "win_len").
%   Note: The window lengths are chosen on a logarithmic scale to
%   facilitate downstream analysis of fitting power law to the fluctuation
%   function.
%
%   [Fval, win_lens] = detrended_fluctuation_analysis(signal, time, fs, cfg)
%
% Input:
%  signal:  N x T, N = number of signals, T = number of time points
%  time:    1 x T, time points associated with the signals, in second
%  fs:      sampling frequency of the signals, in Hz
%  cfg:     struct of configurations
%   cfg.max_win_len:    maximal window length, in second; default: signal
%                       length divided by 8 (corresponding to at least 8
%                       segments for averaging the standard deviations)
%   cfg.min_win_len:    minimal window length, in second; default:
%                       10*sampling period (corresponding to at least 10
%                       time points in each window)
%   cfg.num_win_lens:   number of window length values; default: 12
%   cfg.order_polyfit:  order of polynomial fitted to profile time series;
%                       default: 1
% Output:
%  Fval:    values of fluctuation function
%  win_len: window lengths associated with each of Fval
%
%                                       WU Kejian, YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
[N, T] = size(signal);
if ~isfield(cfg, 'max_win_len'), max_win_len = floor(T/8);  else, max_win_len = floor(cfg.max_win_len*fs);  end
if ~isfield(cfg, 'min_win_len'), min_win_len = 10;  else, min_win_len = floor(cfg.min_win_len*fs);  end
if ~isfield(cfg, 'num_win_lens'), num_win_lens = 12;  else, num_win_lens = cfg.num_win_lens;  end
if ~isfield(cfg, 'order_polyfit'), order_polyfit = 1;  else, order_polyfit = cfg.order_polyfit;  end

% Form window lengths at which fluctuation function is evaluated
win_lens = logspace(log10(min_win_len), log10(max_win_len), num_win_lens);
win_lens = floor(win_lens);
M = length(win_lens);
fprintf("Window lengths of interest: %.4f sec (%d time points) to %.4f sec (%d time points).\n", ...
    min_win_len*fs, min_win_len, max_win_len*fs, max_win_len)

% Pre-allocate arrays to store results
Fval = zeros(N,M);

% Loop over signals
for n = 1:N
    % get signal and compute profile
    x = signal(n,:);
    x = x - mean(x);
    y = cumsum(x);

    % pre-allocate fluctuation function for this signal
    Fval_n = zeros(1,M);

    % Loop over window lengths
    for m = 1:M
        % get window length (in time points) and compute the number of segments
        win_len = win_lens(m);
        n_segm = floor(T/win_len);

        % initialize standard deviations
        std_running_sum = 0;

        % Loop over segments (or time windows)
        for w = 1:n_segm
            % prepare time segment
            idx = (1:win_len) + (w - 1)*win_len;    % time indices for the window
            y_win = y(idx);    % profile time series in the window
            t_win  = (1:win_len)';    % time points in the window

            % fit polynomial to the profile time series
            p = polyfit(t_win, y_win, order_polyfit);
            trend = polyval(p, t_win)';

            % standard deviation of the detrended profile
            std_running_sum = std_running_sum + std(y_win - trend);
        end
        % averaged standard deviation, which is exactly, the value of
        % fluctuation function at this window length
        Fval_n(m) = std_running_sum / n_segm;
    end
    Fval(n,:) = Fval_n;
end

end