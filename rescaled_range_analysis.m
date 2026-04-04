function [RS, win_lens] = rescaled_range_analysis(signal, time, fs, cfg)
%RESCALED_RANGE_ANALYSIS This function performs rescaled range analysis
%on each signal in the "signal".
%
% Rescaled range analysis is a technique to estimate long-range
% correlations in a time series. The steps are described below.
% Given a signal of length T (number of time points),
%   1. Divide the signal into segments. The length of each segment is L.
%   2. For each segment i,
%   - de-mean the segment
%   - culculate cumulative sum of the de-meaned segment, namely the profile
%   - calculate the range R_i and standard deviation S_i of the profile, as
%     well as the ratio (R/S)_i.
%   3. The rescaled range R/S is the average of (R/S)_i across segments.
%   4. Loop over different segment lengths L, resulting in R/S values for
%      each window length ("RS" and "win_lens").
%   Note: The window lengths L are chosen on a logarithmic scale to
%   facilitate downstream analysis of fitting power law to the (R/S).
%
% Key reference(s):
% Dong, J., Jing, B., Ma, X., Liu, H., Mo, X., & Li, H. (2018). Hurst
% exponent analysis of resting-state fMRI signal complexity across the
% adult lifespan. Frontiers in neuroscience, 12, 34.
%
%   [RS, win_lens] = rescaled_range_analysis(signal, time, fs, cfg)
%
% Input:
%  signal:  N x T, N = number of signals, T = number of time points
%  time:    1 x T, time points associated with the signals, in second
%  fs:      sampling frequency of the signals, in Hz
%  cfg:     struct of configurations
%   cfg.max_win_len:    maximal window length, in second; default: signal
%                       length divided by 8 (corresponding to at least 8
%                       segments for averaging the R/S)
%   cfg.min_win_len:    minimal window length, in second; default:
%                       10*sampling period (corresponding to at least 10
%                       time points in each window)
%   cfg.num_win_lens:   number of window length values; default: 12
% Output:
%  RS:          rescaled range values
%  win_lens:    window lengths in number of time points
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
[N, T] = size(signal);
if ~isfield(cfg, 'max_win_len'), max_win_len = floor(T/8);  else, max_win_len = floor(cfg.max_win_len*fs);  end
if ~isfield(cfg, 'min_win_len'), min_win_len = 10;  else, min_win_len = floor(cfg.min_win_len*fs);  end
if ~isfield(cfg, 'num_win_lens'), num_win_lens = 12;  else, num_win_lens = cfg.num_win_lens;  end

% Form window lengths for which rescaled range is evaluated
win_lens = logspace(log10(min_win_len), log10(max_win_len), num_win_lens);
win_lens = floor(win_lens);
M = length(win_lens);
disp('Perform rescaled range analysis on each time series.')
fprintf("Window lengths of interest: %.4f sec (%d time points) to %.4f sec (%d time points).\n", ...
    min_win_len*fs, min_win_len, max_win_len*fs, max_win_len)

% Pre-allocate arrays to store results
RS = zeros(N,M);

% Loop over signals
for n = 1:N
    % get signal
    x = signal(n,:);

    % pre-allocate fluctuation function for this signal
    RS_n = zeros(1,M);

    % Loop over window lengths
    for m = 1:M
        % get window length (in time points) and compute the number of segments
        win_len = win_lens(m);
        n_segm = floor(T/win_len);

        % initialize rescaled range
        RS_running_sum = 0;

        % Loop over segments (or time windows)
        for w = 1:n_segm
            % prepare time segment
            idx = (1:win_len) + (w - 1)*win_len;    % time indices for the window
            x_win = x(idx);

            % de-mean the segment
            x_win = x_win - mean(x_win);

            % calculate signal profile
            y_win = cumsum(x_win);

            % range of profile
            R = range(y_win);

            % standard deviation of de-meaned signal
            S = std(x_win);

            % rescaled range of this segment
            RS_running_sum = RS_running_sum + R / S;
        end
        % average rescaled range, for this window length
        RS_n(m) = RS_running_sum / n_segm;
    end
    RS(n,:) = RS_n;
end

end