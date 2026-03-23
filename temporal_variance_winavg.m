function [variance] = temporal_variance_winavg(signal, time, fs, cfg)
%TEMPORAL_VARIANCE_WINAVG This function calculates the temporal variance of
%multiple signals by averaging variance across multiple windows.
%
%   [variance] = temporal_variance_winavg(signal, ~, ~, cfg)
%
% Input:
%  signal:  N x T, N = number of signals, T = number of time points
%  time:    1 x T, time points associated with the signals, in second
%  fs:      sampling frequency of the signals, in Hz
%  cfg:     struct of configurations
%   cfg.win:    length of each window in second; required
%   cfg.step:   length of moving step in second; required
% Output:
%  variance:    N x 1 variance data
%
%                                                    YANG Hao, 2025-26 Fall
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
if ~isfield(cfg, 'win'), error('Please specify the window length!'); else, win = cfg.win; end
if ~isfield(cfg, 'step'), error('Please specify the step length!'); else, step = cfg.step; end

% Prepare time indices of each window
n_time = length(time);
win_tpnts = round(win*fs);
step_tpnts = round(step*fs);
n_wins = floor((n_time-win_tpnts)/step_tpnts) + 1;
time_idx = (1:win_tpnts) + step_tpnts*(0:(n_wins-1))';
fprintf("Window length turns out to be %.4f s.\n", win_tpnts/fs)
fprintf("Step length turns out to be %.4f s.\n", step_tpnts/fs)

% Loop over time windows
variance = zeros(size(signal,1),1);
for i = 1:n_wins
    sig_temp = signal(:,time_idx(i,:));
    var_temp = var(sig_temp, [], 2);
    variance = variance + var_temp;
end
variance = variance / n_wins;

end