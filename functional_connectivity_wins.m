function [fc_tseries, time_center] = functional_connectivity_wins(signal, time, fs, cfg)
%FUNCTIONAL_CONNECTIVITY_WINS This function calculates the functional
%connectivity matrix in consecutive time windows, returning FC matrix as a
%function of time. This function depends on function
%FUNCTIONAL_CONNECTIVITY.
%
% Functional connectivity is the temporal correlation between two spatially
% remote neurophysiology events. The (i, j)-th element of a functional
% connectivity matrix is the correlation coefficient of a pair of neural
% signals recorded from site i and j. This function separates the given
% signals into multiple slices based on the given time window and step and
% calculate functional connectivity matrix for each slice of signals.
%
%   [fc_tseries, time] = functional_connectivity_dynamic(signal, time, fs, cfg)
%
% Input:
%  signal:  N x T, N = number of signals, T = number of time points
%  time:    1 x T, time points associated with the signals, in second
%  fs:      sampling frequency of the signals, in Hz
%  cfg:     struct of configurations
%   cfg.method: method to calculate correlation, must be one of the
%               "Pearson", "Spearman", "Kendall"; default: 'Pearson'
%   cfg.win:    length of each window in second; required
%   cfg.step:   length of moving step in second; required
% Output:
%  fc_tseries:  N x N x T', time series of functional connectivity matrix,
%               where T' is the number of windows
%  time_center: 1 x T', time points associated with each FC matrix, chosen
%               as the center of each time window
%
%                                                    YANG Hao, 2025-26 Fall
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
if ~isfield(cfg, 'method'), method = 'Pearson'; else, method = cfg.method; end
if ~isfield(cfg, 'win'), error('Please specify the window length!'); else, win = cfg.win; end
if ~isfield(cfg, 'step'), error('Please specify the step length!'); else, step = cfg.step; end

% time indices of each window
n_time = length(time);
win_tpnts = round(win*fs);
step_tpnts = round(step*fs);
n_wins = floor((n_time-win_tpnts)/step_tpnts) + 1;
time_idx = (1:win_tpnts) + step_tpnts*(0:(n_wins-1))';
fprintf("Window turns out to be %.4f s.\n", win_tpnts/fs)
fprintf("Step turns out to be %.4f s.\n", step_tpnts/fs)

% Prepare for looping over time windows
cfgIn = [];
cfgIn.method = method;
fc_tseries = zeros(size(signal,1),size(signal,1),n_wins);
time_center = zeros(1,n_wins);

% Loop over time windows
for i = 1:n_wins
    sig_temp = signal(:,time_idx(i,:));
    time_temp = time(time_idx(i,:));
    fc_tseries(:,:,i) = functional_connectivity(sig_temp, [], [], cfgIn);
    time_center(i) = median(time_temp);
end

end
