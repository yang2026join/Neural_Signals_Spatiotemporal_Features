function [psd, freq] = power_spectrum_winavg(signal, time, fs, cfg)
%POWER_SPECTRUM_WINAVG This function calculates the power spectrum of
%multiple signals by averaging spectrum across multiple windows. This
%function depends on function POWER_SPECTRUM.
%
% Power spectrum, or power spectral density (PSD), is a fundamental signal
% analysis tool that shows how a signal's total energy is spread across
% different frequencies. This function separates the given signals into
% multiple slices based on the given time window and step, calculate
% PSD for each slice of signals, and then average the spectrum across
% slices or windows (essentially, the Welch's method).
%
%   [psd, freq] = power_spectrum_winavg(signal, time, fs, cfg)
%
% Input:
%  signal:  N x T, N = number of signals, T = number of time points
%  time:    1 x T, time points associated with the signals, in second
%  fs:      sampling frequency of the signals, in Hz
%  cfg:     struct of configurations, must include the following fields
%   cfg.win:    length of each window in second; required
%   cfg.step:   length of moving step in second; required
%   cfg.taper:  taper applied to the signals, must be one of the following:
%               'hann', 'hamm', 'blackman', 'bartlett', 'none'; default:
%               'none', i.e., not tapering
% Output:
%  psd:     N x M power density data, M = number of frequency points
%  freq:    frequency points in Hz
%
%                                                    YANG Hao, 2025-26 Fall
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
if ~isfield(cfg, 'taper'), taper = 'none'; else, taper = cfg.taper; end
if ~isfield(cfg, 'win'), error('Please specify the window length!'); else, win = cfg.win; end
if ~isfield(cfg, 'step'), error('Please specify the step length!'); else, step = cfg.step; end
disp("Calculate power spectrum using Welch's method.")

% Prepare time indices of each window
n_time = length(time);
win_tpnts = round(win*fs);
step_tpnts = round(step*fs);
n_wins = floor((n_time-win_tpnts)/step_tpnts) + 1;
time_idx = (1:win_tpnts) + step_tpnts*(0:(n_wins-1))';
fprintf("Window length turns out to be %.4f s.\n", win_tpnts/fs)
fprintf("Step length turns out to be %.4f s.\n", step_tpnts/fs)

% Prepare frequency span
freq = (0:floor(win_tpnts/2)-1)*fs/win_tpnts;

% Prepare for looping over time windows
cfgIn = [];    % configurations passed to the function for a single time window
cfgIn.taper = taper;
psd = zeros(size(signal,1),length(freq));
n_skips = 0;

% Loop over time windows
for i = 1:n_wins
    sig_temp = signal(:,time_idx(i,:));
    time_temp = time(time_idx(i,:));
    if any(diff(time_temp) > 1.01/fs)    % check discontinuous time points
        n_skips = n_skips + 1;    % record and do not calculate
    else
        [psd_temp, ~] = power_spectrum(sig_temp, time_temp, fs, cfgIn);
        psd = psd + psd_temp;
    end
end
psd = psd / (n_wins - n_skips);

end