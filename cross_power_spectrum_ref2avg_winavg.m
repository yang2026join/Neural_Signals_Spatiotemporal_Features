function [cps_ref2avg, freq] = cross_power_spectrum_ref2avg_winavg(signal, time, fs, cfg)
%CROSS_POWER_SPECTRUM_REF2AVG_WINAVG This function calculates the cross
%power spectrum of multiple signals referenced to their spatial averaging.
%
% Cross power spectrum (CPS) measures the relationship between two signals
% in the frequency domain. The magnitude of CPS reveals how the two signals
% share energy at different frequencies, and the phase of CPS is the phase
% difference between the two signals' frequency components. This function
% calculates the cross power spectrum of each signal with the spatially
% averaged signal. CPS is calculated in each time window and then averaged
% across windows.
%
%   [cps_ref2avg, freq] = cross_power_spectrum_refavg_winavg(signal, time, fs, cfg)
%
% Input:
%  signal:  N x T, N = number of signals, T = number of time points
%  time:    1 x T, time points associated with the signals, in second
%  fs:      sampling frequency of the signals, in Hz
%  cfg:     struct of configurations
%   cfg.win:    length of each window in second; required
%   cfg.step:   length of moving step in second; required
%   cfg.taper:  taper applied to the signals, must be one of the following:
%               'hann', 'hamm', 'blackman', 'bartlett', 'none'; default:
%               'none', i.e., not tapering
% Output:
%  cps:     N x N x F cross power spectral density data
%  freq:    frequency points in Hz
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
if ~isfield(cfg, 'taper'), taper = 'none'; else, taper = cfg.taper; end
if ~isfield(cfg, 'win'), error('Please specify the window length!'); else, win = cfg.win; end
if ~isfield(cfg, 'step'), error('Please specify the step length!'); else, step = cfg.step; end
disp("Calculate cross power spectrum with spatial average using Welch's method.")

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
cps_ref2avg = zeros(size(signal,1),length(freq));    % pre-allocate array to store results
n_skips = 0;

% Spatially averaged signal as the reference signal
sig_ref = mean(signal, 1);

% Loop over time windows
for i = 1:n_wins

    % get time slice
    sig_temp = signal(:,time_idx(i,:));
    sig_ref_temp = sig_ref(:,time_idx(i,:));
    time_temp = time(time_idx(i,:));

    % check discontinuous time points
    if any(diff(time_temp) > 1.01/fs)
        n_skips = n_skips + 1;    % record and do not calculate
    else
        % FFT of signals and reference signal
        sig_hat = fft(sig_temp')';    sig_hat = sig_hat(:,1:floor(win_tpnts/2));    % N x F complex
        sig_ref_hat = fft(sig_ref_temp')';    sig_ref_hat = sig_ref_hat(:,1:floor(win_tpnts/2));    % 1 x F complex
        % cross power spectrum of each signal with the reference signal
        cps_temp = sig_hat .* conj(sig_ref_hat);
        cps_ref2avg = cps_ref2avg + cps_temp;
    end
end
cps_ref2avg = cps_ref2avg / win_tpnts^2;    % scale by signal length
cps_ref2avg = cps_ref2avg / (n_wins - n_skips);

end