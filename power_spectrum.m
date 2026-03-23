function [psd, freq] = power_spectrum(signal, time, fs, cfg)
%POWER_SPECTRUM This function calculates the power spectrum of multiple
%signals.
%
% Power spectrum, or power spectral density (PSD), is a fundamental signal
% analysis tool that shows how a signal's total energy is spread across
% different frequencies. This function estimates the PSD of each signal in
% the input 'signal', with multiple preprocessing options offered.
% 
%   [psd, freq] = power_spectrum(signal, time, fs, cfg)
%
% Input:
%  signal:  N x T, N = number of signals, T = number of time points
%  time:    1 x T, time points associated with the signals, in second
%  fs:      sampling frequency of the signals, in Hz
%  cfg:     struct of configurations
%   cfg.taper:  taper applied to the signals, must be one of the following:
%               'hann', 'hamm', 'blackman', 'bartlett', 'none'; default:
%               'none', i.e., not tapering
% Output:
%  psd:     N x M power density data, M = number of frequency points
%  freq:    frequency points, in Hz
%
%                                                    YANG Hao, 2025-26 Fall
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
if ~isfield(cfg, 'taper'), taper = 'none'; else, taper = cfg.taper; end

% Tapering
sig_temp = signal;  n_time = length(time);
switch taper
    case 'none'
        taper = ones(1, n_time);
    case 'hann'
        taper = 0.5 * (1 - cos(2*pi*(0:n_time-1)/(n_time-1)));
    case 'hamm'
        taper = 0.54 - 0.46 * cos(2*pi*(0:n_time-1)/(n_time-1));
    case 'blackman'
        taper = 0.42 - 0.5 * cos(2*pi*(0:n_time-1)/(n_time-1)) + ...
            0.08 * cos(4*pi*(0:n_time-1)/(n_time-1));
    case 'bartlett'
        taper = 1 - abs((0:n_time-1) - (n_time-1)/2) / ((n_time-1)/2);
    otherwise
        error("Taper must be one of the following: 'hann', 'hamm', 'blackman', 'bartlett'.")
end
sig_temp = sig_temp .* taper;
signal = sig_temp;  clear sig_temp;

% Prepare frequency span
n_time = length(time);
freq = (0:floor(n_time/2)-1)*fs/n_time;

% FFT of each signal
sig_hat = fft(signal')';   % Matlab fft returns the Fourier transform of each column.
sig_hat = sig_hat(:,1:floor(n_time/2));

% Power spectrum
psd = (abs(sig_hat) / n_time).^2;

end