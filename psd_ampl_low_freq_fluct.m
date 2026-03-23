function alff = psd_ampl_low_freq_fluct(psd, freq, cfg)
%PSD_AMPL_LOW_FREQ_FLUCT This function calculates the amplitude of
%low-frequency fluctuations of resting-state fMRI signals.
%
% Amplitude of low-frequency fluctuations (ALFF) in resting-state fMRI
% analysis is defined as the averaged magnitude of Fourier transform over a
% low-frequency interval. This function calculates ALFF from given power
% spectra by integrating the square root of spectra. Power spectra data
% should be given as the magnitude of Fourier transform, not in dB.
%
%   alff = psd_ampl_low_freq_fluct(psd, freq, cfg)
%
% Input:
%  psd:     N x F, N = number of spectrums, F = number of frequency points
%  freq:    1 x F, frequency points associated with each spectrum, in Hz
%  cfg:     struct of configurations
%   cfg.freq_range: frequency range of interest, a 1 x 2 array; required
% Output:
%  alff:    N x 1, amplitude of low-frequency fluctuations
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

if ~isfield(cfg, 'freq_range'), error('Please specify a frequency range of interest!');
else, freq_range = cfg.freq_range; end

cfgIn = [];
cfgIn.xrange = freq_range;

alff = area_under_curve(sqrt(psd), freq, cfgIn) / (freq_range(2) - freq_range(1));

end