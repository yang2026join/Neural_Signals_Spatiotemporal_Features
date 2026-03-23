function mean_psd = psd_band_mean_psd(psd, freq, cfg)
%PSD_BAND_MEAN_PSD This function calculates the average spectral
%density of a power spectrum in a given frequency range.
%
% For each power spectrum in "psd" array, calculate the average spectral
% density in a given frequency range by integrating the spectrum over the
% band using trapezoid method and then dividing the resulting total power
% by the band width. Power spectra data should be given as the magnitude of
% Fourier transform, not in dB.
%
%   mean_psd = psd_band_mean_psd(psd, freq, cfg)
%
% Input:
%  psd:     N x F, N = number of spectrums, F = number of frequency points
%  freq:    1 x F, frequency points associated with each spectrum, in Hz
%  cfg:     struct of configurations
%   cfg.freq_range: frequency range of interest, a 1 x 2 array; required
% Output:
%  mean_psd:    N x 1, mean spectral density of each spectrum
%
%                                                    YANG Hao, 2025-26 Fall
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

if ~isfield(cfg, 'freq_range'), error('Please specify a frequency range of interest!');
else, freq_range = cfg.freq_range; end

cfgIn = [];
cfgIn.xrange = freq_range;

mean_psd = area_under_curve(psd, freq, cfgIn) / (freq_range(2) - freq_range(1));

end