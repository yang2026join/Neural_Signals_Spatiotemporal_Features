function falff = psd_frac_ampl_low_freq_fluct(psd, freq, cfg)
%PSD_FRAC_AMPL_LOW_FREQ_FLUCT This function calculates the fraction of
%amplitude of low-frequency fluctuations for resting-state fMRI signals.
%
% Amplitude of low-frequency fluctuations (ALFF) in resting-state fMRI
% analysis is defined as the averaged magnitude of Fourier transform over a
% low-frequency interval. This function calculates the ratio of ALFF over a
% low-frequency interval in ALFF over a wide interval. Power spectra data
% should be given as the magnitude of Fourier transform, not in dB.
%
%   falff = psd_frac_ampl_low_freq_fluct(psd, freq, cfg)
%
% Input:
%  psd:  N x F, N = number of spectrums, F = number of frequency points
%  freq:    1 x F, frequency points associated with each spectrum, in Hz
%  cfg:     struct of configurations
%   cfg.foi:    frequency range of interest, a 1 x 2 array giving the lower
%               and upper bound; required
%   cfg.fwhole: frequency range for the total power, also a 1 x 2 array; if
%               not specified, the total power will be calculated by
%               integrating the whole spectrum
% Output:
%  falff:   N x 1, fraction of ALFF
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

if ~isfield(cfg, 'foi'), error('Please specify a frequency range of interest!'); else, foi = cfg.foi; end
if isfield(cfg, "fwhole"), fwhole = cfg.fwhole; else, fwhole = freq; end

cfgIn = [];
cfgIn.xrange = foi;
power_foi = area_under_curve(sqrt(psd), freq, cfgIn);

cfgIn = [];
cfgIn.xrange = fwhole;
power_fwhole = area_under_curve(sqrt(psd), freq, cfgIn);

falff = power_foi ./ power_fwhole;

end