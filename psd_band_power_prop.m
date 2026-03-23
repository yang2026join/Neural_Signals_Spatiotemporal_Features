function power_prop = psd_band_power_prop(psd, freq, cfg)
%PSD_BAND_POWER_PROP This function calculates the proportion of power over
%a given frequency interval in the total power.
%
% For each power spectrum in "psd", calculate the proportion of power over
% a given band to that over the whole band. Power is calculated by
% integrating the spectrum over a frequency interval with trapezoid method.
% Power spectra data should be given as the magnitude of Fourier transform,
% not in dB.
%
%   power_prop = spectrum_band_power_prop(psd, freq, cfg)
%
% Input:
%  psd:  N x F, N = number of spectrums, F = number of frequency points
%  freq:    1 x F, frequency points associated with each spectrum, in Hz
%  cfg:     struct of configurations
%   cfg.foi:    frequency range of interest, a 1 x 2 array; required
%   cfg.fwhole: frequency range for the total power, a 1 x 2 array;
%               default: the whole frequency range
% Output:
%  power_prop:  N x 1, proportion of power for each spectrum
%
%                                                    YANG Hao, 2025-26 Fall
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

if ~isfield(cfg, 'foi'), error('Please specify a frequency range of interest!'); else, foi = cfg.foi; end
if isfield(cfg, "fwhole"), fwhole = cfg.fwhole; else, fwhole = freq; end

addpath('math_dependencies\')

cfgIn = [];
cfgIn.xrange = foi;
power_foi = area_under_curve(psd, freq, cfgIn);

cfgIn = [];
cfgIn.xrange = fwhole;
power_fwhole = area_under_curve(psd, freq, cfgIn);

power_prop = power_foi ./ power_fwhole;

end