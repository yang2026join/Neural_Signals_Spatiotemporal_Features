function psd_exp = psd_power_law_exponent(psd, freq, cfg)
%PSD_POWER_LAW_EXPONENT This function fits power law to each spectrum and
%returns the exponent in the fitted curve.
%
% Human resting-state signals like EEG/MEG and fMRI typically show a
% power-law or "1/f" structure in power spectrum. A larger power-law
% exponent means the spectrum decays faster, which implies a larger
% proportion of low-frequency power in the signal. This function fits one
% power law to each spectrum and returns the exponent. The power spectra
% data should be given as the magnitude of Fourier transform, not in dB.
%
%   psd_exp = psd_power_law_exponent(psd, freq, cfg)
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

cfgIn = [];
if ~isfield(cfg, 'freq_range'), error('Please specify a frequency range of interest!');
else, cfgIn.xrange = cfg.freq_range; end

psd_exp = fit_power_law(psd, freq, cfgIn);

end