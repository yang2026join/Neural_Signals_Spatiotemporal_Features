function pc_scores = psd_principal_component_scores(psd, freq, cfg)
%PSD_PRINCIPAL_COMPONENT_SCORES This function performs PCA on power spectra
%data at multiple recording sites and returns the PC scores.
%
% This function performs PCA on the power spectra data given in the "psd"
% array, taking each frequency point as a feature and each spectrum as one
% observation. This captures spatial variation in spectra. The output PC
% scores are each spectrum projected to corresponding PC axes. The power
% spectra data are given as the magnitude of Fourier transform, not in dB.
%
%   pc_scores = psd_principal_component_scores(psd, freq, cfg)
%
% Input:
%  psd:     N x F, N = number of spectrums, F = number of frequency points
%  freq:    1 x F, frequency points associated with each spectrum, in Hz
%  cfg:     struct of configurations
%   cfg.freq_range: frequency range of interest, a 1 x 2 array; default:
%                   spectra over the whole band will be passed to PCA.
%   cfg.n_comps:    number of principal components to keep; default: keep
%                   all components
% Output:
%  pc_scores:   N x F', PC scores for each spectrum, F' = number of PCs
%               selected
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

cfgIn = [];
if isfield(cfg, 'freq_range'), cfgIn.xrange = cfg.freq_range; end
if isfield(cfg, 'n_comps'), cfgIn.n_comps = cfg.n_comps; end

disp('Perform PCA on power spectra data.')
fprintf("Frequency range of interest: [%.4f, %.4f] Hz.\n", cfg.freq_range(1), cfg.freq_range(2));

psd_dB = 10*log10(psd);

pc_scores = principal_component_scores(psd_dB, freq, cfgIn);

end