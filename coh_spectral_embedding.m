function [coh_embeds, eigvals] = coh_spectral_embedding(coh, freq, cfg)
%COH_EIGENDECOMP This function calculates the eigenvectors and associated
%eigenvalues of a coherence matrix at a specific frequency.
%
%   [coh_embeds, eigvals] = coh_spectral_embedding(coh, freq, cfg)
%
%Input:
%  coh:     N x N x F coherence matrix for different frequencies
%  freq:    1 x F, frequency points
%  cfg:     struct of configurations
%   cfg.foi:        frequency of interest, a single frequency value
%   cfg.n_comps:    number of embeddings to keep; default: N - 1 (skipping
%                   the first mode)
% Output:
%  coh_embeds:  N x k array, k embeddings selected
%  coh_eigvals: k x 1 array, eigenvalues in descending order
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

if ~isfield(cfg, 'foi'), error('Please specify a frequency of interest!'); else, foi = cfg.foi; end
if ~isfield(cfg, 'n_comps'), n_comps = size(coh,1) - 1; else, n_comps = cfg.n_comps; end
disp('Perform spectral embedding on coherence matrix at a specific frequency.')
fprintf("Number of components to return: %d.\n", n_comps)

freq_idx = dsearchn(freq', foi);
fprintf("Frequency selected: %.4f Hz.\n", freq(freq_idx));

coh_sel = squeeze(coh(:,:,freq_idx));

cfgIn = [];
cfgIn.n_comps = n_comps;
[coh_embeds, eigvals] = spectral_embedding(coh_sel, cfgIn);

end