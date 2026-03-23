function [coh_eigvecs, coh_eigvals] = coh_eigendecomp(coh, freq, cfg)
%COH_EIGENDECOMP This function calculates the eigenvectors and associated
%eigenvalues of a coherence matrix at a specific frequency.
%
%   [coh_eigvecs, coh_eigvals] = coh_eigendecomp(coh, freq, cfg)
%
%Input:
%  coh:     N x N x F coherence matrix for different frequencies
%  freq:    1 x F, frequency points
%  cfg:     struct of configurations
%   cfg.foi:        frequency of interest, a single frequency value
%   cfg.n_modes:    number of eigenmodes selected; default: N (all modes)
% Output:
%  coh_eigvecs: N x k array, eigenvectors of the matrix, where k is the
%               number of modes selected
%  coh_eigvals: k x 1 array, eigenvalues in descending order
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

if ~isfield(cfg, 'foi'), error('Please specify a frequency of interest!'); else, foi = cfg.foi; end
if ~isfield(cfg, 'n_modes'), n_modes = size(coh,1); else, n_modes = cfg.n_modes; end
disp('Perform eigen-decomposition on coherence matrix at a specific frequency.')

freq_idx = dsearchn(freq', foi);
fprintf("Frequency selected: %.4f.\n", freq(freq_idx));
fprintf("Number of eigen-pairs to return: %d.\n", n_modes)

coh_sel = squeeze(coh(:,:,freq_idx));

cfgIn = [];
cfgIn.n_modes = n_modes;
[coh_eigvecs, coh_eigvals] = eigen_decomposition(coh_sel, cfgIn);

end