function [fc_eigvecs, fc_eigvals] = fc_eigendecomp(fc, cfg)
%FC_EIGENDECOMP This function calculates the eigenvectors and associated
%eigenvalues of a functional connectivity matrix.
%
% Eigenmodes of FC represent different modes of brain activity. If the FC
% matrix is calculated by Pearson correlation coefficients, then the
% eigenvectors of FC are exactly the principal components of z-score
% normalized signals, taking each region as a feature and each time point
% as an observation.
%
%   [fc_eigvecs, fc_eigvals] = fc_eigendecomp(fc, cfg)
%
% Input:
%  fc:      N x N FC matrix
%  cfg:     struct of configurations
%   cfg.n_modes:    number of eigenmodes selected; default: N (all modes)
% Output:
%  fc_eigvecs:  N x k array, eigenvectors of the matrix, where k is the
%               number of modes selected
%  fc_eigvals:  k x 1 array, eigenvalues in descending order
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

if ~isfield(cfg, 'n_modes'), n_modes = size(fc,1); else, n_modes = cfg.n_modes; end
disp('Perform eigen-decomposition on FC matrix.')
fprintf("Number of eigen-pairs to return: %d.\n", n_modes)

cfgIn = [];
cfgIn.n_modes = n_modes;
[fc_eigvecs, fc_eigvals] = eigen_decomposition(fc, cfgIn);

end