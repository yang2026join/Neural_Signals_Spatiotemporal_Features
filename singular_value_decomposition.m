function [spat_modes, temp_evol, singular_val, sig_recon] = singular_value_decomposition(signal, ~, ~, cfg)
%SINGULAR_VALUE_DECOMPOSITION This function performs SVD on neural signals
%at multiple sites and returns the spatial modes and temporal evolution, 
%squared singular values associated with each mode, as well as
%reconstructed signals using the first few selected modes.
%
% Singular Value Decomposition (SVD) is a powerful data-driven technique
% to decompose spatiotemporal signals into orthogonal spatial modes,
% associated temporal coefficients, and singular values that are ordered by
% the energy or importance of each mode. By SVD, a z-score normalized time
% series data matrix X (shape: N x T, recording sites x time points) is
% decomposed into X = U S V' (using economic SVD and assuming N < T),
%  U: N x N, each column represents a spatial mode;
%  S: N x N, a diagonal matrix of singular values in descending order
%  V: T x N, temporal evolution of each spatial mode, each row 
%     corresponding to one time point.
% You might realize that the spatial modes obtained in this way are exactly
% the principal components of the signals, taking each recording site as a
% feature and each time point as an observation.
%
%   [spat_modes, temp_evol, sigma_sqr] = eigenmode_signal_svd(signal, ~, ~, cfg)
%   [spat_modes, temp_evol, sigma_sqr, sig_recon] = eigenmode_signal_svd(signal, ~, ~, cfg)
%
% Input:
%  signal:  N x T, N = number of signals, T = number of time points
%  cfg:     struct of configurations
%   cfg.n_modes:    number of eigenmodes to return; default: all eigenmodes
%                   (i.e., minimum of N and T)
%   cfg.norml:      logical array (true or false), indicating whether to
%                   apply normalization to each signal; default: false
% Output:
%  spat_modes:  spatial modes (shape: N x k), i.e., the matrix U, where k
%               is the number of modes selected, i.e., cfg.n_modes;
%  temp_evol:   temporal evolution of each mode (shape: k x T), i.e., the
%               matrix V transposed
%  singular_val:    singular values (shape: k x 1) of each mode
%  sig_recon:   signals reconstructed from the first k modes
%
%                                                    YANG Hao, 2025-26 Fall
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
if ~isfield(cfg, 'n_modes'), k = min(size(signal,1), size(signal,2)); else, k = cfg.n_modes; end
if ~isfield(cfg, 'norml'), norml = false; else, norml = cfg.norml; end

% Option: normalization
if norml
    signal = zscore(signal, [], 2);
end

% Perform SVD
[U, S, V] = svd(signal, "econ");

% Rank truncation
spat_modes = U(:,1:k);
temp_evol = V(:,1:k)';
singular_val = S(1:k,1:k);

% Reconstruct signals if requested
if nargout == 4
    sig_recon = spat_modes * singular_val * temp_evol;
end

end