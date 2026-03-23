function [eigvecs, eigvals] = eigen_decomposition(M, cfg)
%EIGEN_DECOMPOSITION This function performs eigen-decomposition on a square
%matrix and returns the eigenvectors and eigenvalues. This function
%basically depends on Matlab built-in eig function, but allows user to
%select first k modes and format the eigenvectors so that the last
%components are positive.
%
%   [eigvecs, eigvals] = fc_eigendecomp(fc, cfg)
%
% Input:
%  M:   N x N array, N = number of brain regions
%  cfg: struct of configurations
%   cfg.n_modes:    number of eigenmodes to return; default: all eigenmodes
% Output:
%  eigvecs: N x k array, eigenvectors of the matrix, where k is the number
%           of modes selected
%  eigvals: k x 1 array, eigenvalues in descending order

% Get configurations
if ~isfield(cfg, 'n_modes'), k = size(M,1); else, k = cfg.n_modes; end

% Perform eigen-decomposition
[V, L] = eig(M, "vector");

% Sort the modes by decreasing eigenvalues
[L_sorted, idx] = sort(L, "descend");
L = L_sorted;
V = V(:, idx);

% Format the eigenvectors so that the last components are positive
V = V .* sign(V(end,:));

% Return the first k modes
eigvecs = V(:,1:k);
eigvals = L(1:k);

end