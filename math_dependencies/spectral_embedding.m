function [embed, eigval] = spectral_embedding(data, cfg)
%SPECTRAL_EMBEDDING This function performs spectral embedding on a set of
%F-dimensional data points (F: interpreted as the number of features)
%
% Spectral embedding is a classical dimensionalilty reduction technique
% that applies eigendecomposition to graph Laplacian matrix. It seeks a
% low-dimensional representation of the high-dimensional feature data,
% which preserves local and global structures in the original data.
%
%   [embed, eigval] = spectral_embedding(data, cfg)
%
% Input:
%  data:    N x F, N = number of samples, F = number of features.
%           For example, an N-by-N functional connectivity matrix can be
%           interpreted as N samples of connectivity profiles (each row of
%           the FC matrix), each profile consisting of N features (N
%           columns).
%  cfg:     struct of configurations
%   cfg.distance:   distance function for calculating the pairwise distance
%                   between data points, should be one of the following:
%                   'euclidean', 'minkowski', 'chebychev', 'cosine';
%                   default: 'euclidean'
%   cfg.kernel:     kernel function for calculating similarity matrix from
%                   feature data, should be one of the following:
%                   'gaussian': Gaussian function, distance^2/(2*width^2)
%                   'nearest': nearest neighbors are assigned 1, others 0
%                   default: 'gaussian'
%   cfg.width:      width for Gaussian kernel; required when cfg.kernel is
%                   'gaussian'; default: 1 / F
%   cfg.n_neighbors:    number of neighbors k for nearest-neighbor kernel;]
%                       required when cfg.kernel is 'nearest'; default: 5
%   cfg.n_comps:    number of eigenvalues and the associated eigenvectors
%                   (i.e., embeddings) to return, skipping the first
%                   eigenmode; default: all the eigenmodes (i.e., F)
% Output:
%  embed:   embeddings (shape: N x k), eigenvectors of graph Laplacian
%           associated with the lowest-k eigenvalues
%  eigval:  eigenvalues associated with each embedding (shape: k x 1)
%
%                                                    YANG Hao, 2025-26 Fall
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
n_points = size(data,1);
n_features = size(data,2);
if ~isfield(cfg, 'distance'), distance_method = 'euclidean'; else, distance_method = cfg.distance; end
if ~isfield(cfg, 'kernel'), kernel_method = 'gaussian'; else, kernel_method = cfg.kernel; end
if ~isfield(cfg, 'n_comps'), k = n_features - 1; else, k = cfg.n_comps; end

% Calculate pairwise distance
distance = pdist(data, distance_method);
distance_squareForm = squareform(distance);

% Form similarity matrix (affine matrix)
% (using the specified method)
switch kernel_method

    case 'gaussian'
        if ~isfield(cfg, 'width'), width_type = 'median'; else, width_type = cfg.width; end
        if ischar(width_type)
            switch width_type
                case 'scaling'
                    sigma = max(distance) / sqrt(2*n_points);
                case 'median'
                    sigma = median(distance);
                otherwise
                    error("Please specify a width or use 'scaling' or 'median'!")
            end
        elseif isnumeric(width_type)
            sigma = width_type;
        else
            error('cfg.width must be a char or double')
        end
        A = exp(-distance_squareForm.^2) / (2*sigma^2);

    case 'nearest'
        if ~isfield(cfg, 'n_neighbors'), n_neighbors = 5; else, n_neighbors = cfg.n_neighbors; end
        A = zeros(n_points);
        for i = 1:size(A,1)
            [~, idx] = mink(distance_squareForm(i,:), n_neighbors);
            A(i, idx) = 1;
        end

    otherwise
        error('Unsupported kernel type!')

end

% Calculate graph Laplacian matrix and normalize
d = sum(A, 2);    % degree vector
L = diag(d) - A;    % unnormalized graph Laplacian
D_inv_sqrt = diag(1 ./ sqrt(d));    % D^(-1/2)
L_hat = D_inv_sqrt * L * D_inv_sqrt;    % normalized Laplacian

% Eigen-decomposition of the normalized graph Laplacian
[eigvec, eigval] = eig(L_hat, 'vector');
[eigval_sorted, idx] = sort(eigval, 'ascend');    % sort by ascending order
eigvec_sorted = eigvec(:,idx);

% Skip the first trivial eigenmode
eigval_sorted = eigval_sorted(2:end);
eigvec_sorted = eigvec_sorted(:,2:end);

% Sign alignment
eigvec_sorted = eigvec_sorted .* sign(eigvec_sorted(end,:));

% Return the results
eigval = eigval_sorted(1:k);
embed = eigvec_sorted(:,1:k);

end