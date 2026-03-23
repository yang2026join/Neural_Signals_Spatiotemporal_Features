function [embed, eigval] = diffusion_embedding(data, cfg)
%DIFFUSION_EMBEDDING This function performs diffusion embedding on a set of
%F-dimensional data points (F: interpreted as the number of features)
%
% Diffusion embedding is an embedding method based on constructing a
% diffusion process on the data points to find meaningful structures in the
% dataset. The resulting embeddings are essentially eigenvectors of a
% Markov matrix constructed from the affinity matrix of the data points.
% Coifman & Lafon, Appl. Comput. Harmon. Anal. 21 (2006) 5–30
%
%   [embed, eigval] = diffusion_embedding(data, cfg)
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
%   cfg.n_comps:    number of eigenvalues and the associated eigenvectors
%                   (i.e., embeddings) to return, skipping the first
%                   eigenmode; default: all the eigenmodes (i.e., F)
%   cfg.alpha:      exponent for the normalization of affinity matrix, must
%                   between 0 and 1; generally, alpha corrects for
%                   non-uniform density in the data points; a larger alpha
%                   applies more correction to dense regions, while alpha =
%                   0 means to ignore difference in density across regions
%                   and dense regions will look similar in the affinity
%                   matrix; default: 0.5
%   cfg.diff_time:  steps of the constructed diffusion/Markov process;
%                   generally speaking, a larger diffusion time can reveal
%                   coarser structures, say large clusters in the network,
%                   while a smaller diffusion time focuses on finer or more
%                   local structures; default: 0.
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
if ~isfield(cfg, 'alpha'), alpha = 0.5; else, alpha = cfg.alpha; end
if ~isfield(cfg, 'diff_time'), t = 0; else, t = cfg.diff_time; end

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

% Normalize the affinity matrix with alpha
d = sum(A, 2);    % degree vector
D_alpha = diag(d.^(-alpha));
A_hat = D_alpha * A * D_alpha;

% Build a Markov matrix by normalization
% (which is the diffusion operator)
d_hat = sum(A_hat, 2);    % row sums
norm_factor = d_hat .^ (-1);
M = diag(norm_factor) * A_hat;    % Markov matrix: sum of each column is 1

% Eigen-decomposition of the Markov matrix
[eigvec, eigval] = eig(M, 'vector');
[eigval_sorted, idx] = sort(eigval, 'descend');    % sort by descending order
eigvec_sorted = eigvec(:,idx);

% Normalize the eigenvectors and skip the first trivial eigenmode
eigvec_norm = eigvec_sorted ./ eigvec_sorted(:,1);
eigvec_norm = eigvec_norm(:,2:end);
eigval = eigval_sorted(2:end);

% Select the first-k eigenmodes and normalize by diffusion time
if t == 0
    embed = eigvec_norm(:,1:k) .* (eigval(1:k)' ./ (1 - eigval(1:k)) );
else
    embed = eigvec_norm(:,1:k) .* (eigval(1:k)'.^t);
end

% Sign alignment
embed = embed .* sign(embed(end,:));

end