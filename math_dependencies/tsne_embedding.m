function [embeds] = tsne_embedding(data, cfg)
%TSNE_EMBEDDING This function performs t-distributed stochastic embedding
%(t-SNE) on observed features at multiple recording sites.
%
% t-SNE is a non-linear dimensionalilty reduction technique that reduces
% high-dimensional data for 2D or 3D visualization. Student t-distribution
% is used to measure similarities between low-dimensional points.
%
%   [embeds] = tsne_embedding(data, cfg)
%
% Input:
%  data:    N x F, N = number of samples, F = number of features.
%           For example, an N-by-N functional connectivity matrix can be
%           interpreted as N samples of connectivity profiles (each row of
%           the FC matrix), each profile consisting of N features (N
%           columns).
%  cfg:     struct of configurations
%   cfg.distance:   distance function for calculating the pairwise distance
%                   between data points, should be one of the values in
%                   MATLAB function pdist. Typical examples include:
%                   'euclidean', 'minkowski', 'chebychev', 'cosine',
%                   'mahalanobis'; default: 'euclidean'
%   cfg.n_dims:     number of dimensions to return; default: 2
% Output:
%  embeds:  embeddings (shape: N x k)
%
%                                                    YANG Hao, 2025-26 Fall
%            Center for Nonlinear Studies, Hong Kong Baptist University, HK

if ~isfield(cfg, 'distance'), distance_method = 'euclidean'; else, distance_method = cfg.distance; end
if ~isfield(cfg, 'n_dims'), k = 2; else, k = cfg.n_comps; end

embeds = tsne(data, "Distance", distance_method, "NumDimensions", k);

end