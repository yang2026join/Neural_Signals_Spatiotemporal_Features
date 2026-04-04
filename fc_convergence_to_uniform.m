function fc_diversity = fc_convergence_to_uniform(fc, cfg)
%FC_CONVERGENCE_TO_UNIFORM This function constructs a probability
%distribution of FC values by discretizing them into bins and then
%calculates convergence (1 - divergence) of FC distribution to a
%corresponding uniform distribution. This function depends on
%divergence_prob_dist.m which offers multiple divergence measures between
%two probability distributions.
%
% Divergence of the probability distribution of FC values to a uniform
% distribution measures the degree to which FC values are concentrated at
% certain values. A smaller divergence, or a larger convergence, means that
% the FC values are more uniformly distributed, indicating that the brain
% is more balanced between strong and weak connections, so this convergence
% is informally named as FC diversity in literature, e.g., Xu et al. (2021)
% defines FC diversity as such convergence measured by L-1 norm, normalized
% by the maximal divergence.
%
% Key reference(s):
% Xu, L., Feng, J., & Yu, L. (2022). Avalanche criticality in individuals,
% fluid intelligence, and working memory. Human brain mapping, 43(8),
% 2534-2553.
%
%   fc_diversity = fc_divergence_to_uniform(fc, cfg)
%
% Input:
%  fc:      N x N FC matrix
%  cfg:     struct of configurations
%   cfg.num_bins:   number of bins to discretize FC values; default: 30
%   cfg.method: method to calculate divergence, must be one of the
%               following, 'kl', 'js', 'l1_norm', 'max_norm', 'hellinger',
%               referring to divergence_prob_dist.m; default: 'l1_norm'
% Output:
%  fc_diversity:    a scalar between 0 and 1.
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
if ~isfield(cfg, 'num_bins'), num_bins = 30; else, num_bins = cfg.num_bins; end
if ~isfield(cfg, 'method'), method = 'l1_norm'; else, method = cfg.method; end
disp('Calculate FC diversity by convergence to uniform distribution.')
fprintf("Number of bins to discretize FC values: %d.\n", num_bins)
disp(['Method to calculate divergence: ' method])

% Construct probability distribution of FC values
fc_tri = triu(fc,1);  fc_tri = fc_tri(:);  fc_tri = fc_tri(fc_tri ~= 0);
edges = linspace(0, 1, num_bins+1);    % this excludes negative FC values
fc_probs = histcounts(fc_tri, edges);  fc_probs = fc_probs / sum(fc_probs);

% Uniform distribution of the same size
uniform_dist = ones(1,num_bins) * (1 / num_bins);

% Pass to divergence function
cfgIn = [];
cfgIn.method = method;
fc_div = divergence_prob_dist(fc_probs, uniform_dist, cfgIn);

% Normalization factor: divergence of Delta to uniform
delta = zeros(1,num_bins);  delta(1) = 1;
norm_factor = divergence_prob_dist(delta, uniform_dist, cfgIn);

fc_diversity = 1 - fc_div / norm_factor;
end