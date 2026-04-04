function fc_entropy = fc_shannon_entropy(fc, cfg)
%FC_SHANNON_ENTROPY This function constructs a probability
%distribution of FC values by discretizing them into bins and then
%calculates the Shannon entropy. This function depends on shannon_entropy.m
%
% Shannon entropy of FC measures the richness or diversity in connection
% strength. A higher FC entropy means that the FC values are more uniformly
% distributed, indicating that the brain is more balanced between strong
% and weak connections.
%
%   fc_entropy = fc_shannon_entropy(fc, cfg)
%
% Input:
%  fc:      N x N FC matrix
%  cfg:     struct of configurations
%   cfg.num_bins:   number of bins to discretize FC values; default: 30
% Output:
%  fc_entropy:  normalized Shannon entropy of FC, a scalar between 0 and 1.
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
if ~isfield(cfg, 'num_bins'), num_bins = 30; else, num_bins = cfg.num_bins; end
%disp('Calculate Shannon entropy of FC values.')
%fprintf("Number of bins to discretize FC values: %d.\n", num_bins)

% Pass FC values to Shannon entropy function
cfgIn = [];
cfgIn.num_bins = num_bins;
cfgIn.range = [0, 1];    % exclude negative correlations
cfgIn.norml = true;
fc_tri = triu(fc,1);  fc_tri = fc_tri(:);  fc_tri = fc_tri(fc_tri ~= 0);
fc_entropy = shannon_entropy(fc_tri, cfgIn);

end