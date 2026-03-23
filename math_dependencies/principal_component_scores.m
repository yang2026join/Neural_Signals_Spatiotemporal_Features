function [pc_scores, pca_results] = principal_component_scores(data, xspan, cfg)
%PINCIPAL_COMPONENT_SCORES This function performs principal component
%analysis on a set of F-dimensional data points (F: number of features) and
%returns the principal component (PC) scores. This function basically
%depends on Matlab built-in function pca, but allows the user to select the
%range of coordinates and specify a cumulative variance explained.
%
%   [pc_scores, pca_results] = principal_component_scores(data, xspan, cfg)
%
% Input:
%  data:    N x F, N = number of samples, F = number of features.
%           For example, autocorrelation functions of N regions stacked in
%           an N-by-L matrix, where L is the number of time lags.
%  xspan:   1 x F, coordinates associated with each column in "data"
%           For example, the time lags for autocorrelation functions.
%  cfg:     struct of configurations
%   cfg.xrange:     interval of coordinates that will be used in PCA;
%                   default: xspan([1 end]) - all features will be used
%   cfg.cum_var:    required cumulative variance explained; default: 1.00 -
%                   all PCs will be kept
%   cfg.n_comps:    required number of components to keep; if specified,
%                   cumulative variance explained will be ignored.
% Output:
%  pc_scores:   N x F', principal component scores associated with each
%               sample, F' = number of PCs selected
%  pca_results: other PCA results, a struct including the following fields:
%   pca_results.num_pcs:        number of PCs selected, i.e., needed to
%                               exceed the required variance explained.
%   pca_results.cum_var_expl:   cumulative variance explained
%   pca_results.data_reconst:   data reconstructed from selected PCs
%
%                                                    YANG Hao, 2025-26 Fall
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
if isfield(cfg, "cum_var"), cum_var_required = cfg.cum_var; end
if isfield(cfg, "xrange"), xrange = cfg.xrange;  else, xrange = xspan([1 end]); end

% Perform PCA
data2use = data(:, xspan >= xrange(1) & xspan <= xrange(2));
[coeff, score, ~, ~, var_expl, mu] = pca(data2use);

% Find the number of PCs to keep
if isfield(cfg, "cum_var")
    cum_var_expl = cumsum(var_expl);
    num_pcs = find(cum_var_expl >= cum_var_required, 1, 'first');
else
    num_pcs = length(var_expl);
end
if isfield(cfg, "n_comps"),  num_pcs = cfg.n_comps;  end

% Return the results
pc_scores = score(:,1:num_pcs);
pca_results = [];
pca_results.num_pcs = num_pcs;
pca_results.cum_var_expl = cumsum(var_expl(1:num_pcs));
pca_results.data_reconst = score(:,1:num_pcs) * coeff(:,1:num_pcs)' + mu;

disp(['Number of components: ' num2str(num_pcs)])
disp(['Cumulative variance explained: ' num2str(pca_results.cum_var_expl') ' (%)'])

end