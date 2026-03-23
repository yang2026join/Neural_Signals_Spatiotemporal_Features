function fc = functional_connectivity(signal, ~, ~, cfg)
%FUNCTIONAL_CONNECTIVITY This function calculates the functional
%connectivity matrix of neural signals recorded from multiple sites
%(EEG/MEG channels, regions in a brain atlas, etc.)
%
% Functional connectivity is the temporal correlation between two spatially
% remote neurophysiology events. The (i, j)-th element of a functional
% connectivity matrix is the correlation coefficient of a pair of neural
% signals recorded from site i and j.
%
%   fc = functional_connectivity(signal, ~, ~, cfg)
%
% Input:
%  signal:  N x T, N = number of signals, T = number of time points
%  cfg:     struct of configurations
%   cfg.method: method of correlation coefficient, must be one of the
%               'Pearson', 'Spearman', 'Kendall'; default: 'Pearson'
% Output:
%  fc:      N x N functional connectivity data
%
%                                                    YANG Hao, 2025-26 Fall
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

if ~isfield(cfg, 'method'), method = 'Pearson'; else, method = cfg.method; end

fc = corr(signal', 'Type', method);

end

