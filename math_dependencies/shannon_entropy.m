function H = shannon_entropy(data, ~, ~, cfg)
%SHANNON_ENTROPY This function calculates Shannon entropy of each phase
%time series given in input "phase".
%
% a short introduction to [Shannon entropy]...
%
%   H = shannon_entropy(phase, ~, ~, cfg)
%
% Input:
%  data:    an array with finite values in a certain range, typically,
%           a phase time series in [-pi, +pi] (resulting in phase Shannon
%           entropy, phase transfer entropy, etc.), functional connectivity
%           matrix in [0, 1] (resulting in FC diversity, etc.), Kuramoto
%           order parameter in [0, 1] (giving synchrony entropy).
%  cfg:     struct of configurations
%   cfg.range:      range of values to be divided; default: smallest and
%                   largest values of data
%   cfg.num_bins:   number of bins to discretize the range of values;
%                   default: following an equation proposed by Scott (1992)
%   cfg.norml:      logical value (true or false) indicating whether to
%                   apply normalization by the maximal entropy under the
%                   given number of bins; default: false
% Output:
%  H:   Shannon entropy of the data
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
data = data(:);
if isfield(cfg, 'num_bins')
    num_bins = cfg.num_bins;
else
    bin_size = 3.49 * mean(std(data)) * length(data)^(-1/3);
    num_bins = round(2*pi/bin_size);
end
if ~isfield(cfg, 'norml'), norml = false;  else, norml = cfg.norml;  end
if ~isfield(cfg, 'range'), range = [min(data) max(data)];  else, range = cfg.range;  end

% Assign each phase value to a bin
edges = linspace(range(1), range(2), num_bins+1);

% Count occurrences of each bin and calculate frequencies
counts = histcounts(data, edges);
probs = counts / sum(counts);    % sum(counts) equals length(data)

% Exclude bins with no occurrence
probs = probs(probs > 0);

% Compute Shannon entropy
H = -sum(probs .* log2(probs));

% Normalization by the maximum entropy with num_bins possible outcomes
if norml, H = H / log2(num_bins); end

end