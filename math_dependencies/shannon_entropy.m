function H = shannon_entropy(data, cfg)
%SHANNON_ENTROPY This function calculates Shannon entropy of input data.
%
% Shannon entropy measures the uncertainty of a random variable with a
% finite number of outcomes. It is based on the probability distribution of
% the random variable. If the probabilities of just a few outcomes are
% high, while the probabilities of the other outcomes are small, in other
% words, we are more certain about which outcome will happen, then the
% entropy values is small. On the other extreme, if the random variable is
% uniformly distributed, then the entropy value is maximized. In practice,
% a random variable may take continuous values, and it's common to divide
% the range of values into multiple bins and count the frequencies of
% each bin. This function calculates Shannon entropy of data. Note that all
% data points will be used to construct the probability distribution
% regardless of the shape of "data" array.
%
% Key reference(s):
% Shannon, C. E. (1948). A mathematical theory of communication. The Bell
% system technical journal, 27(3), 379-423.
%
%   H = shannon_entropy(data, cfg)
%
% Input:
%  data:    an array with finite values in a certain range, typically,
%           a phase time series in [-pi, +pi] (resulting in phase Shannon
%           entropy, phase transfer entropy, etc.), functional connectivity
%           matrix in [0, 1] (resulting in FC diversity, etc.), Kuramoto
%           order parameter in [0, 1] (giving synchrony entropy).
%  cfg:     struct of configurations
%   cfg.range:      range of values to be divided into bins; default:
%                   smallest and largest values of data
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
data = data(:);    % vectorize the input array
if isfield(cfg, 'num_bins')
    num_bins = cfg.num_bins;
else
    bin_size = 3.49 * mean(std(data)) * length(data)^(-1/3);
    num_bins = round(2*pi/bin_size);
end
if ~isfield(cfg, 'norml'), norml = false;  else, norml = cfg.norml;  end
if ~isfield(cfg, 'range'), range = [min(data) max(data)];  else, range = cfg.range;  end
% verbose
disp('Calculate Shannon entropy for all input data.')
if isfield(cfg, 'range'), fprintf("Range of values: [%.4f, %.4f].\n", range);
else, disp("Range of values: [min, max] of the data.\n");  end
fprintf("Number of bins to divide the range: %d.\n", num_bins);
if norml, disp('Normalization is applied.');  end

% Form edges to bin the values
edges = linspace(range(1), range(2), num_bins+1);

% Count occurrences of each bin and calculate frequencies
probs = histcounts(data, edges);  probs = probs / sum(probs);

% Exclude bins with no occurrence
probs = probs(probs > 0);

% Compute Shannon entropy
H = -sum(probs .* log2(probs));

% Normalization by the maximum entropy if requied
if norml, H = H / log2(num_bins); end

end