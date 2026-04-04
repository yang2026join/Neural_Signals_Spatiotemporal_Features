function shannonEn = shannon_entropy_seqs(seqs, cfg)
%SHANNON_ENTROPY_SEQS This function calculates Shannon entropy of each
%sequence given in input "seqs".
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
% each bin. This function calculates Shannon entropy of each sequence
% (e.g., time series) given in "seqs".
%
% Note: It's important to note that this function applies a fixed set of
% edges to bin each sequence. This is important to eliminating the effects
% of using different binnings.
%
%   shannonEn = shannon_entropy_seqs(seqs, cfg)
%
% Input:
%  seqs:    N x T, N = number of sequences, T = number of data points
%  cfg:     struct of configurations
%   cfg.range:      range of values to be divided into bins; default:
%                   smallest and largest values of data
%   cfg.num_bins:   number of bins to discretize the range of values;
%                   default: following an equation proposed by Scott (1992)
%   cfg.norml:      logical value (true or false) indicating whether to
%                   apply normalization by the maximal entropy under the
%                   given number of bins; default: false
% Output:
%  H:   N x 1 array, Shannon entropy of each phase time series
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
N = size(seqs,1);
seqs_flat = seqs(:);
if isfield(cfg, 'num_bins')
    num_bins = cfg.num_bins;
else
    bin_size = 3.49 * mean(std(seqs_flat)) * length(seqs_flat)^(-1/3);
    num_bins = round(2*pi/bin_size);
end
if ~isfield(cfg, 'norml'), norml = false;  else, norml = cfg.norml;  end
if ~isfield(cfg, 'range'), range = [min(seqs_flat) max(seqs_flat)];  else, range = cfg.range;  end
% verbose
disp('Calculate Shannon entropy for each sequence.')
if isfield(cfg, 'range'), fprintf("Range of values: [%.4f, %.4f].\n", range);
else, disp("Range of values: [min, max] of all the data.");  end
fprintf("Number of bins to divide the range: %d.\n", num_bins);
if norml, disp('Normalization is applied.');  end

% Form edges to bin the values
edges = linspace(range(1), range(2), num_bins+1);

% Loop over phase time series
shannonEn = zeros(N,1);
for i = 1:N

    % Count occurrences of each bin and calculate frequencies
    probs = histcounts(seqs(i,:), edges);  probs = probs / sum(probs);

    % Exclude bins with no occurrence
    probs = probs(probs > 0);
    
    % Compute Shannon entropy
    H = -sum(probs .* log2(probs));

    % Normalization by the maximum entropy if requied
    if norml, H = H / log2(num_bins); end

    % Save the result
    shannonEn(i) = H;
end

end