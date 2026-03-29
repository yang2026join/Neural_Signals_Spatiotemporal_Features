function permEn = permutation_entropy(seq, cfg)
%PERMUTATION_ENTROPY This function calculates permutation entropy of a
%single sequence.
%
% Permutation entropy is an entropy measure of sequences based on the
% ordering patterns in the sequence, which are obtained by sliding a window
% through the sequence and getting the order of values in each window.
% Permutation entropy is then calculated as the Shannon entropy of the
% ordering patterns. This function calculates the permutation entropy of
% a sequence.
%
% Key reference(s):
% Bandt, C., & Pompe, B. (2002). Permutation entropy: a natural complexity
% measure for time series. Physical review letters, 88(17), 174102.
%
%   permEn = permutation_entropy(seq, cfg)
%
% Input:
%  seq:     a sequence given as a 1-D vector
%  cfg:     struct of configurations
%   cfg.win_size:   number of data points in a window; default: 3;
%   cfg.norml:      logical value (true or false) indicating whether to
%                   apply normalization by the maximal entropy under the
%                   given number of bins; default: false
% Output:
%  permEn:  permutation entropy
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Check input 
if ~isvector(seq)
    error('Input seq must be a row vector.');
end
N = length(seq);

% Get configurations
if isfield(cfg, 'win_size'), win_size = cfg.win_size;  else, win_size = 3;  end
if ~isfield(cfg, 'norml'), norml = false;  else, norml = cfg.norml;  end
n_wins = N - win_size + 1;

% Loop over each window and get ordering patterns for each
patterns = zeros(n_wins, win_size);    % pre-allocation
for i = 1:n_wins

    % get segment from the sequence
    slice = seq(i:i+win_size-1);

    % get rank indices
    [~, idx] = sort(slice);
    patterns(i,:) = idx;
end

% Count occurrences of each ordering pattern
unique_patterns = unique(patterns, "rows");
n_unique = size(unique_patterns,1);
counts = zeros(n_unique,1);
for i = 1:n_unique
    counts(i) = sum(ismember(patterns, unique_patterns(i,:), "rows"));
end

% Probability distribution
probs = counts / sum(counts);    % sum(counts) equals n_unique

% Compute Shannon entropy
permEn = -sum(probs .* log2(probs));

% Normalization by the maximum entropy for all possible outcomes
max_num_patterns = factorial(win_size);
if norml, permEn = permEn / log2(max_num_patterns); end

end