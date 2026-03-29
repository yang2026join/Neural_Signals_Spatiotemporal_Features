function permEn = permutation_entropy_seqs(seqs, cfg)
%PHASE_PERMUTATION_SEQS This function calculates permutation entropy of
%each sequence given in input "seqs". This function depends on function
%permutation_entropy.m.
%
% Permutation entropy is an entropy measure of sequences based on the
% ordering patterns in the sequence, which are obtained by sliding a window
% through the sequence and getting the order of values in each window.
% Permutation entropy is then calculated as the Shannon entropy of the
% ordering patterns. This function calculates the permutation entropy of
% each sequence given in the input "seqs".
%
%   permEn = permutation_entropy_seqs(seqs, cfg)
%
% Input:
%  seqs:    N x T, N = number of sequences, T = number of data points
%  cfg:     struct of configurations
%   cfg.win_size:   number of data points in a window; default: 3;
%   cfg.norml:      logical value (true or false) indicating whether to
%                   apply normalization by the maximal entropy under the
%                   given number of bins; default: false
% Output:
%  permEn:  N x 1 array, permutation entropy of each time series
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

% Get configurations
N = size(seqs,1);
if isfield(cfg, 'win_size'), win_size = cfg.win_size;  else, win_size = 3;  end
if isfield(cfg, 'norml'), norml = cfg.norml;  else, norml = cfg.false;  end
disp('Calculate permutation entropy for each sequence.')
fprintf("Window size: %d.\n", win_size);
if norml, disp('Normalization is applied.');  end

% Prepare for looping over time series
cfgIn = [];    % configurations passed to dependent function
cfgIn.win_size = win_size;
cfgIn.norml = norml;
permEn = zeros(N,1);    % pre-allocate array

% Loop over phase time series
for i = 1:N
    permEn(i) = permutation_entropy(seqs(i,:), cfgIn);
end

end