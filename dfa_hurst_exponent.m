function dfa_exp = dfa_hurst_exponent(Fval, win_len, cfg)
%DFA_HURST_EXPONENT This function fits power law to each fluctuation
%function and returns the exponent in the fitted curve.
%
% Detrended fluctuation analysis (DFA) is a technique to estimate
% long-range temporal correlation in time series. This function receives
% the fluctuation functions obtained from detrended_fluctuation_analysis.m
% as input, fits power law to each fluctuation function, and returns the
% exponent.
%
%   dfa_exp = dfa_exponent(Fval, win_len, cfg)
%
% Input:
%  Fval:    N x L, N = number of signals, L = number of window lengths
%  win_len: 1 x L, window length associated with the fluctuation functions
%  cfg:     struct of configurations
%   cfg.win_len_range:  range of window length to fit, a 1 x 2 array;
%                       required
% Output:
%  dfa_exp: N x 1, exponent in the fitted power law
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

if ~isfield(cfg, 'win_len_range'), error('Please specify a range of window length!');
else, win_len_range = cfg.win_len_range; end

disp('Fit power law to each fluctuation function.')
fprintf("Window length to fit power law: [%.4f, %.4f] time points.\n", win_len_range(1), win_len_range(2));

cfgIn = [];
cfgIn.xrange = win_len_range;

dfa_exp = fit_power_law(Fval, win_len, cfgIn);
dfa_exp = -dfa_exp;

end
