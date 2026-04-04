function hurst_exp = rs_hurst_exponent(RS, win_len, cfg)
%RS_HURST_EXPONENT This function fits power law to each rescaled range and
%returns the Hurst exponent in the fitted curve.
%
% Rescaled range analysis is a technique to estimate long-range
% correlations in a time series. This function receives the rescaled ranges
% for different window lengths, obtained from rescaled_range_analysis.m,
% fits power law to each RS, and returns the Hurst exponent.
%
%   rs_hurst_exponent(RS, win_lens, cfg)
%
% Input:
%  RS:          N x L, N = number of signals, L = number of window lengths
%  win_lens:    1 x L, window length associated with the RS values
%  cfg:     struct of configurations
%   cfg.win_len_range:  range of window length to fit, a 1 x 2 array;
%                       required
% Output:
%  hurst_exp:   N x 1, exponent in the fitted power law
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

if ~isfield(cfg, 'win_len_range'), error('Please specify a range of window length!');
else, win_len_range = cfg.win_len_range; end

disp('Fit power law to each rescaled range function.')
fprintf("Window length to fit power law: [%.4f, %.4f] time points.\n", win_len_range(1), win_len_range(2));

cfgIn = [];
cfgIn.xrange = win_len_range;

hurst_exp = fit_power_law(RS, win_len, cfgIn);
hurst_exp = -hurst_exp;

end