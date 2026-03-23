function auc = acf_area_under_curve(acf, lag, cfg)
%ACF_AREA_UNDER_CURVE This function estimates the timescale for each
%autocorrelation function in "acf", calculated as the area under ACF curve
%over the initially decaying lag regime.
%
% Human resting-state neural signals like fMRI, EEG, and ECoG typically
% show decaying autocorrelation functions, reflecting the "memory" of brain
% activity. How fast an autocorrelation function decays, namely, the
% intrinsic neural timescales, measures the capacity of this memory. This
% function quantifies timescale as the area under ACF curve over an initial
% lag interval in which ACF decays.
%
%   auc = acf_area_under_curve(acf, lag, cfg)
%
% Input:
%  acf: N x L, N = number of autocorrelation functions, L = number of lags
%  lag: 1 x L, time lags associated with each autocorrelation curve, in (s)
%  cfg: struct of configurations
%   cfg.max_lag:    largest lag value to integrate the ACF curves; required
% Output:
%  auc: area under curve for each autocorrelation function
%
%                                                    YANG Hao, 2025-26 Fall
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

if ~isfield(cfg, 'max_lag'), error('Please specify a lag interval!'); else, lag_range = [0.0001 cfg.max_lag]; end
fprintf("Calculate area under ACF curve for lag interval [%.4f, %.4f] sec.\n", lag_range(1), lag_range(2));

cfg = [];
cfg.xrange = lag_range;

auc = area_under_curve(acf, lag, cfg);

end