function tau0 = acf_time_constant(acf, lag, cfg)
%TIMESCALE_TIME_CONSTANT This function estimates the timescale for each
%autocorrelation function in "acf" by fitting an exponential decay to each
%of them. Timescale is the time constant in the fitted exponential decay.
%
% Human resting-state neural signals like fMRI, EEG, and ECoG typically
% show decaying autocorrelation functions, reflecting the "memory" of brain
% activity. How fast an autocorrelation function decays, quantified by
% intrinsic neural timescales, measures the capacity of this memory. This
% function assumes exponential decay, AC(tau) = A*exp(tau/tau0), for the
% given autocorrelation data, fits the data to each autocorrelation, and
% obtain tau0 as the timescales.
%
%   tau0 = acf_fit_exponential(acf, lag, cfg)
%
% Input:
%  acf: N x L, N = number of autocorrelation functions, L = number of lags
%  lag: 1 x L, time lags associated with each autocorrelation curve, in (s)
%  cfg: struct of configurations, must contain the following fields
%   cfg.lag_range:  range of lag to fit the exponential function; default:
%                   fit all data points
% Output:
%  tau0:        N x 1, time constant in the fitted autocorrelation curves
%
%                                                    YANG Hao, 2025-26 Fall
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

if ~isfield(cfg, 'lag_range'), lag_range = [0, max(lag)];  else, lag_range = cfg.lag_range; end
fprintf("Lag interval to fit: [%.4f, %.4f] sec.\n", lag_range(1), lag_range(2)); 

cfgIn = [];
cfgIn.xrange = lag_range;

tau0 = fit_exponetial_decay(acf, lag, cfgIn);

end