function acw = acf_autocorrelation_window(acf, lag, cfg)
%ACF_AUTOCORRELATION_WINDOW This function estimates the timescale for each
%autocorrelation function in "acf", calculated as the lag at which
%autocorrelation decays to a given value, i.e., autocorrelation window in
%technical terms.
%
% Human resting-state neural signals like fMRI, EEG, and ECoG typically
% show decaying autocorrelation functions, reflecting the "memory" of brain
% activity. How fast an autocorrelation function decays, quantified by
% intrinsic neural timescales, measures the capacity of this memory. This
% function quantifies timescale as the autocorrelation window, that is, the
% time lag at which ACF decays to a given value (between 0 and 1). To
% improve numerical accuracy, linear interpolation is performed using the
% two time lags at which the ACFs are closet to the target value.
%
%   acw = acf_autocorrelation_window(acf, lag, cfg)
%
% Input:
%  acf: N x L, N = number of autocorrelation functions, L = number of lags
%  lag: 1 x L, time lags associated with each autocorrelation curve, in (s)
%  cfg: struct of configurations
%   cfg.target_ac:  target autocorrelation value; required
% Output:
%  acw: autocorrelation windows for each autocorrelation function
%
%                                                    YANG Hao, 2025-26 Fall
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

if ~isfield(cfg, 'target_ac'), error('Please specify a target AC value!'); else, target_ac = cfg.target_ac; end
fprintf("Find the lag at which autocorrelation value decays to %.4f.\n", target_ac);

acw = zeros(size(acf,1), 1);
for i = 1:size(acf,1)
    acw(i) = find_roots_of_decreasing_functions(acf(i,:), lag, target_ac);
end

function x0 = find_roots_of_decreasing_functions(f, x, target)
x0 = zeros(size(f,1),1);
for j = 1:size(f,1)
    arr = f(j,:);
    idx = find(arr > target, 1, 'last');
    if isempty(idx)    % Target is above all values
        idx1 = length(arr)-1;
        idx2 = length(arr);
    elseif idx == 1    % Target is below all values
        idx1 = 1;
        idx2 = 2;
    else
        idx1 = idx;
        idx2 = idx+1;
    end
    x0(j) = interp1(arr([idx1 idx2]), x([idx1 idx2]), target, 'linear');
end
end

end