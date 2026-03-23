function [coherence] = cps_coherence(cps)
%CPS_COHERENCE This function calculates the coherence of multiple signals
%given their cross power spectrum.
%
% Coherence is a quantity that measures the linear dependency between two
% signals, ranging from 0 to 1. A value of 1 indicates a perfect linear
% relationship, that is, one signal is well-modeled by another signal via
% a linear system. A value of 0 means the signals are either unrelated or
% the relationship is highly nonlinear.
%
%   [coherence] = cps_coherence(cps)
%
%Input:
%  cps: N x N x F complex-valued array, cross power spectrum data, N =
%       number of signals, F = number of frequency points
% Output:
%  coherence:   N x N x F array
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

coherence = zeros(size(cps));

n_freqs = size(coherence,3);
for ifreq = 1:n_freqs
    cps_mag = abs(squeeze(cps(:,:,ifreq)));    % |S_ij|
    psd = diag(squeeze(cps(:,:,ifreq)));    % S_ii

    cps_mag_sqr = cps_mag.^2;
    norml_factor = psd .* psd';

    coherence(:,:,ifreq) = cps_mag_sqr ./ norml_factor;
end

end