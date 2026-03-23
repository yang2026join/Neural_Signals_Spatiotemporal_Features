function fc_grad = fc_spectral_embedding(fc, cfg)
%FC_SPECTRAL_EMBEDDING This function applies spectral embedding algorithm
%to a functional connectivity matrix.
%
% Spectral embedding is a linear dimensionality-reduction technique based
% a graph Laplacian constructed from the given data points. The embeddings
% represent global and local structures in the data. When applied to a
% functional connectivity matrix, the resulting embeddings reflect
% different functional networks in the brain. This function calculates FC
% embeddings with optional thresholding of each row of FC.
%
% Key references:
% Wang, R., Liu, M., Cheng, X., Wu, Y., Hildebrandt, A., & Zhou, C. (2021).
% Segregation, integration, and balance of large-scale resting brain
% networks configure different cognitive abilities. Proceedings of the
% National Academy of Sciences, 118(23), e2022288118.
%
%   fc_grad = fc_spectral_embedding(fc, cfg)
%
% Input:
%  fc:  N x N functional connectivity matrix calculated from Pearson
%       correlation coefficients
%  cfg: struct of configurations
%   cfg.n_comps:    number of components (i.e., gradients) to keep;
%                   default: min(N, 9)
%   cfg.thres:      percentile for thresholding each row of FC - elements
%                   in each row of FC matrix below this threshold will be
%                   set to 0; default: 1 (to exclude a few outliers);
% Output:
%  fc_embed:    N x k FC embeddings, each column one eigenmode.
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

if ~isfield(cfg, 'n_comps'), n_comps = min(size(fc,1),9); else, n_comps = cfg.n_comps; end
if ~isfield(cfg, 'thres'), thres = 1; else, thres = cfg.thres; end
disp('Perform spectral embedding on functional connectivity matrix.')
fprintf("Number of gradients to return: %d.\n", n_comps)
fprintf("Threshold FC by %d-precentile.\n", thres)

% Thresholding
fc_thres = fc;  thres_ = prctile(fc(:), thres);
fc_thres(fc < thres_) = thres_;

cfgIn = [];
cfgIn.n_comps = n_comps;

fc_grad = spectral_embedding(fc_thres, cfgIn);

mean_fc_embed_1 = mean(fc_grad(:,1));
if mean_fc_embed_1 > 1 || mean_fc_embed_1 < -1
    warning('Abnormal large embedding values!')
    fprintf("Mean of FC gradient 1: %.4f.\n", mean_fc_embed_1)
end

end