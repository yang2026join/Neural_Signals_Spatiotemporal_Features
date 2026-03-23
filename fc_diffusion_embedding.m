function fc_grad = fc_diffusion_embedding(fc, cfg)
%FC_DIFFUSION_EMBEDDING This function applies diffusion embedding algorithm
%to a functional connectivity matrix.
%
% Diffusion embedding is a nonlinear dimensionality-reduction technique
% based on a diffusion operator of a graph constructed from given data
% points and generates embeddings that respect local similarities in the
% data. When applied to a functional connectivity matrix, diffusion
% embedding generates FC gradients which describe hierarchies in functional
% networks. This function calculates FC gradients using diffusion embedding
% algorithm, with optional thresholding of elements in the FC and other
% diffusion embedding parameters.
%
% Key references:
% Margulies, D. S., Ghosh, S. S., Goulas, A., Falkiewicz, M., Huntenburg,
% J. M., Langs, G., ... & Smallwood, J. (2016). Situating the default-mode
% network along a principal gradient of macroscale cortical organization.
% Proceedings of the National Academy of Sciences, 113(44), 12574-12579.
%
%   fc_grad = fc_diffusion_embedding(fc, cfg)
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
%   cfg.diff_time:  steps of the constructed diffusion/Markov process;
%                   generally speaking, a larger diffusion time can reveal
%                   coarser structures; default: 3.
% Output:
%  fc_grad: N x k FC gradients, each column one gradient.
%
%                                                  YANG Hao, 2025-26 Spring
%            Centre for Nonlinear Studies, Hong Kong Baptist University, HK

if ~isfield(cfg, 'n_comps'), n_comps = min(size(fc,1),9); else, n_comps = cfg.n_comps; end
if ~isfield(cfg, 'thres'), thres = 1; else, thres = cfg.thres; end
if ~isfield(cfg, 'diff_time'), diff_time = 3; else, diff_time = cfg.diff_time; end
disp('Perform diffusion embedding on functional connectivity matrix.')
fprintf("Number of gradients to return: %d.\n", n_comps)
fprintf("Threshold FC by %d-precentile.\n", thres)

% Thresholding
fc_thres = fc;  thres_ = prctile(fc(:), thres);
fc_thres(fc < thres_) = thres_;

cfgIn = [];
cfgIn.n_comps = n_comps;
cfgIn.diff_time = diff_time;

fc_grad = diffusion_embedding(fc_thres, cfgIn);

mean_fc_grad_1 = mean(fc_grad(:,1));
if mean_fc_grad_1 > 1 || mean_fc_grad_1 < -1
    warning('Abnormal large gradient values!')
    fprintf("Mean of FC gradient 1: %.4f.\n", mean_fc_grad_1)
end

end