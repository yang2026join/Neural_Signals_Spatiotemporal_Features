function divergence = divergence_prob_dist(P, Q, cfg)
%DIVERGENCE_PROB_DIST This function offers multiple divergenve measures of
%two probability distributions.
%
%   divergence = divergence_prob_dist(P, Q)
%
% Input:
%  P, Q:    two probability distributions, must be two vectors of equal
%           length, each vector being non-negative, and elements sum to 1.
%  cfg:     struct of configurations
%   cfg.method: should be one of the following; required
%               'kl':   Kullback Leibler divergence
%               'js':   Jensen Shannon divergence
%               'l1_norm':  L-1 norm of vector (P - Q).
%               'max_norm': max norm of vector (P - Q), so-called total
%                           variation distance
%               'hellinger':    Hellinger Distance
% Output:
%  divergence:  a single scalar

% Check input
if ~isapprox(sum(P), 1)
    error('Probabilities in vector P does not sum to 1!')
end
if ~isapprox(sum(Q), 1)
    error('Probabilities in vector Q does not sum to 1!')
end
if any(P < 0)
    error('Negative entry(ies) in vector P!')
end
if any(Q < 0)
    error('Negative entry(ies) in vector Q!')
end

% Get method
if ~isfield(cfg, 'method'), error('Please specify a divergence measure!'); else, method = cfg.method; end

% Switch for different methods
switch method
    case 'kl'
        mask = P > 0;
        if any(Q(mask) == 0)
            error('Q has 0 entries where P is nonzero, KL divergence is infinite!')
        end
        divergence = sum(P(mask) .* log(P(mask) ./ Q(mask)));

    case 'js'
        M = 0.5 * (P + Q);
        % For terms P>0 where M>0
        mask_P = P > 0;
        Dkl_PM = sum(P(mask_P) .* log(P(mask_P) ./ M(mask_P)));
        % For terms Q>0 where M>0
        mask_Q = Q > 0;
        Dkl_QM = sum(Q(mask_Q) .* log(Q(mask_Q) ./ M(mask_Q)));
        divergence = 0.5 * (Dkl_PM + Dkl_QM);

    case 'l1_norm'
        divergence = sum(abs(P - Q));
        
    case 'max_norm'    % i.e., total variation distance
        divergence = max(abs(P - Q));
        
    case 'hellinger'
        divergence = sqrt(0.5 * sum((sqrt(P) - sqrt(Q)).^2));
        
    otherwise
        error("Method name must be one of the following: 'kl', 'js', 'l1_norm', 'max_norm', 'hellinger'")
end


end