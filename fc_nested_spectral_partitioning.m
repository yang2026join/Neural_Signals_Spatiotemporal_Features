function [fc_rearranged, cluster_size, cluster_number] = fc_nested_spectral_partitioning(fc)
%FC_NESTED_SPECTRAL_PARTITIONING This function performs nested network
%partitioning (NSP) on a functional connectivity matrix. NSP is an
%algorithm that detects hierarchical clusters in FC, correponding to
%whole-brain activation, left- and right-hemisphere, anterior-posterior
%regions, and so forth. This function is adapted from the scripts in Wang
%et al. (2021).
%
% Key reference(s):
% Wang, R., Liu, M., Cheng, X., Wu, Y., Hildebrandt, A., & Zhou, C. (2021).
% Segregation, integration, and balance of large-scale resting brain
% networks configure different cognitive abilities. Proceedings of the
% National Academy of Sciences, 118(23), e2022288118.
%
%   [fc_rearranged, cluster_size] = fc_nested_network_partitioning(fc)
%
% Input:
%  fc:  N x N symmetric functional connectivity matrix
% Output:
%  fc_rearranged:   N x N FC matrix with rearranged rows and columns
%  cluster_size:    N x 1 cell; the k-th element is an array containing the
%                   numbers of regions in each cluster for the k-th mode.
%                   For example, the first array is N, corresponding to
%                   global activation of the whole brain; the second array
%                   is typically [N/2, N/2], corresponding to left- and
%                   right-hemisphere, and the first N/2 and the second N/2
%                   rows and columns of fc_rearranged are assigned to the
%                   two clusters; higher-order modes contain finer division
%                   of N, corresponding to smaller clusters of brain
%                   regions.
%  cluster_number:  "wavenumber" of brain regions for each mode; (1 /
%                   cluster_number) is the rough number of regions in one
%                   cluster

% Exclude negative correlations
fc(fc < 0) = 0;

% Eigen-decomposition of the raw FC
[eigvecs, ~] = eig(fc);
eigvecs = fliplr(eigvecs);    % start from the largest eigenvalue

% Initial partition
H1_1 = find(eigvecs(:,1) < 0);
H1_2 = find(eigvecs(:,1) >= 0);

%==============================================================
% Main Loop starts here

disp('Perform nested spectral decomposition (NSP) on FC matrix.')

N = size(fc,1);
cluster_number = [1];    % The first level has one module and corresponds to the global integration.  
cluster_size = cell(N,1);

% Loop over eigenmodes
for mode = 2:N
    x = find(eigvecs(:,mode) >= 0);
    y = find(eigvecs(:,mode) < 0);
    H = {};
    for j = 1:2*cluster_number(mode-1)
        % assume the number of cluster in j-1 level is 2^(mode-1)
        H{j} = eval(['H',num2str(mode-1),'_',num2str(j)]);
    end

    % length of each cluster in H
    id = cellfun('length',H);

    % delete the cluster with 0 size
    H(id==0) = [];
    id(id==0) = [];

    % save the cluster size
    cluster_size{mode-1} = id;
    cluster_number = [cluster_number, length(H)];    % number of cluster
    
    k = 1; 
    for j = 1:2:2*cluster_number(mode)    % modular number
         Positive_Node = intersect(H{k},x);
         Negtive_Node = intersect(H{k},y);         
         k = k+1;
         eval(['H',num2str(mode),'_',num2str(j+1), '=', 'Positive_Node', ';'])
         eval(['H',num2str(mode),'_',num2str(j), '=', 'Negtive_Node', ';'])
    end  
    for j = 1:2*cluster_number(mode-1)
         eval(['clear',' H',num2str(mode-1),'_',num2str(j),'']);
    end

    Z = [];    % indices for rearranging the rows and columns
    if (cluster_number(end)==N)
        for j = 1:2*cluster_number(mode)
            Z = [Z; eval(['H',num2str(mode),'_',num2str(j)])];
        end
        break;
    end
end
cluster_number(1) = [];
cluster_number = [cluster_number/N, ones(1,N-length(cluster_number))];

% Rearrange FC
fc_rearranged = zeros(N,N);
for i=1:N
    for j=1:N
        fc_rearranged(i,j) = fc(Z(i),Z(j));
    end
end

end