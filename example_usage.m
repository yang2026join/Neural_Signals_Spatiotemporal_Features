
%% Set up and load data

clear;  clc;

add_math_dependencies;
load("data\example_data.mat")

signal = data_parcel;
%time = time;
fs = (length(time) - 1) / (time(end) - time(1));

n_reg = size(signal,1);

% Some information about the data:
% a resting-state fMRI session data from the Human Connectome Project (HCP)
% The data were obtained from HCP preprocessing pipeline output.
% Subject ID: 100307; session ID: 1; phase encoding: left-to-right;
% Parcellation: aggregate from fsLR gray ordinate to 400 cortical regions
% in Schaefer et al. (2018) aligned to the 7 networks in Yeo et al. (2011).
% Bandpass filtering with cut-off frequency 0.008-0.1 Hz was performed
% before aggregating to regions. 

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%-------------------- Power Spectrum
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
%cfg.taper = 'hann';
cfg.win = 86;
cfg.step = 43;

[psd, freq] = power_spectrum_winavg(signal, time, fs, cfg);

figure;
plot(freq, psd, '.-k');
xlabel('Frequency (Hz)');  ylabel('Spectral density');
title('Power Spectrums');

%% PSD-derived features

%--------------------------------------------------------------------------
%-------------------- mean low-frequency psd
%--------------------------------------------------------------------------

cfg = [];
cfg.freq_range = [0.01 0.05];

mean_low_freq_psd = psd_band_mean_psd(psd, freq, cfg);
% Alternative: amplitude of low-frequency fluctuations (ALFF)
%ampl_low_freq_fluct = psd_ampl_low_freq_fluct(psd, freq, cfg);


%--------------------------------------------------------------------------
%-------------------- fraction of low-frequency power
%--------------------------------------------------------------------------

cfg = [];
cfg.foi = [0.01 0.05];
cfg.fwhole = [0.01 0.09];

frac_low_freq_power = psd_band_power_prop(psd, freq, cfg);
% Alternative: fraction of ALFF
%frac_ampl_low_freq_fluct = psd_frac_ampl_low_freq_fluct(psd, freq, cfg);


%--------------------------------------------------------------------------
%-------------------- power-law exponents
%--------------------------------------------------------------------------

cfg = [];
cfg.freq_range = [0.02 0.1];

psd_exp = psd_power_law_exponent(psd, freq, cfg);


%--------------------------------------------------------------------------
%-------------------- PC scores of spectrums
%--------------------------------------------------------------------------

cfg = [];
cfg.n_comps = 3;
%cfg.cum_var = 95;
cfg.freq_range = [0.02 0.1];

psd_pc_scores = psd_principal_component_scores(psd, freq, cfg);

psd_pc1_scores = psd_pc_scores(:,1);
psd_pc2_scores = psd_pc_scores(:,2);
psd_pc3_scores = psd_pc_scores(:,3);

%--------------------------------------------------------------------------
%-------------------- Plot
%--------------------------------------------------------------------------

varNames = {
    'mean_low_freq_psd'
    'frac_low_freq_power'
    'psd_exp'
    'psd_pc1_scores'
    'psd_pc2_scores'
    'psd_pc3_scores'
    };
n_var = length(varNames);

% value vs region
figure;
for i = 1:n_var
    nexttile;
    var = eval(varNames{i});
    plot(1:n_reg, var, '.-');
    title(varNames{i}, 'Interpreter', 'none')
    xlabel('Region ID')
end

% distribution
figure;
for i = 1:n_var
    nexttile;
    var = eval(varNames{i});
    histogram(var, 30);
    title(varNames{i}, 'Interpreter', 'none')
end

% similarity
sim_mat = zeros(n_var);
for i = 1:n_var
    var1 = eval(varNames{i});
    for j = 1:n_var
        var2 = eval(varNames{j});
        sim_mat(i, j) = corr(var1, var2);
    end
end
figure;
imagesc(sim_mat, [0 1]);
colorbar;
axis square;
set(gca, 'XTick', 1:n_var, 'XTickLabel', varNames, ...
         'YTick', 1:n_var, 'YTickLabel', varNames, ...
         'TickLabelInterpreter', 'none');
title('Pairwise similarity (Pearson correlation) between variables');
[rows, cols] = ndgrid(1:n_var, 1:n_var);  rows = rows(:) - 0.4;  cols = cols(:);
text(rows, cols, num2str(sim_mat(:), 3), 'Color', 'k')


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%---------- Functional connectivity
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.method = 'Pearson';
cfg.win = 21;
cfg.step = 21;

fc = functional_connectivity_winavg(signal, time, fs, cfg);

figure;
imagesc(fc);  colorbar; axis square
xlabel('Region ID');  ylabel('Region ID');
title('Functional Connectivity Matrix')

%% FC-derived features

%--------------------------------------------------------------------------
%-------------------- eigenvectors
%--------------------------------------------------------------------------

cfg = [];
cfg.n_modes = 3;

[fc_eigvecs, ~] = fc_eigendecomp(fc, cfg);

fc_eigvec1 = fc_eigvecs(:,1);
fc_eigvec2 = fc_eigvecs(:,2);
fc_eigvec3 = fc_eigvecs(:,3);


%--------------------------------------------------------------------------
%-------------------- spectral embedding
%--------------------------------------------------------------------------

cfg = [];
cfg.n_comps = 3;
%cfg.thres = 5;

fc_grad_spectral = fc_spectral_embedding(fc, cfg);

fc_grad_spectral1 = fc_grad_spectral(:,1);
fc_grad_spectral2 = fc_grad_spectral(:,2);
fc_grad_spectral3 = fc_grad_spectral(:,3);


%--------------------------------------------------------------------------
%-------------------- diffusion embedding
%--------------------------------------------------------------------------

cfg = [];
cfg.n_comps = 3;
%cfg.thres = 5;
cfg.diff_time = 5;

fc_grad_diff = fc_diffusion_embedding(fc, cfg);

fc_grad_diff1 = fc_grad_diff(:,1);
fc_grad_diff2 = fc_grad_diff(:,2);
fc_grad_diff3 = fc_grad_diff(:,3);


%--------------------------------------------------------------------------
%-------------------- Plot
%--------------------------------------------------------------------------

varNames = {
    'fc_eigvec1'
    'fc_eigvec2'
    'fc_eigvec3'
    'fc_grad_spectral1'
    'fc_grad_spectral2'
    'fc_grad_spectral3'
    'fc_grad_diff1'
    'fc_grad_diff2'
    'fc_grad_diff3'
    };
n_var = length(varNames);

% value vs region
figure;
for i = 1:n_var
    nexttile;
    var = eval(varNames{i});
    plot(1:n_reg, var, '.-');
    title(varNames{i}, 'Interpreter', 'none')
    xlabel('Region ID')
end

% distribution
figure;
for i = 1:n_var
    nexttile;
    var = eval(varNames{i});
    histogram(var, 30);
    title(varNames{i}, 'Interpreter', 'none')
end

% similarity
sim_mat = zeros(n_var);
for i = 1:n_var
    var1 = eval(varNames{i});
    for j = 1:n_var
        var2 = eval(varNames{j});
        sim_mat(i, j) = corr(var1, var2);
    end
end
figure;
imagesc(sim_mat, [0 1]);
colorbar;
axis square;
set(gca, 'XTick', 1:n_var, 'XTickLabel', varNames, ...
         'YTick', 1:n_var, 'YTickLabel', varNames, ...
         'TickLabelInterpreter', 'none');
title('Pairwise similarity (Pearson correlation) between variables');
[rows, cols] = ndgrid(1:n_var, 1:n_var);  rows = rows(:) - 0.4;  cols = cols(:);
text(rows, cols, num2str(sim_mat(:), 3), 'Color', 'k')
