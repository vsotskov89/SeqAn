%% clear the workspace , turn off warnings
clear; clc; close all;
warning('off', 'all'); 

%addpath(genpath('/Users/ventalon/ownCloud/Documents/Matlab/Simulations/cnmfe/'))

%% choose data
%nam = 'Substack_corr-02_with_mask2.tif';  
%nam = 'Substack_corr-02.tif';  
%nam = 'Substack_corr-03_with_mask_500frames_crop.tif'; 
nam = 'simulated-data-clean-preprocess.tif'; 
%nam = 'Substack_corr-02_LS.tif'; 
neuron = Sources2D();
nam = neuron.select_data(nam);  %if nam is [], then select data interactively


%% parameters

% -------------------------    ACTIONS TO MAKE    -------------------------  %
actions = struct;
actions.CreateAndSaveMask = false; % to create a mask corresponding to the bundle, save data with the mask, and use the mask to avoid problems close to the bundle limit. To use if the fiberscope data has no mask yet.  
actions.CreateMask = false; % to extract mask corresponding to the bundle (from the image data), and then use it in CNMFE (apply it on Cn and PNR matrix to avoid problems close to the bundle limit). Mandatory for fiberscope data. 
actions.PickNeuronsFromResidual = false; 
actions.LoopsToRun = 1 ;  % number of optimization loops to run after initialisation. Can be 0. 


% -------------------------    COMPUTATION    -------------------------  %
pars_envs = struct('memory_size_to_use', 26, ...   % GB, memory space you allow to use in MATLAB
    'memory_size_per_patch', 2, ...   % GB, space for loading data within one patch
    'patch_dims', [128, 128]);  %GB, patch size

% distribute data and be ready to run source extraction
neuron.getReady(pars_envs); % I put this here now to get the number of frames


% -------------------------      SPATIAL      -------------------------  %
gSig = 8;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
gSiz = 21;          % pixel, neuron diameter
ssub = 1;           % spatial downsampling factor
with_dendrites = false;   % with dendrites or not
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    updateA_search_method = 'dilate';  %#ok<UNRCH>
    updateA_bSiz = 5;
    updateA_dist = neuron.options.dist;
else
    % determine the search locations by selecting a round area
    updateA_search_method = 'ellipse'; % #ok<UNRCH>
    updateA_bSiz = 5;
    updateA_dist = neuron.options.dist;
end
spatial_constraints = struct('connected', true, 'circular', true);  % you can include following constraints: 'circular'
spatial_algorithm = 'hals_thresh';


% -------------------------      TEMPORAL     -------------------------  %
Fs = 100;             % frame rate (Hz)
tsub = 1;           % temporal downsampling factor
deconv_flag = true;     % run deconvolution or not 
window_baseline = 10 * Fs; % 10 secondes
deconv_method = 'bsd'; % method for running deconvolution {'foopsi', 'constrained', 'thresholded', 'bsd'}

% general parameters for deconvolution algo
tauRise = 0.03; % en s
tauDecay = 0.1; % en s
smin = 1; % definition 
    % 1 - pour les methodes habituelles de cnmfe ('foopsi', 'constrained',
    % 'thresholded') : smin defini par les developpeurs. "When the value is
    % negative, the actual threshold is abs(smin)*noise level"
    % 2 - pour bsd
    % si smin = 0 : aucun seuil n'est mis
    % si smin > 0 : un seuil egal à smin * bruit est choisi
    % si smin < 0 : un seuil egal à smin * Pphys.threshold est choisi

% parameters for AR2
Tar2 = struct;
Tar2.tau_r = tauRise * Fs; % en nb de frames
Tar2.tau_d = tauDecay * Fs;

% parameters for bsd
Obsd = struct; % Struct of experimental conditions & decoding options.
dims = neuron.P.mat_data.dims;
Obsd.Time = dims(3); % Number of time frames. 
Obsd.nNeurons = 1; % Number of neurons.
Obsd.dt = 1/Fs; % interval duration. (s)
Obsd.adaptive = 0; % Not adaptive. Will use provided values for parameters, and estimate the unknown ones.

Pbsd = struct; % Struct of generative model properties.
Pbsd.tauRise = tauRise; % Fluorescence rise time (s)
Pbsd.tauDecay = tauDecay; % Fluorescence decay time (s)
Pbsd.th = smin; % Threshold on spike amplitude (see above and in deconvolveCa.m)
Pbsd.b = 0; % Baseline position

deconv_options = struct('type', 'ar2', ... % model of the calcium traces. {'ar1', 'ar2'}, not used in bsd
    'method', deconv_method, ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded', 'bsd'}
    'Tar2', Tar2,... % tau_d et tau_r for AR2
    'smin', smin, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level, not used in bsd
    'optimize_pars', false, ...  % optimize AR coefficients
    'optimize_b', true, ...% optimize the baseline
    'max_tau', 15, ... % maximum decay time (unit: frame);
    'Obsd', Obsd, ... % parameter for bsd algo, see above
    'Pbsd', Pbsd, ... % parameter for bsd algo, see above
    'window_baseline', window_baseline);  % size of time window to compute the baseline (default: 1000)    

nk = 1;             % detrending the slow fluctuation. usually 1 is fine (no detrending)
% when changed, try some integers smaller than total_frame/(Fs*30)
detrend_method = 'spline';  % compute the local minimum as an estimation of trend.


% -------------------------     BACKGROUND    -------------------------  %
bg_model = 'ring';  % model of the background {'ring', 'svd'(default), 'nmf'}
nb = 1;             % number of background sources for each patch (only be used in SVD and NMF model)
ring_radius = 30;  % when the ring model used, it is the radius of the ring used in the background model.
%otherwise, it's just the width of the overlapping area
num_neighbors = []; % number of neighbors for each neuron
bg_ssub = 2;        % downsample background for a faster speed 


% -------------------------      MERGING      -------------------------  %
show_merge = false;  % if true, manually verify the merging step
merge_thr = 0.65;     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
method_dist = 'max';   % method for computing neuron distances {'mean', 'max'}
dmin = 5;       % minimum distances between two neurons. it is used together with merge_thr
dmin_only = 2;  % merge neurons if their distances are smaller than dmin_only.
merge_thr_spatial = [0.8, 0.4, -inf];  % merge components with highly correlated spatial shapes (corr=0.8) and small temporal correlations (corr=0.1)


% -------------------------  INITIALIZATION   -------------------------  %
K = 150;             % maximum number of neurons per patch. when K=[], take as many as possible.
min_corr = 0.6;     % minimum local correlation for a seeding pixel
min_pnr = 6;       % minimum peak-to-noise ratio for a seeding pixel
min_pixel = gSig^2;      % minimum number of nonzero pixels for each neuron
bd = 0;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
frame_range = [];   % when [], uses all frames
save_initialization = false;    % save the initialization procedure as a video.
use_parallel = true;    % use parallel computation for parallel computing
show_init = true;   % show initialization results
choose_params = false; % manually choose parameters
center_psf = true;  % set the value as true when the background fluctuation is large (usually 1p data)
                    % set the value as false when the background fluctuation is small (2p)


% -------------------------  Residual   -------------------------  %
min_corr_res = 0.7;
min_pnr_res = 6;
seed_method_res = 'auto';  % method for initializing neurons from the residual
update_sn = true;


% ----------------------  WITH MANUAL INTERVENTION  --------------------  %
with_manual_intervention = false;


% -------------------------  FINAL RESULTS   -------------------------  %
save_demixed = true;    % save the demixed file or not
kt = 1;                 % frame intervals


%%

% -------------------------    UPDATE PARAMETERS    -------------------------  %
neuron.updateParams('gSig', gSig, ...       % -------- spatial --------
    'gSiz', gSiz, ...
    'ring_radius', ring_radius, ...
    'ssub', ssub, ...
    'search_method', updateA_search_method, ...
    'bSiz', updateA_bSiz, ...
    'dist', updateA_bSiz, ...
    'spatial_constraints', spatial_constraints, ...
    'spatial_algorithm', spatial_algorithm, ...
    'tsub', tsub, ...                       % -------- temporal --------
    'deconv_flag', deconv_flag, ...
    'deconv_options', deconv_options, ...
    'nk', nk, ...
    'detrend_method', detrend_method, ...
    'background_model', bg_model, ...       % -------- background --------
    'nb', nb, ...
    'ring_radius', ring_radius, ...
    'num_neighbors', num_neighbors, ...
    'bg_ssub', bg_ssub, ...
    'merge_thr', merge_thr, ...             % -------- merging ---------
    'dmin', dmin, ...
    'method_dist', method_dist, ...
    'min_corr', min_corr, ...               % ----- initialization -----
    'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, ...
    'bd', bd, ...
    'center_psf', center_psf);
neuron.Fs = Fs;


%% preliminary step, create / save mask for fibersope data. 

if actions.CreateAndSaveMask 
    %if data is from the fiberscope and no mask has been applied to the data yet. 
    % This step creates a mask multiply the data with is and save new data. 

    A = bigread2(nam); 
    % on fait un masque pour enlever les zones hors bundle
    disp('choose fiber bundle limit')
    imagesc(min(A,[],3)); % on plotte le min pour pouvoir eliminer le bundle partout au cours du mouvement
    axis square; 
    circ = drawcircle('Center',[floor(size(A,1)/2), floor(size(A,2)/2) ], 'Radius', floor(size(A,2)/2) );
    pause;
    mask = createMask(circ);
    mask = mask';
    neuron.mask = mask;
    mask = uint16(mask);
    A = A.*mask;
    A = A(floor(circ.Center(1)-circ.Radius):ceil(circ.Center(1)+circ.Radius), floor(circ.Center(2)-circ.Radius):ceil(circ.Center(2)+circ.Radius),:);
    
    tiff_optn = struct; % structure containing tiff options
    tiff_optn.color = false;
    tiff_optn.compress = 'no'; % 'no', 'lzw' (lossless), 'jpeg' (lossy crap), 'adobe' (lossless)
    tiff_optn.message = false; % no output if false
    tiff_optn.append = false;
    tiff_optn.overwrite = true;
    tiff_optn.big = false; % set to true if output > 4GB

    fprintf('Saving stack... ');
    %saveastiff(uint16(Submovie(:,:,1:3000)), 'downsampled_data.tif', tiff_optn);
    %saveastiff(uint16(A), strcat(nam, num2str(kstacks,'%.2d'), '_with_mask.tif'), tiff_optn);
    new_nam = strcat(erase(nam,".tif"), '_with_mask2.tif');
    saveastiff(uint16(A), new_nam, tiff_optn);
    nam = new_nam;
    nam = neuron.select_data(nam); 
    fprintf('Done');
    
    clear A; 
    clear mask;
    clear circle;
end

if actions.CreateMask 
    % to use on fiberscope data with a mask already applied to it. 
    % automatically extract this mask to be used in CNMF-E. 
    
    % methode : on loade les 100 premières images du stack, on prend le max sur l'axe z, et on binarise cette image pour en faire un masque
    % : tous les pixels >0 deviennent egaux a 1. (avec une seule image il
    % restait des pixels egaux a 0 au sein du bundle)

    mask = bigread2(nam, 1, 100);
    mask = (max(mask,[],3)) >0;
    neuron.mask = mask;
    %imagesc(mask)
end


%% initialize neurons from the video data within a selected temporal range
if choose_params
    % change parameters for optimized initialization
    [gSig, gSiz, ring_radius, min_corr, min_pnr] = neuron.set_parameters();
end

[center, Cn, PNR] = neuron.initComponents_parallel(K, frame_range, save_initialization, use_parallel);
neuron.compactSpatial();
if show_init
    figure();
    ax_init= axes();
    imagesc(Cn, [0, 1]); colormap gray;
    hold on;
    plot(center(:, 2), center(:, 1), '.r', 'markersize', 10);
end


% estimate the background components
neuron.update_background_parallel(use_parallel);
neuron_init = neuron.copy();

%  merge neurons and update spatial/temporal components
neuron.merge_neurons_dist_corr(show_merge);
neuron.merge_high_corr(show_merge, merge_thr_spatial);


%% Add loops to improve results. 

remaining_loops = actions.LoopsToRun;

while remaining_loops > 0
    
    remaining_loops = remaining_loops -1;

    % pick neurons from the residual

    if actions.PickNeuronsFromResidual 

        [center_res, Cn_res, PNR_res] =neuron.initComponents_residual_parallel([], save_initialization, use_parallel, min_corr_res, min_pnr_res, seed_method_res);
        if show_init
            axes(ax_init);
            plot(center_res(:, 2), center_res(:, 1), '.g', 'markersize', 10);
        end
        neuron_init_res = neuron.copy();

    end

    % udpate spatial&temporal components, delete false positives and merge neurons
    % update spatial
    if update_sn
        neuron.update_spatial_parallel(use_parallel, true);
        udpate_sn = false;
    else
        neuron.update_spatial_parallel(use_parallel);
    end
    
    % merge neurons based on correlations 
    neuron.merge_high_corr(show_merge, merge_thr_spatial);

    % update temporal
        neuron.update_temporal_parallel(use_parallel);

    % delete bad neurons
        neuron.remove_false_positives();

    % merge neurons based on temporal correlation + distances 
        neuron.merge_neurons_dist_corr(show_merge);
    
     % update background
        neuron.update_background_parallel(use_parallel);


end


%% save and finalise

% save the workspace for future analysis
neuron.orderROIs('snr');
cnmfe_path = neuron.save_workspace();

% show neuron contours
Coor = neuron.show_contours(0.6);

% save neurons shapes
neuron.save_neurons();

% save results
results = neuron.obj2struct(); 
save  /Users/ventalon/ownCloud/Documents/Analyse/2019-07-Simulations_et_PSF/WorkingFolder/results/simulated-data-clean-preprocess_results-2.mat results
















