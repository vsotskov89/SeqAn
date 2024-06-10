 %% clear the workspace , turn off warnings
clear; clc; close all;
warning('off', 'all'); 

%% Parameters

% ------------------------------------- Files ---------------------------------------------------

usingDeepCAD = false;
only10Hz = false;

ses_name = 'sleeppost';

Files = struct;
Files.folder = '/export/home1/RawCalciumData/';
Files.subfolderResults = [Files.folder, 'CNMF_matlab_results/']; % generated files are saved in this subfolder;
Files.vessels = [Files.subfolderResults, ses_name, '_downsampled_vessels.tif' ]; % vessel cropped image that will be saved in the results structure;
Files.subfolderData = [Files.folder, ses_name, '/']; % subfolder containing experimental data files;
Files.basename = 'Substack_cropped-';
Files.matlab_file = '/export/home1/MatlabFiles/Codes_Vova/Pipeline_CNMF_Vova.m'; % to keep track of parameters
Files.results = [ses_name, '_results_high_pnr.mat']; % output files

tiff_optn = struct; % structure containing tiff options
tiff_optn.color = false;
tiff_optn.compress = 'no'; % 'no', 'lzw' (lossless), 'jpeg' (lossy crap), 'adobe' (lossless)
tiff_optn.message = false; % no output if false
tiff_optn.append = false;
tiff_optn.overwrite = true;
tiff_optn.big = false; % set to true if output > 4GB


% -------------------------------- Parameters for CNMFE ---------------------------------------------

% nb : some of these parameters are not used in the second execution as we are only
% running the function update temporal

% ACTIONS TO MAKE  
actions = struct;
actions.CreateMask = true; % to extract mask corresponding to the bundle (from the image data), and then use it in CNMFE (apply it on Cn and PNR matrix to avoid problems close to the bundle limit). Mandatory for fiberscope data. 
actions.PickNeuronsFromResidual = false; 
actions.LoopsToRun = 1 ;  % number of optimization loops to run after initialisation. Can be 0. 
with_manual_intervention = false;


% COMPUTATION   
pars_envs = struct('memory_size_to_use', 120, ...   % GB, memory space you allow to use in MATLAB
    'memory_size_per_patch', 40, ...   % GB, space for loading data within one patch
    'patch_dims', [128, 128]);  %GB, patch size



% SPATIAL     
gSiz = 17;          % pixel, neuron diameter
gSig_fact = 0.28; 
gSig = gSig_fact * gSiz; % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
ssub = 1;           % spatial downsampling factor
min_pixel_fact = 2.25;
min_pixel = min_pixel_fact* gSig^2;      % minimum number of nonzero pixels for each neuron
with_dendrites = false;   % with dendrites or not
updateA_bSiz = 5;
spatial_constraints = struct('connected', true, 'circular', true);  % you can include following constraints: 'circular'
spatial_algorithm = 'hals_thresh';


% TEMPORAL
tsub = 1;           % temporal downsampling factor
deconv_flag = true;     % run deconvolution or not 
deconv_method = 'bsd'; % method for running deconvolution {'foopsi', 'constrained', 'thresholded', 'bsd'}
%GCaMP6f
%tauRise = 0.07; % en s
%tauDecay = 0.8; 
%GCaMP8m
tauRise = 0.005; %0.07; %0.045; % en s
tauDecay = 0.1; %0.8;% !!! changer � 0.8 apres !!!  %s0.3; %0.142; % en s
Fs_fast = 100; % 100;             % acquisition frame rate (Hz) !!!!!!!!!
N_downsample = 15; % downsample factor for first run of CNMFE 
    % 20 pour une manip totale sommeil veille sommeil; 10 pour une manip de 30min
    % test� N_downsample = 15 : out of memory pour CNMFE (downsample
    % stack OK, poids 6.5Go) 
time_window = 30; % seconds, duration of time window for baseline subtraction
smin_slow = 0.5; 
smin_fast = 0.5; 
    % definition 
    % 1 - pour les methodes habituelles de cnmfe ('foopsi', 'constrained',
    % 'thresholded') : smin defini par les developpeurs. "When the value is
    % negative, the actual threshold is abs(smin)*noise level"
    % 2 - pour bsd
    % si smin = 0 : aucun seuil n'est mis
    % si smin > 0 : un seuil egal � smin * bruit est choisi
    % si smin < 0 : un seuil egal � smin * Pphys.threshold est choisi
nk = 1;             % detrending the slow fluctuation. usually 1 is fine (no detrending)
% when changed, try some integers smaller than total_frame/(Fs*30)
detrend_method = 'spline';  % compute the local minimum as an estimation of trend.


%  BACKGROUND 
bg_model = 'ring';  % model of the background {'ring', 'svd'(default), 'nmf'}
nb = 1;             % number of background sources for each patch (only be used in SVD and NMF model)
ring_radius_fact = 1.4;
ring_radius = ring_radius_fact * gSiz;  % when the ring model used, it is the radius of the ring used in the background model.
%otherwise, it's just the width of the overlapping area
num_neighbors = []; % number of neighbors for each neuron
bg_ssub = 2;        % downsample background for a faster speed 


%   MERGING   
show_merge = false;  % if true, manually verify the merging step
merge_thr = 0.5;     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
method_dist = 'max';   % method for computing neuron distances {'mean', 'max'}
dmin = 5;       % minimum distances between two neurons. it is used together with merge_thr
dmin_only = 2;  % merge neurons if their distances are smaller than dmin_only.
merge_thr_spatial = [0.6, 0.5, -inf];  % merge components with highly correlated spatial shapes (corr=0.8) and small temporal correlations (corr=0.1)


%  INITIALIZATION 
K = 300;             % maximum number of neurons per patch. when K=[], take as many as possible.
%K = 10; 
min_corr = 0.7;    % minimum local correlation for a seeding pixel
min_pnr = 55;       % minimum peak-to-noise ratio for a seeding pixel
bd = 0;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
frame_range = [];   % when [], uses all frames
save_initialization = false;    % save the initialization procedure as a video.
use_parallel = true;    % use parallel computation for parallel computing
show_init = true;   % show initialization results
center_psf = true;  % set the value as true when the background fluctuation is large (usually 1p data)
                    % set the value as false when the background fluctuation is small (2p)

%  RESIDUAL  
min_corr_res = 0.5;
min_pnr_res = 5;
seed_method_res = 'auto';  % method for initializing neurons from the residual
update_sn = true;


%  FINAL RESULTS  
save_demixed = true;    % save the demixed file or not
kt = 1;                 % frame intervals
 
%% Create output folder

if ~exist(Files.subfolderResults, 'dir') 
            mkdir(Files.subfolderResults);
else
     fprintf('Warning! Subfolder Results already exists! Data will be erased. Press any key to continue');
     pause; 
end
            
%% copy matlab file

copyfile(Files.matlab_file,Files.subfolderResults); % to keep track of parameters


%% Downsample data

    cd(Files.subfolderData);
    nStacks = length(dir('*.tif'));
    cd(Files.subfolderResults);
    tStacks = zeros(1, nStacks); % number of time points in each stack

    if nStacks >1
        kStacks = 2; % we initialize the number of elements using the length of the second stack because the first one is shorter 
    else
        kStacks = 1;
    end

    nam = [Files.subfolderData Files.basename num2str(kStacks,'%.2d') '.tif'];
    A = bigread2(nam);
    Nx = size(A,1);
    Ny = size(A,2);
    Nt = size(A,3);
    Submovie = zeros(Nx, Ny, floor(Nt*(nStacks)/N_downsample));

    reliquat = [];
    kt0 = 1;

%     for kStacks = 1 : nStacks
%         kStacks
%        nam = [Files.subfolderData Files.basename num2str(kStacks,'%.2d') '.tif'];
%        A = bigread2(nam);
%        tStacks(kStacks) = size(A, 3); % number of frames in the stack;
%        A = cat(3, reliquat, A);
%        ntk = size(A, 3); % number of frames in the stack;
%        if  mod(ntk,N_downsample) == 0 
%            reliquat = [];
%        else
%            reliquat = A(:,:, end - mod(ntk,N_downsample) +1 : end);
%            A = A(:,:, 1 : end - mod(ntk,N_downsample) );
%        end
%        Ntksub = floor(ntk/N_downsample);
%        Submovie(:,:, kt0 : kt0 + Ntksub -1) =  imresize3(A, [size(A, 1) size(A,2) Ntksub]);
%        kt0 = kt0 + Ntksub;
%     end
%     Submovie = Submovie(:,:, 1 : kt0 -1 );

    % cas ou le laser bleu a des defaillances : on elimine les frames
    % correspondantes et on note le numero de ces frames. le code suivant
    % remplace le code precedent. 
    frameCount = 0;
    removedFrames = [];
    
    % a elminer !!!
%     for kStacks = 2 : nStacks
%         kStacks
%        nam = [Files.subfolderData Files.basename num2str(kStacks,'%.2d') '.tif'];
%        A = bigread2(nam);
%        B =squeeze(sum(A,[1,2]));
%        C =find(B<0.8*median(B));
%        if numel(C)>1
%            A(:,:, C) = [];
%            copyfile(nam,[Files.subfolderData Files.basename num2str(kStacks,'%.2d') '_old.tif']);
%            saveastiff(A,nam,tiff_optn); 
%        end
%     end
        
    for kStacks = 1 : nStacks
       nam = [Files.subfolderData Files.basename num2str(kStacks,'%.2d') '.tif'];
       A = bigread2(nam);
       B =squeeze(sum(A,[1,2]));
       C =find(B<0.8*median(B));
       if numel(C)>1
           A(:,:, C) = [];
           removedFrames = [removedFrames ; C + frameCount ];
           copyfile(nam,[Files.subfolderData Files.basename num2str(kStacks,'%.2d') '_old.tif']);
           saveastiff(A,nam,tiff_optn); 
       end
       tStacks(kStacks) = size(A, 3); % number of frames in the stack;
       A = cat(3, reliquat, A);
       ntk = size(A, 3); % number of frames in the stack;
       if  mod(ntk,N_downsample) == 0 
           reliquat = [];
       else
           reliquat = A(:,:, end - mod(ntk,N_downsample) +1 : end);
           A = A(:,:, 1 : end - mod(ntk,N_downsample) );
       end
       Ntksub = floor(ntk/N_downsample);
       Submovie(:,:, kt0 : kt0 + Ntksub -1) =  imresize3(A, [size(A, 1) size(A,2) Ntksub]);
       kt0 = kt0 + Ntksub;
       frameCount = frameCount + tStacks(kStacks);
    end
    Submovie = Submovie(:,:, 1 : kt0 -1 );
    save('removedFrames.mat','removedFrames');
    
    
    if  usingDeepCAD == true
      %  minImage = imread(minImageFile);
      minImage = min(Submovie,[],3);
        figure; imagesc(minImage);
        BW = minImage > 200 ; %BW = imbinarize(minImage); figure; imagesc(BW);
            s = regionprops(BW,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
            if s(1).Orientation < 0
                s(1).Orientation=s(1).Orientation+360;
            end
             h = drawellipse('Center',s(1).Centroid ,'SemiAxes',[s(1).MajorAxisLength/2 s(1).MinorAxisLength/2], 'RotationAngle',s(1).Orientation,'StripeColor','m');

            %modification de l'ellipse a la main quand �a ne marche pas du tout.
%             h.Center = size(minImage')/2;% + [2 2];
%             h.SemiAxes = size(minImage')/2 ;%- [2 4];
%             h.RotationAngle = 0;

            mask = uint16(createMask(h)); figure; imagesc(mask);
            Submovie = uint16(Submovie).*mask;
    end
    
    
    fprintf('Saving downsampled stack... \n');

    nam = [Files.subfolderResults 'downsampled_movie' '.tif'];
    saveastiff(uint16(Submovie), nam , tiff_optn);
    clear Submovie A;
    save([Files.subfolderResults ,'tStacks.mat'],'tStacks')
    fprintf('Done');
    
%pour faire tourner sans calculer le downsampled movie, executer :
%     nam = [Files.subfolderResults 'downsampled_movie' '.tif'];
%  load([Files.subfolderResults 'tStacks.mat']);
% cd(Files.subfolderData);
%    nStacks = length(dir('*.tif'));
% 
%% First execution of CNMFE (10Hz)

% create neuron object
neuron = Sources2D();
neuron.file = nam;
%nam = neuron.select_data(nam); %if nam is [], then select data interactively, but it won't compute pnr
% distribute data and be ready to run source extraction
neuron.getReady(pars_envs); 

% compute missing parameters

Fs_slow = Fs_fast / N_downsample;             % frame rate (Hz) of downsampled data
window_baseline_slow = time_window * Fs_slow; 
updateA_dist = neuron.options.dist; % je sais pas trop a quoi ca sert
Tar2 = struct;
Tar2.tau_r = tauRise * Fs_slow; % en nb de frames
Tar2.tau_d = tauDecay * Fs_slow;

% parameters for bsd
Obsd = struct; % Struct of experimental conditions & decoding options.
dims = neuron.P.mat_data.dims;
Obsd.Time = dims(3); % Number of time frames. 
Obsd.nNeurons = 1; % Number of neurons.
Obsd.dt = 1/Fs_slow; % interval duration. (s)
Obsd.adaptive = 0; % Not adaptive. Will use provided values for parameters, and estimate the unknown ones.

Pbsd = struct; % Struct of generative model properties.
Pbsd.tauRise = tauRise; % Fluorescence rise time (s)
Pbsd.tauDecay = tauDecay; % Fluorescence decay time (s)
Pbsd.th = smin_slow; % Threshold on spike amplitude (see above and in deconvolveCa.m)
Pbsd.b = 0; % Baseline position

fprintf('Running CNMFE... ');

neuron.file_b = [Files.subfolderResults, 'logFile.txt']; 
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
	updateA_search_method = 'dilate';  % #ok<UNRCH>
else
     % determine the search locations by selecting a round area
    updateA_search_method = 'ellipse'; % #ok<UNRCH>
end

deconv_options = struct('type', 'ar2', ... % model of the calcium traces. {'ar1', 'ar2'}, not used in bsd
            'method', deconv_method, ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded', 'bsd'}
            'Tar2', Tar2,... % tau_d et tau_r for AR2
            'smin', smin_slow, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level, not used in bsd
            'optimize_pars', false, ...  % optimize AR coefficients
            'optimize_b', true, ...% optimize the baseline
            'max_tau', 20, ... % maximum decay time (unit: frame);
            'Obsd', Obsd, ... % parameter for bsd algo, see above
            'Pbsd', Pbsd, ... % parameter for bsd algo, see above
            'window_baseline', window_baseline_slow);  % size of time window to compute the baseline (default: 1000)
 
neuron.updateParams('gSig', gSig, ...       % -------- spatial --------
            'gSiz', gSiz, ...
            'ring_radius', ring_radius, ...
            'ssub', ssub, ...
            'search_method', updateA_search_method, ...
            'bSiz', updateA_bSiz, ...
            'dist', updateA_dist, ...
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
neuron.Fs = Fs_slow;
moreParms.actions = actions;
moreParms.with_dendrites = with_dendrites;
moreParms.dmin_only = dmin_only;
moreParms.merge_thr_spatial = merge_thr_spatial;
moreParms.min_corr_res = min_corr_res;
moreParms.min_pnr_res = min_pnr_res ;


% preliminary step, create / save mask for fibersope data. 

if actions.CreateMask 
            % to use on fiberscope data with a mask already applied to it. 
            % automatically extract this mask to be used in CNMF-E. 

            % methode : on loade les 100 premi�res images du stack, on prend le max sur l'axe z, et on binarise cette image pour en faire un masque
            % : tous les pixels >0 deviennent egaux a 1. (avec une seule image il
            % restait des pixels egaux a 0 au sein du bundle)

    mask = bigread2(nam, 1, 100);  
    % !!!!! recently changed (21/5/2021) to use on data where no mask has
    % been created. Check if it works....
    mask = (mean(mask,3)) >100;
    % mask = (max(mask,[],3)) >10;
    neuron.mask = mask;
    %imagesc(mask)
end


[center, Cn, PNR] = neuron.initComponents_parallel(K, frame_range, save_initialization, use_parallel, 0);
neuron.compactSpatial();
if show_init
    figure();
    ax_init= axes();
   %   imagesc(Cn, [0.7, 1]); %colormap gray;
    %imagesc(Cn>0.9, [0, 1]); %colormap gray;
   imagesc(PNR, [0, 150]);% colormap gray;
  %  imagesc((PNR>20).*(Cn>0.7)); 
  %  imagesc(Cn.*PNR, [0, 150]); %colormap gray;
     %   imagesc(PNR); 
    hold on;
    plot(center(:, 2), center(:, 1), '.r', 'markersize', 10);
    %figure; imagesc(PNR>15)
end


% estimate the background components
neuron.update_background_parallel(use_parallel);
neuron_init = neuron.copy();

%  merge neurons and update spatial/temporal components
neuron.merge_neurons_dist_corr(show_merge);
neuron.merge_high_corr(show_merge, merge_thr_spatial);


% Add loops to improve results. Cc

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


% si on veut sauver les donn�es � 10Hz
if only10Hz ==  true
  	results = neuron.obj2struct();
    results.moreParms = moreParms;
    fname = [Files.subfolderResults, 'Results_10Hz.mat'];
    save(fname, 'results', '-v7.3');    
   return
end


%% Second execution of CNMFE (100Hz)

% parameters
window_baseline_fast = time_window * Fs_fast; 
Tar2_fast = struct;
Tar2_fast.tau_r = tauRise * Fs_fast; % en nb de frames
Tar2_fast.tau_d = tauDecay * Fs_fast;

 % parameters for bsd
Obsd_fast = struct; % Struct of experimental conditions & decoding options.
Obsd_fast.nNeurons = 1; % Number of neurons.
Obsd_fast.dt = 1/Fs_fast; % interval duration. (s)
Obsd_fast.adaptive = 0; % Not adaptive. Will use provided values for parameters, and estimate the unknown ones.

Pbsd_fast = struct; % Struct of generative model properties.
Pbsd_fast.tauRise = tauRise; % Fluorescence rise time (s)
Pbsd_fast.tauDecay = tauDecay; % Fluorescence decay time (s)
Pbsd_fast.th = smin_fast; % Threshold on spike amplitude (see above and in deconvolveCa.m)
Pbsd_fast.b = 0; % Baseline position
        
 % Resampling temporal traces
Cs = neuron.C;
Cs_raw = neuron.C_raw;
Cs_prev = neuron.C_prev;
Ss = neuron.S;
nCells = size(Cs, 1);
nTime = sum(tStacks);
nTime_cut = size(Cs, 2) * N_downsample;
Cs_fast = imresize(Cs, [nCells, nTime_cut]);
Cs_fast = padarray(Cs_fast,[0 nTime-nTime_cut],'replicate','post');
Cs_raw_fast = imresize(Cs_raw, [nCells, nTime_cut]);
Cs_raw_fast = padarray(Cs_raw_fast,[0 nTime-nTime_cut],'replicate','post');
Cs_prev_fast = imresize(Cs_prev, [nCells, nTime_cut]);
Cs_prev_fast = padarray(Cs_prev_fast,[0 nTime-nTime_cut],'replicate','post');
Ss_fast = imresize(Ss, [nCells, nTime_cut]);
Ss_fast = padarray(Ss_fast,[0 nTime-nTime_cut],'replicate','post');

tStacks2 = [0 tStacks];

% initialise final temporal traces
C_fast = zeros(nCells, nTime) ;
C_raw_fast = zeros(nCells, nTime) ; % not sure that this one is usefull
C_prev_fast = zeros(nCells, nTime) ; % not sure that this one is usefull
S_fast = zeros(nCells, nTime) ;
    

for kStacks = 1 : nStacks

    % create object neuron_fast
    neuron_fast = Sources2D();

    nam = [Files.subfolderData Files.basename num2str(kStacks,'%.2d') '.tif'];
    %nam = neuron_fast.select_data(nam); 
    neuron_fast.file = nam;
    neuron_fast.getReady(pars_envs);
    
    T_fast = tStacks(kStacks);
    Obsd_fast.Time = T_fast; % Number of time frames
    neuron_fast.Fs = Fs_fast;
	neuron_fast.frame_range = [1 T_fast];
	neuron_fast.options = neuron.options;
	deconv_options_fast = struct('type', 'ar2', ... % model of the calcium traces. {'ar1', 'ar2'}, not used in bsd
            'method', deconv_method, ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded', 'bsd'}
            'Tar2', Tar2_fast,... % tau_d et tau_r for AR2
            'smin', smin_fast, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level, not used in bsd
            'optimize_pars', false, ...  % optimize AR coefficients
            'optimize_b', true, ...% optimize the baseline
            'max_tau', 15, ... % maximum decay time (unit: frame);
            'Obsd', Obsd_fast, ... % parameter for bsd algo, see above
            'Pbsd', Pbsd_fast, ... % parameter for bsd algo, see above
            'window_baseline', window_baseline_fast);  % size of time window to compute the baseline (default: 1000)   
    neuron_fast.updateParams('deconv_options', deconv_options_fast);    
    
  % Filling results from execution of CNMFE in downsampled data in neuron_fast. 
    Tmin = 1+sum(tStacks2(1:kStacks));
    Tmax = sum(tStacks2(2:kStacks+1));
    neuron_fast.C = Cs_fast(:,Tmin: Tmax);
    neuron_fast.C_raw = Cs_raw_fast(:,Tmin: Tmax);
    neuron_fast.C_prev = Cs_prev_fast(:,Tmin: Tmax);
    neuron_fast.S = Ss_fast(:,Tmin: Tmax);
    neuron_fast.A = neuron.A;
    neuron_fast.A_prev = neuron.A_prev;
    neuron_fast.W = neuron.W;
    neuron_fast.b0 = neuron.b0;
    neuron_fast.b = neuron.b;
    neuron_fast.f = neuron.f;
    neuron_fast.P.log_file = neuron.P.log_file;
    neuron_fast.P.log_data = neuron.P.log_data;
    neuron_fast.P.Ymean = neuron.P.Ymean;

    % Run update temporal
    neuron_fast.update_temporal_parallel(use_parallel);
    
    % save data
    C_fast(:,Tmin: Tmax) = neuron_fast.C ;
    C_raw_fast(:,Tmin: Tmax) = neuron_fast.C_raw ; % not sure that this one is usefull
    C_prev_fast(:,Tmin: Tmax) = neuron_fast.C_prev ; % not sure that this one is usefull
    S_fast(:,Tmin: Tmax) = neuron_fast.S ;

end

% attribute full temporal traces to object neuron_fast and update a few
% options

neuron_fast.C = C_fast;
neuron_fast.C_raw = C_raw_fast; % not sure that this one is usefull
neuron_fast.C_prev = C_prev_fast; % not sure that this one is usefull
neuron_fast.S = S_fast;

% recompute C starting from S to avoid problems at borders between stacks
    delta = exp(-Obsd_fast.dt/Pbsd_fast.tauRise - Obsd_fast.dt/Pbsd_fast.tauDecay);
    gamma = exp(-Obsd_fast.dt/Pbsd_fast.tauRise)+exp(-Obsd_fast.dt/Pbsd_fast.tauDecay)-exp(-Obsd_fast.dt/Pbsd_fast.tauRise - Obsd_fast.dt/Pbsd_fast.tauDecay) ;
    eta = (Pbsd_fast.tauRise/Pbsd_fast.tauDecay)^(Pbsd_fast.tauDecay/(Pbsd_fast.tauDecay - Pbsd_fast.tauRise))*( Pbsd_fast.tauDecay/Pbsd_fast.tauRise - 1)./(exp(-Obsd_fast.dt/Pbsd_fast.tauDecay)-exp(-Obsd_fast.dt/Pbsd_fast.tauRise));
    for k = 1:size(neuron_fast.C, 1)
    neuron_fast.C(k, :) = filter(1,[eta, -eta * (gamma+delta)  ,  eta * delta  ],neuron_fast.S(k, :));
    end

T_fast = sum(tStacks);
Obsd_fast.Time = T_fast; % Number of time frames
neuron_fast.frame_range = [1 T_fast];
deconv_options_fast = struct('type', 'ar2', ... % model of the calcium traces. {'ar1', 'ar2'}, not used in bsd
            'method', deconv_method, ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded', 'bsd'}
            'Tar2', Tar2_fast,... % tau_d et tau_r for AR2
            'smin', smin_fast, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level, not used in bsd
            'optimize_pars', false, ...  % optimize AR coefficients
            'optimize_b', true, ...% optimize the baseline
            'max_tau', 15, ... % maximum decay time (unit: frame);
            'Obsd', Obsd_fast, ... % parameter for bsd algo, see above
            'Pbsd', Pbsd_fast, ... % parameter for bsd algo, see above
            'window_baseline', window_baseline_fast);  % size of time window to compute the baseline (default: 1000)   
neuron_fast.updateParams('deconv_options', deconv_options_fast);    

neuron_fast.orderROIs('snr');
neuron_fast.PNR = PNR; 
neuron_fast.Cn = Cn; 

results = neuron_fast.obj2struct();
results.moreParms = moreParms;

%% calculation of raw data

% load results if necessary
%fname = [Files.subfolderResults, 'Results.mat'];
%load(fname);

nCell = size(results.A, 2);
roifn = single(full(results.A));
pixh = size(results.Cn, 1);
pixw = size(results.Cn, 2);
roifn = reshape(roifn, pixh, pixw, size(roifn,2));
nTimeTot = 0;
        
for kStacks = 1 : nStacks
    
    kStacks
    nam = [Files.subfolderData Files.basename num2str(kStacks,'%.2d') '.tif'];
    rawData = bigread2(nam);
    nTime = size(rawData,3);    
    if kStacks == 1
        moyRoi = zeros(nCell, (nTime + 4)* nStacks); % +4 because the first stack is missing 4 frames
    end
    
    for kCell = 1 : nCell
        imageRoi = roifn(:, :, kCell);
        [row,col] = find(imageRoi);
        xmin = min(row); 
        xmax = min(max(row), size(rawData, 1)); 
        ymin = min(col);
        ymax = min(max(col), size(rawData, 2));
        for kt = 1 : nTime 
            A0 = double(rawData(xmin:xmax, ymin:ymax, kt));
            %moyRoi(kCell, kt + nTimeTot ) = mean(mean( (imageRoi(xmin:xmax, ymin:ymax)) .* A0));
            A0 = imageRoi(xmin:xmax, ymin:ymax) .* A0; 
            moyRoi(kCell, kt + nTimeTot ) = mean(A0(:));
        end 
    end
    nTimeTot = nTimeTot + nTime; 

end

moyRoi = moyRoi(:, 1:nTimeTot);
moyRoiNoBaseline = zeros(nCell, nTimeTot);

% calcul de F0

F0 = zeros(1, nCell); 
for kCell = 1 : nCell
    F0(kCell) = prctile(moyRoi(kCell, : ),10);
    [moyRoiNoBaseline(kCell, : ), b0] = remove_baseline_sliding_window(moyRoi(kCell, : ), 10, 3000);
end

figure;
% for kCell =1 : nCell
%     kCell
%     plot(results.C_raw(kCell, :)); 
%     hold on;
%     plot(results.C(kCell, :));
% %    plot(100+ (moyRoiNoBaseline(kCell, :))*max(results.C_raw(kCell, :))/max((moyRoiNoBaseline(kCell, :) )));
% %    plot(200+ (moyRoi(kCell, :)-F0(kCell))*max(results.C_raw(kCell, :))/max((moyRoi(kCell, :)-F0(kCell))));
%     pause;
%     hold off;
% end
 
results.AvgRoi = moyRoi ;
results.AvgRoiNoBaseline = moyRoiNoBaseline ;
results.F0 = F0;
%results.vessels = imread(Files.vessels);
%figure; imagesc(results.vessels)

% kCell = 5;
% output2 = struct; 
% output2.AvgRoi = moyRoi(kCell, :);
% output2.AvgRoiNoBaseline =moyRoiNoBaseline(kCell, :);
% output2.CNMFE_Craw = results.C_raw(kCell, :);
% output2.CNMFE_C = results.C(kCell, :);
% 
% fname = [Files.subfolderResults, 'Output2.mat'];
% save(fname, 'output2')

%% save results
       
fname = [Files.subfolderResults, 'Results.mat'];
save(fname, 'results', '-v7.3'); 
nam = [Files.subfolderResults 'PNR.tif'];
saveastiff(results.PNR, nam , tiff_optn);