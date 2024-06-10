
%% PREPROCESSING

clear all; close all; clc; 
warning('off', 'all'); 


% Script for preprocessing of the fiberscope images. 
% This script allows (in that order) to subtract the autofluorescence, 
% downsample spatially the stacks by a factor 2, apply a gaussian blur, 
% correct from motion correction with NormCorre and apply a mask before 
% going to the CNMFe pipeline. 

% TO BE DONE: 
% 1/ The first section called "REQUIREMENTS" should be completed for each 
% acquisition.

% 2/ The script allows to save two sequences of data in 3 folders:
% /1-Preprocess : corresponds to the folder where the data after motion 
% correction will be saved 
% /2-Cropped : corresponds to the folder where the data after motion
% correction and cropping will be saved.
% /3-Blood_Vessels : corresponds to the folder where the data after motion
% correction of the blood vessels will be saved.

% Please, make sure there is enough free space not to run out of memory. It
% requires to have half the space of the recording, i.e if the recording
% takes 100GB, you need 50GB to run this pipeline.

tic; % Start of the clock

% change compared to version 5 : adaptation for experiments with only green
% channel recorded. 

% change compared with standard version : multiple experiments analysed
% together (same template for registration)

%% PARAMETERS

% This section should be completed for each acquisition.

% 1- Required packages
% addpath(genpath('/Users/Expresso/Documents/Clara/MotionCorrection/2-NoRMCorre-master/'));

% 2- Paths to files and related parameters


nExp = 1; 
InputFiles = struct('basename', repmat({[]}, 1, nExp) );
Files = struct;


Files.autofluorescence = 'F:\22-06-30\LS_11_12\Autofluorescence.tif';  %%'path to autofluorescence image with the two lasers';
%Files.autofluorescenceBlueLaser = 'J:\20-09-24\LS_8_11\140024_Autofluorescence_488-220-10ms.tif';  %%'path to autofluorescence image with just the blue laser';
Files.offset = 'F:\22-06-30\LS_11_12\Offset.tif';  % path to offset image
Files.matlab = 'C:\Users\ventalon\Documents\MATLAB\Preprocess\Preprocessing6_SleepAwSleep_Rotation2.m'; % to keep track of parameters
Files.analysis = 'G:\Souris154148\22-06-30-DeepCadFirst\130014_sCMOS_154148-awake\0-Preliminary\'; %'path to folder where data will be analysed

% Exp1 
InputFiles(1).fluoDataFolder = 'F:\22-06-30\130014_sCMOS_154148-awake\'; %'path to folder with fluorescence raw data';
InputFiles(1).basename = '154148-awake'; %basename of fluorescence raw data (followed by -, a number and .tif)
InputFiles(1).fluorescenceTiming = 'sCMOS_AbsoluteTimingsSec.csv'; 
InputFiles(1).behavior = 'H:\23-08-03\155145_Basler\';


% % Exp2 
% InputFiles(2).fluoDataFolder = 'H:\23-08-03\162222_sCMOS_Souris172123_awake\'; %'path to folder with fluorescence raw data';
% InputFiles(2).basename = 'Souris172123_awake'; %basename of fluorescence raw data (followed by -, a number and .tif)
% InputFiles(2).fluorescenceTiming = 'sCMOS_AbsoluteTimingsSec.csv';
% InputFiles(2).behavior = 'H:\23-08-03\162222_Basler\';

% % Exp3 
% InputFiles(3).fluoDataFolder = 'I:\21-03-02\150021_sCMOS_133656-sleeppost\'; %'path to folder with fluorescence raw data';
% InputFiles(3).basename = '133656-sleeppost'; %basename of fluorescence raw data (followed by -, a number and .tif)
% InputFiles(3).fluorescenceTiming = 'sCMOS_AbsoluteTimingsSec.csv'; 
% InputFiles(3).behavior = 'M:\Souris133656\21-03-08\135948_Basler\';
% 
% % Exp4 
% InputFiles(4).fluoDataFolder = 'I:\21-03-02\153158_sCMOS_133656-sleeppost2\'; %'path to folder with fluorescence raw data';
% InputFiles(4).basename = '133656-sleeppost2'; %basename of fluorescence raw data (followed by -, a number and .tif)
% InputFiles(4).fluorescenceTiming = 'sCMOS_AbsoluteTimingsSec.csv'; 
% InputFiles(4).behavior = 'M:\Souris133656\21-03-08\152525_Basler\';


% 3- Actions to run

nChannels = 2; % number of channels recorded (2 for bicolor imaging, 1 for monocolor imaging)
% if 2 channels : channel1 is the left image and channel2 the right image
channelCorr = 2; % channel used for motion correction
channelNeurons = 2; % used to correct for inhomogeneous illumination and from which crop images will be extracted and used later by cnmfe

actions = struct;
actions.cropImages = false; % to keep only spatial positions that are illuminated at all times
actions.correctTranslations = false; % if awake recording --> True, if anesthetized --> False
actions.correctRotations = false; 
actions.correctBlueIllum = false; %!!! put false for summer 2020 experiments
actions.computeMousePosition = false; 
actions.fftFilt = false;
actions.correctBundleCropped = false;
actions.subtractBundleFootprint = false; % for now it implies that correctBlueIllum = false !!! let's ses

sfact = 1.5; % !!!! pour la nouvelle camera on multiplie par 1.5 toutes les échelles spatiales 

% 4- Computation parameters
% Data parameters
Exposure_Time_Autofluo = 10; % en ms
Exposure_Time_Acquisition = 10; % en ms
% Motion correction parameters
if channelCorr == channelNeurons
    % motion correction is performed on neurons
    gSig = 7; %8; % values for neurons
    gSiz = 15; %19; 
else
    % motion correction is performed on vessels
    gSig = 6; % values for vessels
    gSiz = 12; 
end

boundx = 50; % number of pixels on the sides to not take into account for motion correction, i.e central window.
boundy = 50;
    % !!! before we had 50 but it didnt work with mouse 136278; try new
    % value.
preFactor = 100; % multiply the stack by this factor in order to avoid quantification noise
% Other parameters

gaussFilt = 1.5; % spatial filtering to erase the bundle cores.
gaussFiltLarge = 40 ; % spatial filtering to the correction for inhomogeneous illumination
tempFilt = 1; %temporal averaging before computing movement. Should be an odd integer.  1: no filtering; Integer >1 : sliding window average. Use in case the SNR is low. 
% parameters for rotation correction
%  downsampleFact = 20; % temporal downsample factor for calculating rotations
%  DownSampleFactSave = 20; % temporal downsample factor for saving downsampled channel1 movie (and check that rotations are well corrected)
% all experiments
downsampleFact = 5; % temporal downsample factor for calculating rotations
DownSampleFactSave = 30; % temporal downsample factor for saving downsampled channel1 movie (and check that rotations are well corrected)
rotation_range = 6; % rotations explored between -rotation_range and + rotation_range 
    % !!!! usually rotation range = 4 is OK but for some exp (mouse 136278) it has to be
    % increased to 6 
rotation_step = 1; % 1 should be fine; should not affect the precision of rotation angle calculation
                   % A fit is performed to increase precision


% parameters to find the bundle boundaries automatically
%bundleDiameter = [638 648]; % [638 650]; % [Y X] (inversed coordinates in matlab)
    % !!! put 638 648 for summer 2020 experiments, roi is smaller; [638
    % 650] for later experiment
bundleDiameter = [900 906]; % [Y X] (inversed coordinates in matlab) %!!!  new camera 

if actions.correctBundleCropped == true  
    % dans ce cas on mesure sur imageJ les limites du bundle, qui sont en
    % plus differentes pour les deux voies...
    % il faut que les 2 images soient de meme taille car on souhaite les superposer. On se
    % base sur l'image cote neurones pour determiner cette taille.
    roiC2Large =  [910 1 1548-910+1 611]; % [xmin ymin deltaX deltaY ], channel 2 (right)
    bundleCenterC2 = [298 1231]; % center of the bundle [ycenter xcenter ], channel 1 (left) (matlab coordinate order)
    roiC1Large =  [22 1 roiC2Large(3) roiC2Large(4)]; % [xmin ymin deltaX deltaY ], channel 1 (left) 
    bundleCenterC1 = [294 344]; % center of the bundle [ycenter xcenter ], channel 1 (left) (matlab coordinate order)
    bundleDiameter = [628 644]; % [Y X] (inversed coordinates in matlab)

end

thickness = 10; % thickness of line to make profile and find bundle border

% Parameters for analysing mouse position
ledIntensityThreshold = 110; % adapt to data if necessary
        % if the maximum value of a basler frame is below this threshold we
        % consider that the led is hidden and interpolate data using other frames.
%subsampleFact= 0.5; % not used for now. Used by clara in the previous
%program, see if it is necessary

%% FIRST STEP : complete parameters

% Savings options
data_type = 'uint16'; % important in order not to loose bit information
tiff_optn = struct; % structure containing tiff options
tiff_optn.color = false;
tiff_optn.compress = 'no'; % 'no', 'lzw' (lossless), 'jpeg' (lossy crap), 'adobe' (lossless)
tiff_optn.message = false; % no output if false
tiff_optn.append = false; % check if that is necessarily true; a drawback is that if there are already files in the folder, new files will be concatenated with old files...
tiff_optn.overwrite = true; 
tiff_optn.big = false; % set to true if output > 4GB

% complete folder names
for k = 1:nExp
    InputFiles(k).fluorescenceTiming2 = [ InputFiles(k).basename , '.csv'];
    InputFiles(k).basename = [ InputFiles(k).basename , '-'];
end
Files.channel1 = [Files.analysis , '1-Channel1\']; 
Files.cropped = [Files.analysis , '2-Cropped\']; 
Files.channel2 = [Files.analysis , '3-Channel2\']; 
Files.results = [Files.analysis , '4-Results\'];

% Useful parameters for rigid motion correction (from CNMFE/normcorre
% programs)
if actions.fftFilt == 1
    psf = imread('G:\Filter.tif');
else
    psf = fspecial('gaussian', round(2*gSiz), gSig); % Filter
    ind_nonzero = (psf(:)>=max(psf(:,1)));
    psf = psf-mean(psf(ind_nonzero));
    psf(~ind_nonzero) = 0;   % only use pixels within the center disk
end

cd(Files.analysis )
% make analysis subdirectories if they do not exist
if ~exist(Files.channel1 , 'dir')
       mkdir(Files.channel1 )
end  
if nChannels == 2 
    if ~exist(Files.channel2, 'dir')
       mkdir(Files.channel2)
    end
end
if actions.cropImages == 1
    if ~exist(Files.cropped , 'dir')
       mkdir(Files.cropped )
    end
end
if ~exist(Files.results , 'dir')
       mkdir(Files.results )
end       

copyfile(Files.matlab , Files.results ); % to keep track of parameters and program version

nCrop = 3 * floor(sfact^2)* nChannels; %number of preprocessed files to be concatenated for CNMFE analysis (2nd part : update temporal on 100Hz data)


%% SECOND STEP : calculate relevant images / image parameters for registration

% 1 - compute image to subttract (combining autofluorescence and offset)

Autofluo = single(imread(Files.autofluorescence));
Offset = single(imread(Files.offset));
toSubtract = Exposure_Time_Acquisition / Exposure_Time_Autofluo * (Autofluo - Offset) + Offset; % Compute the image called "To_Subtract" which corresponds to the image to subtract to the subtstacks
%toSubtract =  Offset;  % !!! use this if autofluo is not correct or is missing
                        % only small impact on the results if we subtract only offset
                        % it is important to subtract one image (offset or autofluo or a combination of the 2 to get rid of
                        % camera artefacts on the central lines


%2 - find the positions of channel #1 and channel #2 images on raw images.

kExp = 1;
cd(InputFiles(kExp).fluoDataFolder )
kStack = 1; % The first movie is used to initialize the template for motion correction

       
% new method : automatic

% read raw images
raw = bigread2([InputFiles(kExp).basename , num2str(kStack), '.tif']);

if actions.correctBundleCropped == false  

    raw_avg = squeeze(squeeze(mean(raw,3)));
    %raw_avg = raw_avg - toSubtract + 1;
    %figure; imagesc(raw_avg); colormap(gray);
    [sy, sx] = size(raw_avg);
    sx1 = floor(sx/nChannels);
    roiC1 = zeros(1,4); roiC1Large = zeros(1,4);
    profileX = raw_avg((floor((bundleDiameter(1)-thickness)/2):floor((bundleDiameter(1)+thickness)/2)),:);
    profileX = mean(profileX, 1);
    profileXChannel1 = profileX(1:sx1); 
    [c,lags] = xcorr(profileXChannel1, ones(1, bundleDiameter(2)), floor(bundleDiameter(2)/2));
    [~,I] = max(c);
    roiC1Large(1)= max(lags(I), 1); 
    roiC1(1) = max(ceil(lags(I)/2/sfact), 1); % x0_1

    profileYChannel1 = raw_avg(:, roiC1(1) + floor((bundleDiameter(2)-thickness)/2): roiC1(1) + floor((bundleDiameter(2)+thickness)/2));
    profileYChannel1 = mean(profileYChannel1, 2);
    [c,lags] = xcorr(profileYChannel1, ones( bundleDiameter(1), 1), floor(bundleDiameter(1)/2));
    [~,I] = max(c);
    roiC1Large(2)= max(lags(I), 1); 
    roiC1(2) = max(ceil(lags(I)/2/sfact), 1);

    roiC1Large(3)=bundleDiameter(2);
    roiC1Large(4)=bundleDiameter(1);
    roiC1(3) = floor(bundleDiameter(2)/2/sfact);
    roiC1(4) = floor(bundleDiameter(1)/2/sfact);
    save([Files.channel1 ,'ROI.mat'],'roiC1');  save([Files.channel1 ,'ROILarge.mat'],'roiC1Large')

    if nChannels == 2 

        profileXChannel2 = profileX(sx1:end); 
        [c,lags] = xcorr(profileXChannel2, ones(1, bundleDiameter(2)), floor(bundleDiameter(2)/2));
        %stem(lags,c);
        [~,I] = max(c);
        roiC2 = zeros(1,4); roiC2Large = zeros(1,4);
        roiC2Large(1) = lags(I) + sx1;
        roiC2(1)  = ceil((lags(I) + sx1)/2/sfact);

        profileYChannel2 = raw_avg(:, roiC2Large(1)  + floor((bundleDiameter(2)-thickness)/2): roiC2Large(1)  + floor((bundleDiameter(2)+thickness)/2)); % !!! erreur, ce n'est pas roic2 car roic2 est dans l'autre espace ! enregistrer aussi les coordonnées du bundle dans l'espace initial des images brutes, ça sera utile pour la suite
        profileYChannel2 = mean(profileYChannel2, 2);
        [c,lags] = xcorr(profileYChannel2, ones( bundleDiameter(1), 1), floor(bundleDiameter(1)/2));
        [~,I] = max(c);
        roiC2Large(2) =  max(1, lags(I));
        roiC2(2)   = max(1, ceil(lags(I)/2/sfact));

        roiC2Large(3)=bundleDiameter(2);
        roiC2Large(4)=bundleDiameter(1);
        roiC2(3)  = floor(bundleDiameter(2)/2/sfact);
        roiC2(4)  = floor(bundleDiameter(1)/2/sfact);
        save([Files.channel2 ,'ROI.mat'],'roiC2'); save([Files.channel2 ,'ROILarge.mat'],'roiC2Large')

    end

else
    roiC2 =  [floor(roiC2(1)/2/sfact)  ceil(roiC2(2)/2/sfact) floor(roiC2(3)/2/sfact) floor(roiC2(4)/2/sfact)]; % [xmin ymin deltaX deltaY ], channel 2 (right)
    roiC1 =  [ceil(roiC1(1)/2/sfact)  ceil(roiC1(2)/2/sfact) floor(roiC1(3)/2/sfact) floor(roiC1(4)/2/sfact)]; % [xmin ymin deltaX deltaY ], channel 1 (left) 
    save([Files.channel2 ,'ROI.mat'],'roiC2');
    save([Files.channel1 ,'ROI.mat'],'roiC1');
end

% 3- compute the spatial distribution of the blue illumination to correct for inhomogeneities

% autofluoBlue = single(imread(Files.autofluorescenceBlueLaser));
% blueIllum = autofluoBlue - Offset; % image with which we will compute the spatial distribution of the blue illumination
% blueIllum = imresize(blueIllum,0.5); % rescale autofluo by a factor 2
% figure; imagesc(blueIllum); colormap(gray); 

blueIllum = prctile(raw,50,3); % Computation of F0 (getting rid of activity)
blueIllum = single(blueIllum) - (Offset);
blueIllum = imresize(blueIllum,0.5/sfact);
%imagesc(blueIllum); colormap(gray);

% create a mask on the neurons
if channelNeurons == 1 
    roiN = roiC1;
    neuronFolder = Files.channel1; 
else
    roiN = roiC2;
    neuronFolder =  Files.channel2; 
end

if channelCorr == 1 
    roiCorr = roiC1;
    corrFolder = Files.channel1; 
    if nChannels == 2
        roiNCorr =roiC2;
        nCorrFolder = Files.channel2;
    end
else
    roiCorr = roiC2;
    corrFolder = Files.channel2; 
    nCorrFolder = Files.channel1;
    roiNCorr =roiC1;
end

[xx,yy] = meshgrid(1:roiN(3),1:roiN(4));
if actions.correctBundleCropped == false  
    mask = hypot( (yy - roiN(4)/2)/(roiN(4)/2), (xx - roiN(3)/2)/(roiN(3)/2)) <= 1;
else
    mask = hypot( (yy - floor(bundleCenterC2(1)/2)+roiN(2))/(bundleDiameter(1)/4), (xx - floor(bundleCenterC2(2)/2)+roiN(1))/(bundleDiameter(2)/4)) <= 1;
end
%figure; imagesc(mask);

% put NaN outside of mask and correct for residual offset (if offset image
% was not acquired properly ?)
blueIllum =  blueIllum(roiN(2):roiN(2)+roiN(4)-1,roiN(1) :roiN(1) +roiN(3)-1); 
negativeMask = not(mask); 
residualOffset = mean(blueIllum(negativeMask)); 

if actions.correctBlueIllum == 1 

    blueIllum = (blueIllum-residualOffset)./ mask;
    blueIllum(isinf(blueIllum)) = NaN;

    gaussim = zeros(100,100);
    gaussim(50,50) = 1000;
    gaussim = imgaussfilt(gaussim, gaussFiltLarge);
    gaussim = gaussim./max(gaussim(:));
    blueIllumFilt = nanconv(blueIllum, gaussim).*mask;
    blueIllumFilt = blueIllumFilt / max(blueIllumFilt(:));
    blueIllumFilt(isnan(blueIllumFilt))=0;
    % figure; imagesc( blueIllumFilt);
    
    cd(Files.results);
    saveastiff(cast(blueIllumFilt,'single'),'blueIllumFilt.tif',tiff_optn); 

end

% 4- adapt image that will be subtracted to size and normalisation of the following section

if actions.subtractBundleFootprint == true
       cd(InputFiles(kExp).fluoDataFolder )   
       % first estimation of a "min" image but less noisy
       % save also a small substack to visualise results
       kStack = 1;
       nam = [InputFiles(kExp).basename , num2str(kStack), '.tif'];
       Y = bigread2(nam);
       kStack = 2;
       nam = [InputFiles(kExp).basename , num2str(kStack), '.tif'];
       Y = cat(3, Y(:,:,5:end), bigread2(nam));
       Y1p = prctile(Y,50,3); figure; imagesc(Y1p); %Y1p = prctile(Y,1,3);
       cd(Files.results);
       saveastiff(cast(Y1p,'single'),'Y1p.tif',tiff_optn); 
       %saveastiff(cast(Y1p(6:end,3:end),'single'),'Y1p2.tif',tiff_optn); 
       %saveastiff(cast(Y(6:end,3:end,1:100),'single'),'Y.tif',tiff_optn); 
       
       % then estimation of the cores footprint by the image of fluorescent
       % paper
       kStack = 2;
       FolderNameFluoPaper = 'H:\23-08-18\164509_sCMOS_papier fluo\papier fluo-';
       nam = [FolderNameFluoPaper , num2str(kStack), '.tif'];
       P = bigread2(nam);
       kStack = 3;
       nam = [FolderNameFluoPaper , num2str(kStack), '.tif'];
       P = cat(3, P, bigread2(nam));
       Pmean = mean(P, 3); figure; imagesc(Pmean);
       saveastiff(cast(Pmean,'single'),'Pmean.tif',tiff_optn); 
       saveastiff(cast(Pmean(1:end-5, 1:end-2),'single'),'Pmean2.tif',tiff_optn); 
       
       %shift imageFluoPaper compared to raw images
       % for now I measured them but later I will compute them automatically
       shiftPy = 5;
       shiftPx = 2; 
       
       % all this will be included in the preprocess loop later on. here it
       % is just to test. 
       %compute scaling factor of Pmean compared to Y
       Y1pS = single(Y1p) - single(toSubtract);
       Y1pS = imresize(Y1pS,0.5/sfact);
       Y1pN = Y1pS(roiN(2):roiN(2)+roiN(4)-1,roiN(1) :roiN(1) +roiN(3)-1);
       meanY1P = mean(Y1pN(mask));
       figure; imagesc(Y1pN.*mask);
       
       PmeanS = single(Pmean) - single(toSubtract);
       PmeanS = imresize(PmeanS,0.5/sfact);
       PmeanN = PmeanS(roiN(2)-2:roiN(2)+roiN(4)-1-2,roiN(1)-1 :roiN(1) +roiN(3)-1-1);
       meanPmean = mean(PmeanN(mask));
       figure; imagesc(PmeanN.*mask);
       
       factP = meanY1P/ meanPmean;
       
       Y = single(Y); Pmean = single(Pmean);
       Ycorr = Y(6:end,3:end,1:100) - toSubtract(6:end,3:end) -  factP * (Pmean(1:end-5, 1:end-2) - toSubtract(1:end-5, 1:end-2));
       saveastiff(cast(Ycorr,'single'),'Ycorr.tif',tiff_optn);
       YcorrS = imresize(Ycorr,0.5/sfact);
       YcorrSfilt = imgaussfilt(YcorrS, gaussFilt); 
       saveastiff(cast(YcorrSfilt,'single'),'YcorrSfilt.tif',tiff_optn);
       % on se rend compte d'un pb : a cause de l'absorption par les
       % vaisseaux la partie haute de l'image a moins de lumière que la
       % partie basse et la normalisation simple ne fonctionne pas
       
       % on essaye de normaliser par les images filtrées
       % dans un premier temps on passe dans l'espace final de l'image des neurones a l'echelle
       % des neurones pour plus de simplicité mais a terme il faudra faire
       % ça dans l'espace initial. 
       % pour Y pas de pb on prend les coordonnees qu'on connait
       YS = imresize(Y-toSubtract,0.5/sfact);
       YS = YS(roiN(2):roiN(2)+roiN(4)-1,roiN(1) :roiN(1) +roiN(3)-1, 1:100);
       saveastiff(cast(imgaussfilt(YS, gaussFilt),'single'),'YSfilt.tif',tiff_optn);
       
       PmeanS = imresize(Pmean-toSubtract,0.5/sfact);
       PmeanS = PmeanS(roiN(2)-2:roiN(2)-1+roiN(4)-2,roiN(1)-1 :roiN(1) +roiN(3)-1-1);
       PmeanS = PmeanS./ mask;
       PmeanS(isinf(PmeanS)) = NaN;
       PmeanSfilt = nanconv(PmeanS, gaussim).*mask;
       PmeanSfilt = PmeanSfilt / max(PmeanSfilt(:));
       PmeanSfilt(isnan(PmeanSfilt))=0;
       figure; imagesc(PmeanSfilt);
       figure; imagesc(blueIllumFilt);
       normIm =  blueIllumFilt./PmeanSfilt;
       normIm(isnan(normIm)) = 0;
       figure; imagesc(normIm);
       PmeanCorr = PmeanS.*normIm;
       figure; imagesc(PmeanCorr);
       factP = mean(YS(mask))/mean(PmeanCorr(mask));
       Ycorr2 = YS - factP * PmeanCorr;
       Ycorr2filt = imgaussfilt(Ycorr2, gaussFilt); 
       saveastiff(cast(Ycorr2filt,'single'),'Ycorr2filt.tif',tiff_optn);
       
       figure; imagesc(mask)
%        PmeanFilt = imgaussfilt(Pmean, gaussFiltLarge); 
%        PmeanFilt = PmeanFilt / max(max(PmeanFilt));
%        saveastiff(cast(PmeanFilt,'single'),'PmeanFilt.tif',tiff_optn); 
%        
       
       toSubtract = single(preFactor * prctile(Y,1,3));
       % save toSubtract
else 
    toSubtract = preFactor * (toSubtract + residualOffset); % multiply by Pre_factor in order to reduce the quantification noise 
        % essai. On soustrait l'autofluorescence + l'offset residuel calcule a
        % l'exterieur du bundle et la moitie de la std mesuree sur 
        % !!! for mouse 136278 : toSubtract = preFactor * (toSubtract + residualOffset - 1.5); 
    %toSubtract = preFactor * prctile(raw, 0.1, 3); 
end
toSubtract = imresize(toSubtract,0.5/sfact) ; % rescale autofluo by a factor 2


%% STEP 3 : perform registration

[d1,d2,~] = size(toSubtract);
if channelCorr == channelNeurons
    %!!!! recently changed us_fact -> should accelerate things but not
    %loose in precision
    options = NoRMCorreSetParms('d1',d1-2*boundy,'d2',d2-2*boundx,'init_batch',size(raw,3)-10,'max_shift',40,'correct_bidir',false, 'upd_template', false, 'us_fac', 5, 'use_parallel', true);
else
    options = NoRMCorreSetParms('d1',d1-2*boundy,'d2',d2-2*boundx,'init_batch',200,'max_shift',40,'correct_bidir',false, 'upd_template', false, 'us_fac', 5, 'use_parallel', true);
end

kStackAllExp = 1;
nStacks = zeros(nExp, 1);
corr1_full = [];
corr2_full = []; 
rotation_angle_full = [];
cross_corr_rotation_full  = [];
cross_corr_norotation_full = [];

 
for kExp = 1 : nExp 

    cd(InputFiles(kExp).fluoDataFolder);
    files = dir(strcat(InputFiles(kExp).basename, '*'));
    nStacks(kExp) = length(files);
    
    if actions.correctTranslations == true

        for kStack = 1: nStacks(kExp) %!!!!!
            % first step : read and resize
            nam = [InputFiles(kExp).basename , num2str(kStack), '.tif'];
            Y = bigread2(nam);
            Y = single(Y);
            Ys = imresize(Y,0.5/sfact); % rescale the substack by a factor 2 
            if kStack == 1
                Ys = Ys(:,:,5:end); % remove the first frames because of the problem with the start of the line-scan recording.    
            end    
            Yt = preFactor * Ys - toSubtract;  % multiply by Pre_factor in order to reduce the quantification noise and subtract autofluorescence to the substack
            Yf = imgaussfilt(Yt,gaussFilt); % Preparation for motion correction
            nT = size(Ys,3);
   
            % image1: used for motion correction
            % image2: the other channel if it exists
            image1 = Yf(roiCorr(2):roiCorr(2)+roiCorr(4)-1,roiCorr(1):roiCorr(1)+roiCorr(3)-1,:);
            if nChannels == 2
                image2 = Yf(roiNCorr(2) : roiNCorr(2) + roiNCorr(4) -1, roiNCorr(1) : roiNCorr(1) + roiNCorr(3) -1,:);
            end
            
            if actions.correctBlueIllum == 1 
                % correct for inhomogeneous illumination
                if channelNeurons == channelCorr 
                   image1 = single(image1)./repmat(blueIllumFilt,[1 1 nT]);
                   image1(isinf(image1))=0;
                   image1(isnan(image1))=0;
                else
                   image2 = single(image2)./repmat(blueIllumFilt,[1 1 nT]);
                   image2(isinf(image2))=0;
                   image2(isnan(image2))=0;
                end
            end
            
            image1 = uint16(image1);
            if nChannels == 2
                image2 = uint16(image2);
            end
            
                % run correction of rotation if needed
                if actions.correctRotations == true
                    if actions.fftFilt == 1
                            image1Fourier = fftshift(fft2(image1));
                            psf2 = imresize(psf,[size(image1Fourier,1) size(image1Fourier,2)]);
                            image1Fourier = image1Fourier .* psf2;
                            %Y0 = real(ifft2(image1Fourier));
                            Y0 = real(ifft2(ifftshift(image1Fourier)));
                    else
                            Y0 = imfilter(image1,psf,'symmetric');
                             %Y0 = image1; %!!!!!! test               
                    end
                    % compute reference image if  kExp * kStack == 1
                    if kExp * kStack == 1
                       ref_image = mean(Y0(boundy:end-boundy, boundx:end-boundx, 1:downsampleFact), 3);
                       % ref_image = median(Y0(boundy:end-boundy, boundx:end-boundx, :), 3);    %!!!! test
                        %!!! autres idees a tester : changer gSig et gSiz,
                        % voir si la tf peut fonctionner...                       
                       % ref_image = mean(Y0(bound:end-90,bound:end-90,1:downsampleFact), 3); % !!! souris 136278
                       saveastiff(cast(ref_image,data_type),[Files.results ,'refImageRotation.tif'],tiff_optn); % save the movie as Substack_channel2-*.tif
                      % figure; imagesc(ref_image);
                       center_of_FOV(1) = floor(size(ref_image, 1)/2);
                       center_of_FOV(2) = floor(size(ref_image, 2)/2);
                    end
                    % compute angles of rotation
                    [rotation_angle, cross_corr, cross_corr_norotation] = rotate_angle_stack_v2(ref_image, Y0, downsampleFact, rotation_range, rotation_step);
                    %[rotation_angle, cross_corr, cross_corr_norotation] = rotate_angle_stack(ref_image, Y0, downsampleFact, rotation_range, rotation_step); %!!! changed
                    rotation_angle_full = cat(1,rotation_angle_full,rotation_angle');
                    cross_corr_rotation_full  = cat(1,cross_corr_rotation_full,cross_corr');
                    cross_corr_norotation_full = cat(1,cross_corr_norotation_full, cross_corr_norotation');
                    parfor kz = 1 : size(image1, 3)
                        image1(:,:,kz) = rotate_image_interp(image1(:,:,kz),rotation_angle(kz),[0 0],center_of_FOV);
                    end
                    if nChannels == 2
                         parfor kz = 1 : size(image1, 3)
                            image2(:,:,kz) = rotate_image_interp(image2(:,:,kz),rotation_angle(kz),[0 0],center_of_FOV);
                         end
                    end
                end
                if actions.fftFilt == 1
                            image1Fourier = fftshift(fft2(image1));
                            psf2 = imresize(psf,[size(image1Fourier,1) size(image1Fourier,2)]);
                            image1Fourier = image1Fourier .* psf2;
                            Y0 = real(ifft2(ifftshift(image1Fourier)));
                else
                            Y0 = imfilter(image1,psf,'symmetric');
                             %Y0 = image1; %!!!!!! test               
                end
              
               %  Y0 = image1; %!!!!! test
                % run normcorre 
                if tempFilt > 1
                    Y0 = imboxfilt3(Y0,[1 1 tempFilt]);
                end

                if kExp * kStack == 1
                    template1 = [];
                end
                % find shifts using normcorre
                 [~,shifts1,template1,options] = normcorre_batch(Y0(boundy+1:end-boundy,boundx+1:end-boundx,:),options,template1);         

                % apply shifts 
                  M1 = apply_shifts(image1,shifts1,options,boundy,boundx); % apply shifts to channel 1 dataset
                  tiff_filename_1 = [corrFolder ,'Substack-',num2str(kStackAllExp,'%03d'), '.tif'];
                  saveastiff(cast(M1,data_type),tiff_filename_1,tiff_optn); % save the movie as Substack_channel1-*.tif
                % compute metrics for evaluation of the performance of registration
                    [corr1,~,~] = motion_metrics(M1(boundy:end-boundy,boundx:end-boundx,:),options.max_shift);
                    corr1_full = cat(1,corr1_full,corr1);
                if nChannels == 2
                    M2 = apply_shifts(image2,shifts1,options,boundy,boundx); % apply shifts to channel 2 imaging dataset
                    tiff_filename_2 = [nCorrFolder, 'Substack-',num2str(kStackAllExp,'%03d'), '.tif']; 
                    saveastiff(cast(M2(:,:,1:end),data_type),tiff_filename_2,tiff_optn); % save the movie as Substack_channel2-*.tif
                    [corr2,~,~] = motion_metrics(M2(boundy:end-boundy,boundx:end-boundx,:),options.max_shift);
                    corr2_full = cat(1,corr2_full,corr2);
                end
                shifts_r = squeeze(cat(3,shifts1(:).shifts));

                save([corrFolder ,'shiftsyx',num2str(kStackAllExp,'%03d'), '.txt'], 'shifts_r', '-ascii', '-tabs'); % save the shift computed with Normcorre accross the recording.
               % save([corrFolder ,'shifts1',num2str(kStackAllExp,'%03d'), '.mat'],'shifts1', '-struct' ); % save the shift computed with Normcorre accross the recording.

            disp(['Motion correction done for movie ' num2str(kStackAllExp) ' and the movie has been saved...'])
            kStackAllExp =  kStackAllExp + 1;

        end
        
                %save correlations and rotation angles
                save([corrFolder ,'corr1_full.mat'],'corr1_full')        
                save([Files.results ,'template1.mat'],'template1');
                if nChannels == 2
                        save([corrFolder ,'corr2_full.mat'],'corr2_full')
                end

                if actions.correctRotations == true    
                    save([corrFolder ,'rotation_angle_full.mat'],'rotation_angle_full')
                    save([corrFolder ,'cross_corr_rotation_full.mat'],'cross_corr_rotation_full')
                    save([corrFolder ,'cross_corr_norotation_full.mat'],'cross_corr_norotation_full')
                end
        
    else
        for kStack = 1:length(files)
            nam = [InputFiles.basename , num2str(kStackAllExp), '.tif'];
            Y = bigread2(nam);
            Y = single(Y);
            %Y = single(read_file([Files.basename , num2str(id), '.tif'])); 
            Ys = imresize(Y,0.5/sfact); 
            Yt = preFactor * Ys - toSubtract; 
            Yf = uint16(imgaussfilt(Yt,gaussFilt));
            [d1,d2,T] = size(Yf);
            % !!!! je ne fais pas les vaisseaux, pas utile pour le test en cours
%             image1 = Yf(roiC1(2):roiC1(2)+roiC1(4)-1,roiC1(1):roiC1(1)+roiC1(3)-1,:);
%             tiff_filename_1 = [Files.channel1 ,'Substack_channel1-',num2str(kStackAllExp,'%03d'), '.tif']; 
%             saveastiff(cast(image1,data_type),tiff_filename_1,tiff_optn); % dans cette version on n'a pas enlevé les 4 premières frames
            if nChannels == 2
                image2 = Yf(roiC2(2) :roiC2(2) +roiC2(4) -1,roiC2(1) :roiC2(1) +roiC2(3) -1,:);
                tiff_filename_2 = [Files.channel2 ,'Substack_channel2-',num2str(kStackAllExp,'%03d'), '.tif']; 
                saveastiff(cast(image2,data_type),tiff_filename_2,tiff_optn); % dans cette version on n'a pas enlevé les 4 premières frames
            end
            disp(['Preprocessing done for movie ' num2str(kStackAllExp) 'and the movie has been saved...'])
            kStackAllExp =  kStackAllExp + 1;
        end
    end

end
    
save([corrFolder ,'nStacks.mat'],'nStacks')
nStackAllExp = kStackAllExp -1;

%% compute downsampled version of channel 1 to see if the correction of motion (of rotation for ex) is OK

fprintf('Compute full-session downsampled channel1 movie... \n');

cd(Files.channel1);
nStacksTot = length(dir('*.tif'));
cd(Files.results);
tStacks = zeros(1, nStacksTot); % number of time points in each stack
if nStacksTot >1
        kStacks = 2; % we initialize the number of elements using the length of the second stack because the first one is shorter 
else
        kStacks = 1;
end
nam = [Files.channel1 'Substack-' num2str(kStacks,'%.3d') '.tif'];
A = bigread2(nam);
Nx = size(A,1);
Ny = size(A,2);
Nt = size(A,3);
subMovie = zeros(Nx, Ny, floor(Nt*(nStacksTot)/DownSampleFactSave));

reliquat = [];
kt0 = 1;

for kStacks = 1 : nStacksTot
       kStacks
       nam = [Files.channel1 'Substack-' num2str(kStacks,'%.3d') '.tif'];
       A = bigread2(nam);
       tStacks(kStacks) = size(A, 3); % number of frames in the stack;
       A = cat(3, reliquat, A);
       ntk = size(A, 3); % number of frames in the stack;
       if  mod(ntk,DownSampleFactSave) ==0 
           reliquat = [];
       else
           reliquat = A(:,:, end - mod(ntk,DownSampleFactSave) +1 : end);
           A = A(:,:, 1 : end - mod(ntk,DownSampleFactSave) );
       end
       Ntksub = floor(ntk/DownSampleFactSave);
       subMovie(:,:, kt0 : kt0 + Ntksub -1) =  imresize3(A, [size(A, 1) size(A,2) Ntksub]);
       kt0 = kt0 + Ntksub;
end
subMovie = subMovie(:,:, 1 : kt0 -1 );

fprintf('Saving downsampled stack... \n');

nam = [Files.results 'downsampled_Channel1Movie' '.tif'];
saveastiff(uint16(subMovie), nam , tiff_optn);
downsampled_Channel1Movie_Avg = uint16(squeeze(mean(subMovie, 3))); 
clear subMovie A
fprintf('Done');





%% Mask and crop channel 1 to elminate parts of the field of view that have not been illuminated continuously

% We just saved the stacks after motion correction. However, if you open
% one of them, you observe that the neurons are fixed while the fiber bundle
% is moving to compensate from motion. In order not to detect abberant
% activity due to motion of the bundle at the boudaries, we apply a mask to
% eliminate these bundle movements.

% Creation of the mask based on maximum motion of the FOV accross the
% recording.


if actions.cropImages == 1
       
    % method : create an ellipse mask using an image representing the
    % minimum value of all pixels
   
   % minImageVect = zeros(319,325,84);
    
   % check where the pbs are
%    vectMean = [];
%     for kStacks = 1 : nStacksTot 
%            kStacks
%            nam = [neuronFolder 'Substack-' num2str(kStacks,'%.3d') '.tif'];
%            A = bigread2(nam);
%            vectMean = [vectMean ; squeeze(mean(A,[1 2]))];
%     end
%     figure; plot(vectMean);  
   
     for kStacks = 1 : nStacksTot %/2 
           kStacks
           nam = [neuronFolder 'Substack-' num2str(kStacks,'%.3d') '.tif'];
           A = bigread2(nam);
           if kStacks ==1 
               minImage = min(A,[], 3);
               %meanImage = mean(A, 3);
           else 
               %minImage = min(minImage, min(A,[], 3)); %!!! this line should be uncommented
               %meanImage = meanImage + mean(A, 3);
               B =  min(A,[], 3); % this line and the next 3 should be commented
               if (sum(sum(B))) >  0
                    minImage = min(minImage, B);
               end
           end
      end
        cd(neuronFolder)
 
    % tracer la puissance moyenne de chaque image en fct du temps et
    % regarder si le pb d'intensité d'exc nulle se produit souvent
    % en fonction, décider de comment faire la suite...
        
        
    figure; imagesc(minImage);
     %figure; imagesc(meanImage>100);
        % Create the mask
        % Binarize
        BW = imbinarize(minImage);
        % Extract the maximum area
       % BW = imclearborder(BW);
       % BW = bwareafilt(BW,1);
        % Calculate centroid, orientation and major/minor axis length of the ellipse
        s = regionprops(BW,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
        if s(1).Orientation < 0
            s(1).Orientation=s(1).Orientation+360;
        end
        h = drawellipse('Center',s(1).Centroid ,'SemiAxes',[s(1).MajorAxisLength/2-1 s(1).MinorAxisLength/2-1], 'RotationAngle',s(1).Orientation,'StripeColor','m');
  %      h = drawellipse('Center',s(1).Centroid +[2 -2],'SemiAxes',[s(1).MajorAxisLength/2-4 s(1).MinorAxisLength/2-4], 'RotationAngle',s(1).Orientation,'StripeColor','g');
 %    h = drawellipse('Center',s(1).Centroid+[-4 2],'SemiAxes',[s(1).MajorAxisLength/2-2 s(1).MinorAxisLength/2-1], 'RotationAngle',s(1).Orientation-10+180,'StripeColor','g');
%      h = drawellipse('Center',s(1).Centroid-[2 -2],'SemiAxes',[s(1).MajorAxisLength/2-3 s(1).MinorAxisLength/2-1], 'RotationAngle',s(1).Orientation-20+180,'StripeColor','k');
     %   h = drawellipse('Center',s(1).Centroid-[3 1],'SemiAxes',[s(1).MajorAxisLength/2-2 s(1).MinorAxisLength/2-0], 'RotationAngle',s(1).Orientation,'StripeColor','b');
     % h = drawellipse('Center',s(1).Centroid+[-1.5 1.5],'SemiAxes',[s(1).MajorAxisLength/2-1.5 s(1).MinorAxisLength/2+1.5], 'RotationAngle',s(1).Orientation-4,'StripeColor','m');

        ellipseMask = createMask(h);
        %imshow(ellipseMask)
        %figure; imagesc(minImage.*uint16(ellipseMask));
        % check that there is no null value of minImage within the mask
%         kloop=2;
%         while min(minImage(ellipseMask)) == 0
%             h = drawellipse('Center',s(1).Centroid,'SemiAxes',[s(1).MajorAxisLength/2-kloop s(1).MinorAxisLength/2-kloop], 'RotationAngle',s(1).Orientation,'StripeColor','m');
%             kloop=kloop+1;
%             ellipseMask = createMask(h);
%         end
        
        %find contours of mask
        profile1 = max(ellipseMask,[], 1 ) ;
        %figure; plot(profile1)
        I2  = find(profile1);
        xmin = min(I2);
        xmax = max(I2);
        profile2 = max(ellipseMask,[], 2 ) ;
        %figure; plot(profile1)
        I1  = find(profile2);
        ymin = min(I1);
        ymax = max(I1);
        dx1 = xmax - xmin + 1; % final size of green channel image along x
        dy1 = ymax - ymin + 1;

        %maskold = mask;
        mask = uint16(ellipseMask(ymin:ymax, xmin:xmax));
        
        save([Files.results ,'ROI_Crop.mat'],'xmin', 'ymin', 'dx1', 'dy1')
        nam = [Files.results 'mask.tif'];
        saveastiff(mask,nam ,tiff_optn); 
        nam = [Files.results 'downsampled_Channel1Movie_Crop_Avg' '.tif'];
        saveastiff(downsampled_Channel1Movie_Avg(ymin:ymax, xmin:xmax), nam , tiff_optn); % !! not checked yet
        nam = [Files.results 'minImage.tif'];
        saveastiff(minImage, nam ,tiff_optn); 
        nam = [Files.results 'ellipseMask_fullROI.tif'];
        saveastiff(uint16(ellipseMask), nam ,tiff_optn); 
        % Apply the mask to all the previously motion corrected substacks
        kConcat = 1;
        kStackAllExp = 1;

    for kExp = 1 : nExp
         for kStack = 1:nStacks(kExp)
            nam = ['Substack-', num2str(kStackAllExp,'%03d'), '.tif'];
            kStackAllExp = kStackAllExp+1;
            Y = bigread2(nam);
            Y = Y(ymin:ymax, xmin:xmax, :);
            if mod(kStack,nCrop) == 1
                croppedImage = Y.*mask; % apply the mask
            elseif mod(kStack,nCrop) == 0
                croppedImage = cat(3, croppedImage, Y.*mask) ;
                tiff_filename_crop = [Files.cropped ,'Substack_cropped-',num2str(kConcat,'%02d'), '.tif'];
                saveastiff(croppedImage,tiff_filename_crop,tiff_optn); % save the movie as Substack_cropped-*.tif
                disp(['Cropped and concatenated movie ' num2str(kConcat) ' has been saved...'])
                kConcat = kConcat + 1;
            else 
                croppedImage =   cat(3, croppedImage, Y.*mask) ;
            end
        end

        if  mod(kStack,nCrop) > 0
                tiff_filename_crop = [Files.cropped ,'Substack_cropped-',num2str(kConcat,'%02d'), '.tif'];
                saveastiff(croppedImage, tiff_filename_crop,tiff_optn); % save the movie as Substack_cropped-*.tif
                disp(['Cropped and concatenated movie ' num2str(kConcat) ' has been saved...'])
                kConcat = kConcat + 1;
        end
    end
clear croppedImage;

end


for kExp = 1 : nExp
    copyfile([InputFiles(kExp).fluoDataFolder, InputFiles(kExp).fluorescenceTiming ] , [Files.results, InputFiles(kExp).basename, InputFiles(kExp).fluorescenceTiming ]); % save fluorescenceTiming in analysis folder
    copyfile([InputFiles(kExp).fluoDataFolder, InputFiles(kExp).fluorescenceTiming2 ] , Files.results ); % save fluorescenceTiming in analysis folder
end
disp(['Preprocessing done in ', num2str(toc), 'seconds']); % End of the clock



%% Compute mouse position (in pixels) and save behavior camera timings

if actions.computeMousePosition == 1 

    fprintf('Find mouse position accross time... \n');

    for kExp = 1 : nExp 
        % fist step: load behavior images
        dirBehavior = [ dir([InputFiles(kExp).behavior, 'Basler', '.tif'] ); dir([InputFiles(kExp).behavior, 'Basler-*', '.tif'])];
        imData = bigread2([InputFiles(kExp).behavior, 'Basler.tif']);
        %imData = imresize(imData,subsampleFact);
        if length(dirBehavior)>1
            for i=2:length(dirBehavior)
                %imData = cat(3, imData, imresize(bigread2([InputFiles(kExp).behavior, 'Basler-', num2str(i), '.tif']), subsampleFact));
                imData = cat(3, imData, bigread2([InputFiles(kExp).behavior, 'Basler-', num2str(i), '.tif']));
            end
        end
        nFrames=size(imData,3);
        % second step : extract mouse position. Could find a more elegant method than a loop...
        mousePos1 = zeros(3, nFrames);
        for i = 1:nFrames
           % A = imData(1:end,1:end,i);
           % A = imData(1:end,20:1162,i); %!!! changer aussi 4 lignes en dessous
            A = imData(1:end,20:end,i); %!!! changer aussi 4 lignes en dessous
            if max(max(A))>ledIntensityThreshold
                A = imgaussfilt(A, 2);
               [~,I] = max(A(:));
               [mousePos1(2, i), mousePos1(1, i)] = ind2sub(size(A),I) ; 
               mousePos1(1, i) = 19 + mousePos1(1, i) ; %!!!
               mousePos1(3, i) = 1;
           else
               mousePos1(1, i) = NaN; 
               mousePos1(2, i) = NaN;
               mousePos1(3, i) = 0;
           end
        end
       % mousePos1(mousePos1<10) = NaN;
       % mousePos1(mousePos1>570) = NaN;
       
        mousePos1(1,:) = fillmissing(mousePos1(1,:), 'linear','EndValues','nearest');
        mousePos1(2,:) = fillmissing(mousePos1(2,:), 'linear','EndValues','nearest');
       
         mousePos1x = mousePos1(1,:);
    %    mousePos1x(mousePos1x>1235) = NaN; 
        mouseSpeed1x = gradient(mousePos1x);
        artefacts = mouseSpeed1x' - smooth(mouseSpeed1x', 'rlowess');
  
        figure; 
        ax1 = subplot(2,1,1);
        plot(mousePos1x);
        ax2 = subplot(2,1,2);
        plot(gradient(mousePos1x));  hold on; plot(artefacts); 
        linkaxes([ax1,ax2],'x');
        
        speedTh = 40;
        k1  = 1;
        kseg = 5; 
        while k1 < size(mousePos1, 2)-kseg -1
            %k1 
            k2max = 0;
            if abs(artefacts(k1))>speedTh
                s = sign(artefacts(k1));
                mousePos1x(k1+1) = NaN;
                for k2= 1:2*kseg
                    if s*artefacts(k1+k2)<-speedTh
                        mousePos1x(k1+2 : k1+k2+1) = NaN(1, k2);
                        k2max = k2;
                        %break;
                    end
                end
                if k2max == 0
                        for k2=1:kseg
                            if artefacts(k1+k2)>speedTh
                                mousePos1x(k1+2 : k1+k2+1) = NaN(1, k2);
                            end
                        end
                end
            end
            k1 = k1 + k2max + 1;
        end
        mousePos1x = fillmissing(mousePos1x,'linear','EndValues','nearest');

        figure; 
        ax1 = subplot(2,1,1);
        plot(mousePos1x);
        ax2 = subplot(2,1,2);
        plot(gradient(mousePos1x)); % hold on; plot(artefacts); 
        linkaxes([ax1,ax2],'x');
       
        mousePos1(1,:) =  mousePos1x;
       
        save([Files.results, 'mousePos1.mat' ],'mousePos1');
        %copyfile([InputFiles(kExp).behavior, 'Basler_AbsoluteTimingSec.csv' ] , [Files.results, InputFiles(kExp).basename, 'Basler_AbsoluteTimingSec.csv' ]); % save behavior timings in analysis folder
        mousePos1(3,1)=-0.1; mousePos1(3,2)=1.1;
        figure; 
        %subplot(3,1,1); plot( mousePos1(1,:));subplot(3,1,2); plot( mousePos1(2,:)); subplot(3,1,3); plot(mousePos1(3,:))
        ax(1) = subplot(3,1,1); plot( mousePos1(1,:)); 
        ax(2) = subplot(3,1,2); plot( mousePos1(2,:)); 
        ax(3) = subplot(3,1,3); plot(mousePos1(3,:));
        linkaxes(ax,'x')
    end
    
    clear imData;
    clear A; 
end


%%

% % Preprocessing of the Basler movies
% cd(Files.behaviorDataFolder )
% files = dir('*.tif');
% 
% for i = 1:length(files)
%     basename = files(i).name;
%     movie = read_file(basename);
%     movie_new = imresize(movie,0.5); % rescale the movie by a factor 2 
%     saveastiff(movie_new,['Preprocess_', basename],tiff_optn); % save the movie as Preprocess-Basler*.tif
% end

