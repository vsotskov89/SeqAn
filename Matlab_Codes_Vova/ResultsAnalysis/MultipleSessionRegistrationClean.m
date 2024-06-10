
%% Parameters

clear all; close all; clc; 
warning('off', 'all'); 

% Savings options
data_type = 'uint16'; % important in order not to loose bit information
tiff_optn = struct; % structure containing tiff options
tiff_optn.color = false;
tiff_optn.compress = 'no'; % 'no', 'lzw' (lossless), 'jpeg' (lossy crap), 'adobe' (lossless)
tiff_optn.message = false; % no output if false
tiff_optn.append = false; % check if that is necessarily true; a drawback is that if there are already files in the folder, new files will be concatenated with old files...
tiff_optn.overwrite = true; 
tiff_optn.big = false; % set to true if output > 4GB


files = struct; 

% souris 133656
files.mouseFolder = '/Users/ventalon/Documents/Cathie/Work/Analyse/Souris133656/';
files.session1SubFolder = '21-03-08_135948_awake/';
files.session2SubFolder = '21-03-02_134248_awake/';
files.vessels = 'VesselsCropped.tif';
files.PNR = 'PNR.tif';
files.subfolderResults = 'Results/'; 

gSig = 6; 
gSiz = 12; 

    


%%

files.vessels1 = [files.mouseFolder ,files.session1SubFolder ,files.vessels ] ;
files.vessels2 = [files.mouseFolder ,files.session2SubFolder ,files.vessels ] ;
files.PNR1 = [files.mouseFolder ,files.session1SubFolder ,files.PNR ] ;
files.PNR2 =  [files.mouseFolder ,files.session2SubFolder ,files.PNR ] ;
files.results = [files.mouseFolder ,files.subfolderResults ] ;

vessels1 = read_file(files.vessels1);
vessels2 = read_file(files.vessels2);
PNR1 = read_file(files.PNR1);
PNR2 = read_file(files.PNR2);

% Useful parameters for rigid motion correction
psf = fspecial('gaussian', round(2*gSiz), gSig); % Filter
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;   % only use pixels within the center disk

actionFilterImages = 1; % 1 : filter the images before trying to find tranformation; 0 : do not filter 

figure; 
waitfor(msgbox('When the image appears, draw the ROI that will be used for image registration.'));
imagesc(vessels1); colormap(gray); caxis([min(min(vessels1)) max(max(vessels1))]);       

roi = drawrectangle; % Draw a rectangle for the region that will be used to match sessions;
        pause();
        xmin = round(roi.Position(1));
        ymin = round(roi.Position(2));
        Dx = round(roi.Position(3));
        Dy = round(roi.Position(4));


  
        % recalage avec les vaisseaux
       
        % select different parameters for imregtform registration
        
        %[optimizer,metric] = imregconfig('multimodal');
       
        optimizer = registration.optimizer.OnePlusOneEvolutionary;
        %optimizer = registration.optimizer.RegularStepGradientDescent;
        
        %metric    = registration.metric.MattesMutualInformation;
        metric    = registration.metric.MeanSquares;

            % Tune the properties of the optimizer to get the problem to converge
            % on a global maxima and to allow for more iterations.
            %  optimizer.InitialRadius = 0.009;
            %  optimizer.Epsilon = 1.5e-4;
            % optimizer.GrowthFactor = 1.01;
            % optimizer.MaximumIterations = 300;

       if actionFilterImages == 1
           vessels2FiltCrop = imfilter(vessels2,psf,'symmetric');
           vessels1FiltCrop = imfilter(vessels1,psf,'symmetric'); 
           vessels2FiltCrop = vessels2FiltCrop(ymin:ymin+Dy, xmin:xmin+Dx);
           vessels1FiltCrop = vessels1FiltCrop(ymin:ymin+Dy, xmin:xmin+Dx);
       else
           vessels2FiltCrop = vessels2(ymin:ymin+Dy, xmin:xmin+Dx);
           vessels1FiltCrop = vessels1(ymin:ymin+Dy, xmin:xmin+Dx);
       end
            
       tform = imregtform(vessels2FiltCrop,  vessels1FiltCrop, 'rigid', optimizer, metric);
       vessels2Registered = imwarp(vessels2,tform,'OutputView',imref2d(size(vessels1)));
       PNR2Registered = imwarp(PNR2,tform,'OutputView',imref2d(size(PNR1)));
        
  
        figure; imshowpair(vessels1, vessels2Registered,'Scaling','joint');
        figure; imshowpair(PNR1, PNR2Registered,'Scaling','joint');
        saveastiff(PNR2Registered,[files.results ,'PNR2Registered.tif' ] ,tiff_optn); 
        saveastiff(vessels2Registered,[files.results ,'vessels2Registered.tif' ] ,tiff_optn); 

   
   
