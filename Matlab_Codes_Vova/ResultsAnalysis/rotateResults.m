function [results] = rotateResults(vessels1, results)

% input
    % vessels1 : image of vessels to align results with 
    % results : results structure created by cnmfe and containing spatial
                % information to align with vessels1

% output
    % results : rotated results
    
% Files = struct;
% Files.vessel1 = 'E:\Souris133656\21-03-02\134248_sCMOS_133656-awake\4-Results\VesselsCropped.tif';
% Files.results = 'E:\Souris133656\21-03-08\135948_sCMOS_133656-awake\4-Results\First_Analysis_CroppedMovies\Results.mat';
% vessels1 =  read_file(Files.vessel1);
% load(Files.results);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0- Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gSig = 6; 
gSiz = 12; 
actionFilterImages = 1; % 1 : filter the images before trying to find tranformation; 0 : do not filter 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1- Find geometrical tranformation between sessions and rotate vessels2 and PNR 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vessels2 = results.vessels;
PNR2 = results.PNR;
s2 = size(PNR2);
s1 = size(vessels1);

% Useful parameters for rigid motion correction
psf = fspecial('gaussian', round(2*gSiz), gSig); % Filter
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;   % only use pixels within the center disk

figure; 
waitfor(msgbox('When the image appears, draw the ROI that will be used for image registration.'));
imagesc(vessels1); colormap(gray); caxis([min(min(vessels1)) max(max(vessels1))]);       

roi = drawrectangle; % Draw a rectangle for the region that will be used to match sessions;
        pause();
        xmin = round(roi.Position(1));
        ymin = round(roi.Position(2));
        Dx = round(roi.Position(3));
        Dy = round(roi.Position(4));
         
% select different parameters for imregtform registration 
optimizer = registration.optimizer.OnePlusOneEvolutionary;
%optimizer = registration.optimizer.RegularStepGradientDescent;
%metric    = registration.metric.MattesMutualInformation;
metric    = registration.metric.MeanSquares;

            % Tune the properties of the optimizer to get the problem to converge
            % on a global maxima and to allow for more iterations.
             optimizer.InitialRadius = 0.003;
            %  optimizer.Epsilon = 1.5e-4;
            % optimizer.GrowthFactor = 1.01;
            optimizer.MaximumIterations = 1000;

if actionFilterImages == 1
           vessels2FiltCrop = imfilter(vessels2,psf,'symmetric');
           vessels1FiltCrop = imfilter(vessels1,psf,'symmetric'); 
           vessels2FiltCrop = vessels2FiltCrop(ymin:ymin+Dy, xmin:xmin+Dx);
           vessels1FiltCrop = vessels1FiltCrop(ymin:ymin+Dy, xmin:xmin+Dx);
else
           vessels2FiltCrop = vessels2(ymin:ymin+Dy, xmin:xmin+Dx);
           vessels1FiltCrop = vessels1(ymin:ymin+Dy, xmin:xmin+Dx);
end

% find tranformation between vessels1 and vessels2
tform = imregtform(vessels2FiltCrop,  vessels1FiltCrop, 'rigid', optimizer, metric);

% apply tranformation to vessel2 and PNR2
results.vessels = imwarp(vessels2,tform,'OutputView',imref2d(s1));
results.PNR = imwarp(PNR2,tform,'OutputView',imref2d(s1));

% visual check that tranformation is OK
figure; imshowpair(vessels1, results.vessels,'Scaling','joint');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2- Rotate results.A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A2Reg = sparse(s1(1)*s1(2), size(results.A,2));

for k = 1 : size(results.A,2)
    Ak2 = reshape(results.A(:,k), s2(1), s2(2)) ;
    Ak2Registered = sparse(imwarp(full(Ak2),tform,'OutputView',imref2d(s1)));
    Ak2Registered = reshape(Ak2Registered, s1(1)*s1(2), 1);
    A2Reg(:,k) =  Ak2Registered;
end

results.A = A2Reg;
clear A2Reg;
