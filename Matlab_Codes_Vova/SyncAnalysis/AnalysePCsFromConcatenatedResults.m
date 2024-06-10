
% Dans ce programme on va analyser les PCs à partir des données analysées en manip entière (sommeil et comportement). 
% Pour simplifier on va essayer d'utiliser Master tel qu'il est et analyser les PCs avec. 
% Du coup pour ça : 
% Je crée un dossier (à la main) : G:\Souris154148\22-06-30\130014_sCMOS_154148-awake_FromAllExp
% Autre dossier Awake où je mets les fichiers timings de la caméra sCMOS, ceux de la caméra Basler, le fichier de position, un dossier tmp pour les pdfs. 
% Ensuite j'extrais le fichier results

%MainPath='G:\Souris154148\22-07-07'; nbExp=[1,2,1];
MainPath='G:\Souris154148\22-06-30'; nbExp=[1,1,1];
%MainPath='G:\Souris168679\23-05-10'; nbExp=[0,1,1];

expPaths=GetExpsPaths('NbExp',nbExp,'Path',MainPath);               % Load the paths to the differents file (or ask the user for them if it wasn't done before)
results=LoadExtractedCalciumData(MainPath,expPaths,'',1,'');             % Load Calcium data (extract it from the concatenated file if it wasn't done before)

Params=GetParams; % bien rentrer tous les paramètres y compris ceux des dossiers dans GetParams (un peu redondant avec MainPath mais plus simple pour l'instant)

exp=1;
if ~exist('currentResultsFile','var'); currentResultsFile=[]; results=[] ; end      % Initialize a variable used later to avoid the (long) reloading of results.mat if it does not exist
    
disp(['Starting analysis : experiment ' strrep(Params.Paths{exp},[Params.DataFolder '/'],'')])
Experiment=[];
Experiment.path=[Params.Paths{exp} '/'];                                            % Select one experiment
    
% Call ComputePlaceActivity.m 
disp('Start Computing Place-related Activities')
[Params,Neurons,Experiment,results,currentResultsFile]=ComputePlaceActivity(Params,Experiment,currentResultsFile,results);
    
% Call FindPlaceCells.m
disp('Start Finding Place Cells')
[Experiment,Neurons,Params]=FindPlaceCells(Experiment,Neurons,Params,results);
     
% If saving Matrices is wanted
    if Params.SaveMats
        disp('Saving Matrices')
        save([Experiment.path 'Neurons'], 'Neurons')
        save([Experiment.path 'Params'], 'Params')
        save([Experiment.path 'Experiment'], 'Experiment')
    end
close all;

