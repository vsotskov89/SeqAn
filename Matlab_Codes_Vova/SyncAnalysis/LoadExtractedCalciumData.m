function [ExtractedResults, nameExtract]=LoadExtractedCalciumData(MainPath,expPaths,FreeCAD,nDownsample,tagLowRate)

LED = questdlg('Do you want to load extracted calcium data ?','Extracted data Loading','Yes','No','Yes');
switch LED
    case 'Yes'
        [file,path] =uigetfile([MainPath '/*.mat'],'Select Calcium results file to load');
        ExtractedResults=load(fullfile(path,file));
        ExtractedResults=ExtractedResults.results;
        nameExtract = path;
    case 'No'
        [ExtractedResults,~,nameExtract]=ExtractExpData(expPaths,MainPath,FreeCAD,nDownsample,tagLowRate);
end




