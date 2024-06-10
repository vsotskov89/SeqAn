function expPaths=GetExpsPaths(varargin)
% Use example : expPaths=GetExpsPaths('NbExp',[0,1,1],'Path','/mnt/cortex-data-312/tmpTheo/Data/Souris154147')
% Advice : run it once and save all the path in one file once and for all.


MainPath='/mnt/cortex-data-312/tmpTheo/Data/Souris154147';
NbExp=0;
ExpPhases={'SleepPRE','Behavior', 'SleepPOST'};
i=1;
while i <= length(varargin)
	if ~ischar(varargin{i})
		error(['Parameter ' num2str(i) ' is not a property']);
	end
	switch(lower(varargin{i}))
		case 'nbexp'
			NbExp = varargin{i+1};
			if ~isdvector(NbExp,'#3','>=0')
				error('Incorrect value for ''NbExp'' ; It must be a 1x3 double vector');
			end
			varargin = {varargin{[1:(i-1) (i+2):length(varargin)]}};
        case 'path'
            MainPath = varargin{i+1};
            if ~isastring(MainPath)
                error('Incorrect value for ''Path'' ; It must be a string');
            end
            varargin = {varargin{[1:(i-1) (i+2):length(varargin)]}};
        otherwise
			i = i+2;
	end
end
if isequal(NbExp,0)
    prompt = {'Number of Sleep PRE recordings : ','Number of Behavioral recordings :','Number of Sleep POST recordings :'};
    dlgtitle = 'Number of experiments'; dims = [1 35]; definput = {'0','1','1'}; 
    NbExp = inputdlg(prompt,dlgtitle,dims,definput); NbExp = str2num(cell2mat(NbExp));
end

LoadPaths = questdlg('Do you want to load previous paths ?','Previous Paths Loading','Yes','No','Yes');
switch LoadPaths
    case 'Yes'
        [file,path] =uigetfile([MainPath '/*.mat'],'Select Paths file to load');
        expPaths=load(fullfile(path,file));expPaths=expPaths.expPaths;
    case 'No'
        expPaths=[];
        for phase=1:numel(ExpPhases)
            expPaths.(ExpPhases{phase}).TimingsCMOS={};expPaths.(ExpPhases{phase}).TimingsBasler={};expPaths.(ExpPhases{phase}).Position={};expPaths.(ExpPhases{phase}).ePhy={};
            for nExp=1:NbExp(phase)
                expPaths.(ExpPhases{phase}).TimingsCMOS=[expPaths.(ExpPhases{phase}).TimingsCMOS ; uigetdir(MainPath,['Select sCMOS Timings folder for ' ExpPhases{phase} ' ' num2str(nExp)])];
                expPaths.(ExpPhases{phase}).TimingsBasler=[expPaths.(ExpPhases{phase}).TimingsBasler ; uigetdir(MainPath,['Select Basler Timings folder for ' ExpPhases{phase} ' ' num2str(nExp)])];
                expPaths.(ExpPhases{phase}).Position=[expPaths.(ExpPhases{phase}).Position ; uigetdir(MainPath,['Select Position Data folder for ' ExpPhases{phase} ' ' num2str(nExp)])];
                expPaths.(ExpPhases{phase}).ePhy=[expPaths.(ExpPhases{phase}).ePhy ; uigetdir(MainPath,['Select ePhy Data folder for ' ExpPhases{phase} ' ' num2str(nExp)])];
            end
        end
        expPaths.Concatenated.Calcium={};
        expPaths.Concatenated.Calcium=[expPaths.Concatenated.Calcium ; uigetdir(MainPath,'Select Concatenated Calcium Data folder')];
        SavePaths= questdlg('Do you want to save new paths ?','New Paths Saving','Yes','No','Yes');
        switch SavePaths
            case 'Yes'
                uisave('expPaths',[MainPath '/expPaths'])
        end
end
end