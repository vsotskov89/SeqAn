function [Experiment,Neurons,Params]=FindPlaceCells(Experiment,Neurons,Params,results)

%% Get neurons characteristics and Find PCs
    % Define PFs
    disp('    Defining Neurons putative PFs : in progress')
        Neurons=GetPFs(Experiment,Neurons,Params,'Dir1');
        Neurons=GetPFs(Experiment,Neurons,Params,'Dir2');
    fprintf([repmat('\b',1,12) '%s\n'],'Done')

    % Calculate characteristics
    disp('    Calculating characteristics : in progress')
        Neurons=GetChars(Experiment,Neurons,'Dir1',Params);
        Neurons=GetChars(Experiment,Neurons,'Dir2',Params);
    fprintf([repmat('\b',1,12) '%s\n'],'Done')

    
    % Get PC status
    disp('    Checking PCs status : in progress')
        [Neurons,Experiment]=GetPCstatus(Experiment,Neurons,Params,'Dir1');
        [Neurons,Experiment]=GetPCstatus(Experiment,Neurons,Params,'Dir2');
    fprintf([repmat('\b',1,12) '%s\n'],'Done')
    
% %% Create Summary
%     if Params.PrintPDFs
%         disp('    Printing summary PDFs : neuron 001')
%             plotSummaries(Experiment,Neurons,Params,results);
%         fprintf([repmat('\b',1,11) '%s\n'],'Done')
%     end
%     
%     disp('    Plotting Averages and saving : in progress')
%          plotAvg(Experiment,Neurons,Params,results)
%     fprintf([repmat('\b',1,12) '%s\n'],'Done')