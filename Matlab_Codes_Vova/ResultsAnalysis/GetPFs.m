function [Neurons]=GetPFs(Experiment,Neurons,Params,Direction)

Occupation=Experiment.(Direction).TimeSpent(:,Experiment.(Direction).MvtZone)>0;
for neuron = 1:Experiment.nNeurons
    thrPF=Params.critPCs.thrPFs*mean(Neurons.SummedPlaceActivity(neuron).(Direction));
  
        % method 1 : only 1 PF per neuron 
%        bounds=diff(Neurons.SummedPlaceActivity(neuron).(Direction)>thrPF);
%         PFslims=find(bounds);
%         if bounds(find(bounds,1,'first'))==-1; PFslims=[1 PFslims]; end
%         if bounds(find(bounds,1,'last'))==1; PFslims = [PFslims sum(Experiment.(Direction).MvtZone)]; end
%         PFslims=reshape(PFslims,2,[]);
% 
%         Max=find(Neurons.SummedPlaceActivity(neuron).(Direction)==max(Neurons.SummedPlaceActivity(neuron).(Direction)),1,'last');
%         BinMax=find(PFslims(1,:)<=Max,1,'last');
%         Neurons.PFs(neuron).(Direction)=PFslims(1,BinMax):PFslims(2,BinMax);
% 
%         Neurons.PassagePFs(neuron).(Direction)=find(sum(Occupation(:,Neurons.PFs(neuron).(Direction)),2)>numel(Neurons.PFs(neuron).(Direction))*Params.fractPassagePF);

        
        % method 2 : multiple PFs per neuron
       % nBins = size(Neurons.SummedPlaceActivity(neuron).(Direction),2); 
       % x = 1:1:nBins;
        
       Neurons.PFs(neuron).(Direction) = cell(Params.maxNumberPFs,1); 
       Neurons.PassagePFs(neuron).(Direction) = cell(Params.maxNumberPFs,1);
       ActivitySpace = Neurons.SummedPlaceActivity(neuron).(Direction);

       for kPF = 1 : Params.maxNumberPFs
       
           xmax = findpeaks_chronux(ActivitySpace, thrPF); 
           
           if size(xmax.loc, 1)>0
                % trouver la position du pic maxPos, et sa valeur maxVal
               %[maxVal,imax] = max(Neurons.SummedPlaceActivity(neuron).(Direction)(xmax.loc));
               [maxVal,imax] = max(ActivitySpace(xmax.loc));
               maxPos = xmax.loc(imax);

               % on applique la methode precedente pour trouver les PFs potentiels
               %bounds=diff(Neurons.SummedPlaceActivity(neuron).(Direction)>maxVal *Params.sizePFsThr);
               bounds=diff(ActivitySpace>maxVal *Params.sizePFsThr);
                PFslims=find(bounds);
                if bounds(find(bounds,1,'first'))==-1; PFslims=[1 PFslims]; end
                if bounds(find(bounds,1,'last'))==1; PFslims = [PFslims sum(Experiment.(Direction).MvtZone)]; end
                PFslims=reshape(PFslims,2,[]);

                % parmi ceux la celui qui nous interesse contient maxPos
                % c'est le plus petit PFslims(1,:) qui est supérieur a maxPos
                diffValues = PFslims(1,:) - maxPos;
                diffValues(diffValues > 0) = -inf;
                [~, BinMax] = max(diffValues);
                Neurons.PFs(neuron).(Direction)(kPF) = { [PFslims(1,BinMax):PFslims(2,BinMax)]};
               
                Neurons.PassagePFs(neuron).(Direction)(kPF)={find(sum(Occupation(:,cell2mat(Neurons.PFs(neuron).(Direction)(kPF))),2)>numel(cell2mat(Neurons.PFs(neuron).(Direction)(kPF)))*Params.fractPassagePF)};
         
                xmin = findpeaks_chronux(-ActivitySpace);
                %  trouver le plus grand minimum local inférieur à
                %  PFslims(1,BinMax) et le plus petit minimum local
                %  supérieur à PFslims(2,BinMax)
                xmin.loc = [1; xmin.loc; sum(Experiment.(Direction).MvtZone)];
                diffValues = xmin.loc - PFslims(2,BinMax);
                diffValues(diffValues < 0) = +inf;
                [xSegmentMax, ~] = min(diffValues) ; 
                xSegmentMax  = min(xSegmentMax + PFslims(2,BinMax), sum(Experiment.(Direction).MvtZone)) ;

                diffValues = xmin.loc - PFslims(1,BinMax);
                diffValues(diffValues > 0) = -inf;
                [xSegmentMin, ~] = max(diffValues) ;
                xSegmentMin  =   max(xSegmentMin+PFslims(1,BinMax), 1);
                
                ActivitySpace(xSegmentMin:xSegmentMax)=0;

           else
               break;
           end
       end
        
         Neurons.PFs(neuron).(Direction) = Neurons.PFs(neuron).(Direction)(~cellfun('isempty',Neurons.PFs(neuron).(Direction)));
         Neurons.PassagePFs(neuron).(Direction) = Neurons.PassagePFs(neuron).(Direction)(~cellfun('isempty',Neurons.PassagePFs(neuron).(Direction)));
         
end

end
