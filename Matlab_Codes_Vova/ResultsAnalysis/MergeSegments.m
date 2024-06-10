function Experiment=MergeSegments(Experiment)



segs=[[Experiment.Dir1.Starts Experiment.Dir1.Ends ones(numel(Experiment.Dir1.Starts),1)] ;....
    [Experiment.Dir2.Starts Experiment.Dir2.Ends 2*ones(numel(Experiment.Dir2.Starts),1)] ];
[~,idx]=sort(segs(:,1));
segs=segs(idx,:);

seg=1;
while seg ~=size(segs,1)
    if segs(seg,3)==segs(seg+1,3)
       segs(seg,2)=segs(seg+1,2);
       segs(seg+1,:)=[];
    else
        seg=seg+1;
    end
end

Experiment.Dir1.Starts=segs(segs(:,3)==1,1);
Experiment.Dir1.Ends=segs(segs(:,3)==1,2);
Experiment.Dir2.Starts=segs(segs(:,3)==2,1);
Experiment.Dir2.Ends=segs(segs(:,3)==2,2);

Experiment.Dir1.NbPassage=numel(Experiment.Dir1.Starts);
Experiment.Dir2.NbPassage=numel(Experiment.Dir2.Starts);


Experiment.Dir1.Segments=zeros(size(Experiment.Dir1.Segments  ,1),1);                                     
idxDir=arrayfun(@(B,C)B:C,Experiment.Dir1.Starts,Experiment.Dir1.Ends,'UniformOutput',0);               % Get idx of mouvement periods
for seg = 1 : numel(idxDir)
     Experiment.Dir1.Segments (idxDir{seg}) = 1;                                            % Rebuild segments Mat with new periods
end
Experiment.Dir2.Segments=zeros(size(Experiment.Dir2.Segments  ,1),1);                                     
idxDir=arrayfun(@(B,C)B:C,Experiment.Dir2.Starts,Experiment.Dir2.Ends,'UniformOutput',0);               % Get idx of mouvement periods
for seg = 1 : numel(idxDir)
     Experiment.Dir2.Segments (idxDir{seg}) = 1;                                            % Rebuild segments Mat with new periods
end














