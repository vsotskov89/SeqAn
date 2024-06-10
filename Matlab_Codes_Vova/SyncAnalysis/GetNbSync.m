function [syncEvtsA,nbSyncEvtsA]=GetNbSync(eventsA,eventsB,MaxSep)

syncEvtsA=zeros(1,numel(eventsA));
 for evtA=1:numel(eventsA)
     prevEvtB=find(eventsB<eventsA(evtA),1,'last');
     nextEvtB=find(eventsB>eventsA(evtA),1,'first');
     if isempty(prevEvtB)
         closestEvtB=abs(eventsA(evtA)-eventsB(nextEvtB));
     elseif isempty(nextEvtB)
         closestEvtB=abs(eventsA(evtA)-eventsB(prevEvtB));
     else
         closestEvtB=min(abs(eventsA(evtA)-eventsB(prevEvtB)),abs(eventsA(evtA)-eventsB(nextEvtB)));
     end
     
     if closestEvtB<MaxSep
         syncEvtsA(evtA)=1;
     end
 end
 nbSyncEvtsA=sum(syncEvtsA);
%  disp(sum(syncEvtsA))
%  disp(sum(syncEvtsA)/numel(syncEvtsA))