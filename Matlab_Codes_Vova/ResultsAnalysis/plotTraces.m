function plotTraces(neuron,Neurons,results)

sampleAct=[(1:numel(results.C_raw(neuron,:)))' results.C_raw(neuron,:)'];                                                        % Data to FMA Sample
FiltAct = FilterLFP(sampleAct,'passband',12.5*[0 10]); FiltAct=FiltAct(:,2:end)';
FiltTheta = FilterLFP(sampleAct,'passband',12.5*[0 10]); FiltTheta=FiltTheta(:,2:end)';


risesRaw=results.C_raw(neuron,:); risesRaw(~Neurons.Rises(neuron).Matrix)=nan;
thetaRaw=results.C_raw(neuron,:); thetaRaw(~Neurons.Theta(neuron).Matrix)=nan;

risesFilt=FiltAct; risesFilt(~Neurons.Rises(neuron).Matrix)=nan;
thetaFilt=FiltAct; thetaFilt(~Neurons.Theta(neuron).Matrix)=nan;

f=figure; hold on
t=(1:numel(results.C_raw(neuron,:)))*0.01;
p1 = plot(t,results.C_raw(neuron,:));
    t1 = plot(t,thetaRaw,'LineWidth',2.5,'Visible','off');
    e1 = plot(t,risesRaw,'LineWidth',2.5,'Visible','off');
p2 = plot(t,FiltAct,'Visible','off');
    t2 = plot(t,thetaFilt,'LineWidth',2.5,'Visible','off');
    e2 = plot(t,risesFilt,'LineWidth',2.5,'Visible','off');
% p3 = plot(risesRaw,'LineWidth',2.5,'Visible','off');
%     t3
%     e3

% Controls =uifigure('Name','Controls','Units','Normalized','Position',[0.02 0.1 0.2 0.8]);
Controls =uifigure('Name','Controls','Position',[50 100 300 700]); test=[250 650 250 650]; test2=[.8 .15 .8 .15].*test;
    pan1 = uipanel('Parent',Controls,'Title','Raw Trace','FontSize',12,'Position',[.1 .85 .8 .15].*test);
        cbx1=uicheckbox(pan1, 'Text','Visible','Value', 1,'Position',[.05 .45 .4 .3].*test2);
        cbxE1=uicheckbox(pan1, 'Text','Events','Value', 0,'Position',[.5 .45 .4 .3].*test2,'ValueChangedFcn',@(cbx,evt)EventB(cbx,evt,e1,cbx1));
        cbxT1=uicheckbox(pan1, 'Text','Theta','Value', 0,'Position',[.5 .1 .4 .3].*test2,'ValueChangedFcn',@(cbx,evt)ThetaB(cbx,evt,t1,cbx1));
        cbx1.ValueChangedFcn=@(cbx,evt)VisibleB(cbx,evt,p1,e1,t1,cbxE1,cbxT1);

    pan2 = uipanel('Parent',Controls,'Title','Filtered Trace','FontSize',12,'Position',[.1 .65 .8 .15].*test);
        cbx2=uicheckbox(pan2, 'Text','Visible','Value', 0,'Position',[.05 .45 .4 .3].*test2,'ValueChangedFcn',@(cbx,evt)VisibleB(cbx,evt,p2));
        cbxE2=uicheckbox(pan2, 'Text','Events','Value', 0,'Position',[.5 .45 .4 .3].*test2,'ValueChangedFcn',@(cbx,evt)EventB(cbx,evt,e2,cbx2));
        cbxT2=uicheckbox(pan2, 'Text','Theta','Value', 0,'Position',[.5 .1 .4 .3].*test2,'ValueChangedFcn',@(cbx,evt)ThetaB(cbx,evt,t2,cbx2));
        cbx2.ValueChangedFcn=@(cbx,evt)VisibleB(cbx,evt,p2,e2,t2,cbxE2,cbxT2);
%     pan3 = uipanel('Parent',Controls,'Title','Theta Trace','FontSize',12,'Position',[.1 .45 .8 .15].*test);
%         cbx3=uicheckbox(pan3, 'Text','Visible','Value', 0,'Position',[.05 .45 .4 .3].*test2,'ValueChangedFcn',@(cbx3,evt)VisibleB(cbx3,evt,p3));
        
        
function VisibleB(cbx,~,p,e,t,cE,cT)
    p.Visible=cbx.Value;
    e.Visible=min(cbx.Value,cE.Value);
    t.Visible=min(cbx.Value,cT.Value);
end
function EventB(cbx,~,e,c)
    e.Visible=min(cbx.Value,c.Value);
end
function ThetaB(cbx,~,t,c)
    t.Visible=min(cbx.Value,c.Value);
end
end