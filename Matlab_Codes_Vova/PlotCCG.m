function PlotCCG(tCCG, RipSceCCG, permCCG)
    figure; 
    bar(tCCG,RipSceCCG(:,1,2)); 
    xlabel('Delay (s)'); ylabel('Occurence (#)');title('CCG of Ripples x SCEs') 
    hold on
    plot(tCCG,prctile(permCCG,95),'color','r','LineWidth',2)
end
    
