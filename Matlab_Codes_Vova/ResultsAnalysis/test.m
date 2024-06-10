function test(SuggestList)
% f=figure;
% aa=axes('Parent',f);
% 
% a=triu(ones(7));
% imagesc(a)
% hold on
% for roi=1%:10%size(results.C_raw,1)
%     [~,c]=contour(a,[1,1],'LineColor','r', 'linewidth', 1,'Fill','on','FaceColor','r','ButtonDownFcn',@test2,'Tag','t');
% %     set(c.FacePrims,'ColorType','truecoloralpha');
% %     set(c.FacePrims,'ColorData',uint8([255 0 0 100]'));
%     drawnow
%     set(c.FacePrims,'ColorType','truecoloralpha');
%     set(c.FacePrims,'ColorData',uint8([255 0 0 100]'));
% %     c.FacePrims.ColorType='truecoloralpha';
% %     c.FacePrims.ColorData(4) = 1;
% end
% 
% ff=figure;
% aaa=axes('Parent',ff);
% a=tril(ones(7));
% imagesc(a)
% hold on
% cc=copyobj(findobj(aa,'Type','Contour','Tag','t'),aaa);
f=uifigure;
uit = uitable(f,'Data',SuggestList);
uit.CellSelectionCallback=@(~,evt)test2(evt,SuggestList);
function test2(evt,SuggestList)
    pair = table2array(SuggestList(evt.Indices(1),:));
    disp(pair)
end
end