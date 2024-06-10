function [c_out, b0] = remove_baseline_sliding_window(c_in, p, kwindow)

% remove baseline using a sliding window and computing the baseline as the
% pth percentile on this window, of size kwindow

%kwindow = 1000; % taille de la fenetre glissante en nb de points
%p = 10; %10th percentile

%tic;

Ntime = size(c_in,2);
if kwindow > Ntime
    kwindow = Ntime;
end

kwindow =  10 * (floor(kwindow/10)); % on veut une fen�tre multiple de 10

halfwindow = kwindow/2;
c0 = NaN(1, Ntime);

% avec une sliding window qui se deplace d'un element a la fois
% for ktime = halfwindow+1 : Ntime - halfwindow
%            sort_moy_roi = sort(c_in( ktime - halfwindow : ktime + halfwindow -1)) ; 
%            c0(ktime) = sort_moy_roi(floor(p/100*kwindow)); %10th percentile
% end
% c0( Ntime - halfwindow +1 : end) = c0( Ntime - halfwindow) * ones(1, halfwindow);
% c0( 1: halfwindow) = c0( halfwindow+1 ) * ones(1, halfwindow);

% avec une sliding window qui se deplace d'un dizieme de window a la fois

for ktime = halfwindow + 1 : kwindow/10 : Ntime - halfwindow +1
           %sort_moy_roi = sort(c_in( ktime - halfwindow : ktime + halfwindow -1)) ; 
           %c0(ktime) = sort_moy_roi(floor(p/100*kwindow)); %pth percentile
           c0(ktime) = prctile(c_in( ktime - halfwindow : ktime + halfwindow -1),p);
end
c0 = fillmissing(c0,'linear','EndValues','nearest');
b0 = mean(c0);

c_out = c_in - c0;

% soustraction supplementaire d'une baseline constante prise au 40th
% percentile sur la totalit� de la trace

sort_moy_roi = sort(c_out);
c_out = c_out - sort_moy_roi(floor(40/100*Ntime));

%toc; 
