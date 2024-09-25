root = '/export/home1/MSU_data_Vitya/wvt_default_thr_startev/';
mouse_ses = 'WorkSpace_CC_H01_1D_Features.mat';
load ([root, mouse_ses])
i = 10;
traces = csvread('/export/home1/MSU_data_Vitya/CC_H01_1D_traces.csv',1,1);
scale_fact = 10;
colors = jet(scale_fact);
spike_in_locomotion = find(file_NV(:,i).*velcam');
spike_t_good = find(file_NV(:,i));
trace = traces(:, i);
trace = (trace - min(trace))/(max(trace) - min(trace));
if length(spike_t_good) >= min_spike
    figure('Position', [1 1 Screensize(4) Screensize(4)]);
    hold on
    %axis([ax_xy(1) ax_xy(2) ax_xy(3) ax_xy(4)]);
    title(sprintf('Trajectory of mouse with n = %d (%d) Ca2+ events (in mov, red) of cell #%d', length(spike_t_good),length(spike_in_locomotion), i), 'FontSize', FontSizeTitle);
    xlabel('X coordinate, sm','FontSize', FontSizeLabel);ylabel('Y coordinate, sm','FontSize', FontSizeLabel);
    for col = 1:scale_fact
        idx1 = find(trace >= (col-1)/scale_fact);
        idx2 = find(trace <= col/scale_fact);
        idx = intersect(idx1, idx2);
        hold on, scatter(x_int_sm(idx),y_int_sm(idx), 5*ones(size(idx)), colors(col,:), 'filled')
    end     
    

    scatter(x_vel_b(spike_t_good),y_vel_b(spike_t_good), 40*ones(size(spike_t_good)), 'h', 'k', 'filled');
    scatter(x_vel_b(spike_in_locomotion),y_vel_b(spike_in_locomotion), 60*ones(size(spike_in_locomotion)), 'h', 'r', 'filled');
    set(gca, 'FontSize', FontSizeLabel);
    %saveas(h, sprintf('%s\\Single_spikes\\%s_single_spikes_%d.png',path,FilenameOut,i));
    %delete(h);
end



  %% common fields stat
root = {'/export/home1/MSU_data_Vitya/MAT_3sigma/', '/export/home1/MSU_data_Vitya/wvt_default_thr_startev/'};
FieldStat = [];

for i = 1:2
    filenames = dir([root{i}, '*.mat']);
    FieldStat = [FieldStat, struct('CellIndex',[],'CellNumberAll',[],'CellNumberActivePercent',[],'FiringRate',[],'FiringRateMean',[],'PlaceCellNumberPercent',[],'CellICAll',[],'CellICAllMean',[],'CellICPlaceCells',[],'FieldsNumber',[],'FieldsNumberPerNeuron',[],'FieldsNumberPerNeuronMean',[],'FieldsSquare',[],'FieldsSquareMean',[])];
    
    for mouse = 1:length(filenames)   
        load([root{i}, filenames(mouse).name], 'Cell_IC', 'FieldsIC','test_zone5','MapFieldsIC', 'TimeSession', 'bin_size');
        % creation struct
        FieldStat(i).CellIndex = [FieldStat(i).CellIndex; ones(size(Cell_IC,2),1)*mouse];
        FieldStat(i).CellNumberAll = [FieldStat(i).CellNumberAll; sum(Cell_IC(7,:) > 5)];
        FieldStat(i).CellNumberActivePercent = [FieldStat(i).CellNumberActivePercent; sum(Cell_IC(7,:) > 5)/size(Cell_IC,2)*100];
        FieldStat(i).FiringRate = [FieldStat(i).FiringRate; (Cell_IC(7,:)/TimeSession*60)'];
        FieldStat(i).FiringRateMean = [FieldStat(i).FiringRateMean; mean(Cell_IC(7,:)/TimeSession*60)];
        FieldStat(i).PlaceCellNumberPercent = [FieldStat(i).PlaceCellNumberPercent; sum(Cell_IC(2,:))/FieldStat(i).CellNumberAll(mouse)*100];
        FieldStat(i).CellICAll = [FieldStat(i).CellICAll; Cell_IC(6,~isnan(Cell_IC(6,:)))'];        
        FieldStat(i).CellICAllMean = [FieldStat(i).CellICAllMean; mean(Cell_IC(6,~isnan(Cell_IC(6,:))))];
        FieldStat(i).CellICPlaceCells = [FieldStat(i).CellICPlaceCells; mean(Cell_IC(6,Cell_IC(2,:)==1))];
        FieldStat(i).FieldsNumber = [FieldStat(i).FieldsNumber; size(FieldsIC,2)];
        FieldStat(i).FieldsNumberPerNeuron = [FieldStat(i).FieldsNumberPerNeuron; test_zone5(test_zone5>0)'];
        FieldStat(i).FieldsNumberPerNeuronMean = [FieldStat(i).FieldsNumberPerNeuronMean; mean(test_zone5(test_zone5>0))];
        FieldStat(i).FieldsSquare = [FieldStat(i).FieldsSquare; (sum((reshape(MapFieldsIC, [], size(MapFieldsIC, 3))>0)).*bin_size^2/100)']; % in dm^2
        FieldStat(i).FieldsSquareMean = [FieldStat(i).FieldsSquareMean; mean(sum((reshape(MapFieldsIC, [], size(MapFieldsIC, 3))>0)).*bin_size^2/100)];
        
        clear 'Cell_IC' 'FieldsIC' 'SpikeFieldsReal' 'test_zone5' 'MapFieldsIC' 'x_int_sm' 'y_int_sm' 'TimeSession' 'bin_size';
     end
end

f = fields(FieldStat(1));
for i=1:length(f)
    figure, hold on
    title(f{i})
    histogram(FieldStat(1).(f{i}))
    histogram(FieldStat(2).(f{i}))
end