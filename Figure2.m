%%%%%%%%%%%%%%  plot subplot1 group level fc-ECs %%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close gcf
load('./NB_2back.mat');
groups = {'HC', 'SCZ'};
groups_new = {'HV', 'SCZ'};
%% plot
HCP_ROI = readtable('./HCP_ROI_label.csv');
direction = {'forward', 'backward'};
fc = NB_2back.fc_name(:,2:3);
roi = unique(fc);
roi_new = [];
idx_final = [];
for i = 1:length(roi)
    tmp = roi{i};
    tmp1 = strfind(tmp,'_');
    tmp = tmp(tmp1(end) + 1:end);
    roi_new = [roi_new; {tmp}];
    if ismember(tmp, HCP_ROI.CorticalArea)
        idx = find(ismember(HCP_ROI.CorticalArea, tmp));
        HCP_ROI.CorticalArea_new{idx,1} = roi{i};
        idx_final = [idx_final; idx];
    end
end
HCP_ROI = HCP_ROI(idx_final,:);
HCP_ROI = sortrows(HCP_ROI, 'RegionNumber');

roi = HCP_ROI.CorticalArea_new;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp = unique(HCP_ROI.RegionNumber);
YVarLabels = [];
chk = [];
for i = 1:length(tmp)
    idx = find(ismember(HCP_ROI.RegionNumber, tmp(i)));
    if length(idx) == 1
        tmp_label = HCP_ROI.RegionName(idx(1));
    elseif length(idx) == 2
        tmp_label = [HCP_ROI.RegionName{idx(1)}; {''}];
    elseif length(idx) == 3
        tmp_label = [{''}; HCP_ROI.RegionName{idx(1)}; {''}];
    elseif length(idx) > 3 && mod(length(idx),2) == 1
        tmp_label = [repmat({''}, (length(idx)-1)/2, 1); HCP_ROI.RegionName{idx(1)}; repmat({''}, (length(idx)-1)/2, 1)];
    elseif length(idx) > 3 && mod(length(idx),2) == 0
        tmp_label = [repmat({''}, length(idx)/2 - 1, 1); HCP_ROI.RegionName{idx(1)}; repmat({''}, length(idx)/2, 1)];
    end
    YVarLabels = [YVarLabels; tmp_label];
    clear tmp_label
end

region_color = {
    'navy',[0,0,128]/255;
    'violet',[238,130,238]/255;
    'saddlebrown',[139,69,19]/255;
    'navajowhite',[255,222,173]/255;
    'hotpink',[255,105,180]/255;
    'lightskyblue',[135,206,250]/255;
    'lime',[0,255,0]/255;
    'lightslategray',[119,136,153]/255;
    'purple',[128,0,128]/255;
    'olivedrab',[107,142,35]/255;
    'green',[0,128,0]/255;
    'burlywood',[222,184,135]/255;
    'tomato',[255,99,71]/255;
    'maroon',[128,0,0]/255;
    'darkcyan',[0,139,139]/255;
    'cyan',[0,255,255]/255;
    'orange',[255,165,0]/255;
    'gold',[255,215,0]/255
    'darkslategray',[47,79,79]/255;
    'lightgreen',[144,238,144]/255;
    'royalblue',[65,105,225]/255;
    'red',[255,0,0]/255;
    };

title_loc = [0.35, 0.95; 0.65, 0.95];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hCB = colorbar('north', ...
    'EdgeColor','none', ...
    'TickLength', 0, ...
    'Box', 'off',...
    'FontSize', 16, 'FontWeight', 'bold');
hCB.Position = [0.54 0.03 0.2 0.03];

axisleft = [0.285 0.6 0.1 0.1 + 0.21 0.1 + 0.21*2 0.1 + 0.21*3 0.1 + 0.21*4];
axisbottom = [0.55 0.53 0.15 0.15 0.15 0.15];
axiswidth = [0.38 0.38 0.15 0.15 0.15 0.15];
axisheight = [0.38 0.38 0.3 0.3 0.3 0.3];

set(gcf,'units','normalized','position',[0 0 0.75 1])
set(gcf,'units','pixels')
gcf_position = get(gcf,'position');
figure_width = gcf_position(3);
figure_height = gcf_position(4);


for dc_idx = 1:length(direction)
    for fc_idx = 1:length(fc)
        if strcmp(direction{dc_idx}, 'forward')
            row_idx = find(contains(roi,fc{fc_idx,1}));
            col_idx = find(contains(roi,fc{fc_idx,2}));
            NVC_frequency(col_idx, row_idx) = NB_2back.HC_B_forward_NVC_frequency(fc_idx);
        elseif strcmp(direction{dc_idx}, 'backward')
            row_idx = find(contains(roi,fc{fc_idx,2}));
            col_idx = find(contains(roi,fc{fc_idx,1}));
            NVC_frequency(col_idx, row_idx) = NB_2back.HC_B_backward_NVC_frequency(fc_idx);
        end
    end
end

ph(1) = subplot('position', [axisleft(1), axisbottom(1), ...
    axiswidth(1)*figure_height/figure_width, axisheight(1)]);

imagesc(NVC_frequency);
roi_new = [];
idx_final = [];
for i = 1:length(roi)
    tmp = roi{i};
    tmp1 = strfind(tmp,'_');
    tmp = tmp(tmp1(end) + 1:end);
    roi_new = [roi_new; {tmp}];
    if ismember(tmp, HCP_ROI.CorticalArea)
        idx = find(ismember(HCP_ROI.CorticalArea, tmp));
        HCP_ROI.CorticalArea_new{idx,1} = roi{i};
        idx_final = [idx_final; idx];
    end
end
HCP_ROI = HCP_ROI(idx_final,:);
HCP_ROI = sortrows(HCP_ROI, 'RegionNumber');
% plot line to segment big regions
tmp = unique(HCP_ROI.RegionNumber);
for i = 1:length(tmp)
    idx_tmp = find(ismember(HCP_ROI.RegionNumber, tmp(i)));
    idx(i,1) = idx_tmp(1);
end
idx = [idx; length(roi)];
% axis off;
xline(idx(2:(end-1)) - 0.5,'k-.');%, 'Color', [0, 0, 0]);
yline(idx(2:(end-1)) - 0.5,'k-.');%, 'Color', [0, 0, 0]);
set(gca, 'XTick', {}, 'YTick', {}, 'TickLength',[0 0])
axis square;
color_idx = jet(1200);
colormap([1, 1, 1; ...
    color_idx(1:end,:)]);


%     pos_cb =  get(cb1,'Position'); % get the position of the current figure
pos_axis =  get(gca,'Position'); % get the position of the current axis
n_area = size(HCP_ROI,1);
tmp = unique(HCP_ROI.RegionNumber);
sp = pos_axis(2);
left_start = axisleft(1) - 0.0104;
cube_width = axisheight(1)/20;

% plot color cubes on y axis
for i = length(tmp):-1:1
    idx = find(ismember(HCP_ROI.RegionNumber, tmp(i)));
    annotation('rectangle',...
        [left_start, sp, cube_width, pos_axis(4)/n_area*length(idx)],...
        'LineStyle', 'None', ...
        'FaceColor',region_color{HCP_ROI.RegionNumber(idx(1)),2});
    xl = [left_start - 0.01, left_start - 0.005];
    yl = [sp + pos_axis(4)/n_area*length(idx)/2, sp + pos_axis(4)/n_area*length(idx)/2];
    annotation('line',xl,yl, 'LineWidth', 1.5, ...
        'Color',region_color{HCP_ROI.RegionNumber(idx(1)),2})
    sp = sp + pos_axis(4)/n_area*length(idx);
end

% plot color cubes on x axis
sp = left_start + cube_width*1.132; % start point
for i = 1:length(tmp)
    idx = find(ismember(HCP_ROI.RegionNumber, tmp(i)));
    annotation('rectangle', ...
        [sp, pos_axis(2) + pos_axis(4) + cube_width*0.25, pos_axis(3)*0.92/n_area*length(idx), cube_width*1.7],...
        'LineStyle', 'None', ...
        'FaceColor',region_color{HCP_ROI.RegionNumber(idx(1)),2});
    sp = sp + pos_axis(3)*0.92/n_area*length(idx);
end

% % title
% annotation('textbox', ...
%     [title_loc(1,:), ...
%     0.2, 0.04],...
%     'String', ['group level fc-EC of ' groups_new{1}], ...
%     'LineStyle', 'None', ...
%     'FontSize', 14, 'FontWeight', 'normal');

% add annotations to colors
tmp = unique(HCP_ROI.RegionNumber);
sp = pos_axis(2);
sp1 = pos_axis(2);
for i = length(tmp):-1:1
    idx = find(ismember(HCP_ROI.RegionNumber, tmp(i)));
    annotation('textbox',...
        [0 sp + pos_axis(4)/length(tmp)/1.8 0.255 pos_axis(3)/length(tmp)],...
        'String',HCP_ROI.RegionName(idx(1)), ...
        'HorizontalAlignment', 'right', ...
        'LineStyle', 'none', ...
        'FontSize', 11, 'FontWeight', 'normal');

    xl = [left_start - 0.01 left_start - 0.015];
    yl = [sp1 + pos_axis(4)/n_area*length(idx)/2, sp + pos_axis(4)/length(tmp)/2];
    annotation('line',xl,yl, 'LineWidth', 1.5, ...
        'Color',region_color{HCP_ROI.RegionNumber(idx(1)),2})

    x2 = [left_start - 0.015 left_start - 0.02];
    y2 = [sp + pos_axis(4)/length(tmp)/2, sp + pos_axis(4)/length(tmp)/2];
    annotation('line',x2,y2, 'LineWidth', 1.5, ...
        'Color',region_color{HCP_ROI.RegionNumber(idx(1)),2})

    sp = sp + pos_axis(4)/length(tmp);
    sp1 = sp1 + pos_axis(4)/n_area*length(idx);
end

cb = colorbar('southoutside', ...
    'TickLength', 0.03, ...
    'FontSize', 11, 'FontWeight', 'normal');
cb.Label.String = 'WM-related fc-EC to DLPFC in HV (Hz)';
cb.Label.FontSize = 11;
cb.Label.Rotation = 0;
cb.Ruler.Color = 'k';
cb.Position(1:2) = [0.645, 0.1];
set(cb, 'EdgeColor', 'none');
cb.Position(4) = cb.Position(4) * 1.25;

cmap = colormap;
cmap = struct('b', num2cell(cmap([1, (10:10:120) * 10 + 1], 3) .* 255), ...
    'g', num2cell(cmap([1, (10:10:120) * 10 + 1], 2) .* 255), ...
    'r', num2cell(cmap([1, (10:10:120) * 10 + 1], 1) .* 255), ...
    'val', num2cell([0, 10:10:120])');
cmap = jsonencode(cmap, PrettyPrint=true);

% 将 JSON 字符串写入文件
fileID = fopen(['.\brain_map_FS\cmap'], 'w');
fprintf(fileID, cmap);
fclose(fileID);

%%%%%%%%%%%%%%%%%%% plot subplot2 hierachical clustring %%%%%%%%%%%%%%%%%%%

clear
load('./NB_2back.mat');
HCP_ROI = readtable('./HCP_ROI_label.csv');
direction = {'forward', 'backward'};
fc = NB_2back.fc_name(:,2:3);
roi = unique(fc);
roi_new = [];
idx_final = [];
for i = 1:length(roi)
    tmp = roi{i};
    tmp1 = strfind(tmp,'_');
    tmp = tmp(tmp1(end) + 1:end);
    roi_new = [roi_new; {tmp}];
    if ismember(tmp, HCP_ROI.CorticalArea)
        idx = find(ismember(HCP_ROI.CorticalArea, tmp));
        HCP_ROI.CorticalArea_new{idx,1} = roi{i};
        idx_final = [idx_final; idx];
    end
end
HCP_ROI = HCP_ROI(idx_final,:);
HCP_ROI = sortrows(HCP_ROI, 'RegionNumber');
ROI_color = {
    'Primary Visual Cortex', 'navy',[0,0,128]/255;
    'Early Visual Cortex', 'violet',[238,130,238]/255;
    'Dorsal Stream Visual Cortex', 'saddlebrown',[139,69,19]/255;
    'Ventral Stream Visual Cortex', 'navajowhite',[255,222,173]/255;
    'MT+ Complex and Neighboring Visual Areas', 'hotpink',[255,105,180]/255;
    'Somatosensory and Motor Cortex', 'lightskyblue',[135,206,250]/255;
    'Paracentral Lobular and Mid Cingulate Cortex', 'lime',[0,255,0]/255;
    'Premotor  Cortex', 'lightslategray',[119,136,153]/255;
    'Posterior Opercular Cortex', 'purple',[128,0,128]/255;
    'Early Auditory Cortex', 'olivedrab',[107,142,35]/255;
    'Auditory Association Cortex', 'green',[0,128,0]/255;
    'Insular and Frontal Opercular Cortex', 'burlywood',[222,184,135]/255;
    'Medial Temporal Cortex', 'tomato',[255,99,71]/255;
    'Lateral Temporal Cortex', 'maroon',[128,0,0]/255;
    'Temporo-Parieto-Occipital Junction', 'darkcyan',[0,139,139]/255;
    'Superior Parietal Cortex', 'cyan',[0,255,255]/255;
    'Inferior Parietal Cortex', 'orange',[255,165,0]/255;
    'Posterior Cingulate Cortex', 'red',[255,215,0]/255;
    'Anterior Cingulate and Medial Prefrontal Cortex', 'darkslategray',[47,79,79]/255;
    'Orbital and Polar Frontal Cortex', 'lightgreen',[144,238,144]/255;
    'Inferior Frontal Cortex', 'royalblue',[65,105,225]/255;
    'DorsoLateral Prefrontal Cortex', 'gold',[255,0,0]/255
    };

a = [NB_2back.HC_B_forward_NVC_time', NB_2back.HC_B_forward_NVC_frequency';
    NB_2back.HC_B_backward_NVC_time', NB_2back.HC_B_backward_NVC_frequency']; % heterogeneity matrix
a(find(a(:,1) > 1000),1) = 1000;

fc = [NB_2back.fc_name(:,2), NB_2back.fc_name(:,3);
    NB_2back.fc_name(:,3), NB_2back.fc_name(:,2)];

idx = find(sum(a,2) == 0);
a(idx,:) = [];
fc(idx,:) = [];


aa(:,1) = zscore(a(:,1));
aa(:,2) = zscore(a(:,2));
b = pdist(aa,'euclidean'); % only frequency or plus time
Z = linkage(b,'average');

% rescale the distance
Z(:,3) = Z(:,3)/(max(Z(:,3))-min(Z(:,3)))*10;


axisw = 0.5;
axish = 0.75;


axisl = 0.6;
axisb = 0.5;

ph(2)=subplot('position', [axisl, axisb, axisw, axish]);



[h,T,perm] = dendrogram(Z,length(aa),'colorthreshold',5);

set(h,'LineWidth',1.5)

%--------------------------------------------------------------------------
%Changing the line colours
lineColours = cell2mat(get(h,'Color'));
colourList = unique(lineColours, 'rows');

myColours = [0.4660 0.6740 0.1880;
    0 0.4470 0.7410;
    0.6350 0.0780 0.1840;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    [14, 102, 85]./255;
    [236, 64, 122]./255];
%// Replace each colour (colour by colour). Start from 2 because the first colour are the "unclustered" black lines
for colour = 2:size(colourList,1)
    %// Find which lines match this colour
    idx = ismember(lineColours, colourList(colour,:), 'rows');
    %// Replace the colour for those lines
    lineColours(idx, :) = repmat(myColours(colour-1,:),sum(idx),1);
end
%// Apply the new colours to the chart's line objects (line by line)
for line = 1:size(h,1)
    set(h(line), 'Color', lineColours(line,:));
end


feet_points = [];
for i=1:size(h)
    xdata = get(h(i),'XData');
    xdata = xdata([1,3]); % only take the feet points of this fork
    color = get(h(i),'Color');
    for j = 1:length(xdata)
        if mod(xdata(j),1) == 0
            feet_points = [feet_points; xdata(j), color];
        end
    end
end
idx_rep = find(hist(feet_points(:,1),unique(feet_points(:,1)))>1)';
for i = 1:length(idx_rep)
    idx_tmp = find(feet_points(:,1) == idx_rep(i));
    feet_points(idx_tmp(2:end),:) = [];
end
[~,idx_sort] = sort(feet_points(:,1), 'ascend');
feet_points = feet_points(idx_sort,:);
turn_point = [];
for i = 2:length(feet_points)
    if ~isequal(feet_points(i,2:4), feet_points(i-1, 2:4))
        turn_point = [turn_point; i-1];
    end
end
idx_cluster_new = [[1; turn_point - 1], [turn_point; length(feet_points)]];

% relabel the new big cluster to the orignal order and get the index
label = get(gca, 'XTickLabel');
label = str2num(label);
for i = 1:length(idx_cluster_new)
    label(idx_cluster_new(i, 1):idx_cluster_new(i, 2),2) = i;
end
[~,idx_sort] = sort(label(:,1),'ascend');
T_idx = label(idx_sort,2);

comp_color = feet_points(idx_cluster_new(:,2),2:end);
%--------------------------------------------------------------------------

%Get x and y ranges
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
minx = xlim(1);
maxx = xlim(2);
miny = ylim(1);
maxy = ylim(2);
xrange = maxx-minx;
yrange = maxy-miny;

%Reshape into a polar plot
for i=1:size(h)
    xdata = get(h(i),'XData');
    ydata = get(h(i),'YData');
    %Rescale xdata to go from pi/12 to 2pi - pi/12
    xdata = (((xdata-minx)/xrange)*(pi*11/6))+(pi/12);
    %Rescale ydata to go from 1 to 0, cutting off lines
    %which drop below the axis limitregion_name
    ydata = max(ydata,miny);
    ydata = 1-((ydata-miny)/yrange);
    %To make horizontal lines look more circular,
    %insert ten points into the middle of the line before
    %polar transform
    newxdata = [xdata(1), linspace(xdata(2),xdata(3),10), xdata(4)];
    newydata = [ydata(1), repmat(ydata(2),1,10), ydata(4)];
    %Transform to polar coordinates
    [xdata,ydata]=pol2cart(newxdata,newydata);
    %Reset line positions to new polar positions
    set(h(i),'XData',xdata);
    set(h(i),'YData',ydata);
end


% Plot ratio of brain regions in each cluster
warning off
for i = 1:length(idx_cluster_new)
    weight_MX_tmp  = cell2table(cell(0,4), 'VariableNames', {'region_name', ...
        'importance', 'afferent', 'efferent'});
    fc_tmp = fc(find(T_idx == i),:);
    fc_tmp = reshape(fc_tmp, size(fc_tmp,1)*size(fc_tmp,2), 1);
    for j = 1:length(fc_tmp)
        idx_tmp = find(ismember(HCP_ROI.CorticalArea_new,fc_tmp{j,1}));
        fc_tmp{j,1} = HCP_ROI.RegionName{idx_tmp};
    end
    region_name = unique(fc_tmp);
    for j = 1:length(region_name)
        region_ratio_tmp(j,1) = length(find(ismember(fc_tmp, region_name{j})))/length(fc_tmp);
        weight_MX_tmp.region_name{j} = region_name{j};
        weight_MX_tmp.importance(j) = region_ratio_tmp(j,1);
    end
    [region_ratio_tmp,idx_tmp] = sort(region_ratio_tmp, 'descend');
    region_name = region_name(idx_tmp);

    t1 = (((idx_cluster_new(i,1) -minx)/xrange)*(pi*11/6))+(pi*1/12);
    t2 = (((idx_cluster_new(i,2) -minx)/xrange)*(pi*11/6))+(pi*1/12);

    t_in = t1:0.001:t2;
    %     t_out = t2:-0.001:t1;

    region_seg_tmp = round(region_ratio_tmp*length(t_in),0);
    s_tmp = 0;
    for j = 1:length(region_seg_tmp)
        region_seg_new_tmp(j, 1) = s_tmp + 1;
        region_seg_new_tmp(j, 2) = s_tmp + region_seg_tmp(j);
        s_tmp = s_tmp + region_seg_tmp(j);
    end
    region_seg_new_tmp(2:end, 1) = region_seg_new_tmp(2:end, 1) - 1;

    % minor modification
    if sum(region_seg_tmp) ~= length(t_in)
        idx = find(region_seg_new_tmp(:,2) > length(t_in));
        if ~isempty(idx)
            region_name(idx(1):end) = [];
            region_seg_new_tmp(idx(1):end,:) = [];
            region_seg_new_tmp(idx(1) - 1,2) = length(t_in);
        elseif isempty(idx)
            region_seg_new_tmp(length(region_seg_new_tmp),2) = length(t_in);
        end
    end


    if i == 6
        plot(1.26*cos(t_in) , ...
            1.26*sin(t_in), ...
            'Color', comp_color(i,:), 'LineWidth', 1.5); 
        % 画最大culster的圆弧
    end


    r1 = 1.15;                                                      % Outer radius
    r = 1.05;                                                     % Inner radius
    for k=1:length(region_name)
        si = t_in(region_seg_new_tmp(k,1)); % Start index of polygon to be filled
        t_in_tmp = t_in(region_seg_new_tmp(k,1):region_seg_new_tmp(k,2));
        t_out_tmp = t_in(region_seg_new_tmp(k,2):-1:region_seg_new_tmp(k,1));
        x=[r*cos(t_in_tmp) r1*cos(t_out_tmp) r*cos(si)];
        y=[r*sin(t_in_tmp) r1*sin(t_out_tmp) r*sin(si)];
        %         plot(x,y,'w-');
        hold on
        idx_roi_color = find(ismember(ROI_color(:,1), region_name{k}));
        if ~isempty(idx_roi_color)
            fill(x,y, ROI_color{idx_roi_color, 3},'edgecolor', ROI_color{idx_roi_color, 3});
            %         elseif isempty(idx_roi_color)
            %             fill(x,y, [17 17 17]/255,'LineColor', [17 17 17]/255);

        end

        % plot ratio of afferent and efferent connections
        %         if i == 1 && k < 20 % only plot the 1st component
        fc_tmp_ae = fc(find(T_idx == i),:);
        for g = 1:size(fc_tmp_ae,1)
            for h = 1:size(fc_tmp_ae,2)
                idx_tmp = find(ismember(HCP_ROI.CorticalArea_new,fc_tmp_ae{g,h}));
                fc_tmp_ae{g,h} = HCP_ROI.RegionName{idx_tmp};
            end
        end

        idx_efferent = length(find(ismember(fc_tmp_ae(:,1), region_name{k})));
        idx_afferent = length(find(ismember(fc_tmp_ae(:,2), region_name{k})));
        r3 = 1.21;                                                      % Outer radius
        r2 = 1.18;                                                     % Inner radius
        idx_roi = find(ismember(weight_MX_tmp.region_name, region_name{k}));
        weight_MX_tmp.afferent(idx_roi) = idx_afferent/(idx_efferent + idx_afferent);
        weight_MX_tmp.efferent(idx_roi) = idx_efferent/(idx_efferent + idx_afferent);
        if idx_efferent > 0 && idx_afferent > 0
            turn_point = round(length(t_in_tmp) * idx_afferent/(idx_efferent + idx_afferent), 0);
            % afferent part
            si_a = t_in_tmp(1); % Start index of polygon to be filled
            t_in_tmp_a = t_in_tmp(1:turn_point);
            t_out_tmp_a = t_in_tmp(turn_point:-1:1);
            x_a = [r2*cos(t_in_tmp_a) r3*cos(t_out_tmp_a) r2*cos(si_a)];
            y_a = [r2*sin(t_in_tmp_a) r3*sin(t_out_tmp_a) r2*sin(si_a)];
            hold on
            fill(x_a,y_a, [231,76,60]/255, 'edgecolor', [231,76,60]/255);

            if i == 6 && k == 1
                plot(1.31*cos(t_in_tmp(1:turn_point)) , 1.31*sin(t_in_tmp(1:turn_point)), ...
                    'Color', ROI_color{idx_roi_color, 3}, 'LineWidth', 1.5);
                % 画最大culster的DLPFC的afferent圆弧
            end

            % efferent part
            si_e = t_in_tmp(turn_point); % Start index of polygon to be filled
            t_in_tmp_e = t_in_tmp(turn_point:end);
            t_out_tmp_e = t_in_tmp(end:-1:turn_point);
            x_e = [r2*cos(t_in_tmp_e) r3*cos(t_out_tmp_e) r2*cos(si_e)];
            y_e = [r2*sin(t_in_tmp_e) r3*sin(t_out_tmp_e) r2*sin(si_e)];
            hold on
            fill(x_e,y_e, [36, 113,163]/255, 'edgecolor', [36, 113,163]/255);
        elseif idx_efferent == 0 && idx_afferent > 0
            hold on
            x_a = [r2*cos(t_in_tmp) r3*cos(t_out_tmp) r2*cos(si)];
            y_a = [r2*sin(t_in_tmp) r3*sin(t_out_tmp) r2*sin(si)];
            fill(x_a,y_a, [231,76,60]/255, 'edgecolor', [231,76,60]/255);
        elseif idx_efferent > 0 && idx_afferent == 0
            hold on
            x_e = [r2*cos(t_in_tmp) r3*cos(t_out_tmp) r2*cos(si_e)];
            y_e = [r2*sin(t_in_tmp) r3*sin(t_out_tmp) r2*sin(si_e)];
            fill(x_e,y_e, [36, 113,163]/255, 'edgecolor', [36, 113,163]/255);
        end
        %         end
    end
    eval(['weight_MX.cluster', num2str(i), ' = weight_MX_tmp;']);
    clear *_tmp
end
warning on

%Add and label gridlines
hold on
lineh(1) = polar([0,0],[0,1],'-');
lineh(2) = polar(linspace(0,2*pi,50),ones(1,50),'--');
lineh(3) = polar(linspace(0,2*pi,50),ones(1,50)*0.8,'--');
lineh(4) = polar(linspace(0,2*pi,50),ones(1,50)*0.6,'--');
lineh(5) = polar(linspace(0,2*pi,50),ones(1,50)*0.4,'--');
lineh(6) = polar(linspace(0,2*pi,50),ones(1,50)*0.2,'--');
set(lineh,'Color',[0.5,0.5,0.5],'LineWidth',1.5);
for i=1:5
    [x,y]=pol2cart(0,i/5);
    %     str = sprintf('%2.1f',((maxy-miny)*(5-i)/5)+miny);
    str = sprintf('%2.0f',25-5*i);
    text(x,y,str,'VerticalAlignment','bottom');
    if x == 1
        text(x+0.1,y,['Rescaled' newline 'Distance'],...
            'VerticalAlignment','middle',...
            'FontSize',11, ...
            'FontWeight', 'normal');
    end
end

%Prettier
set(ph(2),'XLim',[-1.5,2.5],'YLim',[-1.3,2.7],'Visible','off');
view(3)
axis fill
daspect([1,1,100]);
zoom(5);
view(2);


%%%%%%%%%%%%%%%%%%% plot subplot3 hierachical clustring %%%%%%%%%%%%%%%%%%%
axisl = 0.15;
axisb = 0.03;
axisw = 0.35;
axish = 0.45;

ph(3)=subplot('position', [axisl, axisb, axisw, axish]);
comp_num = unique(T_idx);
RT = 427.2082189359591;
RT_std = 188.8560094481999;
rectangle('Position',[RT - RT_std, 0,  2*RT_std, 120], ... % be carefule of the criteria here [x1 y1 delta_x y2]
    'FaceColor', [[128,128,128]./255, 0.3], ...
    'EdgeColor', [0, 0, 0, 0]);
hold on


comp_name = {'Pre-Response Alpha/Beta', ...
    'Pre-Response High Gamma',...
    'In-Response Low Gamma',...
    'Post-Response Low Gamma', ...
    'Pre-Response Low Gamma', ...
    'In-Response Alpha/Beta', ...
    'Post-Response Alpha/Beta', ...
    'Post-Response High Gamma'};
comp_num_reorder = [2, 8, 4, 7, 1, 3, 5, 6];

for i = 1:length(comp_num_reorder)
    idx = find(T_idx == comp_num_reorder(i));
    weight_MX.ec_cluster_num(comp_num_reorder(i)) = length(idx);
    weight_MX.ec_cluster_name{comp_num_reorder(i)} = comp_name{i};
    chk(i,:) = mean(a(idx,:));
    scatter(a(idx,1), a(idx,2), 42, comp_color(comp_num_reorder(i),:), ...
        'filled', ...
        'MarkerEdgeColor',[1 1 1],...
        'LineWidth',1);
    hold on
    display(num2str(length(idx)))
end


plot([RT, RT],[0,120], 'Color', [[0,0,0]./255, 1], 'LineWidth', 2);
set(gca, 'FontSize', 11, 'FontWeight', 'normal'); % 这将设置刻度标签的字体大小为12号字体

hold on
text(RT * 0.95, 120*1.07,'RT', 'Color', [0, 0, 0], 'FontSize', 11, 'FontWeight', 'normal');

legend(comp_name, 'FontSize',11, 'location', 'southoutside', ...
    'box', 'off', 'NumColumns',2, 'FontWeight', 'normal');
xlabel('Time (ms)', 'Color', [0, 0, 0], 'FontSize',11, 'FontWeight', 'normal');
ylabel('Frequency (Hz)', 'Color', [0, 0, 0], 'FontSize',11, 'FontWeight', 'normal');
set(ph(2), 'FontWeight', 'normal', 'TickLength',[0, 0], 'box', 'on');

box on

%%%%%%%%%%%%%%%%%%% plot subplot4 fc-ec to DLPFC %%%%%%%%%%%%%%%%%%%

axisl = 0.565;
axisb = 0.15;
axisw = 0.2;
axish = 0.3;

ph(4)=subplot('position', [axisl, axisb, axisw, axish]);

y = imread('.\brain_map_FS/fc_EC_2_DLPFC_HC_lateral.png', 'BackgroundColor', [1 1 1]);
y(:,[1:100, 1100:end],:) = [];
y([1:70,780:end],:,:) = [];
imshow(y);

axisl = 0.775;
axisb = 0.15;
axisw = 0.2;
axish = 0.3;

ph(5)=subplot('position', [axisl, axisb, axisw, axish]);

y = imread('.\brain_map_FS/fc_EC_2_DLPFC_HC_medial.png', 'BackgroundColor', [1 1 1]);
y(:,[1:100, 1100:end],:) = [];
y([1:50,760:end],:,:) = [];
imshow(y);



%%%%%%%%%%%%%%%%%% add annotations %%%%%%%%%%%%%%%%%%%%%%

annotation('rectangle', [0.295, 0.549, 0.266, 0.032], 'Color', 'red', 'LineWidth',1.5);
annotation('line',[0.561 0.61],[0.565 0.565], 'LineWidth', 1.5, 'Color',[255,0,0]/255)
annotation('line',[0.61 0.61],[0.565 0.78], 'LineWidth', 1.5, 'Color',[255,0,0]/255)
annotation('arrow',[0.61 0.625],[0.78 0.78], 'LineWidth', 1.5, 'Color',[255,0,0]/255)
annotation('arrow',[0.61 0.61],[0.565 0.42], 'LineWidth', 1.5, 'Color',[255,0,0]/255)

annotation('line',[0.675 0.675],[0.59 0.525], 'LineWidth', 1.5, 'Color', [0.466, 0.674, 0.188])
annotation('line',[0.675 0.44],[0.525 0.525], 'LineWidth', 1.5, 'Color', [0.466, 0.674, 0.188])
annotation('arrow',[0.44 0.44],[0.525 0.479], 'LineWidth', 1.5, 'Color', [0.466, 0.674, 0.188])





annotation('rectangle',...
    [0.9 0.55 0.02 0.01],...
    'LineStyle', 'None', ...
    'FaceColor',[231,76,60]/255);
annotation('rectangle',...
    [0.9 0.53 0.02 0.01],...
    'LineStyle', 'None', ...
    'FaceColor',[36, 113,163]/255);
annotation('textbox',...
    [0.92 0.555 0.07 0.015],...
    'String','Afferents', ...
    'LineStyle', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontWeight', 'normal', ...
    'FontSize', 11);
annotation('textbox',...
    [0.92 0.535 0.07 0.015],...
    'String','Efferents', ...
    'LineStyle', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontWeight', 'normal', ...
    'FontSize', 11);

annotation('textbox',...
    [0.08, 0.88, 0.1, 0.1],...
    'String','(A)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);

annotation('textbox',...
    [0.57, 0.88, 0.1, 0.1],...
    'String','(B)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);

annotation('textbox',...
    [0.08, 0.42, 0.1, 0.1],...
    'String','(C)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);

annotation('textbox',...
    [0.57, 0.42, 0.1, 0.1],...
    'String','(D)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);





set(gcf,'color','w');
print(gcf,'./Figure2.png','-dpng','-r300');
close gcf

