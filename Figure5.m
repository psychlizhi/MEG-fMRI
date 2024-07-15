clear
%%%%%%%%%%%%%%%%%%%  plot gene expression map %%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf,'units','normalized','position',[0 0 0.6 0.8])
axisleft =   [0.02, 0.02, 0.24, 0.47, 0.3, 0.75, 0.75, 0.75];
axisbottom = [0.71,  0.48,  0.16, 0.16, 0.46, 0.76, 0.43, 0.1];
axiswidth =  [0.22,  0.22,  0.22, 0.22, 0.38, 0.23, 0.23, 0.23];
axisheight = [0.22,  0.22,  0.22, 0.22, 0.52, 0.22, 0.22, 0.22];


load('./Restricted_Access/MEG_EC_NVC_idv.mat');

%%%%%%%%%%%%% use EC as sensitivity test for individual fwEC %%%%%%%%%%%%%%
test_EC = 1; % 0 = use fwEC; 1 = use EC
if test_EC == 1
    load('./Restricted_Access/NB_2back.mat')
    ec_name = [NB_2back.fc_name(:, 2:3); NB_2back.fc_name(:, 3:-1:2)];
    meg_fc_name = MEG_EC_NVC_idv.dlpfc_whole.ec_fmri_name;
    idx_ec = [];
    for i = 1:size(meg_fc_name,1)
        tmp = find(ismember(ec_name(:,1), meg_fc_name(i,1)) & ismember(ec_name(:,2), meg_fc_name(i,2)));
        if ~isempty(tmp)
            idx_ec = [idx_ec; tmp];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bev = readtable('./Restricted_Access/bev_final.csv');
HC = bev.subjects_MEG(find(ismember(bev.group, 'HC')));
SCZ = bev.subjects_MEG(find(ismember(bev.group, 'SCZ')));
HC_MRI = bev.subjects_MRI(find(ismember(bev.group, 'HC')));
SCZ_MRI = bev.subjects_MRI(find(ismember(bev.group, 'SCZ')));


% bev.MEG_RT_2back = bev.MEG_RT_2back - bev.MEG_RT_0back;
% bev.MRI_RT_2back = bev.MRI_RT_2back - bev.MRI_RT_0back;
% bev.MEG_ACC_2back = bev.MEG_ACC_2back - bev.MEG_ACC_0back;
% bev.MRI_ACC_2back = bev.MRI_ACC_2back - bev.MRI_ACC_0back;
bev.MEG_RT_2back = (bev.MEG_RT_2back + bev.MRI_RT_2back)/2;
bev.MEG_ACC_2back = (bev.MEG_ACC_2back + bev.MRI_ACC_2back)/2;


roi = unique(MEG_EC_NVC_idv.dlpfc_whole.ec_fmri_name(:,1));
Model_Residuals_NVC = [];
Model_Residuals_NVC = MEG_EC_NVC_idv.dlpfc_whole.Model_Residuals_NVC;


for i = 1:length(HC)
    HC_idx_bev(i,1) = find(ismember(bev.subjects_MEG, HC{i}));
    HC_idx_meg(i,1) = find(ismember(MEG_EC_NVC_idv.dlpfc_whole.subjects_MEG, HC{i}));
    if test_EC == 1
        HC_idx_ec(i,1) = find(ismember(NB_2back.subjects, HC_MRI{i}));
    end
end
for i = 1:length(SCZ)
    SCZ_idx_bev(i,1) = find(ismember(bev.subjects_MEG, SCZ{i}));
    SCZ_idx_meg(i,1) = find(ismember(MEG_EC_NVC_idv.dlpfc_whole.subjects_MEG, SCZ{i}));
    if test_EC == 1
        SCZ_idx_ec(i,1) = find(ismember(NB_2back.subjects, SCZ_MRI{i}));
    end
end


prs_var = table2array(bev(HC_idx_bev,14:23));

bev_other = table2array(bev(HC_idx_bev,[9, 36,39,42,43]));
bev_other_name = bev.Properties.VariableNames([9, 36,39,42,43]);

bev_cov = table2array(bev(HC_idx_bev,[6,7,8]));

% tmp_HC = Model_Residuals_NVC(:,HC_idx_meg)';
% for i = 1:size(tmp_HC,1)
%     tmp_HC(i, find(isnan(tmp_HC(i,:)))) = nanmean(tmp_HC(i,:));
% end

tmp_HC = Model_Residuals_NVC(:,HC_idx_meg)';
for i = 1:size(tmp_HC,2)
    idx = find(isnan(tmp_HC(:, i)));
    if ~isempty(idx)
        tmp_HC(idx, i) = nanmean(tmp_HC(:, i));
    end
end


if test_EC == 1
    EC = [NB_2back.B_forward, NB_2back.B_backward];
    tmp_HC = EC(HC_idx_ec,idx_ec);
end


idx_ex = find(isnan(mean([prs_var, tmp_HC],2)));
prs_var(idx_ex, :) = [];
bev_other(idx_ex, :) = [];
tmp_HC(idx_ex, :) = [];
bev_cov(idx_ex, :) = [];

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(tmp_HC,prs_var,10);


[r1,p1] = partialcorr([XS(:,1), prs_var(:,2), bev_other], bev_cov, 'rows','pairwise');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
roi = unique(MEG_EC_NVC_idv.dlpfc_whole.ec_fmri_name(:,1));

use_all_HC = 0; % 0 = use HC have no prs for case control analysis
% only use the HC have no prs so two analysis have no shared
% participants
if use_all_HC == 0
    HC_new = bev.subjects_MEG(find(ismember(bev.group, 'HC') & isnan(bev.p0_001)));
    HC_idx_meg_idx = [];
    for i_chk = 1:length(HC_idx_meg)
        if ~isempty(find(ismember(HC_new, MEG_EC_NVC_idv.dlpfc_whole.subjects_MEG(HC_idx_meg(i_chk)))))
            HC_idx_meg_idx = [HC_idx_meg_idx; i_chk];
        end
    end
    HC_idx_meg = HC_idx_meg(HC_idx_meg_idx);
    HC_idx_bev = HC_idx_bev(HC_idx_meg_idx);
end

Model_Residuals_NVC = [];
for i = 1:length(roi)
    idx_roi = find(ismember(MEG_EC_NVC_idv.dlpfc_whole.ec_fmri_name(:,1), roi{i}));
    
    XL_roi(i,1) = mean(XL(idx_roi,1),'all','omitnan');
    
    Model_Residuals_NVC_SCZ(i,1) = mean(MEG_EC_NVC_idv.dlpfc_whole.Model_Residuals_NVC(idx_roi,SCZ_idx_meg),'all','omitnan');
    Model_Residuals_NVC_HC(i,1) = mean(MEG_EC_NVC_idv.dlpfc_whole.Model_Residuals_NVC(idx_roi,HC_idx_meg),'all','omitnan');
    HC = MEG_EC_NVC_idv.dlpfc_whole.Model_Residuals_NVC(idx_roi,HC_idx_meg);
    SCZ = MEG_EC_NVC_idv.dlpfc_whole.Model_Residuals_NVC(idx_roi,SCZ_idx_meg);
    
    
    bev_cov_all = table2array(bev([HC_idx_bev; SCZ_idx_bev], 6:8));
    tmp_HC = nanmean(HC)';
    tmp_SCZ = nanmean(SCZ)';
    
    
    tmp_bev = array2table([[tmp_HC; tmp_SCZ], ...
        [ones(length(tmp_HC'),1); zeros(length(tmp_SCZ'),1)], ...
        bev_cov_all(:,1),bev_cov_all(:,2),bev_cov_all(:,3)], ...
        'VariableNames', {'RES', 'Group', 'age', 'gender', 'education'});
    tmp_bev.Group = categorical(tmp_bev.Group);
    tmp_bev.gender = categorical(tmp_bev.gender);

    md = fitlm(tmp_bev,'RES~Group + age + gender + education');

    T(i, 1) = md.Coefficients.tStat(2);
    clear idx_roi

end


%%%%%%%%%%%%%%%%%%%  plot gene expression map %%%%%%%%%%%%%%%%%%%%%%%%%



fig_names = {'T_group_diff_lateral.png', 'T_group_diff_medial.png', 'XL_roi_lateral.png', 'XL_roi_medial.png'};
for i = 1:length(fig_names)
    ph(i)=subplot('position', [axisleft(i), axisbottom(i), axiswidth(i), axisheight(i)]);
    y = imread(['.\brain_map_FS\', fig_names{i}], 'BackgroundColor', [1 1 1]);
    if i == 1 || i == 3
        y(:,[1:100, 1100:end],:) = [];
        y([1:70,780:end],:,:) = [];
    elseif i == 2 || i == 4
        y(:,[1:100, 1100:end],:) = [];
        y([1:50,760:end],:,:) = [];
    end
    imshow(y);
end

ph(5)=subplot('position', [axisleft(5), axisbottom(5), axiswidth(5), axisheight(5)]);

Y = zscore(XL_roi);
X = zscore(T);
[r, p] = corr(X, Y, "rows","complete");

% First the first plot.
coefficients = polyfit(X, Y, 1);
% Make new x coordinates
x = linspace(min(X), max(X), 100);
y = polyval(coefficients, x);
s = sprintf(['h' num2str(i) ' = plot(%s, ''%s'', ''%s'', ''%s'', ''%s'', %s)'], ...
    'x, y', '-', 'Color', '#D95319', 'LineWidth', num2str(2));
eval(s);
%                     h = plot(x, y, '-', 'Color', colors{i_layer}, 'LineWidth', 2);
hold on;

% Plot scatter
scatter(X, Y, ...
    50, ...
    [0.8500 0.3250 0.0980], ...
    'LineWidth',0.2);
hold on;
scatter(X, Y, ...
    50, 'filled', ...
    'MarkerFaceColor', '#D95319', ...
    'MarkerFaceAlpha', 0.4, ...
    'LineWidth',1.5);
hold on;

xlabel('PRS-weighted fc-ECs acrosscortical parcels (  )', ...
    'Color', [0, 0, 0], 'FontSize',11, 'FontWeight', 'normal');
ylabel('Group difference between HC and SCZ on fc-EC (  )', ...
    'Color', [0, 0, 0], 'FontSize',11, 'FontWeight', 'normal');

hold on;
dummyh= line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
tmp_p = num2str(p, '%.e');
legend(dummyh, {['{\it r} = ' num2str(round(r,3)) ', ' '{\itp} = ' tmp_p]}, ...
    'FontSize', 11, ...
    'Location', 'southeast', 'Box', 'off', 'FontWeight', 'normal');
        

set(gca,  'YTick', [-2, -1, 0, 1, 2], 'YTickLabel', {'-2', '-1', '0', '1', '2'}, ...
    'xColor', 'k', 'yColor', 'k', 'ylim', [-2.5, 2.2], 'xlim', [-2.5, 2.5], ...
    'TickLength', [0,0], 'FontWeight', 'normal', 'FontSize',11);

tmp1 = xlim;
tmp2 = ylim;
hold on

% increase the distance between label and axis
% xh = get(gca,'xlabel'); % handle to the label object
% p = get(xh,'position'); % get the current position property
% p(2) = 0.3 + p(2);       % double the distance,
% set(xh,'position',p)   % set the new position

% yh = get(gca,'ylabel'); % handle to the label object
% p = get(yh,'position'); % get the current position property
% p(1) = -0.01 + p(1);       % double the distance,
% set(yh,'position',p)   % set the new position

box on

%%%%%%%%%%%%%%%%%% plot correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
XS_one = XS(:,1);

bev_cand = [prs_var(:,2), bev_other(:,[1,3])];

labels = {'(B)', '(C)', '(D)'};
subplot_idx = [4, 8, 12];
colors = {'#0072BD', '#A2142F', '#7E2F8E'};
colors_vals = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.4940 0.1840 0.5560]];
round_num = [11,3,3];
xlabels = {'PRS for SCZ', 'IQ', 'RT'};
for i = 1:3
    
    ph(i + 5)=subplot('position', [axisleft(i + 5), axisbottom(i + 5), axiswidth(i + 5), axisheight(i + 5)]);

    Y = zscore(bev_cand(:,i));
    X = zscore(XS_one);
    
    [r, p] = corr(X, Y, "rows","complete");
    
    [r, p] = partialcorr([X, Y], bev_cov, "rows","complete");
    r = r(1,2:end);
    p = p(1,2:end);
    
    % First the first plot.
    coefficients = polyfit(X, Y, 1);
    % Make new x coordinates
    x = linspace(min(X), max(X), 100);
    y = polyval(coefficients, x);
    s = sprintf(['h' num2str(i) ' = plot(%s, ''%s'', ''%s'', ''%s'', ''%s'', %s)'], ...
        'x, y', '-', 'Color', colors{i}, 'LineWidth', num2str(2));
    eval(s);
    %                     h = plot(x, y, '-', 'Color', colors{i_layer}, 'LineWidth', 2);
    hold on;
    
    % Plot scatter
    scatter(X, Y, ...
        50, colors_vals(i,:), ...
        'LineWidth',0.2);
    hold on;
    scatter(X, Y, ...
        50, 'filled', ...
        'MarkerFaceColor',colors{i}, ...
        'MarkerFaceAlpha', 0.4, ...
        'LineWidth',1);
    hold on;
   
    xlabel({'PRS-weighted fc-ECs', ...
        'across participants (  )'},...
        'Color', [0, 0, 0], 'FontSize',11, 'FontWeight', 'normal');

    ylabel([xlabels{i} ' (  )'], ...
        'Color', [0, 0, 0], 'FontSize',11, 'FontWeight', 'normal');
    
    hold on;
    dummyh= line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
    
    
%     legend(dummyh, {['{\it r} = ' num2str(round(r,3)) ', ' '{\it p} = ' num2str(round(p, round_num(i)))]}, ...
%         'FontSize', 12, ...
%         'Location', 'southeast', 'Box', 'off', 'FontWeight', 'normal');
    
    if p > 0.05
        legend(dummyh, {['{\it r} = ' num2str(round(r,3)) ', ' '{\it p} = ' num2str(round(p, round_num(i)))]}, ...
            'FontSize', 11, ...
            'Location', 'southeast', 'Box', 'off', 'FontWeight', 'normal');
        hold on
    elseif p <= 0.05 && p >= 0.001
        legend(dummyh, {['{\it r} = ' num2str(round(r,3)) ', ' '{\itp} =' num2str(round(p, 3))]}, ...
            'FontSize', 11, ...
            'Location', 'southeast', 'Box', 'off', 'FontWeight', 'normal');
        hold on
    elseif p < 0.001
        tmp_p = num2str(p, '%.e');
        legend(dummyh, {['{\it r} = ' num2str(round(r,3)) ', ' '{\itp} = ' tmp_p]}, ...
            'FontSize', 11, ...
            'Location', 'southeast', 'Box', 'off', 'FontWeight', 'normal');
        hold on
        clear tmp_p*
    end
    
    if test_EC == 1
        set(gca,  'YTick', [-2, 0, 2], 'YTickLabel', {'-2', '0', '2'}, ...
            'xColor', 'k', 'yColor', 'k', 'ylim', [-4, 3], 'TickLength', ...
            [0,0], 'FontWeight', 'normal', 'FontSize', 11);
    elseif test_EC == 0
        set(gca,  'YTick', [-2, 0, 2], 'YTickLabel', {'-2', '0', '2'}, ...
            'XTick', [-3, 0, 3], 'XTickLabel', {'-3', '0', '3'}, ...
            'xColor', 'k', 'yColor', 'k', 'ylim', [-4, 3], 'xlim', [-3, 3],'TickLength', ...
            [0,0], 'FontWeight', 'normal', 'FontSize', 11);
    end

    box on

    tmp1 = xlim;
    tmp2 = ylim;
    hold on

end

%%%%%%%%%%%%%%%%%%% plot just a colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
row_cmap = 1000;  %define the color scale
color_map=ones(row_cmap,3);  %define the color matrix
bar_range = [-1,1]; % define the range of value reflected by the colorbar
ratio = round((0 - min(bar_range)) / (max(bar_range) - min(bar_range)),2);
color_neg = [0, 51, 255]./255; % blue, fix B change R and G
color_pos = [255, 87, 51]./255; % orange, fix R change B and G
color_map(1:row_cmap,:) = [
    [color_neg(1):(1-color_neg(1))/(row_cmap*ratio - 1):1;
    color_neg(2):(1-color_neg(2))/(row_cmap*ratio - 1):1;
    ones(1, int16(row_cmap*ratio))]';
    flip([ones(1, int16(row_cmap*(1-ratio)));
    color_pos(2):(1-color_pos(2))/(row_cmap*(1-ratio) - 1):1;
    color_pos(3):(1-color_pos(3))/(row_cmap*(1-ratio) - 1):1;]')
    ];
ph(9) = subplot('position', [0.01, 0.01, 0.001, 0.001]);
colormap(gca, color_map);
% hCB = colorbar('north');
hCB = colorbar('north', ...
    'YTick',[],...
    'EdgeColor','none', ...
    'TickLength', 0, ...
    'Box', 'off',...
    'FontSize', 16, 'FontWeight', 'normal');
hCB.Position = [0.34 0.08 0.25 0.03];

annotation(gcf,'textbox',[0.3 0.081 0.05 0.03],...
    'LineStyle','none',...
    'String','-2.5',...
    'FontSize', 11, 'FontWeight', 'normal');
annotation(gcf,'textbox',[0.59 0.081 0.05 0.03],...
    'LineStyle','none',...
    'String','2.5',...
    'FontSize', 11, 'FontWeight', 'normal');
annotation(gcf,'textbox',[0.45 0.082 0.05 0.03],...
    'LineStyle','none',...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');
ph(9).Visible = 'off';

annotation('textbox',...
    [0.253 0.98 0.02 0.02],...
    'String','(A)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);
annotation('textbox',...
    [0.703 0.98 0.02 0.02],...
    'String','(B)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);
annotation('textbox',...
    [0.703 0.65 0.02 0.02],...
    'String','(C)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);
annotation('textbox',...
    [0.703 0.32 0.02 0.02],...
    'String','(D)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);

annotation(gcf,'textbox',[0.277 0.935 0.05 0.03],...
    'LineStyle','none',...
    'Rotation', 90, ...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');

annotation(gcf,'textbox',[0.727 0.918 0.05 0.03],...
    'LineStyle','none',...
    'Rotation', 90, ...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');

annotation(gcf,'textbox',[0.727 0.538 0.05 0.03],...
    'LineStyle','none',...
    'Rotation', 90, ...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');

annotation(gcf,'textbox',[0.727 0.212 0.05 0.03],...
    'LineStyle','none',...
    'Rotation', 90, ...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');

annotation(gcf,'textbox',[0.635 0.401 0.05 0.03],...
    'LineStyle','none',...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');

annotation(gcf,'textbox',[0.922 0.675 0.05 0.03],...
    'LineStyle','none',...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');

annotation(gcf,'textbox',[0.922 0.344 0.05 0.03],...
    'LineStyle','none',...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');

annotation(gcf,'textbox',[0.922 0.015 0.05 0.03],...
    'LineStyle','none',...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');

set(gcf,'color','w');
print(gcf,'./Figure5.png','-dpng','-r300');
close gcf
