clear
%%%%%%%%%%%%%%%%%%%  plot gene expression map %%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf,'units','normalized','position',[0 0 0.5 1])

axisleft1 = [0.05, 0.05, 0.05, 0.21, 0.37, 0.37, 0.37];
axisbottom1 = [0.6, 0.73, 0.86, 0.73, 0.6, 0.73, 0.86];
axiswidth1 = [0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12];
axisheight1 = [0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12];

figs = {'co_exp_layer1_z_lateral.png', 'co_exp_layer2_z_lateral.png', ...
    'co_exp_layer3_z_lateral.png', 'fc_EC_2_DLPFC_HC_z_lateral.png', ...
    'co_exp_layer4_z_lateral.png', 'co_exp_layer5_z_lateral.png', ...
    'co_exp_layer6_z_lateral.png'};

for i = 1:length(figs)
    ph(i)=subplot('position', [axisleft1(i), axisbottom1(i), axiswidth1(i), axisheight1(i)]);

    y = imread(['.\brain_map_FS/' figs{i}], 'BackgroundColor', [1 1 1]);
    y(:,[1:100, 1100:end],:) = [];
    y([1:70,780:end],:,:) = [];
    imshow(y);
end

%%%%%%%%%%%%%%%%%%%%%% subplot 2 the laminal difference %%%%%%%%%%%%%%%%%%%
clear

axisleft = [0.63, 0.5];
axisbottom = [0.57, 0.055];
axiswidth = [0.35, 0.48];
axisheight = [0.41, 0.41];

load('./CoExp_group.mat');
fc_var_name = {'fc_var_group_diff', 'fc_var_EC', 'fc_var_NVC_frequency'};
gene_sets_name = {'Dopamine', 'GABA','Glutamine'};

%% just take the set with at least one < 0.05
for i = 1:length(fc_var_name)
    for j = 1:length(gene_sets_name)
        eval(['tmp = CoExp_All.', fc_var_name{i}, '.', gene_sets_name{j}, ';'])
        tmp1 = [tmp.afferent_p_Layer1, tmp.afferent_p_Layer2, tmp.afferent_p_Layer3, ...
            tmp.afferent_p_Layer4, tmp.afferent_p_Layer5, tmp.afferent_p_Layer6];
        idx = find(any(tmp1 < 0.05, 2));
        idx = setdiff(1:length(tmp1), idx);
        idx_tmp = find(ismember(tmp.Annotations, 'Sensory perception of sweet, bitter, and umami (glutamate) taste'));
        if ~isempty(idx_tmp)
            idx = unique([idx_tmp, idx]);
        end
        eval(['CoExp_All.', fc_var_name{i}, '.', gene_sets_name{j}, '(idx,:) = [];'])
        clear tmp*
    end
end

%--------------------------------------------------------------------------
% Plot both healthy controls and patients with schizophrenia
ph(8) = subplot('position', [axisleft(1), axisbottom(1), axiswidth(1), axisheight(1)]);

SUPP = [];
DEEP = [];
G_names_supp = [];
G_names_deep = [];

clear layer_mean layer_error
for i = 1:length(gene_sets_name)
    eval(['tmp = CoExp_All.fc_var_NVC_frequency.', gene_sets_name{i}, ';'])
    supp = [tmp.afferent_r_Layer3; tmp.afferent_r_Layer2];
    deep = [tmp.afferent_r_Layer5; tmp.afferent_r_Layer6];
    [~,P,~,~] = ttest2(supp, deep);
%     supp = (tmp.afferent_r_Layer3 + tmp.afferent_r_Layer2)/2;
%     deep = (tmp.afferent_r_Layer5 + tmp.afferent_r_Layer6)/2;
%     [~,P,~,~] = ttest(supp, deep);
    p_diff(i,1) = P;
    N(i,1) = length(supp)/2;
    layer_mean(i,1) = nanmean(supp);
    layer_mean(i,2) = nanmean(deep);
    layer_error(i,1) = nanstd(supp)/sqrt(length(supp));
    layer_error(i,2) = nanstd(deep)/sqrt(length(deep));
    
    G_names_supp = [G_names_supp; repmat({gene_sets_name{i}}, length(supp), 1)];
    G_names_deep = [G_names_deep; repmat({gene_sets_name{i}}, length(deep), 1)];
    SUPP = [SUPP; supp];
    DEEP = [DEEP; deep];
    
end


bc = bar(layer_mean, 'FaceAlpha', 0.5, 'EdgeColor','flat', 'LineWidth', 1);
bc_FaceColor = get(bc, 'FaceColor');
hold on
% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(layer_mean);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, layer_mean(:,i), [], layer_error(:,i), ...
        'k', 'linestyle', 'none', 'LineWidth', 1.5);
end
hold off


legend({'Layer 2/3', 'Layer 5/6'}, ...
    'FontSize', 11, 'location', 'northeast', ...
    'box', 'off', 'NumColumns',1);
xlabel('Neurotransmitter-related human gene sets', 'Color', [0, 0, 0], 'FontSize', 11, 'FontWeight', 'normal');
ylabel({'Correlation of co-expression with', ...
    'WM-related fc-EC to DPLFC'}, 'Color', [0, 0, 0], 'FontSize', 11, 'FontWeight', 'normal');
set(gca, 'TickLength', [0,0], ...
    'ylim', [0.1, 0.35], ...
    'xticklabel', ...
    {'Dopamine', 'GABA', 'Glutamine'}, ...
    'xColor', [0, 0, 0], 'yColor', [0, 0, 0], 'FontWeight', 'normal', 'FontSize', 11);
for i = 1:length(p_diff)
    if p_diff(i) < 0.05
        text(i - 0.35,0.285,['{\it p} = ' num2str(round(p_diff(i),3))], ...
            'FontSize', 11, 'FontWeight', 'normal', 'Color', 'r');
    else
        text(i - 0.13,0.285,'{\it n.s}', 'FontSize', 11, 'FontWeight', 'normal', 'Color', 'b');
    end
    x = [i-0.25, i + 0.25];
    y = [0.27, 0.27];
    hold on
    if p_diff(i) < 0.05
        plot(x, y, 'k', 'LineWidth', 2, 'HandleVisibility','off', 'Color', 'r');
    else
        plot(x, y, 'k', 'LineWidth', 2, 'HandleVisibility','off', 'Color', 'b');
    end
    hold on
end
box on
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%% plot subplot 3 the GABAergic co-expression %%%%%%%%%%%%%%
clearvars -except  axisleft axisbottom axiswidth axisheight bc_*
load('./CoExp_group.mat');
% Select the gene set with nominally significant correlation with NVC or EC
% to display and eamine the layer effect.
% fc_var_name = {'fc_var_group_diff_PostR_HG', ...
%     'fc_var_group_diff', 'fc_var_EC', 'fc_var_NVC_frequency'};
fc_var_name = {'fc_var_NVC_frequency'};

gene_sets_name = {'GABA'};

for i = 1:length(fc_var_name)
    for j = 1:length(gene_sets_name)
        eval(['tmp = CoExp_All.', fc_var_name{i}, '.', gene_sets_name{j}, ';'])
        tmp1 = [tmp.afferent_p_Layer1, tmp.afferent_p_Layer2, tmp.afferent_p_Layer3, ...
            tmp.afferent_p_Layer4, tmp.afferent_p_Layer5, tmp.afferent_p_Layer6];
        idx = find(any(tmp1 < 0.05, 2));
        idx = setdiff(1:length(tmp1), idx);
        idx_tmp = find(ismember(tmp.Annotations, 'Sensory perception of sweet, bitter, and umami (glutamate) taste'));
        if ~isempty(idx_tmp)
            idx = unique([idx_tmp, idx]);
        end
        eval(['CoExp_All.', fc_var_name{i}, '.', gene_sets_name{j}, '(idx,:) = [];'])
        clear tmp*
    end
end
%--------------------------------------------------------------------------
colors = {[0.6350 0.0780 0.1840];
    [0.9290 0.6940 0.1250];
    [0.8500 0.3250 0.0980];
    [0.4940 0.1840 0.5560];
    [0 0.4470 0.7410];
    [0.4660 0.6740 0.1880]};

colors = colors([2,5], :);

for i = 1:length(fc_var_name)
    for j = 1:length(gene_sets_name)
        ph(9)=subplot('position', [axisleft(2), axisbottom(2), axiswidth(2), axisheight(2)]);
        eval(['tmp = CoExp_All.', fc_var_name{i}, '.', gene_sets_name{j}, ';'])
        tmp1 = [tmp.afferent_p_Layer6, tmp.afferent_p_Layer5, tmp.afferent_p_Layer4, ...
            tmp.afferent_p_Layer3, tmp.afferent_p_Layer2, tmp.afferent_p_Layer1];
        [tmp1,idx_sort] = sortrows(tmp1,6,'descend');
        tmp1 = -log10(tmp1);

        tmp1 = [mean(tmp1(:, 2:3), 2), mean(tmp1(:, 5:6), 2)];

%         tb = barh(tmp1, 'FaceAlpha', 0.5, 'EdgeColor','flat');

        hold on
        w = 0.4;
        tb1 = barh((1:length(tmp1))-w/2,tmp1(:,1),0.3, 'FaceColor', bc_FaceColor{2}, 'FaceAlpha', 0.5, 'EdgeColor','flat');
        tb2 = barh((1:length(tmp1))+w/2,tmp1(:,2),0.3, 'FaceColor', bc_FaceColor{1}, 'FaceAlpha', 0.5, 'EdgeColor','flat');

        
%         for k = 1:length(bc_FaceColor)
%             tb(k).FaceColor = bc_FaceColor{k};
%         end

        annotation_label = tmp.Annotations(idx_sort);
        annotation_label = strrep(annotation_label, '_', ' ');
        set(gca, 'TickLength', [0,0], 'ylim', [0.25,17.75], 'YTick', 1:length(annotation_label), ...
            'yticklabel',annotation_label, 'XColor', 'k', 'YColor', 'k', ...
            'FontSize', 11, 'FontWeight', 'normal');
        hold on
        if max(tmp1, [], "all") > -log10(0.001)
            plot([-log10(0.05), -log10(0.05)],ylim,'k--');
            plot([-log10(0.001), -log10(0.001)],ylim,'k--');
            text(-log10(0.05 * 1.5),max(ylim) + 0.028 * max(ylim),'-log10({\itp} = 0.05)', 'FontSize', 11, 'FontWeight', 'normal');
            text(-log10(0.001 * 1.5),max(ylim) + 0.028 * max(ylim),'-log10({\itp} = 0.001)', 'FontSize', 11, 'FontWeight', 'normal');
        elseif max(tmp1, [], "all") > -log10(0.05)
            plot([-log10(0.05), -log10(0.05)],ylim,'k--');
            text(-log10(0.05 * 1.5),max(ylim) + 0.028 * max(ylim),'-log10({\itp} = 0.05)', 'FontSize', 11, 'FontWeight', 'normal');
        end

        %         legend(tb([6:-1:1]), {'Layer 1', 'Layer 2', 'Layer 3', 'Layer 4', 'Layer 5', 'Layer 6'}, ...
        %             'FontSize', 11, 'location', 'southeast', ...
        %             'box', 'off', 'NumColumns',1);

        legend([tb2, tb1], {'Layer 2/3', 'Layer 5/6'}, ...
            'FontSize', 11, 'location', 'southeast', ...
            'box', 'off', 'NumColumns',1);

        xlabel('-log10({\itp})', 'Color', [0, 0, 0], 'FontSize', 11, 'FontWeight', 'normal');

        box on

        clear tmp* tb
    end
end







%%%%%%%%%%%%%% add title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
annotation('textbox',...
    [0.02 0.98 0.02 0.02],...
    'String','(A)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);
annotation('textbox',...
    [0.545 0.98 0.02 0.02],...
    'String','(B)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);
annotation('textbox',...
    [0.02 0.47 0.02 0.02],...
    'String','(C)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);


annotation('arrow',[0.21 0.17],[0.788, 0.788], 'LineWidth', 1.5, 'Color',[0,0,0]/255);
annotation('line',[0.195 0.195],[0.788 0.918], 'LineWidth', 1.5, 'Color',[0,0,0]/255);
annotation('line',[0.195 0.195],[0.788 0.658], 'LineWidth', 1.5, 'Color',[0,0,0]/255);
annotation('arrow',[0.195 0.17],[0.918 0.918], 'LineWidth', 1.5, 'Color',[0,0,0]/255);
annotation('arrow',[0.195 0.17],[0.658 0.658], 'LineWidth', 1.5, 'Color',[0,0,0]/255);

annotation('arrow',[0.33 0.37],[0.788, 0.788], 'LineWidth', 1.5, 'Color',[0,0,0]/255);
annotation('line',[0.345 0.345],[0.788 0.918], 'LineWidth', 1.5, 'Color',[0,0,0]/255);
annotation('line',[0.345 0.345],[0.788 0.658], 'LineWidth', 1.5, 'Color',[0,0,0]/255);
annotation('arrow',[0.345 0.37],[0.918 0.918], 'LineWidth', 1.5, 'Color',[0,0,0]/255);
annotation('arrow',[0.345 0.37],[0.658 0.658], 'LineWidth', 1.5, 'Color',[0,0,0]/255);


annotation('textbox',...
    [0.08 0.51 0.1 0.115],...
    'String','Layer 3', ...
    'LineStyle', 'none', ...
    'FontWeight', 'normal', ...
    'FontSize', 11);
annotation('textbox',...
    [0.08 0.64 0.1 0.115],...
    'String','Layer 2', ...
    'LineStyle', 'none', ...
    'FontWeight', 'normal', ...
    'FontSize', 11);
annotation('textbox',...
    [0.08 0.77 0.1 0.115],...
    'String','Layer 1', ...
    'LineStyle', 'none', ...
    'FontWeight', 'normal', ...
    'FontSize', 11);
annotation('textbox',...
    [0.40 0.51 0.1 0.115],...
    'String','Layer 6', ...
    'LineStyle', 'none', ...
    'FontWeight', 'normal', ...
    'FontSize', 11);
annotation('textbox',...
    [0.40 0.64 0.1 0.115],...
    'String','Layer 5', ...
    'LineStyle', 'none', ...
    'FontWeight', 'normal', ...
    'FontSize', 11);
annotation('textbox',...
    [0.40 0.77 0.1 0.115],...
    'String','Layer 4', ...
    'LineStyle', 'none', ...
    'FontWeight', 'normal', ...
    'FontSize', 11);
annotation('textbox',...
    [0.2 0.635 0.15 0.115],...
    'String','fc-EC to DLPFC', ...
    'LineStyle', 'none', ...
    'FontWeight', 'normal', ...
    'FontSize', 11);

% plot a colorbar for standardized fc-EC to DLPFC in HV
hold on
row_cmap = 1000;  %define the color scale
color_map=ones(row_cmap,3);  %define the color matrix
bar_range = [-1,1]; % define the range of value reflected by the colorbar
ratio = round((0 - min(bar_range)) / (max(bar_range) - min(bar_range)),2);
color_neg = [153, 51, 255]./255; % blue, fix B change R and G
color_pos = [255, 204, 0]./255; % orange, fix R change B and G
color_map(1:row_cmap,:) = [
    [color_neg(1):(1-color_neg(1))/(row_cmap*ratio - 1):1;
    color_neg(2):(1-color_neg(2))/(row_cmap*ratio - 1):1;
    ones(1, int16(row_cmap*ratio))]';
    flip([ones(1, int16(row_cmap*(1-ratio)));
    color_pos(2):(1-color_pos(2))/(row_cmap*(1-ratio) - 1):1;
    color_pos(3):(1-color_pos(3))/(row_cmap*(1-ratio) - 1):1;]')
    ];
ph(11) = subplot('position', [0.3, 0.65, 0.001, 0.001]);
colormap(gca, color_map);
% hCB = colorbar('north');
hCB = colorbar('north', ...
    'YTick',[],...
    'EdgeColor','none', ...
    'TickLength', 0, ...
    'Box', 'off',...
    'FontSize', 16, 'FontWeight', 'normal');
hCB.Position = [0.09 0.55 0.12 0.02];

annotation(gcf,'textbox',[0.045 0.545 0.05 0.03],...
    'LineStyle','none',...
    'String','-1.7',...
    'FontSize', 11, 'FontWeight', 'normal');
annotation(gcf,'textbox',[0.21 0.545 0.05 0.03],...
    'LineStyle','none',...
    'String','2.6',...
    'FontSize', 11, 'FontWeight', 'normal');
annotation(gcf,'textbox',[0.137 0.547 0.05 0.03],...
    'LineStyle','none',...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');
ph(11).Visible = 'off';
annotation('textbox',...
    [0.075 0.398 0.2 0.155],...
    'String','WM-related fc-EC', ...
    'LineStyle', 'none', ...
    'FontWeight', 'normal', ...
    'FontSize', 11);

% plot a colorbar for co-expression
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
ph(12) = subplot('position', [0.35, 0.65, 0.001, 0.001]);
colormap(gca, color_map);
% hCB = colorbar('north');
hCB = colorbar('north', ...
    'YTick',[],...
    'EdgeColor','none', ...
    'TickLength', 0, ...
    'Box', 'off',...
    'FontSize', 16, 'FontWeight', 'normal');
hCB.Position = [0.34 0.55 0.12 0.02];

annotation(gcf,'textbox',[0.298 0.545 0.05 0.03],...
    'LineStyle','none',...
    'String','-0.9',...
    'FontSize', 11, 'FontWeight', 'normal');
annotation(gcf,'textbox',[0.46 0.545 0.05 0.03],...
    'LineStyle','none',...
    'String','0.9',...
    'FontSize', 11, 'FontWeight', 'normal');
annotation(gcf,'textbox',[0.385 0.547 0.05 0.03],...
    'LineStyle','none',...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');
ph(12).Visible = 'off';
annotation('textbox',...
    [0.315 0.398 0.2 0.155],...
    'String','Gene co-expression', ...
    'LineStyle', 'none', ...
    'FontWeight', 'normal', ...
    'FontSize', 11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf,'color','w');
print(gcf,'./Figure3.png','-dpng','-r300');
close gcf