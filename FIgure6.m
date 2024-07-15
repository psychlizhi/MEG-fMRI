%%%%%%%%%%%%%%%%%%%%%%% plot layer effect %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%plot subplot 1 ~ 3 %%%%%%%%%%%%%%%%%%%%%%
clear
set(gcf,'units','normalized','position',[0 0 0.7 1])


fc_EC_vs_EC = 1; % 1 = use fc-EC, 0 = use EC


load('./Restricted_Access/CoExp_idv_corr.mat');
bev = readtable('./Restricted_Access/bev_final.csv');
HC = bev.subjects_MEG(find(ismember(bev.group, 'HC')));
SCZ = bev.subjects_MEG(find(ismember(bev.group, 'SCZ')));
HC_idx = find(ismember(CoExp_All.subjects_MEG, HC));
SCZ_idx = find(ismember(CoExp_All.subjects_MEG, SCZ));

for i = 1:length(HC_idx)
    HC_idx_bev(i,1) = find(ismember(bev.subjects_MEG, CoExp_All.subjects_MEG(HC_idx(i))));
end
for i = 1:length(SCZ_idx)
    SCZ_idx_bev(i,1) = find(ismember(bev.subjects_MEG, CoExp_All.subjects_MEG(SCZ_idx(i))));
end



annots = CoExp_All.GABA.Annotations;
idx_annots_new = find(ismember(annots, {'GABA-Transaminase Deficiency'...
    'GABA receptor activity'}));
annots = annots(idx_annots_new);
annot_names = [annots, {'SCZ GWAS'}];
annots = [annots, {'SCZ_GWAS_2022'}];
t = tiledlayout(3,3);

for i_annot = 1:length(annots)
    
    nexttile()
    
    
    if fc_EC_vs_EC == 1
        if i_annot ~= 3
            tmp = squeeze(CoExp_All.GABA.m_NVC_who.r(idx_annots_new(i_annot), :, :));
        else
            tmp = squeeze(CoExp_All.SCZ_GWAS_2022.m_NVC_who.r(1, :, :));
        end
    elseif fc_EC_vs_EC == 0
        if i_annot ~= 3
            tmp = squeeze(CoExp_All.GABA.m_EC_to_DLPFC.r(idx_annots_new(i_annot), :, :));
        else
            tmp = squeeze(CoExp_All.SCZ_GWAS_2022.m_EC_to_DLPFC.r(1, :, :));
        end
    end
    
    supp = (tmp(3,:)' + tmp(2,:)')./2;
    deep = (tmp(5,:)' + tmp(6,:)')./2;
    layer_effect = supp - deep;
    
    [~,p_hc] = ttest(supp(HC_idx), deep(HC_idx));
    NVC_frq_mean(1,1) = nanmean(supp(HC_idx));
    NVC_frq_mean(1,2) = nanmean(deep(HC_idx));
    NVC_frq_error(1,1) = nanstd(supp(HC_idx))/sqrt(length(supp(HC_idx)));
    NVC_frq_error(1,2) = nanstd(deep(HC_idx))/sqrt(length(deep(HC_idx)));
    [~,p_scz] = ttest(supp(SCZ_idx), deep(SCZ_idx));
    NVC_frq_mean(2,1) = nanmean(supp(SCZ_idx));
    NVC_frq_mean(2,2) = nanmean(deep(SCZ_idx));
    NVC_frq_error(2,1) = nanstd(supp(SCZ_idx))/sqrt(length(supp(SCZ_idx)));
    NVC_frq_error(2,2) = nanstd(deep(SCZ_idx))/sqrt(length(deep(SCZ_idx)));
    
    %     [~,p_group] = ttest2(supp(HC_idx) - deep(HC_idx), ...
    %         supp(SCZ_idx) - deep(SCZ_idx));
    
    
    bev_cov = table2array(bev(:, 6:8));
    bev_cov = bev_cov([HC_idx_bev; SCZ_idx_bev], :);
    tmp_HC = layer_effect(HC_idx);
    tmp_SCZ = layer_effect(SCZ_idx);
    
    
    tmp_bev = array2table([[tmp_HC; tmp_SCZ], ...
        [ones(length(tmp_HC'),1); zeros(length(tmp_SCZ'),1)], ...
        bev_cov(:,1),bev_cov(:,2),bev_cov(:,3)], ...
        'VariableNames', {'RES', 'Group', 'age', 'gender', 'education'});
    tmp_bev.Group = categorical(tmp_bev.Group);
    tmp_bev.gender = categorical(tmp_bev.gender);
    
    md = fitlm(tmp_bev,'RES~Group + age + gender + education');
    
    p_group = md.Coefficients.pValue(2);
    
    %%%%%%%%% plot
    bc = bar(NVC_frq_mean, 'FaceAlpha', 0.5, 'EdgeColor','none');
    hold on
    
    ngroups = size(NVC_frq_mean, 1);
    nbars = size(NVC_frq_mean, 2);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        if mean(NVC_frq_mean(:,i)) < 0
            errorbar(x, NVC_frq_mean(:,i), NVC_frq_error(:,i), [], ...
                'k', 'linestyle', 'none', 'LineWidth', 1.5);
        elseif mean(NVC_frq_mean(:,i)) > 0
            errorbar(x, NVC_frq_mean(:,i), [], NVC_frq_error(:,i), ...
                'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    y_range = ylim;
    if i_annot ~= 3
        if abs(y_range(2)) > abs(y_range(1))
            y_range(2) = y_range(2)*1.3;
            y_loc = (y_range(2) -  y_range(1)) * 3.2/4;
        elseif abs(y_range(1)) > abs(y_range(2))
            y_range(1) = y_range(1)*1.3;
            y_loc = (y_range(1) -  y_range(2)) * 3.2/4;
        end
    else
        if abs(y_range(2)) > abs(y_range(1))
            y_range(2) = y_range(2)*1.5;
            y_loc = (y_range(2) -  y_range(1)) * 2.8/4;
        elseif abs(y_range(1)) > abs(y_range(2))
            y_range(1) = y_range(1)*1.5;
            y_loc = (y_range(1) -  y_range(2)) * 2.8/4;
        end
    end
    
    set(gca, 'TickLength', [0,0], ...
        'ylim', y_range, ...
        'xticklabel', {'HV', 'SCZ'}, ...
        'xColor', [0, 0, 0], 'yColor', [0, 0, 0], 'FontWeight', 'normal', 'FontSize',10);
    
    %     title([labels{i_annot} ' ' annots_name{i_annot}], 'FontSize',10, 'Color', 'k');
    %     title(labels{i_annot}, 'FontSize',10, 'Color', 'k');
    %     s = sprintf(['h' num2str(i_annot) '.TitleHorizontalAlignment = ''%s'';'], ...
    %         'left');
    %     eval(s);
    %     hold on
    
    %     xlabel(['Gene Co-Expression on ' annots_name{i_annot}], 'Color', [0, 0, 0], 'FontSize',10, 'FontWeight', 'normal');
    xlabel(annot_names{i_annot}, 'Color', [0, 0, 0], 'FontSize',10, 'FontWeight', 'normal');
    
    if fc_EC_vs_EC == 1
        ylabel({'Correlation between co-expression', ...
            'and fc-ECs to DPLFC in WM'}, 'Color', [0, 0, 0], 'FontSize',10, 'FontWeight', 'normal');
    elseif fc_EC_vs_EC == 0
        ylabel({'Correlation between co-expression', ...
            'and ECs to DPLFC in WM'}, 'Color', [0, 0, 0], 'FontSize',10, 'FontWeight', 'normal');
    end
    
    % make laminar effect in HC
    %     text(1 - 0.25,y_loc*1.07,['{\it p} = ' num2str(round(p_hc,round_num(i_annot)))], 'FontSize',10, 'FontWeight', 'normal', 'Color', 'k');
    %     hold on
    %     x = [0.85 1.15];
    %     y = [y_loc, y_loc];
    %     plot(x, y, 'k', 'LineWidth', 2, 'HandleVisibility','off', 'Color', 'k');
    %     hold on
    if p_hc > 0.05
        text(1 - 0.1,y_loc*1.1,'{\it n.s}', 'FontSize',10, 'FontWeight', 'normal', 'Color', 'b');
        hold on
        x = [0.75 1.25];
        y = [y_loc, y_loc];
        plot(x, y, 'b', 'LineWidth', 2, 'HandleVisibility','off');
        hold on
    elseif p_hc <= 0.05 && p_hc >= 0.001
        text(1 - 0.25,y_loc*1.1, ['{\itp} = ' num2str(round(p_hc, 3))], 'FontSize',10, 'FontWeight', 'normal', 'Color', 'r');
        hold on
        x = [0.75 1.25];
        y = [y_loc, y_loc];
        plot(x, y, 'r', 'LineWidth', 2, 'HandleVisibility','off');
        hold on
    elseif p_hc < 0.001 && p_hc >= 1e-08
        tmp_p = num2str(p_hc, '%.e');
        text(1 - 0.25,y_loc*1.1,['{\itp} = ' tmp_p], 'FontSize',10, 'FontWeight', 'normal', 'Color', 'r');
        hold on
        x = [0.75 1.25];
        y = [y_loc, y_loc];
        plot(x, y, 'r', 'LineWidth', 2, 'HandleVisibility','off');
        hold on
        clear tmp_p*
    elseif p_hc < 1e-08
        text(1 - 0.25,y_loc*1.1, '{\itp} < 1e-08', 'FontSize',10, 'FontWeight', 'normal', 'Color', 'r');
        hold on
        x = [0.75 1.25];
        y = [y_loc, y_loc];
        plot(x, y, 'r', 'LineWidth', 2, 'HandleVisibility','off');
        hold on
    end
    
    % make laminar effect in SCZ
    %     text(2 - 0.25,y_loc*1.07,['{\it p} = ' num2str(round(p_scz,3))], 'FontSize',10, 'FontWeight', 'normal', 'Color', 'k');
    %     hold on
    %     x = [1.85 2.15];
    %     y = [y_loc, y_loc];
    %     plot(x, y, 'k', 'LineWidth', 2, 'HandleVisibility','off', 'Color', 'k');
    %     hold on
    if p_scz > 0.05
        text(2 - 0.1,y_loc*1.1,'{\it n.s}', 'FontSize',10, 'FontWeight', 'normal', 'Color', 'b');
        hold on
        x = [1.75 2.25];
        y = [y_loc, y_loc];
        plot(x, y, 'b', 'LineWidth', 2, 'HandleVisibility','off');
        hold on
    elseif p_scz <= 0.05 && p_scz >= 0.001
        text(2 - 0.25,y_loc*1.1, ['{\itp} = ' num2str(round(p_scz, 3))], 'FontSize',10, 'FontWeight', 'normal', 'Color', 'r');
        hold on
        x = [1.75 2.25];
        y = [y_loc, y_loc];
        plot(x, y, 'r', 'LineWidth', 2, 'HandleVisibility','off');
        hold on
    elseif p_scz < 0.001 && p_scz >= 1e-08
        tmp_p = num2str(p_scz, '%.e');
        text(2 - 0.25,y_loc*1.1,['{\itp} = ' tmp_p], 'FontSize',10, 'FontWeight', 'normal', 'Color', 'r');
        hold on
        x = [1.75 2.25];
        y = [y_loc, y_loc];
        plot(x, y, 'r', 'LineWidth', 2, 'HandleVisibility','off');
        hold on
        clear tmp_p*
    elseif p_scz < 1e-08
        text(2 - 0.25,y_loc*1.1, '{\itp} < 1e-08', 'FontSize',10, 'FontWeight', 'normal', 'Color', 'r');
        hold on
        x = [1.75 2.25];
        y = [y_loc, y_loc];
        plot(x, y, 'r', 'LineWidth', 2, 'HandleVisibility','off');
        hold on
    end
    
    %%%%%%%%%
    
    if i_annot == 3
        legend({'Layer 2/3', 'Layer 5/6'}, ...
            'FontSize',10, 'location', 'northeastoutside', ...
            'box', 'off', 'NumColumns',1);
    end
    box on
    
end


%%%%%%%%%%%%%%% plot subplots 4 ~ 6%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except annots bev CoExp_All fc_EC_vs_EC HC HC_idx HC_idx_bev idx_annots_new SCZ SCZ_idx SCZ_idx_bev idx_annots_new annot_names
for i_annot = 1:length(annots)
    
    nexttile()
    
    
    if fc_EC_vs_EC == 1
        if i_annot ~= 3
            tmp = squeeze(CoExp_All.GABA.m_NVC_who.r(idx_annots_new(i_annot), :, :));
        else
            tmp = squeeze(CoExp_All.SCZ_GWAS_2022.m_NVC_who.r(1, :, :));
        end
    elseif fc_EC_vs_EC == 0
        if i_annot ~= 3
            tmp = squeeze(CoExp_All.GABA.m_EC_to_DLPFC.r(idx_annots_new(i_annot), :, :));
        else
            tmp = squeeze(CoExp_All.SCZ_GWAS_2022.m_EC_to_DLPFC.r(1, :, :));
        end
    end
    
    supp = (tmp(3,:)' + tmp(2,:)')./2;
    deep = (tmp(5,:)' + tmp(6,:)')./2;
    layer_effect = supp - deep;
    
    NVC_frq_mean(1,1) = nanmean(layer_effect(HC_idx));
    NVC_frq_error(1,1) = nanstd(layer_effect(HC_idx))/sqrt(length(layer_effect(HC_idx)));
    NVC_frq_mean(1,2) = nanmean(layer_effect(SCZ_idx));
    NVC_frq_error(1,2) = nanstd(layer_effect(SCZ_idx))/sqrt(length(layer_effect(SCZ_idx)));
    
    
    bev_cov = table2array(bev(:, 6:8));
    bev_cov = bev_cov([HC_idx_bev; SCZ_idx_bev], :);
    tmp_HC = layer_effect(HC_idx);
    tmp_SCZ = layer_effect(SCZ_idx);
    
    
    tmp_bev = array2table([[tmp_HC; tmp_SCZ], ...
        [ones(length(tmp_HC'),1); zeros(length(tmp_SCZ'),1)], ...
        bev_cov(:,1),bev_cov(:,2),bev_cov(:,3)], ...
        'VariableNames', {'RES', 'Group', 'age', 'gender', 'education'});
    tmp_bev.Group = categorical(tmp_bev.Group);
    tmp_bev.gender = categorical(tmp_bev.gender);
    
    md = fitlm(tmp_bev,'RES~Group + age + gender + education');
    
    p_group = md.Coefficients.pValue(2);
    
    %%%%%%%% plot
    bc = bar(diag(NVC_frq_mean,0), 'stacked', 'FaceAlpha', 0.5, 'EdgeColor','none');
    bc(1).FaceColor = [0.4660 0.6740 0.1880];
    bc(2).FaceColor = [0.6350 0.0780 0.1840];
    hold on
%     h1=bar(1, NVC_frq_mean(1), 'FaceAlpha', 0, 'EdgeColor','#0072BD', 'LineWidth', 1);
%     hold on
%     h2=bar(2, NVC_frq_mean(2), 'FaceAlpha', 0, 'EdgeColor','#D95319', 'LineWidth', 1);
%     hold on
    
    
    if mean(NVC_frq_mean) < 0
        errorbar([1,2], NVC_frq_mean, NVC_frq_error, [], ...
            'k', 'linestyle', 'none', 'LineWidth', 1.5);
    elseif mean(NVC_frq_mean) > 0
        errorbar([1,2], NVC_frq_mean, [], NVC_frq_error, ...
            'k', 'linestyle', 'none', 'LineWidth', 1.5);
    end
    
    y_range = ylim;
    if i_annot ~= 3
        if abs(y_range(2)) > abs(y_range(1))
            y_range(2) = y_range(2)*1.5;
            y_loc = (y_range(2) -  y_range(1)) * 4/5;
        elseif abs(y_range(1)) > abs(y_range(2))
            y_range(1) = y_range(1)*1.5;
            y_loc = (y_range(1) -  y_range(2)) * 4/5;
        end
    else
        if abs(y_range(2)) > abs(y_range(1))
            y_range(2) = y_range(2)*1.7;
            y_loc = (y_range(2) -  y_range(1)) * 4/5;
        elseif abs(y_range(1)) > abs(y_range(2))
            y_range(1) = y_range(1)*1.7;
            y_loc = (y_range(1) -  y_range(2)) * 4/5;
        end
    end
    
    set(gca, 'TickLength', [0,0], ...
        'ylim', y_range, ...
        'xticklabel', {'HV', 'SCZ'}, ...
        'xColor', [0, 0, 0], 'yColor', [0, 0, 0], 'FontWeight', 'normal', 'FontSize',10);
    
    xlabel(annot_names{i_annot}, 'Color', [0, 0, 0], 'FontSize',10, 'FontWeight', 'normal');
    
    if fc_EC_vs_EC == 1
        ylabel({'Laminar difference in correlation of', ...
            'fc-EC-transcriptomic mapping'}, ...
            'Color', [0, 0, 0], 'FontSize',10, 'FontWeight', 'normal');
    elseif fc_EC_vs_EC == 0
        ylabel({'Laminar difference in correlation of', ...
            'EC-transcriptomic mapping'}, ...
            'Color', [0, 0, 0], 'FontSize',10, 'FontWeight', 'normal');
    end
    
    if p_group > 0.05
        text(1.5 - 0.16,y_loc*1.1,'{\it n.s}', 'FontSize',10, 'FontWeight', 'normal', 'Color', 'b');
        hold on
        x = [1 2];
        y = [y_loc, y_loc];
        plot(x, y, 'b', 'LineWidth', 2, 'HandleVisibility','off');
        hold on
    elseif p_group <= 0.05 && p_group >= 0.001
        text(1.5 - 0.4,y_loc*1.1, ['{\itp} = ' num2str(round(p_group, 3))], 'FontSize',10, 'FontWeight', 'normal', 'Color', 'r');
        hold on
        x = [1 2];
        y = [y_loc, y_loc];
        plot(x, y, 'r', 'LineWidth', 2, 'HandleVisibility','off');
        hold on
    elseif p_group < 0.001
        tmp_p = num2str(p_group, '%.e');
        text(1.5 - 0.4,y_loc*1.1,['{\itp} = ' tmp_p], 'FontSize',10, 'FontWeight', 'normal', 'Color', 'r');
        hold on
        x = [1 2];
        y = [y_loc, y_loc];
        plot(x, y, 'r', 'LineWidth', 2, 'HandleVisibility','off');
        hold on
        clear tmp_p*
    end
    %%%%%%%%
    
    if i_annot == 3
        legend({'HV', 'SCZ'}, ...
            'FontSize',10, 'location', 'northeastoutside', ...
            'box', 'off', 'NumColumns',1);
    end
    box on
    
end







%%%%%%%%%%%%%  plot correlation %%%%%%%%%%%%%%%


bev_idx = 15;
bev_names = bev.Properties.VariableNames(bev_idx);

colors = {'#0072BD', '#D95319', '#7E2F8E'};
colors_vals = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.4940 0.1840 0.5560]];

for i_annot = 1:length(annots)
    clear tmp* supp deep layer_effect

    if fc_EC_vs_EC == 1
        if i_annot ~= 3
            tmp = squeeze(CoExp_All.GABA.m_NVC_who.r(idx_annots_new(i_annot), :, :));
        else
            tmp = squeeze(CoExp_All.SCZ_GWAS_2022.m_NVC_who.r(1, :, :));
        end
    elseif fc_EC_vs_EC == 0
        if i_annot ~= 3
            tmp = squeeze(CoExp_All.GABA.m_EC_to_DLPFC.r(idx_annots_new(i_annot), :, :));
        else
            tmp = squeeze(CoExp_All.SCZ_GWAS_2022.m_EC_to_DLPFC.r(1, :, :));
        end
    end

    tmp = tmp(:, HC_idx);
    supp = (tmp(3,:)' + tmp(2,:)')./2;
    deep = (tmp(5,:)' +  tmp(6,:)')./2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    layer_effect = supp - deep;
%     layer_effect = deep - supp;
    bev_dat = table2array(bev(HC_idx_bev, bev_idx));
    bev_cov = table2array(bev(HC_idx_bev, 6:8));
    idx_ex = find(isnan(mean([layer_effect, bev_dat],2)));
    layer_effect(idx_ex,:) = [];
    bev_dat(idx_ex, :) = [];
    bev_cov(idx_ex, :) = [];


%     [R,P] = corr(layer_effect, bev_dat, "rows","pairwise");

    %%%%%%%%%%%%%%   partial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [rho,pval] = partialcorr([layer_effect, bev_dat],...
        bev_cov,"rows","complete");
    R = rho(1,2:end);
    P = pval(1, 2:end);

    %%%%%%%%%%%%%%%%% plot

    eval(['h', num2str(i_annot + 3), ' = nexttile();'])
%     title(labels{i}, 'FontSize',10, 'Color', 'k');
%     s = sprintf(['h' num2str(i) '.TitleHorizontalAlignment = ''%s'';'], ...
%         'left');
%     eval(s);
%     hold on

    X = zscore(bev_dat);
    Y = zscore(layer_effect);
%     [r, p] = corr(X, Y, "rows","complete");
%     [rho,pval] = partialcorr([X, Y],bev_cov,"rows","complete");
%     r = rho(1,2:end);
%     p = pval(1, 2:end);
    
    % First the first plot.
    coefficients = polyfit(X, Y, 1);
    % Make new x coordinates
    x = linspace(min(X), max(X), 100);
    y = polyval(coefficients, x);
    s = sprintf(['h' num2str(i_annot) ' = plot(%s, ''%s'', ''%s'', ''%s'', ''%s'', %s)'], ...
        'x, y', '-', 'Color', colors{i_annot}, 'LineWidth', num2str(2));
    eval(s);
    %                     h = plot(x, y, '-', 'Color', colors{i_layer}, 'LineWidth', 2);
    hold on;

    % Plot scatter
    scatter(X, Y, ...
        50, colors_vals(i_annot,:), ...
        'LineWidth',0.2);
    hold on;
    scatter(X, Y, ...
        50, 'filled', ...
        'MarkerFaceColor',colors{i_annot}, ...
        'MarkerFaceAlpha', 0.4, ...
        'LineWidth',1);
    hold on;

    ylabel('PRS for SCZ at {\itp} < 0.01 (  )', ...
        'Color', [0, 0, 0], 'FontSize',10, 'FontWeight', 'normal');
%     xlabel(['Layer Effect on ' annots_name{i_annot} ' (\sigma)'], ...
%         'Color', [0, 0, 0], 'FontSize',10, 'FontWeight', 'normal');
    if fc_EC_vs_EC == 1
        xlabel({'Laminar difference in correlation of', ...
            'fc-EC-transcriptomic mapping (  )'}, ...
            'Color', [0, 0, 0], 'FontSize',10, 'FontWeight', 'normal');
    elseif fc_EC_vs_EC == 0
        xlabel({'Laminar difference in correlation of', ...
            'EC-transcriptomic mapping (  )'}, ...
            'Color', [0, 0, 0], 'FontSize',10, 'FontWeight', 'normal');
    end

    hold on;
    dummyh= line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
    if P > 0.5
        legend(dummyh, {['\color[rgb]{0,0,1} ' '{\it r} = ' num2str(round(R,3)) ', ' '{\it p} = ' num2str(round(P, 3))]}, ...
            'FontSize',10, ...
            'Location', 'southeast', 'Box', 'off', 'FontWeight', 'normal');
    else
        legend(dummyh, {['\color[rgb]{1,0,0} ' '{\it r} = ' num2str(round(R,3)) ', ' '{\it p} = ' num2str(round(P, 3))]}, ...
            'FontSize',10, ...
            'Location', 'southeast', 'Box', 'off', 'FontWeight', 'normal', 'Color', 'r');
    end

    set(gca, 'xColor', 'k', 'yColor', 'k', 'ylim', [-3.9, 3], 'TickLength', ...
        [0,0], 'FontWeight', 'normal', 'FontSize',10);
    box on

end
clear annots

annotation('textbox',...
    [0.1 0.94 0.15 0.02],...
    'String','(A)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'FontSize',10);
annotation('textbox',...
    [0.385 0.94 0.15 0.02],...
    'String','(B)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'FontSize',10);
annotation('textbox',...
    [0.665 0.94 0.15 0.02],...
    'String','(C)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'FontSize',10);
annotation('textbox',...
    [0.1 0.645 0.15 0.02],...
    'String','(D)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'FontSize',10);
annotation('textbox',...
    [0.385 0.645 0.15 0.02],...
    'String','(E)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'FontSize',10);
annotation('textbox',...
    [0.665 0.645 0.15 0.02],...
    'String','(F)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'FontSize',10);
annotation('textbox',...
    [0.1 0.345 0.15 0.02],...
    'String','(G)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'FontSize',10);
annotation('textbox',...
    [0.385 0.345 0.15 0.02],...
    'String','(H)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'FontSize',10);
annotation('textbox',...
    [0.665 0.345 0.15 0.02],...
    'String','(I)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'FontSize',10);

annotation(gcf,'textbox',[0.097 0.297 0.05 0.03],...
    'LineStyle','none',...
    'Rotation', 90, ...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');

annotation(gcf,'textbox',[0.381 0.297 0.05 0.03],...
    'LineStyle','none',...
    'Rotation', 90, ...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');

annotation(gcf,'textbox',[0.665 0.297 0.05 0.03],...
    'LineStyle','none',...
    'Rotation', 90, ...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');

annotation(gcf,'textbox',[0.287 0.037 0.05 0.03],...
    'LineStyle','none',...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');

annotation(gcf,'textbox',[0.571 0.037 0.05 0.03],...
    'LineStyle','none',...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');

annotation(gcf,'textbox',[0.855 0.037 0.05 0.03],...
    'LineStyle','none',...
    'String','$\it{z}$','Interpreter','latex',...
    'FontSize', 14, 'FontWeight', 'bold');

set(gcf,'color','w');
print(gcf,'./Figure6.png','-dpng','-r300');
close gcf
