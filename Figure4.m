clear
%%%%%%%%%%%%%%%%%%%  plot gene expression map %%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf,'units','normalized','position',[0 0 0.55 0.75])

axisleft1 = [0.08, 0.36];
axisbottom1 = [0.65, 0.65];
axiswidth1 = [0.24, 0.24];
axisheight1 = [0.36, 0.36];

figs = {'fc_EC_2_DLPFC_SCZ_lateral.png', 'fc_EC_2_DLPFC_SCZ_medial.png'};

h(1)=subplot('position', [axisleft1(1), axisbottom1(1), axiswidth1(1), axisheight1(1)]);
y = imread(['.\brain_map_FS/' figs{1}], 'BackgroundColor', [1 1 1]);
y(:,[1:100, 1100:end],:) = [];
y([1:70,780:end],:,:) = [];
imshow(y);

h(2)=subplot('position', [axisleft1(2), axisbottom1(2), axiswidth1(2), axisheight1(2)]);
y = imread(['.\brain_map_FS/' figs{2}], 'BackgroundColor', [1 1 1]);
y(:,[1:100, 1100:end],:) = [];
y([1:50,760:end],:,:) = [];
imshow(y);


%%%%%%%%%%%%%%%%%%%%%%%% subplot size %%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
axisleft = [0.72, 0.11, 0.45];
axisbottom = [0.56, 0.075, 0.075];
axiswidth = [0.25, 0.25, 0.52];
axisheight = [0.41, 0.41, 0.41];

%%%%%%%%%%%%%%%% subplot 3  group difference on gamma oscillation %%%%%%%%%
clearvars -except  axisleft axisbottom axiswidth axisheight
bev = readtable('./Restricted_Access/bev_final.csv');
prs_sczs = bev.Properties.VariableNames(14:23);
prs_sczs = prs_sczs([10, 9, 8, 1:7]);
load('./Restricted_Access/NB_2back.mat');
HC = bev.subjects_MEG(find(ismember(bev.group, 'HC')));
SCZ = bev.subjects_MEG(find(ismember(bev.group, 'SCZ')));
HC_idx = find(ismember(NB_2back.subjects, HC));
SCZ_idx = find(ismember(NB_2back.subjects, SCZ));
bd_range = [1 120; 1 4; 4 8; 8 12; 12 30; 30 70; 70, 120];
x_label = {'Whole Band\newline(1~120Hz)', '   Delta\newline (1~4Hz)', ...
    '    Theta\newline (4~8Hz)', '    Alpha\newline (8~12Hz)', '     Beta\newline(12~30Hz)', ...
    'Low Gamma\newline(30~70Hz)', 'High Gamma\newline(70~120Hz)'};
bd_range = [1 120; 1 4; 4 8; 8 12; 12 30; 30 120];
x_label = {'Whole Band\newline(1~120Hz)', '   Delta\newline (1~4Hz)', ...
    '    Theta\newline (4~8Hz)', '    Alpha\newline (8~12Hz)', '     Beta\newline(12~30Hz)', ...
    '    Gamma\newline(30~120Hz)'};
HC_NVC_frq = [NB_2back.HC_B_forward_NVC_frequency, NB_2back.HC_B_backward_NVC_frequency];
SCZ_NVC_frq = [NB_2back.SCZ_B_forward_NVC_frequency, NB_2back.SCZ_B_backward_NVC_frequency];


%%%%%%%%%%%%%%%%%%%%% only the afferent to DLPFC %%%%%%%%%%%%%%%%%%%%%%%%%%

idx_DLPFC = 1; % 1 = only the afferent to DLPFC; 1 = all the EC
load('./NB_2back.mat');
par.DLPFC = {'GM_067R_8Av'
    'GM_070R_8BL'
    'GM_071R_9p'
    'GM_073R_8c'
    'GM_083R_p9-46v'
    'GM_084R_46'
    'GM_085R_a9-46v'
    'GM_086R_9-46d'
    'GM_087R_9a'
    'GM_097R_i6-8'
    'GM_098R_s6-8'};
par.ec_fmri_name_all = [NB_2back.fc_name(:,2:3);NB_2back.fc_name(:,3:-1:2)];
par.idx_dlpfc = find(ismember(par.ec_fmri_name_all(:,2), par.DLPFC)); % afferent

if idx_DLPFC == 1
    bd_range = [12 30; 30 120];
    x_label = {'Beta (12~30Hz)', 'Gamma (30~120Hz)'};
    y_label = {'beta', 'gamma'};
    titles = {'(A)', '(B)', '(C)', '(D)'};
    HC_NVC_frq = HC_NVC_frq(par.idx_dlpfc);
    SCZ_NVC_frq = SCZ_NVC_frq(par.idx_dlpfc);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 2:2


    h(3)=subplot('position', [axisleft(1), axisbottom(1), axiswidth(1), axisheight(1)]);


    eval(['tmp1 = HC_NVC_frq(find(HC_NVC_frq> ', num2str(bd_range(i,1)),...
        ' & HC_NVC_frq < ', num2str(bd_range(i,2)), '));'])
    eval(['tmp2 = SCZ_NVC_frq(find(SCZ_NVC_frq> ', num2str(bd_range(i,1)),...
        ' & SCZ_NVC_frq < ', num2str(bd_range(i,2)), '));'])
    [~,p] = ttest2(tmp1, tmp2);
    NVC_frq_mean(1,1) = nanmean(tmp1);
    NVC_frq_mean(1,2) = nanmean(tmp2);
    NVC_frq_error(1,1) = nanstd(tmp1)/sqrt(length(tmp1));
    NVC_frq_error(1,2) = nanstd(tmp2)/sqrt(length(tmp2));
    
    data(i, 1:2) = NVC_frq_mean;
    data(i, 3) = p;
    data(i, 4:5) = [length(tmp1), length(tmp2)];
    
    bc = bar(diag(NVC_frq_mean,0), 'stacked', 'FaceAlpha', 0.5, 'EdgeColor','none');
    bc(1).FaceColor = [0.4660 0.6740 0.1880];
    bc(2).FaceColor = [0.6350 0.0780 0.1840];
    hold on
    h1=bar(1, NVC_frq_mean(1), 'FaceAlpha', 0, 'EdgeColor',[0.4660 0.6740 0.1880], 'LineWidth', 1);
    hold on
    h2=bar(2, NVC_frq_mean(2), 'FaceAlpha', 0, 'EdgeColor',[0.6350 0.0780 0.1840], 'LineWidth', 1);
    hold on
    
    errorbar(1:2, NVC_frq_mean, [], NVC_frq_error, ...
        'k', 'linestyle', 'none', 'LineWidth', 1.5);
    hold off
    
    if i == length(x_label)
        legend({'HV', 'SCZ'}, ...
            'FontSize',11, 'location', 'northeast', ...
            'box', 'off', 'NumColumns',1);
    end
    xlabel(x_label{i}, 'Color', [0, 0, 0], 'FontSize',11, 'FontWeight', 'normal');
    ylabel({'WM-related fc-EC to DLPFC (Hz)'}, 'Color', [0, 0, 0], 'FontSize',11, 'FontWeight', 'normal');

    
    set(gca, 'TickLength', [0,0], ...
        'ylim', bd_range(i,:), ...
        'xticklabel', '', ...
        'xColor', [0, 0, 0], 'yColor', [0, 0, 0], 'FontWeight', 'normal', 'FontSize',11);
    
    loc_lab = (bd_range(i,2) - bd_range(i,1)) * 0.77 + bd_range(i,1);
    

    if p > 0.05
        text(1.3,loc_lab,'{\it n.s}', 'FontSize',11, 'FontWeight', 'normal', 'Color', 'b');
        hold on
        loc_lab = (bd_range(i,2) - bd_range(i,1)) * 0.72 + bd_range(i,1);
        x = [1 2];
        y = [loc_lab, loc_lab];
        plot(x, y, 'b', 'LineWidth', 2, 'HandleVisibility','off');
        hold on
    elseif p <= 0.05 && p >= 0.001
        text(1.25,loc_lab, ['{\itp} =' num2str(round(p, 3))], 'FontSize',11, 'FontWeight', 'normal', 'Color', 'r');
        hold on
        loc_lab = (bd_range(i,2) - bd_range(i,1)) * 0.72 + bd_range(i,1);
        x = [1 2];
        y = [loc_lab, loc_lab];
        plot(x, y, 'r', 'LineWidth', 2, 'HandleVisibility','off');
        hold on
    elseif p < 0.001
        tmp_p = num2str(p, '%.e');
        text(1.05,loc_lab,['{\itp} = ' tmp_p], 'FontSize',11, 'FontWeight', 'normal', 'Color', 'r');
        hold on
        loc_lab = (bd_range(i,2) - bd_range(i,1)) * 0.72 + bd_range(i,1);
        x = [1 2];
        y = [loc_lab, loc_lab];
        plot(x, y, 'r', 'LineWidth', 2, 'HandleVisibility','off');
        hold on
        clear tmp_p*
    end
    
    hold on
    box on
end



%%%%%%%%%%% subplot 5  group difference on gamma oscillation in PRS %%%%%%%

% bd_range = [12 30; 45 90];


for i_prs = 1:length(prs_sczs)
    load(['./NVC_PRS_grouplevel/' ...
        prs_sczs{i_prs} '_lprs_vs_hprs.mat']);
    eval(['HC_lprs_NVC_frq = [', prs_sczs{i_prs}, '.HC_lprs_B_forward_NVC_frequency, ', ...
        prs_sczs{i_prs}, '.HC_lprs_B_backward_NVC_frequency];'])
    eval(['HC_mprs_NVC_frq = [', prs_sczs{i_prs}, '.HC_mprs_B_forward_NVC_frequency, ', ...
        prs_sczs{i_prs}, '.HC_mprs_B_backward_NVC_frequency];'])
    eval(['HC_hprs_NVC_frq = [', prs_sczs{i_prs}, '.HC_hprs_B_forward_NVC_frequency, ', ...
        prs_sczs{i_prs}, '.HC_hprs_B_backward_NVC_frequency];'])

    if idx_DLPFC == 1
        HC_lprs_NVC_frq = HC_lprs_NVC_frq(par.idx_dlpfc);
        HC_mprs_NVC_frq = HC_mprs_NVC_frq(par.idx_dlpfc);
        HC_hprs_NVC_frq = HC_hprs_NVC_frq(par.idx_dlpfc);
    end

    for i = 1:length(bd_range)
        eval(['tmp1 = HC_lprs_NVC_frq(find(HC_lprs_NVC_frq> ', num2str(bd_range(i,1)),...
            ' & HC_lprs_NVC_frq < ', num2str(bd_range(i,2)), '));'])
        eval(['tmp2 = HC_mprs_NVC_frq(find(HC_mprs_NVC_frq> ', num2str(bd_range(i,1)),...
            ' & HC_mprs_NVC_frq < ', num2str(bd_range(i,2)), '));'])
        eval(['tmp3 = HC_hprs_NVC_frq(find(HC_hprs_NVC_frq> ', num2str(bd_range(i,1)),...
            ' & HC_hprs_NVC_frq < ', num2str(bd_range(i,2)), '));'])
        [~,p32] = ttest2(tmp3, tmp2);
        [~,p31] = ttest2(tmp3, tmp1);
        [~,p21] = ttest2(tmp2, tmp1);

        NVC_frq_low(i,i_prs) = nanmean(tmp1);
        NVC_frq_low_error(i,i_prs) = nanstd(tmp1)/sqrt(length(tmp1));
        NVC_frq_low_N(i,i_prs) = length(tmp1);

        NVC_frq_med(i,i_prs) = nanmean(tmp2);
        NVC_frq_med_error(i,i_prs) = nanstd(tmp2)/sqrt(length(tmp2));
        NVC_frq_med_N(i,i_prs) = length(tmp2);

        NVC_frq_high(i,i_prs) = nanmean(tmp3);
        NVC_frq_high_error(i,i_prs) = nanstd(tmp3)/sqrt(length(tmp3));
        NVC_frq_high_N(i,i_prs) = length(tmp3);

        NVC_frq_p32(i,i_prs) = p32;
        NVC_frq_p31(i,i_prs) = p31;
        NVC_frq_p21(i,i_prs) = p21;
    end
end

for i = 2:2
    h(5)=subplot('position', [axisleft(3), axisbottom(3), axiswidth(3), axisheight(3)]);
    errorbar(1:10,NVC_frq_low(i,:),NVC_frq_low_error(i,:), 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880]);
    hold on

    errorbar(1:10,NVC_frq_high(i,:),NVC_frq_high_error(i,:), 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840]);
    hold on

    if i == 2
        legend({'LPRS', 'HPRS'}, ...
            'FontSize',11, 'location', 'northeast', ...
            'box', 'off', 'NumColumns',1);
    end

    xlabel('PRS for Schizophrenia at different {\itp} thresholds', 'Color', [0, 0, 0], 'FontSize',11, 'FontWeight', 'normal');

    if i == 1
    ylabel('WM-related fc-EC to DLPFC (Hz)', 'Color', [0, 0, 0], 'FontSize',11, 'FontWeight', 'normal');
    elseif i ==2
    ylabel('WM-related fc-EC to DLPFC (Hz)', 'Color', [0, 0, 0], 'FontSize',11, 'FontWeight', 'normal');
    end

    
    set(gca, 'TickLength', [0,0], ...
        'ylim', bd_range(i,:), ...
        'xlim', [0.5,10.5], ...
        'xtick', [1:10], ...
        'xticklabel', {'5e-08', '1e-06', '1e-04', '0.001', '0.01', '0.05', ...
        '0.1', '0.2', '0.5', '1'}, ...
        'xColor', [0, 0, 0], 'yColor', [0, 0, 0], 'FontWeight', 'normal', 'FontSize',11);


    loc_lab = (bd_range(i,2) - bd_range(i,1)) * 0.68 + bd_range(i,1);
    for j = 1:10
        if NVC_frq_p31(i,j) < 0.05 && NVC_frq_p31(i,j) >= 0.001
            text(j - 0.07,loc_lab,'*', 'FontSize',11, 'FontWeight', 'normal', 'Color', 'r');
            hold on
        elseif NVC_frq_p31(i,j) < 0.001 && NVC_frq_p31(i,j) >= 1e-08
            text(j - 0.1,loc_lab,'**', 'FontSize',11, 'FontWeight', 'normal', 'Color', 'r');
            hold on
        elseif NVC_frq_p31(i,j) < 1e-08
            text(j - 0.15,loc_lab,'***', 'FontSize',11, 'FontWeight', 'normal', 'Color', 'r');
            hold on
        else
            text(j - 0.3,loc_lab + 0.2,'{\it n.s}', 'FontSize',11, 'FontWeight', 'normal', 'Color', 'b');
            hold on
        end
    end


    hold on
end




%%%%%%%%%%%%%%%%%%%%%% subplot 3 group difference and co-expression %%%%%%%
clearvars -except  axisleft axisbottom axiswidth axisheight
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

h(4)=subplot('position', [axisleft(2), axisbottom(2), axiswidth(2), axisheight(2)]);

fc_var_name = {'fc_var_group_diff'};
clear p_diff layer_mean

SUPP = [];
DEEP = [];
G_names_supp = [];
G_names_deep = [];


clear layer_mean layer_error
for i = 1:length(fc_var_name)
    eval(['tmp = CoExp_All.', fc_var_name{i}, '.GABA;'])
    supp = [tmp.afferent_r_Layer3; tmp.afferent_r_Layer2];
    deep = [tmp.afferent_r_Layer5; tmp.afferent_r_Layer6];
    [~,P,~,~] = ttest2(supp, deep);
    p_diff(i,1) = P;
    N(i,1) = length(supp);
    layer_mean(i,1) = nanmean(supp);
    layer_mean(i,2) = nanmean(deep);
    layer_error(i,1) = nanstd(supp)/sqrt(length(supp));
    layer_error(i,2) = nanstd(deep)/sqrt(length(deep));
    
    G_names_supp = [G_names_supp; repmat({gene_sets_name{i}}, length(supp), 1)];
    G_names_deep = [G_names_deep; repmat({gene_sets_name{i}}, length(deep), 1)];
    SUPP = [SUPP; supp];
    DEEP = [DEEP; deep];
end



bc = bar(diag(layer_mean,0), 'stacked', 'FaceAlpha', 0.5, 'EdgeColor','none');
hold on
h1=bar(1, layer_mean(1), 'FaceAlpha', 0, 'EdgeColor','#0072BD', 'LineWidth', 1);
hold on
h2=bar(2, layer_mean(2), 'FaceAlpha', 0, 'EdgeColor','#D95319', 'LineWidth', 1);
hold on
errorbar(1:2, layer_mean, [], layer_error, ...
    'k', 'linestyle', 'none', 'LineWidth', 1.5);
hold off

legend({'Layer 2/3', 'Layer 5/6'}, ...
    'FontSize', 11, 'location', 'northeast', ...
    'box', 'off', 'NumColumns',1);
xlabel('GABA-related gene sets', 'Color', [0, 0, 0], 'FontSize', 11, 'FontWeight', 'normal');
ylabel({'Correlation of co-expression with' ...
    'differences between HV and' ...
    'SCZ in WM-related fc-EC'}, 'Color', [0, 0, 0], 'FontSize', 11, 'FontWeight', 'normal');
set(gca, 'TickLength', [0,0], ...
    'ylim', [0.1, 0.35], ...
    'xticklabel', '', ...
    'xColor', [0, 0, 0], 'yColor', [0, 0, 0], 'FontWeight', 'normal', 'FontSize', 11);

text(1 + 0.02,0.285,['{\it p} = ' num2str(round(p_diff(i),3))], 'FontSize', 11, 'FontWeight', 'normal', 'Color', 'r');
hold on
x = [1-0.25,2+0.25];
y = [0.27, 0.27];
plot(x, y, 'k', 'LineWidth', 2, 'HandleVisibility','off', 'Color', 'r');
hold on
box on


annotation('textbox',...
    [0.025, 0.90, 0.1, 0.1],...
    'String','(A)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);

annotation('textbox',...
    [0.65, 0.90, 0.1, 0.1],...
    'String','(B)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);

annotation('textbox',...
    [0.025, 0.42, 0.1, 0.1],...
    'String','(C)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);

annotation('textbox',...
    [0.38, 0.42, 0.1, 0.1],...
    'String','(D)', ...
    'LineStyle', 'none', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);





% plot just a colorbar
hold on
% row_cmap = 1000;  %define the color scale
% color_map=ones(row_cmap,3);  %define the color matrix
% bar_range = [-2.3,2.3]; % define the range of value reflected by the colorbar
% ratio = round((0 - min(bar_range)) / (max(bar_range) - min(bar_range)),2);
% color_neg = [0, 85, 255]./255; % blue, fix B change R and G
% color_pos = [255, 85, 0]./255; % orange, fix R change B and G
% color_map(1:row_cmap,:) = [
%     [color_neg(1):(1-color_neg(1))/(row_cmap*ratio - 1):1;
%     color_neg(2):(1-color_neg(2))/(row_cmap*ratio - 1):1;
%     ones(1, int16(row_cmap*ratio))]';
%     flip([ones(1, int16(row_cmap*(1-ratio)));
%     color_pos(2):(1-color_pos(2))/(row_cmap*(1-ratio) - 1):1;
%     color_pos(3):(1-color_pos(3))/(row_cmap*(1-ratio) - 1):1;]')
%     ];

ph(11) = subplot('position', [0.01, 0.98, 0.01, 0.01]);
imagesc(floor(rand(120, 120) * 121),'AlphaData', 0);
% colormap(color_map);
colormap('jet');
hCB = colorbar('southoutside', ...
    'Ticks', [0,20,40,60,80,100,120],...
    'TickLabels', cellfun(@num2str, num2cell([0,20,40,60,80,100,120]), 'UniformOutput', false),...
    'TickLength', 0.03, ...
    'EdgeColor', 'none', ...
    'FontSize', 11, 'FontWeight', 'normal');
hCB.Label.Rotation = 0;
hCB.Ruler.Color = 'k';
hCB.Position = [0.2 0.65 0.275 0.03];

% annotation(gcf,'textbox',[0.215 0.65 0.05 0.03],...
%     'LineStyle','none',...
%     'String','0',...
%     'FontSize', 11, 'FontWeight', 'normal');
% annotation(gcf,'textbox',[0.44 0.65 0.05 0.03],...
%     'LineStyle','none',...
%     'String','120',...
%     'FontSize', 11, 'FontWeight', 'normal');
annotation(gcf,'textbox',[0.19 0.6 0.3 0.02],...
    'LineStyle','none',...
    'String','WM-related fc-EC to DLPFC in SCZ (Hz)',...
    'FontSize', 11, 'FontWeight', 'normal', 'HorizontalAlignment', 'center');
ph(11).Visible = 'off';


set(gcf,'color','w');
print(gcf,'./Figure4.png','-dpng','-r300');
close gcf
