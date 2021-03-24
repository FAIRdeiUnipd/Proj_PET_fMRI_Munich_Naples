% Analysis code to support the findings of 
% "Variability of regional glucose metabolism and the topology of
% functional networks in the human brain"
% by Palombit et al., Neuroimage, 2021
%
% Code written by Palombit A., Silvestri E., Volpi T.

clear
close all
clc

% Set path of folder with Code and Data
folder_path = fullfile(pwd,'..');

Figures_path    = fullfile(folder_path,'Results');   
Utilities_path  = fullfile(pwd,'Utilities'); 
Code_path       = fullfile(folder_path,'Code'); 

addpath(genpath(fullfile(Utilities_path)))
addpath(genpath(fullfile(Code_path)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dataset A
load(fullfile(folder_path, 'Data', 'dataset_A.mat'))
% Dataset B
load(fullfile(folder_path, 'Data', 'dataset_B.mat'))
% group
load(fullfile(folder_path, 'Data', 'dataset_merged.mat'))

% GL atlas networks/parcels info
load(fullfile(folder_path, 'Data', 'networks_info_GL_subcort'))
ROI_TO_REMOVE = [303 328 334 335]; % as seen to provide NaN values over FC matrices 
VALID_ROIs = [1:size(net_assignment_Gordon_table,1)]'; 
VALID_ROIs(ROI_TO_REMOVE) = [];
net_separation = net_assignment_Gordon_table.net_separation(VALID_ROIs);
labels_ord = net_assignment_Gordon_table.labels_ord;
node_hemi = [labels_ord>161; 1;0;1; 0;1; 0;1; 0;1; 0;1; 0;1; 0;1; 0;1; 0];
net_size = [39     8    38     8    24    40    23     4     5    32    24    41    45    16];
net_extent = zeros(14, 1);
for net = 1 : 14
    DD = dsqrform(Euclidean_distance_group(find(net_separation==net),find(net_separation==net)));
    net_extent(net) = median(DD(:));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dataset comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 18F-FDG PET data agreement

% Suppl. Fig. 1-A
figure,
subplot(141), 
imagesc(SUVR_group_dataset_A, [1 1.8]),
colormap hot
colorbar, title('SUVR Dataset A')
rr=find(gradient(net_separation));
rr=[1 ; rr(1:2:end); 333]; 
central_pos=[];
for ii=2:length(rr)
    central_pos=[central_pos (rr(ii)+rr(ii-1))/2]; 
end
central_pos=round(central_pos)'; 
atick=central_pos; set(gca,'XTick',atick); 
set(gca,'XTickLabel', net_names_Gordon,'TickLength',[0 0]);
set(gca,'YTick',atick); set(gca,'YTickLabel', net_names_Gordon,'TickLength',[0 0]);

subplot(142)
imagesc(SUVR_group_dataset_B, [1 1.8])
colormap hot
colorbar, title('SUVR Dataset B')
atick=central_pos; set(gca,'XTick',atick); set(gca,'XTickLabel', net_names_Gordon,'TickLength',[0 0]);
set(gca,'YTick',atick); set(gca,'YTickLabel', net_names_Gordon,'TickLength',[0 0]);

subplot(143), imagesc(SUVR_group, [1 1.8]), colormap hot
colorbar, title('SUVR MERGED')
atick=central_pos; set(gca,'XTick',atick); set(gca,'XTickLabel', net_names_Gordon,'TickLength',[0 0]);
set(gca,'YTick',atick); set(gca,'YTickLabel', net_names_Gordon,'TickLength',[0 0]);

set(gcf, 'Position', [100 100 1800 1000])
print('-dpng', '-r300', fullfile(Figures_path,'SF_1A_SUVR_dataset_A_B_merged.png'))
close

% Suppl. Fig. 1-B : scatter plot - SUVR comparison 
figure,
subplot(2,2,[1 3]), C = hsv(14); hold on
for ii = 1:14
    plot(SUVR_group_dataset_B(net_separation==ii)', SUVR_group_dataset_A(net_separation==ii),'o','color',C(ii,:),'linewidth',3), 
end
legend(net_names_Gordon, 'Location', 'southeast'), xlabel('Dataset B'), ylabel('Dataset A'), 
title('Node SUVR in the two datasets'), grid, set(gca,'DataAspectRatio', [1,1,1]);
xlim([0.85 1.8])
ylim([0.85 1.8])
subplot 222, hist(SUVR_group_dataset_A,20), xlim([0.85 1.8]), title('Dataset A')
subplot 224, hist(SUVR_group_dataset_B,20), xlim([0.85 1.8]), title('Dataset B')
set(gcf, 'Position', [100 100 1800 1000])
print('-dpng', '-r300', fullfile(Figures_path,'SF_1B_SUVR_dataset_A_B_merged_distributions.png'))
close

[~,~,~,~,stats] = regress(SUVR_group_dataset_B, [ones(length(SUVR_group_dataset_A),1), SUVR_group_dataset_A]);
R2 = stats(1);
disp(['SUVR dataset A vs. B - linear least squares  R^2 = ' num2str(R2)])


%% fMRI data agreement

% FC comparison
edges_dataset_A = dsqrform(FC_group_dataset_A)';
edges_dataset_B = dsqrform(FC_group_dataset_B)';
[~,~,~,~,stats] = regress(edges_dataset_A, [ones(length(edges_dataset_A),1), edges_dataset_B]);
R2 = stats(1);
disp(['FC dataset A vs. B - linear least squares    R^2 = ' num2str(R2)])

% check scaling effect
e = ([ones(length(edges_dataset_A),1), edges_dataset_A])\edges_dataset_B; 
f = ([ones(length(edges_dataset_A),1), e(2)*edges_dataset_A+e(1)])\edges_dataset_B;

% Suppl. Fig 2
figure,
subplot(221)
imagesc(FC_group_dataset_A, [-.2 .8])
colormap jet
colorbar, title('FC - Dataset A')
atick=central_pos; set(gca,'XTick',atick,'XTickLabelRotation',90); set(gca,'XTickLabel', net_names_Gordon,'TickLength',[0 0]);
set(gca,'YTick',atick); set(gca,'YTickLabel', net_names_Gordon,'TickLength',[0 0]); set(gca,'DataAspectRatio', [1,1,1]);

subplot(222)
imagesc(FC_group_dataset_B, [-.2 .8])
colormap jet
colorbar, title('FC - Dataset B')
atick=central_pos; set(gca,'XTick',atick,'XTickLabelRotation',90); set(gca,'XTickLabel', net_names_Gordon,'TickLength',[0 0]);
set(gca,'YTick',atick); set(gca,'YTickLabel', net_names_Gordon,'TickLength',[0 0]); set(gca,'DataAspectRatio', [1,1,1]);

subplot(223)
imagesc(FC_group, [-.2 .8])
colormap jet
colorbar, title('FC - Merged')
atick=central_pos; set(gca,'XTick',atick,'XTickLabelRotation',90); set(gca,'XTickLabel', net_names_Gordon,'TickLength',[0 0]);
set(gca,'YTick',atick); set(gca,'YTickLabel', net_names_Gordon,'TickLength',[0 0]); set(gca,'DataAspectRatio', [1,1,1]);

subplot(224)
scatter(edges_dataset_A,edges_dataset_B,'o'), 
xlabel('Dataset A'), ylabel('Dataset B')
xlim([-0.4 1]), ylim([-0.4 1])
grid on
lsline, axis square
set(gcf, 'Position', [100 100 1800 1000])
print('-dpng', '-r300', fullfile(Figures_path,'SF_2_FC_dataset_A_B_merged.png'))
close

% FC graph metrics comparison (computed using BCT external toolbox)

disp(' ')
disp('=============      FC graph metrics comparison  ===================')

% Dataset A
CC      = FC_group_dataset_A;
CC(CC<=0) = 0; 
dens_thr = 1:1:30; 
y       = prctile(dsqrform(CC), 100-dens_thr); % regular THRESHOLDING based on SPARSITY
CCb     = CC > y(20);  % 80th prctile
DEG_A   = degrees_und(CCb); % degree of binarized FC matrix
STR_A   = strengths_und(CC); % strength
EC_A    = eigenvector_centrality_und(CC); % eigenvector centrality
[EBC, BC_A] = edge_betweenness_bin(CCb); % BC 
LE_A    = efficiency_bin(CCb, 1); % weighted local efficiency

% Dataset B
CC      = FC_group_dataset_B; 
CC(CC<=0) = 0; 
y       = prctile(dsqrform(CC), 100-dens_thr); % regular THRESHOLDING based on SPARSITY
CCb     = CC > y(20);  % 80th percentile
DEG_B   = degrees_und(CCb); % degree of binarized FC matrix
STR_B   = strengths_und(CC); % strength
EC_B    = eigenvector_centrality_und(CC); % eigenvector centrality
[~, BC_B] = edge_betweenness_bin(CCb); % BC 
LE_B    = efficiency_bin(CCb, 1); % weighted local efficiency

% assess consistency with Pearson's, Spearman's correlation, linear regression
[r, p] = corr(DEG_B', DEG_A','type','Spearman');
disp(['Node Degree:             R=' num2str(r,3), ', p=' num2str(p,3)])
[~,~,~,~,stats] = regress(DEG_B', [ones(length(DEG_B),1), DEG_A']); 
R2 = stats(1);
disp(['Node Degree:             linear regression R2=' num2str(R2,3)])

[r, p] = corr(STR_B', STR_A','type','Spearman');
disp(['Node Strength:           R=' num2str(r,3), ', p=' num2str(p,3)])
[~,~,~,~,stats] = regress(STR_B', [ones(length(STR_B),1), STR_A']); 
R2 = stats(1);
disp(['Node Strength:           linear regression R2=' num2str(R2,3)])

[r, p] = corr(EC_B, EC_A,'type','Spearman');
disp(['Eigenvector Centrality:  R=' num2str(r,3), ', p=' num2str(p,3)])
[~,~,~,~,stats] = regress(EC_B, [ones(length(EC_B),1), EC_A]); 
R2 = stats(1);
disp(['Eigenvector Centrality:  linear regression R2=' num2str(R2,3)])

[r, p] = corr(BC_B, BC_A,'type','Spearman');
disp(['Betweeness Centrality:   R=' num2str(r,3), ', p=' num2str(p,3)])
[~,~,~,~,stats] = regress(BC_B, [ones(length(BC_B),1), BC_A]); 
R2 = stats(1);
disp(['Betweeness Centrality:   linear regression R2=' num2str(R2,3)])

[r, p] = corr(LE_B, LE_A,'type','Spearman');
disp(['Local Efficiency (BIN):  R=' num2str(r,3), ', p=' num2str(p,3)])
[~,~,~,~,stats] = regress(LE_B, [ones(length(LE_B),1), LE_A]); 
R2 = stats(1);
disp(['Local Efficiency (BIN):  linear regression R2=' num2str(R2,3)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Metabolism of Resting State Networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MAP = colormap(hsv);
network_rgb_map = MAP([1 (1:12)*5 end],:);

% Fig. 2-A - boxplot of SUVR distribution in the GL RSNs
figure
boxplot(SUVR_group, net_separation,'labels',net_names_Gordon), grid 
set(gcf, 'Position', [100 100 1800 1000])
xlabel('RSNs')
ylabel('SUVR [adim.]')
print('-dpng', '-r300', fullfile(Figures_path,'Fig_2A_boxplot_SUVR_RSNs.png'))
close all

% Kruskal-Wallis test (Nonparametric ANOVA) for differences in SUVR medians between RSNs
KW_RSNs = kruskalwallis(SUVR_group, net_separation);
close all

% after removing None and Subcortical RSNs
KW_RSNs_woNoneSub = kruskalwallis(SUVR_group(net_separation<13), net_separation(net_separation<13));
close all

% Post-hoc Wilcoxon rank-sum tests with FDR correction
post_hoc = zeros(14, 14);
for G1 = 1 : 14
    for G2 = 1 : 14
        post_hoc(G1, G2) = ranksum(SUVR_group(find(net_separation==G1)), SUVR_group(find(net_separation==G2)),'tail','right');
    end
end
disp(' ')
disp('==============  Statistical significant differences in median SUVR between RSNs =====================')
[h, crit_p] = fdr_bh(post_hoc(:),.05,'pdep','yes');

%% examine possible confounds: network size, distance between ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=====================================================================')

mdn_net = [];
var_net = [];
for ii = 1 : 14    
    mdn_net = [mdn_net; median(SUVR_group(net_separation==ii))]; % network SUVR mean among the group  
    var_net = [var_net; mad(SUVR_group(net_separation==ii))];    % network SUVR mad among the group 
end

% MAD SUVR
[r, p] = corr(var_net, net_size','type','Spearman');
disp(['Correlation between SUVR (mad) and net size:                         R=' num2str(r,2) ', p=' num2str(p,3)])
[r, p] = corr(var_net(1:12), net_size(1:12)','type','Spearman');
disp(['Correlation between SUVR (mad) and net size (no SUBCORT):            R=' num2str(r,2) ', p=' num2str(p,3)])

[r, p] = corr(var_net, net_extent,'type','Spearman');
disp(['Correlation between SUVR (mad) and net extent:                       R=' num2str(r,2) ', p=' num2str(p,3)])
[r, p] = corr(var_net(1:12), net_extent(1:12),'type','Spearman');
disp(['Correlation between SUVR (mad) and net extent (no SUBCORT):          R=' num2str(r,2) ', p=' num2str(p,3)])

[r, p] = corr(var_net, mdn_net,'type','Spearman');
disp(['Correlation between SUVR (mad) and SUVR (median):                    R=' num2str(r,2) ', p= ' num2str(p,3)])
[r, p] = corr(var_net([1:12]), mdn_net([1:12]),'type','Spearman');
disp(['Correlation between SUVR (mad) and SUVR (median) ((no SUBCORT)):     R=' num2str(r,2) ', p= ' num2str(p,3)])


% Median SUVR
[r, p] = corr(mdn_net, net_size','type','Spearman');
disp(['Correlation between SUVR (median) and net size:                      R=' num2str(r,2) ', p= ' num2str(p,3)])
[r, p] = corr(mdn_net([1:12]), net_size([1:12])','type','Spearman');
disp(['Correlation between SUVR (median) and net size (no SUBCORT):         R=' num2str(r,2) ', p= ' num2str(p,3)])
[r, p] = corr(mdn_net, net_extent,'type','Spearman');
disp(['Correlation between SUVR (median) and net extent:                    R=' num2str(r,2) ', p= ' num2str(p,3)])
[r, p] = corr(mdn_net([1:12]), net_extent([1:12]),'type','Spearman');
disp(['Correlation between SUVR (median) and net extent (no SUBCORT):       R=' num2str(r,2) ', p= ' num2str(p,3)])


%% Intrinsic RSNs vs Extrinsic RSNs (as in Hacker et al. 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

extr = [10, 7, 1, 2, 3, 4]; % Extrinsic: DAN, VAN, VIS, SMN (=SMN+AUD)
intr = [11, 12];            % Intrinsic: FPC (Fronto Parietal Control), LAN (Language Netw - not available), DMN

SUVR_extr = []; 
for ii = 1 : length(extr) 
    SUVR_extr = [SUVR_extr; SUVR_group(net_separation==extr(ii))]; 
end
SUVR_intr = []; 
for ii = 1 : length(intr) 
    SUVR_intr = [SUVR_intr; SUVR_group(net_separation==intr(ii))]; 
end


disp(' ')
disp('============= SUVR in intrinsic vs extrinsic RSNs =================')

% Wilcoxon rank-sum test
[p, h] = ranksum(SUVR_extr, SUVR_intr);
disp(['Ranksum between SUVR in intrinsic and extrinsic networks: H= ' num2str(h) ', p=' num2str(p,3)]);

% Suppl. Fig. 3-B
boxplot([SUVR_extr; SUVR_intr],[ones(length(SUVR_extr),1); 2*ones(length(SUVR_intr),1)]), ylabel('SUVR [adim.]','FontSize',14)
set(gcf, 'Position', [100 100 1800 1000]), grid on
xticks(1:2), xticklabels({'Extrinsic', 'Intrinsic'})
print('-dpng', '-r300', fullfile(Figures_path,'SF_3B_boxplot_SUVR_intrinsic_vs_extrinsic.png'))
close


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Relationship between graph FC metrics and metabolism (on the merged dataset)
% Calculate nodal graph metrics
% BCT is required, plus Versatility package
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CC          = FC_group;
CC(CC<=0)   = 0; 
dens_thr    = [1:1:30]; 
y           = prctile(dsqrform(CC), 100-dens_thr); % regular THRESHOLDING based on SPARSITY
CCb         = CC > y(20);  % 80th prctile

load(fullfile(folder_path, 'Data', 'net_metrics.mat')) %% Loading results instead of computing graph measures 

% bin_deg         = degrees_und(CCb); % degree of binarized FC matrix
% wei_strength    = strengths_und(CC); % strength
% clust_coeff_wei = clustering_coef_wu(CC); % weighted clustering coeff
% Loc_eff_wei     = efficiency_wei(CC, 2); % weighted local efficiency
% [~, BC]         = edge_betweenness_bin(CCb); % BC 
% [~, BC_w]       = edge_betweenness_wei(CC); % BC 
% Eig_centr_wei   = eigenvector_centrality_und(CC); % eigenvector centrality
% Part_bin        = participation_coef(CCb, net_separation); % participation coefficient
% 
net_metrics = [bin_deg', wei_strength', clust_coeff_wei, Loc_eff_wei, BC, BC_w, Eig_centr_wei, Part_bin];
net_metrics_names = {'degree', 'strength', 'CC', 'LEw', 'BC', 'BCw', 'EC', 'Participation'}';

% GLOBAL CORRELATIONS (= ACROSS ALL ROIS) BETWEEN SUVR AND GRAPH METRICS
YY = zeros(size(net_metrics,2),1); 
PP = YY;
for yy=1:size(YY,1)
    [a, b] = corr(net_metrics(find(net_metrics(:,yy)>0),yy), SUVR_group(find(net_metrics(:,yy)>0)), 'type', 'Spearman');
    YY(yy)=a; 
    PP(yy)=b;
end
correlation_results = table(net_metrics_names, YY, PP); 

% NETWORK-WISE CORRELATIONS for all graph metrics
YY_subnets = zeros(14, size(net_metrics,2)); 
PP_subnets = zeros(14, size(net_metrics,2));
for ii = 1 : 14
    [A, B] = corr(net_metrics(find(net_separation==ii),:), SUVR_group(find(net_separation==ii)));
    YY_subnets(ii,:)=A;
    PP_subnets(ii,:)=B;
end

net_size_node = zeros(length(net_separation),1);
for nn = 1 : length(net_separation)
    net_size_node(nn) = net_size(net_separation(nn));
end
[R, p] = corr(YY_subnets(:,1),net_size');
disp(['Correlation between graph and SUVR correlation and network size: R=' num2str(R,3), ', p=' num2str(p,3)]);

net_extent = zeros(14, 1);
for net = 1 : 14
    DD=dsqrform(Euclidean_distance_group(find(net_separation==net),find(net_separation==net)));
    net_extent(net) = median(DD(:));
end
[R, p] = corr(YY_subnets(:,2), net_extent);
disp(['Correlation between graph and SUVR correlation and network extent: R=' num2str(R,3), ', p=' num2str(p,3)]);

% Figure 3-B, Suppl. Fig. 4-A,B - Scatterplots
chosen_metrics = [2 1 7]; % STR, DEG, EC 
figures_name = {'Fig_3B', 'SF_4A', 'SF_4B'};
for met = 1 : length(chosen_metrics)
    figure, C = jet(14); hold on
    for ii=1:14
        plot(net_metrics(net_separation==ii, chosen_metrics(met))', SUVR_group(net_separation==ii),'o','color',C(ii,:),'linewidth',3), 
    end
    legend(net_names_Gordon), xlabel([net_metrics_names{chosen_metrics(met)}]), ylabel('SUVR'), 
    title(['Node ' net_metrics_names{chosen_metrics(met)} ' vs SUVR, omitting zeros: corr = ' num2str(YY(chosen_metrics(met)), '%.2f')]), hold off
    
    % add regression line and CI 75% to each plot
    x1      = net_metrics(:, chosen_metrics(met)); 
    y       = SUVR_group; 
    X       = [ones(size(x1)) x1]; 
    [b, bint] = regress(y,X); 
    xval    = min(x1):0.01:max(x1);
    yhat    = b(1)+b(2)*xval; 
    ylow    = bint(1,1)+bint(2,1)*xval;
    yupp    = bint(1,2)+bint(2,2)*xval;
    hold on;
    plot(xval,ylow,'r-.');
    plot(xval,yupp,'r-.');
    plot(xval,yhat,'b','linewidth',2);
    set(gcf, 'Position', [100 100 1800 1000])
    print('-dpng', '-r300', fullfile(Figures_path,[figures_name{met} '_SUVR_vs_' net_metrics_names{chosen_metrics(met)} '_scatter.png']))
    close
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% relationship between STR and SUVR for within-RSN vs. between-RSN/
% short-distance vs. long-distance edges

CC           = FC_group; 
CC(CC<=0)    = 0;
wei_strength =(sum(CC))';

dsqr_Euclidean_distance_group = dsqrform(Euclidean_distance_group);

% Long-range edges (> 75th percentile): LR
LR = Euclidean_distance_group > prctile(dsqr_Euclidean_distance_group(dsqr_Euclidean_distance_group>0), 75); 
wei_strength_LR75 = (sum(CC.*LR))';

% Short-range edges (< 25th percentile): SR
SR = Euclidean_distance_group < prctile(dsqr_Euclidean_distance_group(dsqr_Euclidean_distance_group>0), 25); 
wei_strength_SR25 = (sum(CC.*SR))'; 

disp(' ')
disp('=====================================================================')

[R,p]=corr(wei_strength_LR75(find(wei_strength_LR75>0)), SUVR_group(find(wei_strength_LR75>0)), 'type', 'Spearman');
disp(['Correlation between long range distance (>75th) and group SUVR:      R=' num2str(R,3), ', p= ' num2str(p,3)])
[R,p]=corr(wei_strength_SR25(find(wei_strength_SR25>0)), SUVR_group(find(wei_strength_SR25>0)), 'type', 'Spearman');
disp(['Correlation between short range distance (<25th) and group SUVR:     R=' num2str(R,3), ', p= ' num2str(p,3)])


net_size_associated = zeros(size(wei_strength)); 
for node = 1 : length(net_separation)
    net_size_associated(node) = length(find(net_separation==net_separation(node))); 
end
% Within-network edges: IN
intranet_edges = zeros(size(CC)); 
for net = 1 : length(unique(net_separation))
    intranet_edges(find(net_separation==net),find(net_separation==net)) = 1; 
end
% Between-network edges: OUT
internet_edges = ones(size(CC)); 
internet_edges = internet_edges-intranet_edges;
wei_strength_IN = strengths_und(CC.*intranet_edges);
wei_strength_OUT = strengths_und(CC.*internet_edges);

% Intersection of IN/OUT with LR/SR
wei_strength_SR25_IN	  = (sum(CC.*SR.*intranet_edges))'; 
wei_strength_SR25_OUT	= (sum(CC.*SR.*internet_edges))'; 
wei_strength_LR75_IN	  = (sum(CC.*LR.*intranet_edges))'; 
wei_strength_LR75_OUT	= (sum(CC.*LR.*internet_edges))'; 

[r1, p1] = corr(SUVR_group, wei_strength_SR25_IN, 'type', 'Spearman');
disp(['Correlation between group SUVR and weighted strength in short range intrinsic connections:   R=' num2str(r1,3), ', p= ' num2str(p1,3)])
[r2, p2] = corr(SUVR_group, wei_strength_SR25_OUT, 'type', 'Spearman');
disp(['Correlation between group SUVR and weighted strength in short range estrinsic connections:   R=' num2str(r2,3), ', p= ' num2str(p2,3)])
[r3, p3] = corr(SUVR_group, wei_strength_LR75_IN, 'type', 'Spearman');
disp(['Correlation between group SUVR and weighted strength in long range intrinsic connections:    R=' num2str(r3,3), ', p= ' num2str(p3,3)])
[r4, p4] = corr(SUVR_group, wei_strength_LR75_OUT, 'type', 'Spearman');
disp(['Correlation between group SUVR and weighted strength in long range estrinsic connections:    R=' num2str(r4,3), ', p= ' num2str(p4,3)])


% Figure 4-A: SR,LR,IN,OUT edges
figure, 
subplot(221) 
imagesc(FC_group, [-0.2 0.6]), colorbar, colormap jet
title('Group FC')
atick=central_pos; set(gca,'XTick',atick); set(gca,'XTickLabel', net_names_Gordon,'XTickLabelRotation',90,'TickLength',[0 0]);
set(gca,'YTick',atick); set(gca,'YTickLabel', net_names_Gordon,'TickLength',[0 0]); set(gca,'DataAspectRatio', [1,1,1]);

subplot(222)
imagesc(CC.*SR, [-0.2 0.6]), title('Short range edges')
colorbar
colormap jet
atick=central_pos; set(gca,'XTick',atick); set(gca,'XTickLabel', net_names_Gordon,'XTickLabelRotation',90, 'TickLength',[0 0]);
set(gca,'YTick',atick); set(gca,'YTickLabel', net_names_Gordon,'TickLength',[0 0]); set(gca,'DataAspectRatio', [1,1,1]);

subplot(224)
imagesc(CC.*LR, [-0.2 0.6]), title('Long range edges')
colorbar, colormap jet
atick=central_pos; set(gca,'XTick',atick); set(gca,'XTickLabel', net_names_Gordon,'XTickLabelRotation',90,'TickLength',[0 0]);
set(gca,'YTick',atick); set(gca,'YTickLabel', net_names_Gordon,'TickLength',[0 0]); set(gca,'DataAspectRatio', [1,1,1]);

subplot(223)
imagesc(internet_edges), title('Within vs. Between network edges')
colormap jet
atick=central_pos; set(gca,'XTick',atick); set(gca,'XTickLabel', net_names_Gordon,'XTickLabelRotation',90,'TickLength',[0 0]);
set(gca,'YTick',atick); set(gca,'YTickLabel', net_names_Gordon,'TickLength',[0 0]); set(gca,'DataAspectRatio', [1,1,1]);
set(gcf, 'Position', [100 100 1800 1000])
print('-dpng', '-r300', fullfile(Figures_path,'Fig_4A_FC_short_long_range_within_between_RSNs.png'))
close all

% Figure 4-B: scatter plots SUVR vs STR SR (IN/OUT) and LR (IN/OUT)
figure,
subplot(121)
plot(wei_strength_SR25_IN(find(wei_strength_SR25_IN>0)), ...
    SUVR_group(find(wei_strength_SR25_IN>0)),'ogr', ...
    wei_strength_SR25_OUT(find(wei_strength_SR25_OUT>0)), ...
    SUVR_group(find(wei_strength_SR25_OUT >0)),'or')
legend(['Within-RSN (\rho = ' num2str(round(r1,2)) ')'],...
    ['Between-RSN (\rho = ' num2str(round(r2,2)) ')'], 'Location', 'southeast'), 
xlabel('Strength'), ylabel('SUVR'), lsline
title('Short range edges only')

subplot(122)
plot(wei_strength_LR75_IN(find(wei_strength_LR75_IN>0)), ...
    SUVR_group(find(wei_strength_LR75_IN>0)),'ogr', ...
    wei_strength_LR75_OUT(find(wei_strength_LR75_OUT>0)), ...
    SUVR_group(find(wei_strength_LR75_OUT >0)),'or')
legend(['Within-RSN (\rho = ' num2str(round(r3,2)) ')'],...
    ['Between-RSN (\rho = ' num2str(round(r4,2)) ')'], 'Location', 'southeast'), 
xlabel('Strength'), ylabel('SUVR'), lsline
title('Long range edges only')
set(gcf, 'Position', [100 100 1800 1000])
print('-dpng', '-r300', fullfile(Figures_path,'Fig_4B_scatter_SUVR_vs_short_long_range_within_between_RSNs.png'))
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% How does the SUVR- graph metrics relationship change after regressing out distance? --> Spoiler: No change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CC          = FC_group; 
CC(CC<=0)   = 0; 
dens_thr    = 1:1:30; 
y           = prctile(dsqrform(CC), 100-dens_thr); % regular THRESHOLDING based on SPARSITY
CCb         = CC > y(20);  

X = dsqrform(Euclidean_distance_group)'; 
Y = dsqrform(CC)';
b = robustfit(X((X>0)&(Y>0)), Y((X>0)&(Y>0)));
Y_regr = Y-(b(1)+b(2)*X); 
CC_dist_regr = sqrform(Y_regr');

CC          = CC_dist_regr; % updating CC 
CC(CC<=0)   = 0; 
dens_thr    = 1:1:30; 
y           = prctile(dsqrform(CC), 100-dens_thr); % regular THRESHOLDING based on SPARSITY
CCb         = CC > y(20);  % 80th prctile

load(fullfile(folder_path,'Data','net_metrics_distREG.mat')) % loading results instead of computing graph measures 

% bin_deg_distREG                 = degrees_und(CCb); % degree of binarized FC matrix
% wei_strength_distREG            = strengths_und(CC); % strength
% clust_coeff_wei_distREG         = clustering_coef_wu(CC); % weighted clustering coeff
% Loc_eff_wei_distREG             = efficiency_wei(CC, 2); % weighted local efficiency
% [EBC_distREG, BC_distREG]       = edge_betweenness_bin(CCb); % BC 
% [EBC_w_distREG, BC_w_distREG]   = edge_betweenness_wei(CC); % BC 
% Eig_centr_wei_distREG           = eigenvector_centrality_und(CC); % eigenvector centrality
% Part_bin_distREG                = participation_coef(CCb, net_separation); % participation coefficient

net_metrics_distREG = [bin_deg_distREG', wei_strength_distREG', clust_coeff_wei_distREG, Loc_eff_wei_distREG, BC_distREG, BC_w_distREG, Eig_centr_wei_distREG, Part_bin_distREG];

%%%%%% GLOBAL CORRELATIONS (= ACROSS ALL ROIS) BETWEEN SUVR AND GRAPH METRICS
YY = zeros(size(net_metrics_distREG,2),1); 
PP = YY;
for yy=1:size(YY,1)
    [a, b] = corr(net_metrics_distREG(find(net_metrics_distREG(:,yy)>0),yy), SUVR_group(find(net_metrics_distREG(:,yy)>0)), 'type', 'Spearman');
    YY(yy)=a; 
    PP(yy)=b;
end
global_correlation_results_distance_regressed = table(net_metrics_names, YY, PP); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Selection of the most representative graph measure: stepwise regression analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('  ')
disp('============ Selection of hubness measure: Stepwise Analysis   ==================================')

chosen_graph_metrics = [1:4 6:8];
net_metrics_norm = (net_metrics - repmat(median(net_metrics, 1, 'omitnan'),[size(net_metrics,1),1]))./(repmat(mad(net_metrics,1,1),[size(net_metrics,1),1]));
[b,se,pval,inmodel] = stepwisefit(net_metrics_norm(:,chosen_graph_metrics), SUVR_group);

list_of_chosen_feat = chosen_graph_metrics(inmodel);
list_of_chosen_feat_txt = [] ;
for kk = 1: length(list_of_chosen_feat)
    list_of_chosen_feat_txt = [list_of_chosen_feat_txt ' -' net_metrics_names{list_of_chosen_feat(kk)}]; 
end
disp(['Selected Features: ' list_of_chosen_feat_txt])
disp(['STR - p= ' num2str(pval(2),3)])
disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Metabolic consumption of central network nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEG = bin_deg; % degree of binarized FC matrix
STR = wei_strength; % strength
EC = Eig_centr_wei; % eigenvector centrality
PAR = Part_bin; % participation coefficient

selection_metric = DEG;
nperm = 5000;
net_metrics     = [DEG', STR, EC, BC, PAR];
unique_bin_deg  = unique(selection_metric(find(selection_metric>0))); 

p_perm          = zeros(95,2); p_perm_group = zeros(95,2);
hub_nodes       = zeros(95,1); 
x_axis          = zeros(95,1); 
suv_relation    = zeros(95,1); 
suv_relation_p  = zeros(95,1);

for deg = 1 : 95
    
    x_axis(deg) = prctile(unique_bin_deg, deg);
    hubs        = find(selection_metric > prctile(unique_bin_deg, deg))';
    
    if length(hubs)>2 && not(isempty(hubs))
        hub_nodes(deg) = length(hubs)/length(selection_metric);
        
        % select random subsample of nodes with all degrees to compute p-values
        PERMS = zeros(nperm, 1);
        for pp = 1 : nperm
            indxs       = sort(datasample(1:size(selection_metric,2), length(hubs), 'Replace', false));
            PERMS(pp)   = corr(net_metrics(indxs, 2), SUVR_group(indxs), 'type', 'Spearman');
        end
        
        [p_perm(deg,1), p_perm(deg,2)] = mult_comp_perm_corr(net_metrics(hubs, 2), SUVR_group(hubs), nperm, 1, 0.05, 'rank', 0);
        [YY, PP]            = corr(net_metrics(hubs, 2), SUVR_group(hubs), 'type', 'Spearman');
        suv_relation(deg)   = YY; 
        suv_relation_p(deg) = PP;
        
        p_perm_group(deg,1) = length(find(PERMS<YY))/nperm;
        p_perm_group(deg,2) = length(find(PERMS>YY))/nperm;
        
    end
end
x_axis_prct = 1:95;

disp(' ')
disp('======== Statically significant correlation between SUVR and strength (FDR mul comp corr) ==============')
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(suv_relation_p,.05,'dep','yes');

% Fig 5 - Corr SUVR-Strength (Hubs)
figure,
subplot(211)
plot(x_axis_prct, suv_relation,'-b', x_axis_prct(find(p_perm(:,1)<crit_p)), suv_relation(p_perm(:,1)<crit_p),'*r'), 
title('Association STR-SUVR with increasing DEG (FDR)'), 
xlabel('Distribution percentiles (HUB selection threshold over DEG)'), 
ylabel('\rho SUVR-STR'), grid
hold on
plot(x_axis_prct(find(p_perm_group(:,1)<0.05)), suv_relation(find(p_perm_group(:,1)<0.05)),'-m', x_axis_prct(find(p_perm_group(:,2)<0.05)), suv_relation(find(p_perm_group(:,2)<0.05)),'-gr','LineWidth',2)
hold off
legend('\rho SUVR-STR', '\rho significance', '\rho < \rho random set', '\rho > \rho random set', 'Location', 'northwest')

subplot(212) 
plot(x_axis_prct, hub_nodes,'-b')
title('Number of nodes selected with increasing DEG')
xlabel('Percentiles')
ylabel('Number of nodes')
grid on
set(gcf, 'Position', [100 100 1800 1000])
print('-dpng', '-r300', fullfile(Figures_path, 'Fig_5_STR_SUVR_association_at_increasing_DEG_percentiles.png'))
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STR-SUVR relationship in provincial and connector hubs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')

% Build HUB selection plot by increasing degree (0-95th percentile)
CC          = FC_group; 
CC(CC<=0)   = 0; 
dens_thr    = 1:1:30; 
y           = prctile(dsqrform(CC), 100-dens_thr); % regular THRESHOLDING based on SPARSITY
CCb         = CC > y(20);  % select sparsity level consistent with Buckner 2009 (r>0.25 approx.)

hub_by_degree   = find(bin_deg > prctile(bin_deg, 90));

% Fig. 6-A 
tmp         = zeros(347,347);
tmp(hub_by_degree,:) = 1;
tmp(:,hub_by_degree) = 1;
figure1 = figure('Colormap',[0 0 0.3;0 0 1;0 1 0]);
imagesc(tmp.*CCb +CCb)
colorbar('Ticks',[0.35 1 1.7],'TickLabels',{'Absent Conn','nonHUB','HUB'});
set(gca,'XTick',atick); 
set(gca,'XTickLabel', net_names_Gordon, 'XTickLabelRotation',90,'TickLength',[0 0]);
set(gca,'YTick',atick); set(gca,'YTickLabel', net_names_Gordon,'TickLength',[0 0])
axis square
set(gcf, 'Position', [100 100 1800 1000])
print('-dpng', '-r300', fullfile(Figures_path,'Fig_6A_connector_provincial_HUBs.png'))
close all

connector_hubs  = intersect(hub_by_degree, find(Part_bin > prctile(Part_bin(hub_by_degree), 70)));
provincial_hubs = intersect(hub_by_degree, find(Part_bin < prctile(Part_bin(hub_by_degree), 30)));

disp('======= Connectors vs Provincial Hubs SUVR comparison ====== ')
hubs = hub_by_degree';
disp(['Number of high-DEG nodes: ' num2str(length(hubs))]) % N of high-DEG nodes
disp(['Number of hubs in R hemiphere: ' num2str(sum(node_hemi(hubs)))]) % n Hubs by hemisphere
disp(['Number of hubs in L hemiphere: ' num2str(sum(not(node_hemi(hubs))))]) % n Hubs by hemisphere
list_of_hubs = unique(net_separation(hubs));
list_of_hubs_txt = [];
for jj = 1: length(list_of_hubs)
    list_of_hubs_txt = [list_of_hubs_txt ' ' net_names_Gordon{list_of_hubs(jj)}];
end
disp(['Networks enriched in Hubs: ' list_of_hubs_txt])


allN = 1:length(net_separation); 
allN(hubs) = [];
P = ranksum(SUVR_group(hubs), SUVR_group(allN));
disp(['Wilcoxon rank-sum test - SUVR in HUBs vs. other nodes: p= ' num2str(P)])
net_metrics = [bin_deg', wei_strength, clust_coeff_wei, Loc_eff_wei, BC, BC_w, Eig_centr_wei, Part_bin];

[YY_hubs, PP_hubs] = corr(net_metrics(hubs,:), SUVR_group(hubs));
[YY_connector_hubs, PP_connector_hubs] = corr(net_metrics(connector_hubs,:), SUVR_group(connector_hubs));
[YY_provincial_hubs, PP_provincial_hubs] = corr(net_metrics(provincial_hubs,:), SUVR_group(provincial_hubs));

P = ranksum(SUVR_group(provincial_hubs),SUVR_group(connector_hubs));
relD = 100*(median(SUVR_group(provincial_hubs))-median(SUVR_group(connector_hubs)))/((median(SUVR_group(provincial_hubs))+median(SUVR_group(connector_hubs)))/2);
disp(['Wilcoxon rank-sum test - SUVR in connector vs provincial HUBs: p= ' num2str(P)])
disp(['Wilcoxon rank-sum test - SUVR in connector vs provincial HUBs: relative difference = ' num2str(relD,3) '%'])

 % Figure 7
figure
subplot(231)
boxplot([SUVR_group(provincial_hubs),SUVR_group(connector_hubs)])
xticks(1:2), xticklabels({'Provincial HUBs', 'Connector HUBs'})
ylabel('SUVR')

subplot(232)
plot(net_metrics(hubs,1), SUVR_group(hubs),'.b', net_metrics(connector_hubs,1), ...
    SUVR_group(connector_hubs),'or', net_metrics(provincial_hubs,1), SUVR_group(provincial_hubs),'*gr')
title(['SUVR-DEG relation over HUB nodes, corr = ' num2str(YY_hubs(1))]), 
legend('All HUBs','Connector HUBs','Provincial HUBs', 'Location', 'northwest'), 
xlabel('DEG'), ylabel('SUVR'), grid

subplot(233)
plot(net_metrics(hubs,5), SUVR_group(hubs),'.b', net_metrics(connector_hubs,5), ...
    SUVR_group(connector_hubs),'or', net_metrics(provincial_hubs,6), SUVR_group(provincial_hubs),'*gr')
title(['SUVR-BC relation over HUB nodes, corr = ' num2str(YY_hubs(5))]), 
legend('All HUBs','Connector HUBs','Provincial HUBs', 'Location', 'northwest'), 
xlabel('BC'), ylabel('SUVR'), grid
xlim([0 4000])
subplot(235) 
plot(net_metrics(hubs,2), SUVR_group(hubs),'.b', net_metrics(connector_hubs,2), ...
    SUVR_group(connector_hubs),'or', net_metrics(provincial_hubs,2), SUVR_group(provincial_hubs),'*gr')
title(['SUVR-STRENGTH relation over HUB nodes, corr = ' num2str(YY_hubs(2))]), 
legend('All HUBs','Connector HUBs','Provincial HUBs', 'Location', 'northwest'), 
xlabel('STR'), ylabel('SUVR'), grid
subplot(236)
plot(net_metrics(hubs,8), SUVR_group(hubs),'.b', net_metrics(connector_hubs,8), ...
    SUVR_group(connector_hubs),'or', net_metrics(provincial_hubs,8), SUVR_group(provincial_hubs),'*gr')
title(['SUVR-PAR relation over HUB nodes, corr = ' num2str(YY_hubs(8))]), 
legend('All HUBs','Connector HUBs','Provincial HUBs', 'Location', 'northwest'), 
xlabel('PAR'), ylabel('SUVR'), grid
set(gcf, 'Position', [100 100 1800 1000])
print('-dpng', '-r300', fullfile(Figures_path,'Fig_7_SUVR_graph_metrics_association_connector_provincial_HUBs.png'))
close all

save(fullfile(Figures_path,'Results.mat'),'correlation_results','global_correlation_results_distance_regressed')
