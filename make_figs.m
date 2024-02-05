%% Overview
%
% This code is for making figures in
% "Personalized pain biomarker based on densely sampled data"

% Dependencies
%   https://github.com/cocoanlab/CocoanCore
%   https://github.com/canlab/CanlabCore
%   https://github.com/spm/spm12

% Tested on
%   MATLAB version 2021a
%   macOS Moneterey
%
% Details in README.md


%% Prerequisites

basedir = '/Volumes/cocoanlab01/projects/MPC/MyPainConnectome/sync/writing/codes/temp';  % add your path of github repository
cd(basedir)
addpath(genpath(pwd))

% make sure every path to dependencies are added.
% addpath(genpath('~/PATH/TO/COCOANLAB/'))
% addpath(genpath('~/PATH/TO/CANLAB/'))
% addpath(genpath('~/PATH/TO/SPM12/'))
addpath(genpath('functions'));


%% Fig.1b

% violin plot of model test performances
% Model -> Test set: 1. One -> One; 2. Pop -> One; 3. One -> Pop; 4. Pop -> Pop

load(fullfile(basedir,'data/fig1b_data.mat'))

% bootstrap test to calculate mean correlation and p value
for i=1:4
    r_dat = r_all{i};
    [out{i}, boot_vals] = boot_mean_wani(r_dat, 10000, 'noverbose');
    boot_p(i) = out{i}.bootP;
    boot_mean(i) = out{i}.bootmean;
end

% violin plots
boxplot_wani_2016_ye(r_all,'violin','dots','nobox','dot_size',50,'data_dotcolor',violin_col) %r
set(gcf,'position',[163   377    374   266])
set(gca,'tickdir','out', 'ticklength', [.02 .02], 'fontsize',22,'box','off','linewidth',1.7,'color','w')
xlim([.50 4.5])


% convert correlation r -> z
% perform two sample t-test to see if each pair of test performance are significantly different
for i=1:4
    r_z{i} = reformat_r_new(r_all{i},'r2z');
end


for i=1:4
    for j=1:4
        mean_diff(i,j) = mean(r_all{i}) - mean(r_all{j});
        [~,p(i,j),~,stats{i,j}] =  ttest2(r_z{i},r_z{j});
    end
end

%% Fig. 1c
mask = which('gray_matter_mask.nii');
imgdir = fullfile(basedir,'data/pattern_sim_bf_conj.nii');
conj_img = fmri_data(imgdir,mask);

col = [251,115,63; 102,189,249]./255;
conj_img = region(conj_img,'unique_mask_values');

interval = [-22 -14 -8 9 15 44 65];

close all;

[out, o2] = brain_activations_display(conj_img,'all2','all2_xyz',[-2 4 -36 50 interval],'region_color',col);

%% Fig.1d

% Model explained variance
load('fig1d_data.mat')

figure;
for i=1:size(ensemble_result,1)
    hold on
    line([0 ensemble_result.r2(i)],[i i],'color',ensemble_col(i,:),'linewidth',2.5)
    scatter(ensemble_result.r2(i),i,120,'filled','markerfacecolor',ensemble_col(i,:))
    text(ensemble_result.r2(i)+.03,i,sprintf('%0.3f',ensemble_result.r2(i)),'fontsize',16)
end
hold off
set(gca,'tickdir','out', 'ticklength', [.02 .02], 'fontsize',18,'box','off','linewidth',1.2,'color','w')
yticklabels(ensemble_result.name); ylim([0.5 7.5])
set(gcf,'position', [442   662   395   241])


% pie chart

r2_unique = ensemble_result.r2(7) - ensemble_result.r2(4:6); % I unique, E unique, P unique 
r2_overlap(1) = ensemble_result.r2(5) + ensemble_result.r2(4) - ensemble_result.r2(1) -ensemble_result.r2(7); % I+E
r2_overlap(2) = ensemble_result.r2(5) + ensemble_result.r2(6) - ensemble_result.r2(3) -ensemble_result.r2(7); % P+E
r2_overlap(3) = ensemble_result.r2(6) + ensemble_result.r2(4) - ensemble_result.r2(2) -ensemble_result.r2(7); % I+P
r2_all = ensemble_result.r2(7) - sum(ensemble_result.r2(4:6)) + sum(ensemble_result.r2(1:3)); %I+P+E

% order = I unique, I+P, I+E, I+P+E, P+E, E unique, P unique, un explained 
r2_comb = [r2_unique(1),r2_overlap(3), r2_overlap(1), r2_all, r2_overlap(2), r2_unique(2), r2_unique(3)];
r2_comb(end+1) = 1 - ensemble_result.r2(7); %unexplained variance

if any(r2_comb < 0) % if value is negative, treat it as zero
    zero_idx = find(r2_comb < 0);
    r2_comb(zero_idx) = 0;
end


figure;
ax = gca();
h = pie(ax,r2_comb,{'','','','','','','',''})

ax.Colormap = pie_col;

ax.Colormap(8,:) = repmat(120,1,3)./255;
for i=1:8
ax.Children(2*i).EdgeAlpha = 0;
end
