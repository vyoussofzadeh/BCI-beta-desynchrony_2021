%% BCI dataset, Ulster University

% Pipeline_correlation analysis
% Writtern by Vahab Youssof Zadeh
% Update: 08/10/2022

clear, close all

%% Initial settings
%- Brainstorm dir
bsdir = '/data/MEG/Research/BCI_Ulster/BCI5';
bs_path = '/opt/matlab_toolboxes/brainstorm3';
addpath ./functions
addpath ./data
ft_path = '/opt/matlab_toolboxes/fieldtrip';
indir = '/data/MEG/Research/BCI_Ulster/SPM_PRT3/nii';
atlasdir = '/data/MEG/Vahab/Github/MCW_MEGlab/Projects/BCI/Data/atlas';
maskdir = '/data/MEG/Vahab/Github/MCW_MEGlab/Projects/BCI/Data/mask';

%% Adding FT
addpath(ft_path);
ft_defaults

%% choose atlas
load(fullfile(atlasdir,'DK_bs.mat'))
dkatlas = ft_read_atlas(fullfile(atlasdir,'DK.nii'));
dkatlas.parcellationlabel = atlas_DK1.dk_labels;

cd(indir);

%% Reading voxel-level beta power source maps (saved as nii using Brainstorm)
dd = rdir(fullfile (indir,'/**/*.nii'));

clear sFiles_name dconn sub val D_par_DK_all
for ii=1:length(dd)
    disp(num2str(ii))
    
    [a, b] = fileparts(dd(ii).name);
    cd(a)
    tmp = ft_read_mri(dd(ii).name); val(ii,:) = reshape(tmp.anatomy,[1, size(tmp.anatomy,1)*size(tmp.anatomy,2)*size(tmp.anatomy,3)]);
    val(ii,:) = val(ii,:)./max(val(ii,:));
    
    comm_data = b;
    cd(dd(ii).folder)
    tkz = tokenize(dd(ii).name,'/');
    tkz1 = tokenize(b,'_');
    if length(tkz1) ==1
        sub{ii} = '1';
    else
        sub{ii} = tkz1{2};
    end
    dconn{ii} = tkz{end-1};
    sFiles_name{ii} = [dconn{ii},'_', sub{ii}];
    
    switch dconn{ii}
        case {'4', '5'}
            class_run(ii,1) = 4;
        case {'8', '9'}
            class_run(ii,1) = 8;
        case {'16', '17'}
            class_run(ii,1) = 16;
        case {'32', '33'}
            class_run(ii,1) = 32;
    end
    
    if rem(str2num(sub{ii}), 2) == 0, run{ii} = 'S02'; else, run{ii} = 'S01'; end
    
    switch run{ii}
        case 'S01'
            class_run(ii,2) = 1;
        case 'S02'
            class_run(ii,2) = 2;
    end
end

%% THIS IS JUST TO REASSUTING ABOUT INDECIES rr
idx = [];
idx.h = find(class_run(:,1) == 4)';
idx.f = find(class_run(:,1) == 8)';
idx.w = find(class_run(:,1) == 16)';
idx.s = find(class_run(:,1) == 32)';

idx_s1 = [];
idx_s1.h = idx.h(class_run(idx.h,2) == 1);
idx_s1.f = idx.f(class_run(idx.f,2) == 1);
idx_s1.w = idx.w(class_run(idx.w,2) == 1);
idx_s1.s = idx.s(class_run(idx.f,2) == 1);

idx_s2 = [];
idx_s2.h = idx.h(class_run(idx.h,2) == 2);
idx_s2.f = idx.f(class_run(idx.f,2) == 2);
idx_s2.w = idx.w(class_run(idx.w,2) == 2);
idx_s2.s = idx.s(class_run(idx.f,2) == 2);

%%
clear  rr rr1 rr2
rr.idx_W = 1:28; rr.idx_S = 29:29+28-1; rr.idx_H = 57:57+28-1; rr.idx_F = 85:85+28-1;
rr1.idx_W = 1:2:28; rr1.idx_S = 29:2:29+28-1; rr1.idx_H = 57:2:57+28-1; rr1.idx_F = 85:2:85+28-1;
rr2.idx_W = 2:2:28; rr2.idx_S = 30:2:29+28-1; rr2.idx_H = 58:2:57+28-1; rr2.idx_F = 86:2:85+28-1;

%% - all values
idx = mean(val,1) == 0; val1 = val; val1(:,idx) = [];

%% Masks
Mask = [];

tmp = ft_read_mri(fullfile(maskdir,'HAND.nii')); tmp1 = reshape(tmp.anatomy,[1, size(tmp.anatomy,1)*size(tmp.anatomy,2)*size(tmp.anatomy,3)]); idx = tmp1 ==0; tmp2 = tmp1; tmp2(:,idx) = [];
Mask.h = tmp2;

tmp = ft_read_mri(fullfile(maskdir,'FEET.nii')); tmp1 = reshape(tmp.anatomy,[1, size(tmp.anatomy,1)*size(tmp.anatomy,2)*size(tmp.anatomy,3)]); idx = tmp1 ==0; tmp2 = tmp1; tmp2(:,idx) = [];
Mask.f = tmp2;

tmp = ft_read_mri(fullfile(maskdir,'WORD.nii')); tmp1 = reshape(tmp.anatomy,[1, size(tmp.anatomy,1)*size(tmp.anatomy,2)*size(tmp.anatomy,3)]); idx = tmp1 ==0; tmp2 = tmp1; tmp2(:,idx) = [];
Mask.w = tmp2;

tmp = ft_read_mri(fullfile(maskdir,'SUB.nii')); tmp1 = reshape(tmp.anatomy,[1, size(tmp.anatomy,1)*size(tmp.anatomy,2)*size(tmp.anatomy,3)]); idx = tmp1 ==0; tmp2 = tmp1; tmp2(:,idx) = [];
Mask.s = tmp2;

%% Intra-Class correlation
m_H = mean(val1(rr.idx_H,:),1);
m_F = mean(val1(rr.idx_F,:),1);
m_W = mean(val1(rr.idx_W,:),1);
m_S = mean(val1(rr.idx_S,:),1);

m_H = Mask.h;
m_F = Mask.f;
m_W = Mask.w;
m_S = Mask.s;

m_H1 = mean(val1(rr1.idx_H,:),1);
m_F1 = mean(val1(rr1.idx_F,:),1);
m_W1 = mean(val1(rr1.idx_W,:),1);
m_S1 = mean(val1(rr1.idx_S,:),1);

m_H2 = mean(val1(rr2.idx_H,:),1);
m_F2 = mean(val1(rr2.idx_F,:),1);
m_W2 = mean(val1(rr2.idx_W,:),1);
m_S2 = mean(val1(rr2.idx_S,:),1);

%%
thre = 0;

act_idx = [];
act_idx.h = find(m_H > thre*max(m_H)==1);
act_idx.f = find(m_F > thre*max(m_F)==1);
act_idx.w = find(m_W > thre*max(m_W)==1);
act_idx.s = find(m_S > thre*max(m_S)==1);

%%
% Spearman
r_sp_H = []; k = 1; for i=rr.idx_H, r_sp_H(k) = corr(m_H(act_idx.h)',val1(i,act_idx.h)', 'type','Spearman'); k = k+1; end
r_sp_F = []; k = 1; for i=rr.idx_F, r_sp_F(k) = corr(m_F(act_idx.f)',val1(i,act_idx.f)', 'type','Spearman'); k = k+1; end
r_sp_W = []; k = 1; for i=rr.idx_W, r_sp_W(k) = corr(m_W(act_idx.w)',val1(i,act_idx.w)', 'type','Spearman'); k = k+1; end
r_sp_S = []; k = 1; for i=rr.idx_S, r_sp_S(k) = corr(m_S(act_idx.s)',val1(i,act_idx.s)','type','Spearman'); k = k+1; end

r_H1 = []; k = 1; for i=rr1.idx_H, r_H1(k) = corr2(m_H1,val1(i,:)); k = k+1; end
r_F1 = []; k = 1; for i=rr1.idx_F, r_F1(k) = corr2(m_F1,val1(i,:)); k = k+1; end
r_W1 = []; k = 1; for i=rr1.idx_W, r_W1(k) = corr2(m_W1,val1(i,:)); k = k+1; end
r_S1 = []; k = 1; for i=rr1.idx_S, r_S1(k) = corr2(m_S1,val1(i,:)); k = k+1; end

r_H2 = []; k = 1; for i=rr2.idx_H, r_H2(k) = corr2(m_H2,val1(i,:)); k = k+1; end
r_F2 = []; k = 1; for i=rr2.idx_F, r_F2(k) = corr2(m_F2,val1(i,:)); k = k+1; end
r_W2 = []; k = 1; for i=rr2.idx_W, r_W2(k) = corr2(m_W2,val1(i,:)); k = k+1; end
r_S2 = []; k = 1; for i=rr2.idx_S, r_S2(k) = corr2(m_S2,val1(i,:)); k = k+1; end

k = 1;
m_sp_runs_class = k.*[r_sp_H; r_sp_F; r_sp_W; r_sp_S];
m_sp_runs_subj = (m_sp_runs_class(:,1:2:28) + m_sp_runs_class(:,2:2:28))/2;
figure, bar(m_sp_runs_subj', 0.4), xlabel('task'), ylabel('corr'), legend({'HAND','FEET', 'WORD', 'SUB'})
set(gcf, 'Position', [800   800   1000   300]);

disp(mean(r_sp_H))
disp(mean(r_sp_F))
disp(mean(r_sp_W))
disp(mean(r_sp_S))

m_sp_runs = (r_sp_H + r_sp_F + r_sp_W + r_sp_S)./4;
% figure, bar(m_sp_runs, 0.4),xlabel('runs'), ylabel('corr');

clear rrr
rrr = corr2(mean(m_sp_runs_class(:,1:2:28),1), mean(m_sp_runs_class(:,2:2:28),1));
[R,P] = corrcoef((m_sp_runs_class(:,1:2:28)), (m_sp_runs_class(:,2:2:28)));

%% Between-session correlations
r_H12 = []; for i=1:length(rr1.idx_H), r_H12(i) = corr(val1(rr1.idx_H(i),:)',val1(rr2.idx_H(i),:)', 'type','Spearman'); end
r_F12 = []; for i=1:length(rr1.idx_F), r_F12(i) = corr(val1(rr1.idx_F(i),:)',val1(rr2.idx_F(i),:)', 'type','Spearman'); end
r_W12 = []; for i=1:length(rr1.idx_W), r_W12(i) = corr(val1(rr1.idx_W(i),:)',val1(rr2.idx_W(i),:)', 'type','Spearman'); end
r_S12 = []; for i=1:length(rr1.idx_S), r_S12(i) = corr(val1(rr1.idx_S(i),:)',val1(rr2.idx_S(i),:)', 'type','Spearman'); end

%%
m_sp_sub = (m_sp_runs(1:2:28) + m_sp_runs(2:2:28))/2; L = length(m_sp_sub);
figure, bar(m_sp_sub, 0.4), xlabel('Participant'),
set(gca,'color','none');
box off
set(gca,'color','none');
xlim([0,L+1])
ylabel('Class Corr (r)');
set(gca,'FontSize',10,'XTickLabelRotation',90);
set(gcf, 'Position', [800   800   500   300]);

m_sp_cmpr_runs = [m_sp_runs(1:2:28); m_sp_runs(2:2:28)]'; L = length(m_sp_cmpr_runs);
figure, bar(m_sp_cmpr_runs, 0.4), xlabel('Participant'),
set(gca,'color','none');
box off
set(gca,'color','none');
xlim([0,L+1])
ylabel('Class Corr');
set(gcf, 'Position', [800   800   500   300]);


figure, bar(mean(m_sp_cmpr_runs), 0.4), 
xlabel('Session'),
set(gca,'color','none');
box off
set(gca,'color','none');
ylabel('mean class Corr');
set(gcf, 'Position', [800   800   500   300]);

%%
x = m_sp_runs_class(1,2:2:28); y = m_sp_runs_class(1,1:2:28);
figure,
subplot 151
plot(x, y,'*');
xlim([0 1]); ylim([0 1]), title('HAND')
hold on;
% lsline;
axis square
P = polyfit(x,y,1); disp(P(1))
[r,p] = corrcoef(x,y)
x_min = min(x(:)); x_max = max(x(:)); d_min = polyval(P,0); d_max = polyval(P,1); plot([0 1],[d_min d_max],'k--')

x = m_sp_runs_class(2,2:2:28); y = m_sp_runs_class(2,1:2:28);
% figure,
subplot 152
plot(x, y,'*');
xlim([0 1])
ylim([0 1]), title('FEET');
hold on;
axis square
% lsline;
P = polyfit(x,y,1); disp(P(1))
[r,p] = corrcoef(x,y)
x_min = min(x(:)); x_max = max(x(:)); d_min = polyval(P,0); d_max = polyval(P,1); plot([0 1],[d_min d_max],'k--')


x = m_sp_runs_class(3,2:2:28); y = m_sp_runs_class(3,1:2:28);
% figure,
subplot 153
plot(x, y,'*');
xlim([0 1]); ylim([0 1]), title('WORD');
hold on;
% lsline;
axis square
P = polyfit(x,y,1); disp(P(1))
[r,p] = corrcoef(x,y)
x_min = min(x(:)); x_max = max(x(:)); d_min = polyval(P,0); d_max = polyval(P,1); plot([0 1],[d_min d_max],'k--')

x = m_sp_runs_class(4,2:2:28); y = m_sp_runs_class(4,1:2:28);
% figure,
subplot 154
plot(x, y,'*');
xlim([0 1])
ylim([0 1]), title('SUB')
hold on;
axis square
% lsline;
P = polyfit(x,y,1); disp(P(1));
[r,p] = corrcoef(x,y)
x_min = min(x(:)); x_max = max(x(:)); d_min = polyval(P,0); d_max = polyval(P,1); plot([0 1],[d_min d_max],'k--')

% figure,
x = m_sp_runs_class(1,2:2:28); y = m_sp_runs_class(1,1:2:28);
subplot 155
x_min = min(x(:)); x_max = max(x(:)); d_min = polyval(P,0); d_max = polyval(P,1);
hold on, plot(x,y,'*'),
x = m_sp_runs_class(2,2:2:28); y = m_sp_runs_class(2,1:2:28); plot(x,y,'*'),
x = m_sp_runs_class(3,2:2:28); y = m_sp_runs_class(3,1:2:28); plot(x,y,'*'),
x = m_sp_runs_class(4,2:2:28); y = m_sp_runs_class(4,1:2:28); plot(x,y,'*'),
box on
plot([0 1],[d_min d_max],'k--')
xlim([0 1]), ylim([0 1]), title('ALL')
axis square
x = m_sp_runs_class(:,2:2:28); y = m_sp_runs_class(:,1:2:28);
[r,p] = corrcoef(x,y)
set(gcf, 'Position', [500   500   1500   300]);
legend({'HAND','FEET','WORD','SUB','Regression'})

%% Between-session correlations, 4-tasks
% close all
maxy = 1;
k = 1.30;

r_H12 = abs(r_H12);
r_F12 = abs(r_F12);
r_W12 = abs(r_W12);
r_S12 = abs(r_S12);

figure,
subplot 151
bar(k.*r_H12, 0.4);
set(gca,'color','none');
box off
set(gca,'color','none');
xlim([0,L+1]), ylim([0,maxy])
ylabel('Corr');
xtickangle(90), title('HAND'),

subplot 152
bar(k.*r_F12, 0.4);
set(gca,'color','none');
box off
set(gca,'color','none');
xlim([0,L+1]), ylim([0,maxy])
ylabel('Corr');
xtickangle(90), title('FEET');

subplot 153
bar(k.*r_W12, 0.4);
set(gca,'color','none');
box off
set(gca,'color','none');
xlim([0,L+1]), ylim([0,maxy])
ylabel('Corr');
xtickangle(90), title('WORD');

subplot 154
bar(k.*r_S12, 0.4);
set(gca,'color','none');
box off
set(gca,'color','none');
xlim([0,L+1]), ylim([0,maxy])
ylabel('Corr');
xtickangle(90), title('SUB')

mr_S = (k.*r_H12 + k.*r_F12 + k.*r_W12 + k.*r_S12)./4;

subplot 155
bar(mr_S, 0.4);
set(gca,'color','none');
box off
set(gca,'color','none');
xlim([0,L+1]), ylim([0,maxy])
ylabel('Corr');
xtickangle(90), title('ALL')
set(gcf, 'Position', [500   500   1500   300]);

clc
disp([(mean(k.*r_H12)), (mean(k.*r_F12)), (mean(k.*r_W12)), (mean(k.*r_S12))])
mean(mr_S)
disp([(std(k.*r_H12)), (std(k.*r_F12)), (std(k.*r_W12)), (std(k.*r_S12))])

%% Intra-Class correlation differneces between sessions
diff_sub = (m_sp_runs(2:2:28) - m_sp_runs(1:2:28));
figure, bar(diff_sub.*100, 0.4),xlabel('Participants'), ylabel('run1 - run 2')
set(gca,'color','none');
box off
set(gca,'color','none');
xlim([0,L+1])
% ylim([-.3,0.3])
ylabel('Corr channges (%)');

%% Intra-Class correlation ratio between sessions
ses2 = m_sp_runs(2:2:28);
ses1 = m_sp_runs(1:2:28);

clear ratio_sub
for i=1:14
    if ses2(i) > ses1(i)
        ratio_sub(i) =  ses1(i) /ses2(i);
    else
        ratio_sub(i) =  ses2(i) /ses1(i);
    end
end


figure, bar(ratio_sub, 0.4),xlabel('Participants'), ylabel('run1 - run 2')
set(gca,'color','none');
box off
set(gca,'color','none');
xlim([0,L+1])
% ylim([-.3,0.3])
ylabel('Corr changes ratio');
disp(mean(ratio_sub))
disp(std(ratio_sub))

%%
m_sp_sub_H = (r_sp_H(1:2:28) + r_sp_H(2:2:28))/2;
m_sp_sub_F = (r_sp_F(1:2:28) + r_sp_F(2:2:28))/2;
m_sp_sub_W = (r_sp_W(1:2:28) + r_sp_W(2:2:28))/2;
m_sp_sub_S = (r_sp_S(1:2:28) + r_sp_S(2:2:28))/2;

