%% Sentence recognition task
clear; clc, close('all'); warning off

%% Initial settings
set(0,'DefaultFigureWindowStyle','docked');
% set(0,'DefaultFigureWindowStyle' , 'normal')
% set(gcf,'units','points','position',[500,500,500,500]);

cd '/data/MEG/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));
rmpath('./Failedattemps');

%- Input dir
% indir = '/mnt/file1/binder/binder_data/spendl';
indir = '/mnt/file1/binder/binder_data/spendl/Clozedata';
indir_second = '/data/MEG/Projects/spendl/Raw data';

%- Output dir
outdir = '/data/MEG/Projects/spendl/ft_process';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
% [allpath, atlas] = vy_init(cfg_init);

%%
bsdir = '/data/MEG/Research/BCI';
cd(bsdir)

%%
disp('1: 4')
disp('2: 8');
disp('3: 16');
disp('4: 32');
dc = input('Select data condition:');

%% Time-Freq Res (TFR) analysis
switch dc
    case 1
        stag = '4'; % hand
    case 2
        stag = '8'; % Feet
    case 3
        stag = '16'; % Math
    case 4
        stag = '32'; % Word
end

%%
clear ProtocolSubjects
protocol = fullfile(bsdir,'data','protocol.mat');
load(protocol);
d = rdir([bsdir,['/*',stag,'*IC_data.mat']]);
% d = rdir([bsdir,'/*IC_data.mat']);
L_data =length(d);

Subj = ProtocolSubjects.Subject;
L = length(Subj);
% d = rdir([indir,['/**/tSSS/*',tag,'*_raw.fif']]);
subj = [];
if ~isempty(d)
    for i=1:length(d)
        [pathstr, name] = fileparts(d(i).name);
        datafolder{i} = pathstr;
        datafile{i} = d(i).name;
        Index = strfind(datafile{i}, 'cloze');
        clz{i}= datafile{i}(Index(end)+5);
        Index1 = strfind(datafile{i}, '_');
        subj{i} = datafile{i}(Index1(2)+1:Index1(3)-1);
        %     subj_all{i} = Subj(i).Name;
    end
    ft_data = [];
    ft_data.clz = clz;
    ft_data.subj = subj;
    ft_data.datafile = datafile;
    ft_data.datafolder = datafolder;
    disp(datafile')
else
    disp('no data available to load');
end

%%
set(0,'DefaultFigureWindowStyle', 'normal')
bs_path = '/opt/matlab_toolboxes/brainstorm3';
addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
pause,
db_reload_database('current',1)

%%
d3 = rdir([bsdir,'/data/**/channel_vectorview306_acc1.mat']);
clear datafolder datafile subj_all run_all
d = d3;
for i=1:length(d)
    sStudies = bst_get('Study', i);
    [pathstr, name] = fileparts(sStudies.FileName);
    Index = strfind(pathstr, '/');
    subj{i} =  pathstr(1:Index-1);
end

bs_data = [];
bs_data.subj = subj;

%%
not_imported = [];
for ii=1:L_data
    [a, b] = fileparts(ft_data.datafile{ii});
    if exist(fullfile(bsdir,'data',ft_data.subj{ii},b), 'file') ~= 7
        ss = find(strcmp(ft_data.subj, ft_data.subj{ii}), 1 );
        for i=1:length(ProtocolStudies.Study)
            sSubject = bst_get('Subject', i, 1);
            if strcmp(sSubject.Name,ft_data.subj{ii});ss = i; break,
            end
        end
        if ~isempty(ss)
            import_data(ft_data.datafile{ii}, [], 'FT-TIMELOCK', [], ss, [], []);
            delete(ft_data.datafile{ii});
        else
            not_imported = [not_imported,ii];
        end
    else
        delete(ft_data.datafile{ii});
    end
    disp([num2str(ii),':',ft_data.datafile{ii}])
end

%% Updating channel files
clear d3
% d2 = rdir(fullfile(bsdir,['*',stag,'*BS_channels.mat']));
d3 = rdir(fullfile(bsdir,'data','S*','**', 'channel_vectorview306_acc1.mat'));

% % ----
clear datafolder datafile subj_all run_all clz
d = d3;
k=1;
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    if ~contains(d(i).name,'@raw')
    datafile{k} = d(i).name; k=1+k;
    end
end
BS_chan = [];
BS_chan.datafile = datafile;

%% Est. noise cov
Options = bst_noisecov();
Options.Baseline = [-0.3, -0.001];
Options.DataTimeWindow = [-0.3,1];

for ii=1:length(BS_chan.datafile)
    [a, ~] = fileparts(BS_chan.datafile{ii}); [c,d] = fileparts(a);
    cd(a)
    if exist('noisecov_full.mat', 'file') ~= 2
        for j=1:10000
            sStudies = bst_get('Study', j);
            if strcmp(sStudies.Name,d), break, end
        end
        sDatas = [sStudies.Data];
        iDataStudies = j.*ones(1,length(sDatas));
        db_reload_studies(j, 1);
        bst_noisecov(j, iDataStudies, 1:length(sDatas), Options, 0);
    end
end

%% Co-reg
close all,

% 
if exist(fullfile(bsdir,[stag,'checked_coreg.mat']), 'file') == 2
    load(fullfile(bsdir,[stag,'checked_coreg.mat']));
else
    checked_coreg = [];
end

if size(checked_coreg,1)>0
    for i=1:size(checked_coreg,1)
        checked_coreg_cell{i}=checked_coreg(i,:);
    end
else
    checked_coreg_cell = [];
end

%-
for ii=1:length(BS_chan.datafile)
    
    [a, ~] = fileparts(BS_chan.datafile{ii});
    [c,d] = fileparts(a);
    
    if isempty(find(strcmp(d(1:size(checked_coreg,2)),checked_coreg_cell)==1, 1))
        
        %     if exist(fullfile(subj_all{ii},b),'channel_vectorview306_acc1.mat'), 'file') ~= 2
        disp(a)
        ChannelFile = BS_chan.datafile{ii};
        
        % 1) refine using headshape
        isConfirm = 0; isWarning = 1;
        ChannelMat = [];
        channel_align_auto(ChannelFile, ChannelMat, isWarning, isConfirm)
        close all
        
        % 2) check
        isEdit = 0; Modality = 'MEG';
        SurfaceType = [];
        channel_align_manual( ChannelFile, Modality, isEdit, SurfaceType);
        
        disp('1: Yes; 2: No, edit; 3: save and stop');
        %         in = input('Coreg-check OK?');
        in = 1;
        close all;
        
        % 3) edit
        switch in
            case 2
                isEdit = 1; Modality = 'MEG';
                SurfaceType = [];
                channel_align_manual( ChannelFile, Modality, isEdit, SurfaceType);
                pause,
            case 3
                break
        end
        disp('done')
        checked_coreg = [checked_coreg;d(1:size(checked_coreg,2))];
        save(fullfile(bsdir,[stag,'checked_coreg.mat']),'checked_coreg');
    end
end

%% Est. head model
OPTIONS = [];
OPTIONS.comment = 'Overlapping spheres';
OPTIONS.MEGMethod =  'os_meg';
OPTIONS.EEGMethod  ='';
OPTIONS.ECOGMethod = '';
OPTIONS.SEEGMethod = '';
OPTIONS.SaveFile = 1;

for ii=1:length(BS_chan.datafile)
    [a, ~] = fileparts(BS_chan.datafile{ii}); [c,~] = fileparts(a);[e,f] = fileparts(c);
    OPTIONS.HeadModelFile =  a;
    OPTIONS.HeadModelType  = 'surface';
    A = load(BS_chan.datafile{ii});
    OPTIONS.Channel = A.Channel;
    OPTIONS.CortexFile = fullfile(f,'tess_cortex_concat_15000V.mat');
    OPTIONS.HeadFile =  fullfile(f,'tess_head_mask.mat');
    cd(a)
    if ~~exist(fullfile(bsdir,'anat',OPTIONS.CortexFile),'file') && ~exist('headmodel_surf_os_meg.mat','file')
        bst_headmodeler(OPTIONS);
    end
    disp(a)
end

%% deleting old source models,
d = rdir(fullfile(bsdir,'data','*',['*', stag,'*/results_dics*.mat']));

k=1; del_list = [];
for i=1:length(d)
    tmp = load(d(i).name);
    if contains(tmp.Comment,'20Hz') || ~contains(tmp.Comment,'tv')
        %             if ~contains(tmp.Comment,'tv')
        %     if contains(tmp.Comment,'dics')
        del_list{k} = d(i).name; 
        delete(del_list{k});
        k=k+1;
    end
end
disp(del_list')

%% Source modelling, DICS-BF
db_reload_database('current',1)
%
for ii=1:length(BS_chan.datafile)
    [a, ~] = fileparts(BS_chan.datafile{ii});[c,d] = fileparts(a); [e,f] = fileparts(c);
    tkz = tokenize(d,'_');
    ee = [tkz{1},'_',tkz{2},'_',tkz{3}(1:3)];
    cd(a);
    clear dd1
    dd1 = rdir('./results_dics*.mat');

    if isempty(dd1)
%         db_reload_studies(ii, 1);
        dd = rdir('./*_IC*.mat');
        sFiles1 = [];
        for jj=1:length(dd)
            sFiles1{jj} = fullfile(f,d,dd(jj).name);
        end
        bst_process('CallProcess', 'process_ft_sourceanalysis_tvDICS_BF', sFiles1, [], ...
            'method',     'dics', ...  % DICS beamformer
            'sensortype', 'MEG');  % MEG
        disp('done');
        clc,
    end
end

%% Detect redundant files, remove files manually using BS GUI
clc
cloze = [];
redunant_data = [];
redunant_passed = [];
for ii=1:length(BS_chan.subj_all)
    [a, ~] = fileparts(BS_chan.datafile{ii});[c,d] = fileparts(a); [e,f] = fileparts(c);
    tkz = tokenize(d,'_');
    dd = rdir(['./cloze*',tkz{3},'*/results_dics*.mat']);
    ee = [[tkz{1},'_',tkz{2},'_',tkz{3}(1:3)]];
    sFiles1 = [];
    for jj=1:length(dd)
        sFiles1{jj} = fullfile(f,dd(jj).name);
        tkz = tokenize(dd(jj).name,'_'); tmp = tkz{1};
        cloze(jj) = str2double(tmp(end));
    end
    if length(cloze) == length(unique(cloze))
        %     disp(sFiles1');
        redunant_passed = [redunant_passed;ee];
    else
        redunant_data = [redunant_data;ee];
    end
end
disp('Redundancy check, passed')
disp(redunant_passed)
disp('=====')
disp('not passed, check files in BS')
disp(redunant_data)

%% Deleting redundant saved files, if so!
cd(bsdir);
dd = rdir(fullfile (bsdir,'data','/*/@intra/results*.mat'));
dd1 = rdir(fullfile (bsdir,'data','/Group*/@intra/results*.mat'));

d = []; d2 = [];
for ii=1:length(dd), d{ii}=dd(ii).name; disp(dd(ii).name); end
for ii=1:length(dd1), d2{ii}=dd1(ii).name; disp(dd1(ii).name); end
dd3 = setdiff(d,d2);

subj = []; dconn = [];
for ii=1:length(dd3)
    [a, ~] = fileparts(dd3{ii});
    cd(a)
    tmp = load(dd3{ii});
    comm_data = tmp.Comment;
    tkz = tokenize(comm_data,'_');
    subj{ii}= tkz{2};
    dconn{ii} = tkz{3};
    sFiles_name{ii} = [subj{ii},'_',dconn{ii}];
end

[idx,b] = unique(sFiles_name);
sFiles_name2 = [];
for m=1:length(b)
    sFiles_name2{m} = sFiles_name{b(m)};
end

cd(fullfile(bsdir,'data'));
for k=1:length(dd3)
    if ~(ismember(k,b))
        delete(dd3{k})
        disp(dd3{k})
    end
end

%% 
disp('1: overwrite intra-avg across clozes');
disp('2: not overwrite intra-avg across clozes')
ow = input(':');

subj_all1 = unique(subj_all);

if ow == 1
    clc
    close all,
    for ii = 1:length(subj_all1)
        cd(fullfile(bsdir,'data',subj_all1{ii}))
        subj = []; dconn = []; sFiles_name = [];
        dd = rdir('./@intra/results*.mat');
        for i=1:length(dd)
            tmp = load(dd(i).name);
            comm_data = tmp.Comment;
            tkz1 = tokenize(comm_data,'_');
            subj{i}= tkz1{2};
            dconn{i} = tkz1{3}(1:3);
            sFiles_name{i} = [subj{i},'_',dconn{i}];
        end        
        for i=1:length(sFiles_name)
            if contains(sFiles_name{i},stag(1:3))
                disp(dd(i).name)
                delete(dd(i).name)
            end
        end      
    end
    db_reload_database('current',1)
end

%% Intra-subject averaging
cd(bsdir);
dd = rdir(fullfile (bsdir,'data','/*/@intra/results*.mat'));
dd1 = rdir(fullfile (bsdir,'data','/Group*/@intra/results*.mat'));

d = []; d2 = [];
for ii=1:length(dd), d{ii}=dd(ii).name; disp(dd(ii).name); end
for ii=1:length(dd1), d2{ii}=dd1(ii).name; disp(dd1(ii).name); end
dd3 = setdiff(d,d2);

subj = []; dconn = []; sFiles_name_inter = [];
for ii=1:length(dd3)
    [a, ~] = fileparts(dd3{ii});
    cd(a)
    tmp = load(dd3{ii});
    comm_data = tmp.Comment;
    tkz = tokenize(comm_data,'_');
    subj{ii}= tkz{2};
    dconn{ii} = tkz{3};
    sFiles_name_inter{ii} = [subj{ii},'_',dconn{ii}];
end


clc
close all,
for ii = 1:length(subj_all1)
    
    cd(fullfile(bsdir,'data',subj_all1{ii}))
    dd = rdir(['./cloze*',stag,'*/results_dics*.mat']);
    sFiles1 = []; sFiles_name = [];
    for jj=1:length(dd)
        sFiles1{jj} = fullfile(subj_all1{ii},dd(jj).name);
        [aa,bb] = fileparts(dd(jj).name);
        sFiles_name{jj} = aa;
    end
    
    if length(sFiles1)==5
    else
        warning(['check results of:',subj_all1{ii}])
    end
    if isempty(find(contains(sFiles_name_inter,[subj_all1{ii},'_',stag(1:3)])==1, 1))
        % Process: Average: Everything
        bst_process('CallProcess', 'process_average', sFiles1, [], ...
            'avgtype',         1, ...  % Everything
            'avg_func',        1, ...  % Arithmetic average:  mean(x)
            'weighted',        0, ...
            'scalenormalized', 0);
    end

end
disp('intra-subject source averaging was completed!');

%% 
disp('1: overwrite inter-avg (deleting old estimates)');
disp('2: not overwrite inter-avg across clozes')
ow = input(':');

if ow ==1
    dd = rdir(fullfile (bsdir,'data','/Group*/*/results*.mat'));
    
    k = 1; sFiles = [];
    for ii=1:length(dd)
        [a, ~] = fileparts(dd(ii).name);
        cd(a)
        tmp = load(dd(ii).name);
        comm_data = tmp.Comment;
        tkz = tokenize(comm_data,'_');
        if contains(comm_data,stag) && ~contains(comm_data,'files')
            sFiles{k} = dd(ii).name; k=k+1;
        end
    end
    
    clear subj dconn sFiles_name
    for ii=1:length(sFiles)
        tmp = load(sFiles{ii});
        comm_data = tmp.Comment;
        tkz = tokenize(comm_data,'_');
        subj{ii}= tkz{2};
        dconn{ii} = tkz{3};
        sFiles_name{ii} = [subj{ii},'_',dconn{ii}];
        delete(sFiles{ii})
    end
    db_reload_database('current',1)
end

%% 
cd(bsdir);
dd = rdir(fullfile (bsdir,'data','/*/@intra/results*.mat'));
dd1 = rdir(fullfile (bsdir,'data','/Group*/@intra/results*.mat'));

d = []; d2 = [];
for ii=1:length(dd), d{ii}=dd(ii).name; disp(dd(ii).name); end
for ii=1:length(dd1), d2{ii}=dd1(ii).name; disp(dd1(ii).name); end
dd3 = setdiff(d,d2);

subj = []; dconn = [];
for ii=1:length(dd3)
    [a, ~] = fileparts(dd3{ii});
    cd(a)
    tmp = load(dd3{ii});
    comm_data = tmp.Comment;
    tkz = tokenize(comm_data,'_');
    subj{ii}= tkz{2};
    dconn{ii} = tkz{3};
    sFiles_name{ii} = [subj{ii},'_',dconn{ii}];
end

[~,b] = unique(sFiles_name);
sFiles_name2 = [];
for m=1:length(b)
    sFiles_name2{m} = sFiles_name{b(m)};
end

cd(fullfile(bsdir,'data'));
for k=1:length(dd3)
    if ~(ismember(k,b))
        delete(dd3{k})
        disp(dd3{k})
    end
end

% db_reload_database('current',1)

%%
dd1 = rdir(fullfile (bsdir,'data','/Group*/@intra/results*.mat'));
for ii=1:length(dd1), d2{ii}=dd1(ii).name; disp(dd1(ii).name); end

subj = []; dconn = []; sFiles_name_check = []; k=1;
for ii=1:length(dd1)
    [a, ~] = fileparts(d2{ii});
    cd(a)
    tmp = load(d2{ii});
    comm_data = tmp.Comment;
    tkz = tokenize(comm_data,'_');
    if length(tkz)>2
        subj{k}= tkz{2};
        dconn{k} = tkz{3};
        sFiles_name_check{k} = [subj{k},'_',dconn{k}];
        k=k+1;
    end
end

%% Project on default anatomy (for group mapping)
idx = find(contains(sFiles_name,stag)==1);
destSurfFile = '@default_subject/tess_cortex_crtex_brvisa_15000V.mat';

sFiles_name1 = [];
if length(idx)>1
    for i=1:length(idx)
        sFiles{i} = dd3{idx(i)};
        tkz = tokenize(sFiles{i},'/');
        sFiles_name1{i}= [tkz{end-2},'_', stag];
        if isempty(find(contains(sFiles_name_check, sFiles_name1(i))==1, 1))
            sFiles1{i} = fullfile(tkz{end-2}, tkz{end-1}, tkz{end});
            bst_project_sources ({sFiles1{i}}, destSurfFile);
        end
    end
end

%% Inter-subject (group) averaging
clc
dd = rdir(fullfile (bsdir,'data','/Group*/*/results*.mat'));

k = 1; sFiles = [];
for ii=1:length(dd)
    [a, ~] = fileparts(dd(ii).name);
    cd(a)
    tmp = load(dd(ii).name);
    comm_data = tmp.Comment;
    tkz = tokenize(comm_data,'_');
    if contains(comm_data,stag(1:3)) && ~contains(comm_data,'files')
        sFiles{k} = dd(ii).name; k=k+1;
    end
end

clear subj dconn sFiles_name
k=1;
for ii=1:length(sFiles)
    tmp = load(sFiles{ii});
    comm_data = tmp.Comment;
    if contains(comm_data, [stag,'_IC'])
        tkz = tokenize(comm_data,'_');
        subj{k}= tkz{2};
        dconn{k} = tkz{3};
        sFiles_name{k} = [subj{ii},'_',dconn{ii}];
        comment{k} = comm_data;
        k=k+1;
    end
end
disp(comment')

[idx,a] = unique(sFiles_name);
idx_del = setdiff(1:length(sFiles),a);

if idx_del>0
    for i=1:length(idx_del)
%         delete(sFiles{idx_del(i)})
%         disp(sFiles{idx_del(i)})
        disp(comment{idx_del(i)})
    end
end


%%
clc
dd = rdir(fullfile (bsdir,'data','/Group*/*/results*.mat'));

k = 1; sFiles = [];
for ii=1:length(dd)
    [a, ~] = fileparts(dd(ii).name);
    cd(a)
    tmp = load(dd(ii).name);
    comm_data = tmp.Comment;
    tkz = tokenize(comm_data,'_');
    if contains(comm_data,stag) && ~contains(comm_data,'files')
        sFiles{k} = dd(ii).name; k=k+1;
    end
end
disp(sFiles')

% Process: Average: Everything
bst_process('CallProcess', 'process_average', sFiles, [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'scalenormalized', 0);
