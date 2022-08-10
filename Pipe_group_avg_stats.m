%% BCI dataset, Ulster University

% Pipeline_group_avg_stat pipeline
% Writtern by Vahab Youssof Zadeh
% Update: 02/23/2022

clear, clc, close('all'); warning off

%% Initial settings
%- Brainstorm dir
bsdir = '/data/MEG/Research/BCI_Ulster/BCI5';
bs_path = '/opt/matlab_toolboxes/brainstorm3';
addpath ./functions
ft_path = '/opt/matlab_toolboxes/fieldtrip';

%% Open BS (and select BCI database)
addpath(bs_path);
brainstorm
disp('slected BCI-DB from BS, then press enter');
pause,
% db_reload_database('current',1)

%% Adding FT
addpath(ft_path);
ft_defaults

%%
d = rdir([bsdir,'/**/','channel_vectorview306_acc1.mat']);
clear datafolder datafile subj_all run_all datafl datafdr
j = 1;
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    %     disp(datafile{i})
    %     pause
    idx = strfind(datafile{i}, 'Sub');
    idx1 = strfind(datafile{i}, 'IM');
    idx2 = strfind(datafile{i}, '/');
    idx3 = strfind(datafile{i}, '_');
    tkz = tokenize(datafile{i},'/');
    
    if ~(isempty(tkz) || isempty(idx1) || isempty(idx2) || isempty(idx3))
        subj_all{j} = datafile{i}(idx+7:idx+8);
        run_all{j} = datafile{i}(idx1+4:idx1+5);
        datafdr{j} = pathstr;
        datafl{j} = d(i).name;
        j = j+1;
    end
end
BS_chan = [];
BS_chan.run_all = run_all;
BS_chan.subj_all = subj_all;
BS_chan.datafile = datafl;
BS_chan.datafolder = datafdr;

%%
cd(bsdir);
dd = rdir(fullfile (bsdir,'data','Sub*/**/results*.mat'));

subj = unique(BS_chan.subj_all);

clear sFiles_name
subj = []; dconn = [];
kk=1;
for ii=1:length(dd)
    [a, ~] = fileparts(dd(ii).name);
    cd(a)
    tmp = load(dd(ii).name);
    if ~contains(dd(ii).name, '@intra')
        comm_data = tmp.Comment;
        cd(dd(kk).folder)
        tkz = tokenize(dd(kk).name,'/');
        tkz1 = tokenize(tmp.DataFile,'/'); subj{kk} = tkz1{1};
        tkz2 = tokenize(tkz1{2},'_'); run{kk} = tkz2{3};
        tkz3 = tokenize(tkz1{end},'_'); dconn{kk} = tkz3{2};
        sFiles_name{kk} = [subj{kk},'_',run{kk}, '_', dconn{kk}];
        kk=kk+1;
    end
end
disp(sFiles_name')

%% Inter-subject (group) averaging
tag = ['_4';'_5'];
sFiles_4 = do_find_sfiles(sFiles_name,dd, tag);
% Process: Average: Everything
bst_process('CallProcess', 'process_average', sFiles_4, [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'Comment',       'avg_04', ...
    'scalenormalized', 0);

%- 
tag = ['_8';'_9'];
sFiles_8 = do_find_sfiles(sFiles_name,dd, tag);
% Process: Average: Everything
bst_process('CallProcess', 'process_average', sFiles_8, [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'Comment',       'avg_08', ...
    'scalenormalized', 0);

%- 
tag = ['_16';'_17'];
sFiles_16 = do_find_sfiles(sFiles_name,dd, tag);

% Process: Average: Everything
bst_process('CallProcess', 'process_average', sFiles_16, [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'Comment',       'avg_16', ...
    'scalenormalized', 0);

%- 
tag = ['_32';'_33'];
sFiles_32 = do_find_sfiles(sFiles_name,dd, tag);

% Process: Average: Everything
bst_process('CallProcess', 'process_average', sFiles_32, [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'Comment',       'avg_32', ...
    'scalenormalized', 0);

%%
% Session 1 
bst_process('CallProcess', 'process_average', sFiles_4(1:2:28), [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'Comment',       'avg_04_run1', ...
    'scalenormalized', 0);

% Session 2 
bst_process('CallProcess', 'process_average', sFiles_4(2:2:28), [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'Comment',       'avg_04_run2', ...
    'scalenormalized', 0);

% Session 1  
bst_process('CallProcess', 'process_average', sFiles_8(1:2:28), [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'Comment',       'avg_08_run1', ...
    'scalenormalized', 0);

% Session 2 
bst_process('CallProcess', 'process_average', sFiles_8(2:2:28), [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'Comment',       'avg_08_run2', ...
    'scalenormalized', 0);

% Session 1   
bst_process('CallProcess', 'process_average', sFiles_16(1:2:28), [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'Comment',       'avg_16_run1', ...
    'scalenormalized', 0);

% Session 2   
bst_process('CallProcess', 'process_average', sFiles_16(2:2:28), [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'Comment',       'avg_16_run2', ...
    'scalenormalized', 0);

% Session 1 
% Process: Average: Everything
bst_process('CallProcess', 'process_average', sFiles_32(1:2:28), [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'Comment',       'avg_32_run1', ...
    'scalenormalized', 0);

% Session 2 
% Process: Average: Everything
bst_process('CallProcess', 'process_average', sFiles_32(2:2:28), [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'Comment',       'avg_32_run2', ...
    'scalenormalized', 0);

%%
Run_BSparcellation_BCI

%%
Run_BSstats_BCI

%% GII
cd(bsdir); clear sFiles
tag = ['_4';'_5'];stag = '4';
sFiles = do_find_sfiles(sFiles_name,dd, tag);

cd(bsdir); clear sFiles
tag = ['_8';'_9'];stag = '8';
sFiles = do_find_sfiles(sFiles_name,dd, tag);

cd(bsdir); clear sFiles
tag = ['_16';'_17'];stag = '16';
sFiles = do_find_sfiles(sFiles_name,dd, tag);

cd(bsdir); clear sFiles
tag = ['_32';'_33'];stag = '32';
sFiles = do_find_sfiles(sFiles_name,dd, tag);

%% Process: Export to SPM12 (surface) - for SPM stat analysis
% Process: Export to SPM12 (surface) - for SPM stat analysis
RawFiles = {...
    ['/data/MEG/Research/BCI_Ulster/SPM_PRT3/gii/',stag]};

if exist(RawFiles{1}, 'file') == 0, mkdir(RawFiles{1}); end
cd(RawFiles{1})
% Process: Export to SPM12 (surface)
bst_process('CallProcess', 'process_export_spmsurf', sFiles, [], ...
    'outputdir',  {RawFiles{1}, 'GIFTI'}, ...
    'filetag',    ['SPM export, ', stag], ...
    'timewindow', [1, 1], ...
    'isabs',      1);

%% Process: Export to SPM8/SPM12 (volume) - for machine-learning analysis
RawFiles = {...
    ['/data/MEG/Research/BCI/SPM_PRT3/nii/',stag]};

if exist(RawFiles{1}, 'file') == 0, mkdir(RawFiles{1}); end
cd(RawFiles{1})

for i=1:length(sFiles)
    bst_process('CallProcess', 'process_export_spmvol', sFiles{i}, [], ...
        'outputdir',      {RawFiles{1}, 'Nifti1'}, ...
        'filetag',        ['SPM export, ', stag], ...
        'isconcat',       1, ...
        'timewindow',     [1, 1], ...
        'timedownsample', 3, ...
        'timemethod',     1, ...  % Average time (3D volume)
        'voldownsample',  2, ...
        'isabs',          1, ...
        'iscut',          1);
    
end

cd(RawFiles{1})
