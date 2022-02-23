%% BCI dataset, Ulster University

% preprocessing pipeline
% Writtern by Vahab Youssof Zadeh
% Update: 02/23/2022

clear; clc, close('all'); warning off

%% Initial settings
%- Brainstorm dir
bsdir = '/data/MEG/Research/BCI_Ulster/BCI5';

addpath ./functions


ft_path = '/data/MEG/Vahab/Github/fieldtrip';
addpath(ft_path);
ft_defaults

%%
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
        stag = '4'; % Hand
    case 2
        stag = '8'; % Feet
    case 3
        stag = '16'; % Math
    case 4
        stag = '32'; % Word
end


%% Open BS (and select BCI database)
bs_path = '/opt/matlab_toolboxes/brainstorm3';
addpath(bs_path);
brainstorm
disp('slected BCI-DB from BS, then press enter');
pause,
db_reload_database('current',1)

%%
clear ProtocolSubjects
protocol = fullfile(bsdir,'data','protocol.mat');
load(protocol);

Subj = ProtocolSubjects.Subject;
L = length(Subj);

%% listing data files
d = rdir(fullfile(bsdir,'data','S*','**', 'channel_vectorview306_acc1.mat'));

% % ----
clear datafolder datafile subj_all run_all clz
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

%% Source modelling, DICS-BF
for ii=1:length(BS_chan.datafile)
    if ~contains(BS_chan.datafile{ii},'intra')
        [a, ~] = fileparts(BS_chan.datafile{ii});[c,d] = fileparts(a); [e,f] = fileparts(c);
        tkz = tokenize(d,'_');
        ee = [tkz{1},'_',tkz{2},'_',tkz{3}(1:3)];
        cd(a);
        clear dd1
        dd1 = rdir('./results_dics*.mat');
        
        
        stag1 = stag;
        switch ee
            case {'P002_IM_S01', 'P002_IM_S02', 'P004_IM_S01'}
                switch stag
                    case '4'
                        stag1 = '5';
                    case '8'
                        stag1 = '9';
                    case '16'
                        stag1 = '17';
                    case '32'
                        stag1 = '33';
                end
        end
        
        ex = [];
        if ~isempty(dd1)
            for i=1:length(dd1)
                tmp = load(dd1(i).name);
                tkz = tokenize(tmp.DataFile,'/');
                ex(i) = double(contains(tkz{3},['_',stag1, '_']));
            end
        end
        
        if sum(ex) == 0
            dd = rdir(['./',['data_',stag1, '*.mat']]);
            sFiles1 = [];
            for jj=1:length(dd)
                sFiles1{jj} = fullfile(f,d,dd(jj).name);
            end
            pause,
            bst_process('CallProcess', 'process_ft_sourceanalysis_tvDICS_BF_BCI', sFiles1, [], ...
                'method',     'dics', ...  % DICS beamformer
                'sensortype', 'MEG');  % MEG
            disp('done');
            clc,
        end
    end
end
