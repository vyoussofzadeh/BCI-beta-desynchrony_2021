%% BCI dataset, Ulster University

% preprocessing pipeline
% Writtern by Vahab Youssof Zadeh
% Update: 02/23/2022

clear, clc, close('all'); warning off

%% Initial settings
%- Brainstorm dir
bsdir = '/data/MEG/Research/BCI_Ulster/BCI5';
bs_path = '/opt/matlab_toolboxes/brainstorm3';
addpath ./functions
ft_path = '/opt/matlab_toolboxes/fieldtrip';

%
clc
disp('1: HAND')
disp('2: FEET');
disp('3: WORD');
disp('4: SUB');
disp('5: Cog_vs_Motor');
disp('6: Motor_vs_Cog');
disp('7: FEET_vs_HAND');
disp('8: HAND_vs_FEET');
disp('9: WORD_vs_SUB');
disp('10: S01_vs_S02');
disp('11: Sess1');
disp('12: Sess2');

cd(bsdir)

atlas_sel = input('Select data:');
switch atlas_sel
    case 1
        FileName = 'Group_analysis/@intra/results_average_201103_2106.mat';
        HAND = load(FileName);
        ImageGridAmp = HAND.ImageGridAmp; savetag = 'HAND';
        %         ImageGridAmp1 = HAND_stat.tmap(:,1); savetag = 'HAND_stats';
    case 2
        FileName = 'Group_analysis/@intra/results_average_201103_2105.mat';
        FEET = load(FileName);
        ImageGridAmp = FEET.ImageGridAmp; savetag = 'FEET';
    case 3
        FileName = 'Group_analysis/@intra/results_average_201103_2104.mat';
        WORD = load(FileName);
        ImageGridAmp = WORD.ImageGridAmp; savetag = 'WORD';
    case 4
        FileName = 'Group_analysis/@intra/results_average_201103_2103.mat';
        SUB = load(FileName);
        ImageGridAmp = SUB.ImageGridAmp; savetag = 'SUB';
    case 5
        FileName = 'Group_analysis/@intra/results_abs_201105_1820.mat';
        Cog_vs_Motor = load(FileName);
        ImageGridAmp = Cog_vs_Motor.ImageGridAmp; savetag = 'Cog_vs_Motor';
    case 6
        FileName = 'Group_analysis/@intra/results_abs_201105_1821.mat';
        Motor_vs_Cog = load(FileName);
        ImageGridAmp = Motor_vs_Cog.ImageGridAmp; savetag = 'Motor_vs_Cog';
    case 7
        FileName = 'Group_analysis/@intra/results_201105_2004.mat';
        FEET_vs_HAND = load(FileName);
        ImageGridAmp = FEET_vs_HAND.ImageGridAmp; savetag = 'FEET_vs_HAND';
    case 8
        FileName = 'Group_analysis/@intra/results_abs_200716_2323.mat';
        HAND_vs_FEET = load(FileName);
        ImageGridAmp = HAND_vs_FEET.ImageGridAmp; savetag = 'HAND_vs_FEET';
    case 9
        FileName = 'Group_analysis/@intra//results_abs_200716_2326.mat';
        HAND_vs_FEET = load(FileName);
        ImageGridAmp = HAND_vs_FEET.ImageGridAmp; savetag = 'WORD_vs_SUB';
    case 10
        FileName = 'Group_analysis/@intra/presults_max_220123_1135.mat';
        R1_vs_R2 = load(FileName);
        ImageGridAmp = R1_vs_R2.tmap(:,1); savetag = 'run1_vs_run2';
    case 11
        FileName = 'Group_analysis/@intra/presults_fdr_220122_2040.mat';
        R1 = load(FileName);
        ImageGridAmp = R1.tmap(:,1); savetag = 'run1';
    case 12
        FileName = 'Group_analysis/@intra/presults_fdr_220122_2037.mat';
        R2 = load(FileName);
        ImageGridAmp = R2.tmap(:,1); savetag = 'run2';
end

%%
atlas = load('scout_Desikan-Killiany_68_template.mat'); atag = 'DK';

[parcelval,rois] = do_sourceparcell_surface(atlas,ImageGridAmp);
[roiid,m,idx] = do_barplot(parcelval,rois, 15);

%% Laterality
m_left = mean(parcelval(1:2:end));
m_right = mean(parcelval(2:2:end));
LI = (m_left - m_right)./ (m_left + m_right)

