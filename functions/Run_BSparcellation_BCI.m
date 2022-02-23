%% Source parcellation
atlas = load('scout_Desikan-Killiany_68_template.mat'); atag = 'DK';

cd(bsdir);
Pa_04 = do_sourceparcel(sFiles_4,atlas,bsdir);
Pa_08 = do_sourceparcel(sFiles_8,atlas,bsdir);
Pa_16 = do_sourceparcel(sFiles_16,atlas,bsdir);
Pa_32 = do_sourceparcel(sFiles_32,atlas,bsdir);

clc,[roiid_04,m_04,idx_04] = do_barplot(mean(Pa_04.Parc,1),Pa_04.rois, 15);
clc,[roiid_08,m_08,idx_08] = do_barplot(mean(Pa_08.Parc,1),Pa_08.rois, 15);
clc,[roiid_16,m_16,idx_16] = do_barplot(mean(Pa_16.Parc,1),Pa_16.rois, 15);
clc,[roiid_32,m_32,idx_32] = do_barplot(mean(Pa_32.Parc,1),Pa_32.rois, 15);

%% Laterality analysis
clc
L = length(Pa_04.rois);
for i=1:L
    roiid{i} = [num2str(i),': ',Pa_04.rois{i}];
end
disp(roiid')
AngSmg   = [15,16,63,64];
Front    = [3,4,5,6,11,12,25,26,29,30,33,34,37,38,39,40,41,42,49,50,53,54,55,56,57,58];
LatFront = [5,6,11,12,37,38,39,40,41,42,55,56,57,58];
LatTemp  = [1,2,17,18,31,32,61,62,65,66,67,68];
PeriSyl  = [15,16,37,38,41,42,61,62,63,64];
Tanaka   = [37,38,41,42,61,62,63,64];
Temp     = [1,2,9,10,13,14,17,18,19,20,27,28,31,32,35,36,61,62,65,66,67,68];
Whole    = 1:68;

SensoriMotor = [33,34,49,50,5,6,55,56,57,58];
Salience = [];
Language   = [37,38,41,42,61,62,63,64]; % Tanaka ROIs

%
datain  = Pa_04; roisel = SensoriMotor; clc, close all, 
lat_04 = do_laterality(datain.Parc(:,roisel,:));

datain  = Pa_08; roisel = SensoriMotor; clc, close all, 
lat_08 = do_laterality(datain.Parc(:,roisel,:));

datain  = Pa_16; roisel = SensoriMotor; clc, close all, 
lat_16 = do_laterality(datain.Parc(:,roisel,:), datain.rois(roisel), length(roisel));

datain  = Pa_32; roisel = SensoriMotor; clc, close all, 
lat_32 = do_laterality(datain.Parc(:,roisel,:), datain.rois(roisel), length(roisel));
