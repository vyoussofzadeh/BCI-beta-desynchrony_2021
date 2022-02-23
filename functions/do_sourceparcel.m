function Pa = do_sourceparcel(sFiles_in,atlas,bsdir)

Pa = [];
for iii=1:length(sFiles_in)
    d = [];
    clc
    disp([num2str(iii),'/', num2str(length(sFiles_in))]);
    tmp  = load(fullfile(bsdir,'data',sFiles_in{iii}));
%     tkz = tokenize(tmp.Comment,'/');
%     sub_all{iii} = tkz{1};
    %     tmp = load('200628_1728_sources_Str.mat'); tag = 'Str';
    ImageGridAmp = tmp.ImageGridAmp;
    [pow_parcel,rois] = do_sourceparcell_surface(atlas,ImageGridAmp);
    Parc(iii,:) = pow_parcel;
    pvoxel(iii,:) = ImageGridAmp;
end
Pa.Parc = Parc;
Pa.pvoxel = pvoxel;
Pa.rois = rois;
% Pa_all.submath = sub_all;

end