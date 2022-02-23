function [pow_parcel,rois] = do_sourceparcell_surface(atlas,ImageGridAmp)

rois = [];
pow_parcel = [];
for i=1:length(atlas.Scouts)
    pow_parcel(i) = mean(ImageGridAmp(atlas.Scouts(i).Vertices));
    rois{i} = atlas.Scouts(i).Label;
end

end