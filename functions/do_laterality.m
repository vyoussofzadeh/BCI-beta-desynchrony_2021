function lat = do_laterality(pow_parcel)

% datain = [mean(pow_parcel,1)];
% 
% figure,
% bar(datain);
% L = length(rois);
% set(gca,'Xtick', 1:L,'XtickLabel',1:L);
% % set(gca,'Xtick', 1:L,'XtickLabel',rois);
% box off
% set(gca,'color','none');
% xlim([0,L+1])
% xlabel('ROI');
% ylabel('Source power');
% set(gcf, 'Position', [500   300   500   300]);
% set(gca,'FontSize',10,'XTickLabelRotation',90);
% grid
% 
% for i=1:L
%     roiid{i} = [num2str(i),': ',rois{i}];
% end
% disp(roiid')
% 
% [m,idx]  = sort(datain,'descend');
% 
% 
% roiid(idx(1:nroi))'
% zpow_parcel = zscore(datain);
% zpow_parcel(idx(1:nroi))'


%% laterality
for i=1:size(pow_parcel,1)
    m_left = mean(pow_parcel(i,1:2:end),2);
    m_right = mean(pow_parcel(i,2:2:end),2);
    lat(i) = 100*(m_left-m_right)/(m_left+m_right);    
end

figure,
bar(lat);
L = length(lat);
set(gca,'Xtick', 1:L,'XtickLabel',1:L);
box off
set(gca,'color','none');
xlim([0,L+1])
xlabel('Subjects');