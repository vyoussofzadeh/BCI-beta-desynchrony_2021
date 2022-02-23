%%
function sFiles1 = vy_find_sfiles2(sFiles_name,dd3, tag)

clear idx idx1
for i=1:size(tag,1)
    idx{i} = find(contains(sFiles_name,tag(i,:))==1);
end
idx1 = sort([idx{1}, idx{2}]');

clear sFiles sFiles1
for i=1:length(idx1)
    sFiles{i} = dd3(idx1(i)).name;
    tkz = tokenize(sFiles{i},'/');
    sFiles1{i} = fullfile(tkz{end-2}, tkz{end-1}, tkz{end});
end

end