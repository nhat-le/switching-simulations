function cols = paperaesthetics

% redcol = [178,24,43]/255;%[178,24,43]/255; %[5,113,176]/255
cols.redcol = [0.7 0 0];
cols.bluecol = [33,102,172]/255; %[202,0,32]/255
% bluecol = [0 0 0.7];

cols.colors = brewermap(9, 'Set1');
cols.colors = cols.colors([2,1,5,6,4,3,7,8,9],:);
% colors(4,:) = [0,0,0];

end