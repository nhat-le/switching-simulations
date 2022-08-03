% crawl behav sessions to determine the frame rates for each session
master_folder = '/Users/minhnhatle/Dropbox (MIT)/Nhat/Rigbox/e54';
folders = dir(fullfile(master_folder, '*/*/*Timeline.mat'));

for i = 1:numel(folders)
    load(fullfile(folders(i).folder, folders(i).name), 'Timeline');

    framerate = max(Timeline.rawDAQData(:,2)) / max(Timeline.rawDAQTimestamps);
    fprintf('Session %s: frame rate = %.2f\n', folders(i).name, framerate);

end