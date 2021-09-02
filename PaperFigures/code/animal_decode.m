load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/code/svm_Mdl.mat');
load('/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/expdata/expfit_params.mat');

% Decode for all animal sessions
lapseFlat = reshape(explapses_all, [], 1);
effFlat = reshape(expeff_all, [], 1);
offsetFlat = reshape(expoffsets_all, [], 1);
slopesFlat = reshape(expslopes_all, [], 1);

features_flat = [effFlat lapseFlat slopesFlat -offsetFlat];


statesFlat = Mdl.predict(features_flat);
% eliminate nan's (nan's decode to 1)
statesFlat(isnan(lapseFlat) | isnan(effFlat) | isnan(offsetFlat) | isnan(slopesFlat)) = nan;


% Unflatten
states = reshape(statesFlat, size(expeff_all, 1), []);


%%
plot(nanmean(states))


%% 
for i = 1:15
    subplot(3,5,i)
    plot(states(i,:), 'o')
    ylim([0, 5])
    
end