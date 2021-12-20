function out = pathsetup(project)
% Set up common paths for several projects



switch project
    case 'matchingsim'
        
        out.datapath = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/processed_data';
        out.codepath = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/code';
        out.figpath = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations/PaperFigures/figs';
        
        out.expdatapath = fullfile(out.datapath, 'expdata');
        out.blockhmmfitpath = fullfile(out.datapath, 'blockhmmfit');
        out.simdatapath = fullfile(out.datapath, 'simdata');
        out.svmdatapath = fullfile(out.datapath, 'svm');
        out.svmconfigpath = fullfile(out.datapath, 'svm/configs');
        out.svmmodelpath = fullfile(out.datapath, 'svm/models');

        out.blockhmm_codepath = fullfile(out.codepath, 'blockhmm');
        out.characterize_codepath = fullfile(out.codepath, 'characterization');
        out.decoding_codepath = fullfile(out.codepath, 'decoding');
        out.expfit_codepath = fullfile(out.codepath, 'expfit');
        out.schematic_codepath = fullfile(out.codepath, 'schematic');
        
    case 'tca'
        out.datapath = '/Users/minhnhatle/Documents/ExternalCode/tca/data';
        out.codepath = '/Users/minhnhatle/Documents/ExternalCode/tca/src/matlab';
        out.rawdatapath = '/Volumes/GoogleDrive/Other computers/ImagingDESKTOP-AR620FK/processed/raw/extracted';
        out.tcamatpath = '/Volumes/GoogleDrive/Other computers/ImagingDESKTOP-AR620FK/processed/tca-factors';
        
        
    otherwise
        error('Unrecognized project. Must be tca or matchingsim')
end


end