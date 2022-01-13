function out = pathsetup(project)
% Set up common paths for several projects



switch project
    case 'matchingsim'
        out.rigboxpath = '/Users/minhnhatle/Dropbox (MIT)/Nhat/Rigbox';
        out.rootpath = '/Users/minhnhatle/Dropbox (MIT)/Sur/MatchingSimulations';
        out.datapath = fullfile(out.rootpath, 'processed_data');
        out.codepath = fullfile(out.rootpath, 'PaperFigures/code');
        out.figpath = fullfile(out.rootpath, '/PaperFigures/figs');
        
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
        out.surfaceimgpath = '/Volumes/GoogleDrive/Other computers/ImagingDESKTOP-AR620FK/processed/session_surface_imgs';
        
    case 'wftoolbox'
        out.locanmfpath = '/Users/minhnhatle/Documents/ExternalCode/locaNMF-preprocess';
    otherwise
        error('Unrecognized project. Must be tca or matchingsim')
end


end