% Perform a linear regression to uncover opto effects
% on brain regions, power, choice/outcome

expfitdate = '040122';

load(fullfile('optodata/', expfitdate, 'modecounts.mat'));

NSTATES = numel(modecounts(1).count_opto{1,1,1}); 


%% Count and aggregate
stateID = [4];
mastertbl = table;
for powerID = 1:2
    for periodID = 1:2
        tbl = build_table(modecounts, stateID, powerID, periodID);
        mastertbl = cat(1, mastertbl, tbl);
    end
end
mastertbl.fracs = max(mastertbl.fracs, 1e-2); %for stability
mastertbl.fracs = min(mastertbl.fracs, 1 - 1e-2); %for stability



%% Filter out conditions with low counts
opts = struct;
opts(1).power = 'high';
opts(1).period = 'outcome';
opts(1).animal = 'f26';

opts(2).power = 'high';
opts(2).period = 'outcome';
opts(2).animal = 'f27';

opts(3).power = 'high';
opts(3).period = 'outcome';
opts(3).animal = 'f29';

opts(4).power = 'high';
opts(4).period = 'choice';
opts(4).animal = 'f32';
tbl = filter_tbl(mastertbl, opts);


%% 
tbl.y_trans = log(tbl.fracs ./ (1 - tbl.fracs));
tbl.CVH = tbl.opto & strcmp(tbl.period, 'choice') & strcmp(tbl.regions, 'visual') & ...
    strcmp(tbl.power, 'high');
tbl.CMH = tbl.opto & strcmp(tbl.period, 'choice') & strcmp(tbl.regions, 'motor') & ...
    strcmp(tbl.power, 'high');
tbl.CRH = tbl.opto & strcmp(tbl.period, 'choice') & strcmp(tbl.regions, 'rsc') & ...
    strcmp(tbl.power, 'high');
tbl.CFH = tbl.opto & strcmp(tbl.period, 'choice') & strcmp(tbl.regions, 'frontal') & ...
    strcmp(tbl.power, 'high');

tbl.CVL = tbl.opto & strcmp(tbl.period, 'choice') & strcmp(tbl.regions, 'visual') & ...
    strcmp(tbl.power, 'low');
tbl.CML = tbl.opto & strcmp(tbl.period, 'choice') & strcmp(tbl.regions, 'motor') & ...
    strcmp(tbl.power, 'low');
tbl.CRL = tbl.opto & strcmp(tbl.period, 'choice') & strcmp(tbl.regions, 'rsc') & ...
    strcmp(tbl.power, 'low');
tbl.CFL = tbl.opto & strcmp(tbl.period, 'choice') & strcmp(tbl.regions, 'frontal') & ...
    strcmp(tbl.power, 'low');

tbl.OVH = tbl.opto & strcmp(tbl.period, 'outcome') & strcmp(tbl.regions, 'visual') & ...
    strcmp(tbl.power, 'high');
tbl.OMH = tbl.opto & strcmp(tbl.period, 'outcome') & strcmp(tbl.regions, 'motor') & ...
    strcmp(tbl.power, 'high');
tbl.ORH = tbl.opto & strcmp(tbl.period, 'outcome') & strcmp(tbl.regions, 'rsc') & ...
    strcmp(tbl.power, 'high');
tbl.OFH = tbl.opto & strcmp(tbl.period, 'outcome') & strcmp(tbl.regions, 'frontal') & ...
    strcmp(tbl.power, 'high');

tbl.OVL = tbl.opto & strcmp(tbl.period, 'outcome') & strcmp(tbl.regions, 'visual') & ...
    strcmp(tbl.power, 'low');
tbl.OML = tbl.opto & strcmp(tbl.period, 'outcome') & strcmp(tbl.regions, 'motor') & ...
    strcmp(tbl.power, 'low');
tbl.ORL = tbl.opto & strcmp(tbl.period, 'outcome') & strcmp(tbl.regions, 'rsc') & ...
    strcmp(tbl.power, 'low');
tbl.OFL = tbl.opto & strcmp(tbl.period, 'outcome') & strcmp(tbl.regions, 'frontal') & ...
    strcmp(tbl.power, 'low');

% ------------------

tbl.CV = tbl.opto & strcmp(tbl.period, 'choice') & strcmp(tbl.regions, 'visual');
tbl.CM = tbl.opto & strcmp(tbl.period, 'choice') & strcmp(tbl.regions, 'motor');
tbl.CR = tbl.opto & strcmp(tbl.period, 'choice') & strcmp(tbl.regions, 'rsc');
tbl.CF = tbl.opto & strcmp(tbl.period, 'choice') & strcmp(tbl.regions, 'frontal');

tbl.OV = tbl.opto & strcmp(tbl.period, 'outcome') & strcmp(tbl.regions, 'visual');
tbl.OM = tbl.opto & strcmp(tbl.period, 'outcome') & strcmp(tbl.regions, 'motor');
tbl.OR = tbl.opto & strcmp(tbl.period, 'outcome') & strcmp(tbl.regions, 'rsc');
tbl.OF = tbl.opto & strcmp(tbl.period, 'outcome') & strcmp(tbl.regions, 'frontal');


% Performs a regression
X = double(table2array(tbl(:,8:23)));
y = double(tbl.y_trans);

% mdl = fitlm(tbl, ['y_trans ~ CVH + CMH + CRH + CFH + CVL + CML + CRL + CFL' ...
%     '+ OVH + OMH + ORH + OFH + OVL + OML + ORL + OFL']);

mdl = fitlm(tbl, ['y_trans ~ CVL + CML + CRL + CFL' ...
    '+ OVL + OML + ORL + OFL']);

% mdl = fitlm(tbl, ['y_trans ~ CV + CM + CR + CF' ...
%     '+ OV + OM + OR + OF + OV + OM + OR + OF']);
% mdl = fitlm(X, y);



pvals = mdl.Coefficients.pValue(2:end);
level = 0.05;
bhline = (1:numel(pvals)) * level / numel(pvals);
figure;
plot(sort(pvals), '.', 'MarkerSize', 30);
hold on;
plot(bhline, 'k--', 'LineWidth', 2);
mymakeaxis('x_label', 'Hypothesis #', 'y_label', 'Sorted p-values', ...
    'font_size', 25, 'xticks', 1:numel(pvals))

%%
% tbl


function tbl = filter_tbl(rawtbl, conditions)
% rawtbl: a table with column names
% conditions: a struct where each element specifies an 'and' condition
% to be removed from the table
% for e.g. if conditions(1).animal = 'f26'
% conditions(1).power = 'low'
% then rows where animal == 'f26' and power == 'low' will be filtered out.
% Returns the filtered table, tbl.
tbl = rawtbl;
vars = fieldnames(conditions);
for i = 1:numel(conditions)
    idxfilt = ones(size(tbl, 1), 1);
    
    for j = 1:numel(vars)
        varname = vars{j};
        cond = conditions(i).(varname);
   
        if ischar(class(cond))
            idx_cond = strcmp(tbl.(varname), cond);
            
            
        elseif isa(class(cond), 'double')
            idx_cond = tbl.(varname) == cond;
  
        end
        
        idxfilt = (idxfilt & idx_cond);
        
    end
    
    tbl = tbl(~idxfilt, :);
    
end
    
    
    
    
    
end












function tbl = build_table(modecounts, stateID, idxpower, idxperiod)
% Build a table to observe effects of opto on brain areas, power,
% choice/outcome epoch
% Inputs: stateID: 1-6, HMM mode, if this is an array, will add
% all the contributions from each state in the array
% power_criterion: int, 1/2 corresponding to low/high
% period_criterion: int, 1/2 corresponding to choice/outcome
% Returns: a table with the following columns:
% animal, power, period, opto, regions, fracs

fracs = [];
regions = {};
power = {};
period = {};
opto = [];
animal = {};

regionsOrder = modecounts(1).area_criteria;


for animalID = 1:numel(modecounts)
    counts_opto = helper.block_counter(modecounts(animalID).count_opto, idxperiod, idxpower);
    opto_states_frac = counts_opto ./ sum(counts_opto, 2);
    opto_state_frac = sum(opto_states_frac(:, stateID), 2);

    counts_no_opto = helper.block_counter(modecounts(animalID).count_no_opto, idxperiod, idxpower);
    no_opto_states_frac = counts_no_opto ./ sum(counts_no_opto, 2);
    no_opto_state_frac = sum(no_opto_states_frac(:, stateID), 2);
    
    
    fracs = [fracs; opto_state_frac; no_opto_state_frac];
    for i = 1:numel(opto_state_frac) * 2
        power{end+1} = modecounts(1).power_criteria{idxpower};
        period{end+1} = modecounts(1).period_criteria{idxperiod};
        animal{end+1} = modecounts(animalID).animal;
    end
    
    regions = cat(2, regions, regionsOrder);
    regions = cat(2, regions, regionsOrder);
    
    opto = [opto; ones(numel(opto_state_frac), 1); zeros(numel(opto_state_frac), 1)];
    
end

animal = animal';
power = power';
period = period';
regions = regions';

tbl = table(animal, power, period, opto, regions, fracs);
end

