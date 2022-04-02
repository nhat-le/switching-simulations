% For visualizing the HMM compositions according to the different conditions
load opto_hmm_info_030422.mat

animal = 'f26'; %f26/f27/f29/f32
area = 'Motor'; %Visual/Motor/Rsc/Frontal (note case-sensitive)
power = 'Low'; %Low/High
period = 'Outcome'; %Choice/Outcome


animal = lower(animal);


id = find(contains(animalModeInfo.animals, animal));
areas = animalinfo(id).areas;
powers = animalinfo(id).power;
periods = animalinfo(id).period;
tbl = table(string(areas)', string(powers)', string(periods)');

sess_filt = contains(animalinfo(id).areas , area) & contains(animalinfo(id).period, period) &...
    contains(animalinfo(id).power , power) & ~contains(animalinfo(id).areas, '+');

sess_filt = find(sess_filt);
opto_all = {};
zopto_all = {};
zno_opto_all = {};

for i = 1:numel(sess_filt)
    sess = sess_filt(i);
    zsingle = animalinfo(id).zclassified(sess);
    compsingle = animalinfo(id).composition(sess,:);
    optosingle = animalinfo(id).opto(sess);
    blocklens_single = animalinfo(id).block_lens(sess);
    
    zarr = zsingle{1};
    opto_arr = optosingle{1};
    blocklens = blocklens_single{1}';
    
    valid = ~isnan(opto_arr) & (blocklens >= 15);
    
    zarr = zarr(valid);
    opto_arr = opto_arr(valid);
    
    
    %split into opto/non-opto
    z_opto = zarr(logical(opto_arr));
    z_nonopto = zarr(~logical(opto_arr));
    
    z_all{i} = zarr;
    zopto_all{i} = z_opto;
    zno_opto_all{i} = z_nonopto;
    
end

%%
K = 6;
opto_frac = get_frac_array(zopto_all, K);
nonopto_frac = get_frac_array(zno_opto_all, K);



%% Plot the fractions and associated session
figure('Position', [440,423,736,375], 'Name', 'Fraction');
subplot(121)
h = bar(1:size(opto_frac, 1), opto_frac,'stacked');
cols = paperaesthetics;
colors = cols.colors;

for i = 1:K
    h(i).FaceColor = colors(i,:);
    h(i).ShowBaseLine = 'off';
end
ylim([0 1])

labels = {};
for i = 1:numel(sess_filt)
    labels{i} = num2str(sess_filt(i));
end


mymakeaxis('x_label', 'Session #', 'y_label', 'Fraction', 'xticks', 1:numel(sess_filt),...
    'xticklabels', labels, 'xytitle', 'Opto')


subplot(122)
h = bar(1:size(nonopto_frac, 1), nonopto_frac,'stacked');
cols = paperaesthetics;
colors = cols.colors;

for i = 1:K
    h(i).FaceColor = colors(i,:);
    h(i).ShowBaseLine = 'off';
end
ylim([0 1])
mymakeaxis('x_label', 'Session #', 'y_label', 'Fraction', 'xticks', 1:numel(sess_filt),...
    'xticklabels', labels, 'xytitle', 'No-opto')



function frac = get_frac_array(arr, N)

frac_cell = cellfun(@(x) parse_fractions(x, N), arr, 'UniformOutput', false);
frac = cell2mat(frac_cell');

assert(abs(min(sum(frac, 2)) - 1) < eps);

end



function frac = parse_fractions(arr, N)
frac = [];

for i = 1:N
    frac(i) = sum(arr == i) / numel(arr);
end

end
