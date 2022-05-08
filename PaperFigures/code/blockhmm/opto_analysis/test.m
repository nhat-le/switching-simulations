load shift.mat

shift.Shift = nominal(shift.Shift);
shift.Operator = nominal(shift.Operator);

lme = fitlme(shift, 'QCDev ~ Shift + (1|Operator)', ...
    'FitMethod', 'REML', 'DummyVarCoding', 'effects');

%%
load repeatedmeas.mat