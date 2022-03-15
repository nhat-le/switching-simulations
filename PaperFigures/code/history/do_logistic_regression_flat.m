function B_all = do_logistic_regression_flat(choicearr, outcomearr, N)
% choicearr: array of choices, size Nblocks x T
% outcomearr: array of outcomes, size Nblocks x T
% N: how many trials back to go when doing logistic regression

[Xmat, y] = make_Xy_flat(choicearr * 2 - 1, outcomearr * 2 - 1, N);

% logistic regression
B = mnrfit(Xmat, (y > 0) + 1);
B_all = B(2:end); %first coef is intercept

end