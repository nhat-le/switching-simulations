function B_all = do_logistic_regression(block, N)
[choicelst1, outcomelst1] = unfold_block(block);
[Xmat1, y1] = make_Xy(choicelst1, outcomelst1, N, true);

% logistic regression
B = mnrfit(Xmat1, (y1 > 0) + 1);
B_all = B(2:end); %first coef is intercept
end