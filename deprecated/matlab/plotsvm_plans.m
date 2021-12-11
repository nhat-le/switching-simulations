load svmresults.mat

idx = 5;
idy = 5;

coef0_single = squeeze(coefs0(idx, idy, :, :));
coef1_single = squeeze(coefs1(idx, idy, :, :));
coef2_single = squeeze(coefs2(idx, idy, :, :));
coef3_single = squeeze(coefs3(idx, idy, :, :));

figure;
subplot(221)
imagesc(coef0_single)
caxis([-2 2])
colorbar

subplot(222)
imagesc(coef1_single)
caxis([-2 2])
colorbar

subplot(223)
imagesc(coef2_single)
caxis([-2 2])
colorbar

subplot(224)
imagesc(coef3_single)
caxis([-2 2])
colormap redblue
colorbar