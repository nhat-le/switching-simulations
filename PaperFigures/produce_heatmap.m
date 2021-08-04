function im = produce_heatmap(arr, xlst, ylst, clim, plotname, varargin)

%2d interpolation
if length(varargin) == 2
    [xgrid, ygrid] = meshgrid(xlst, ylst);

    xlstnew = linspace(min(xlst), max(xlst), varargin{1});
    ylstnew = linspace(min(ylst), max(ylst), varargin{2});
    [xgridnew, ygridnew] = meshgrid(xlstnew, ylstnew);
    arr = interp2(xgrid, ygrid, arr, xgridnew, ygridnew);
end


% figure;
im = imshow(arr, 'InitialMagnification', 'fit', 'XData', [min(xlst), max(xlst)],...
    'YData', [min(ylst), max(ylst)]);

% Make colorbar
colormap hot

c = colorbar;
c.Label.String = plotname;
set(c, 'Position', [0.847 0.3433 0.04 0.5805], 'FontSize', 12, 'FontAngle', 'italic');
axis xy

if ~isempty(clim)
    caxis(clim)
end
% caxis([0.5, 1.0])
%xlabel('$K$', 'interpreter', 'latex')
%ylabel('$I$', 'interpreter', 'latex')

set(gca, 'YDir', 'normal');
set(gca, 'XDir', 'normal');
set(gca, 'Visible', 'on')
axis square

mymakeaxis('x_label', '$\epsilon$', 'y_label', '$\gamma$', 'interpreter', 'latex',...
    'offsetRatio', 0, 'font_size', 20, 'yticks', 0:0.5:1.5)
yticks(0:0.2:1.6);
yticklabels(0:0.2:1.6);

% Resize the figure
set(gcf, 'Position', [89 105 550 420]);
set(gca, 'Position', [0 0.03 0.8 0.9]);
set(c, 'Position', [0.8 0.2859 0.04 0.643], 'FontSize', 20);

end

