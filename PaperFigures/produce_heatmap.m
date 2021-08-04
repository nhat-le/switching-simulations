function im = produce_heatmap(arr, xlst, ylst, varargin)

% Parse the varargin and default values
eval(evalargs(varargin));
if ~exist('newfig', 'var')
    newfig = 1;
end

if ~exist('x_label', 'var')
    x_label = '';
end

if ~exist('y_label', 'var')
    y_label = '';
end

if ~exist('ytickvalues', 'var')
    ytickvalues = 0:0.5:1.5;
end

if ~exist('legendname', 'var')
    legendname = '';
end

if ~exist('clim', 'var')
    clim = [0, 1];
end

if ~exist('font_size', 'var')
    font_size = 20;
end




%2d interpolation
if exist('xgridsize', 'var')
    [xgrid, ygrid] = meshgrid(xlst, ylst);

    xlstnew = linspace(min(xlst), max(xlst), xgridsize);
    ylstnew = linspace(min(ylst), max(ylst), ygridsize);
    [xgridnew, ygridnew] = meshgrid(xlstnew, ylstnew);
    arr = interp2(xgrid, ygrid, arr, xgridnew, ygridnew);
end

if newfig
    figure;
end
im = imshow(arr, 'InitialMagnification', 'fit', 'XData', [min(xlst), max(xlst)],...
    'YData', [min(ylst), max(ylst)]);

% Make colorbar
colormap hot

c = colorbar;
c.Label.String = legendname;
if newfig
    set(c, 'Position', [0.847 0.3433 0.04 0.5805], 'FontSize', 12, 'FontAngle', 'italic');
end

% if exist('cbar_pos', 'var')
%     set(c, 'Position', cbar_pos, 'FontSize', 12, 'FontAngle', 'italic');
% end


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

mymakeaxis('x_label', x_label, 'y_label', y_label, 'interpreter', 'latex',...
    'offsetRatio', 0, 'font_size', font_size, 'yticks', ytickvalues)
yticks(0:0.2:1.6);
yticklabels(0:0.2:1.6);

% Resize the figure
if newfig
    set(gcf, 'Position', [89 105 550 420]);
    set(gca, 'Position', [0 0.03 0.8 0.9]);
    set(c, 'Position', [0.8 0.2859 0.04 0.643], 'FontSize', 20);
end

end
