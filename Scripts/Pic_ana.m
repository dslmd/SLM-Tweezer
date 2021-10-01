pic = imread("Camera.bmp");
a = pic(500:1100,2300:2900);
%% Surface Plot
% surf(a);
% shading interp
%% Heatmap Plot: Ax is used to hide the labels
cmap = colormap('hot');
heatmap(a,'Colormap',cmap);
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));