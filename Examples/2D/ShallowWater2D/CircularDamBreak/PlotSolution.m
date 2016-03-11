clear all;
close all;

fname = input('Enter name of file to plot: ','s');

% Load topography file
topo = load('topography_00000.dat');

% find out max grid sizes
imax = max(topo(:,1)) + 1;
jmax = max(topo(:,2)) + 1;

%spatial coordinates
xcoord = reshape(topo(:,3),imax,jmax);
ycoord = reshape(topo(:,4),imax,jmax);

% bottom topography
b = reshape(topo(:,5),imax,jmax);
% figure(1);
% surf(xcoord,ycoord,b);
% hold on;

% read solution
data = load(fname);
h = reshape(data(:,5),imax,jmax);
u = reshape(data(:,6),imax,jmax)./h;
v = reshape(data(:,7),imax,jmax)./h;


% plot
figure(1);
set(1, 'Position', [100, 100, 1300, 500]);
subplot(1,2,1);
surf(xcoord,ycoord,h+b);
colorbar;
subplot(1,2,2);
contour(xcoord,ycoord,h+b,30);
colorbar;
% figure(3);
% quiver(xcoord,ycoord,u,v);