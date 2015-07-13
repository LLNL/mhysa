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
figure(1);
surf(xcoord,ycoord,b);

% read solution
data = load(fname);
h = reshape(data(:,5),imax,jmax);

% plot
figure(2);
contourf(xcoord,ycoord,h+b,30);
grid on;
colorbar;