clear all;
close all;

fname_indices = input('Enter range of file indices to plot: ');

% Load initial solution file
data = load('op_00000.dat');

% find out max grid sizes
imax = max(data(:,1)) + 1;
jmax = max(data(:,2)) + 1;

%spatial coordinates
x = reshape(data(:,3),imax,jmax);
y = reshape(data(:,4),imax,jmax);

count = 0;
for i = fname_indices
    fname = strcat('op_',sprintf('%05d',i),'.dat');
    % read solution
    data = load(fname);
    h = reshape(data(:,5),imax,jmax);
    u = reshape(data(:,6),imax,jmax)./h;
    v = reshape(data(:,7),imax,jmax)./h;

    % plot
    scrsz = get(0,'screensize');
    width = 0.9*scrsz(3);
    figure(1);
    set(1,'Position',[0 0 width 0.8*width]);
    surf(x,y,h);
    view(100,60);
    colorbar;
    
    filename = ['fig_',sprintf('%05d',count),'.png'];
    print(1,filename,'-dpng');
    
    count = count+1;
end