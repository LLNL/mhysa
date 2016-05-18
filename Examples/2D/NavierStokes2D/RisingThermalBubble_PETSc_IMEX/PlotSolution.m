% Script to plot the solution of the rising thermal bubble
% problem. 

% Compile and run aux/PostProcess.c after simulation is
% completed and choose text output.

clear all;
close all;

flag = 1;
if (exist('op.dat','file'))
    flag = 0;
elseif (exist('op_00000.dat','file'))
    flag = 1;
else
  fprintf('Error: No solution files exist.\n');
  return;
end

% get physical domain dimensions from any one solution file
if (flag)
    data = load('op_00000.dat');
else
    data = load('op.dat');
end
imax    = max(data(:,1)) + 1;
jmax    = max(data(:,2)) + 1;
xcoord  = reshape(data(:,3),imax,jmax);
ycoord  = reshape(data(:,4),imax,jmax);
xlen    = max(xcoord(:,1)) - min(xcoord(:,1));
zlen    = max(ycoord(1,:)) - min(ycoord(1,:));

% Get screen size
scrsz = get(0,'ScreenSize');
width = max(min(scrsz(3),1000),640);
height = 0.9 * width;
% open figure window
scrsz = get(0,'ScreenSize');
figSolution = figure('Position',[1 scrsz(4)/2 width height]);

maxfiles = 1000;
if (flag)
    for n = 1:maxfiles
        filename = ['op_',sprintf('%05d',n-1),'.dat'];
        if (~exist(filename,'file'))
            fprintf('No more files found.\n');
            break;
        end
        fprintf('Plotting %s.\n',filename);
        % read in the solution
        data    = load(filename);
        rho     = reshape(data(:, 5),imax,jmax);
        uvel    = reshape(data(:, 6),imax,jmax);
        vvel    = reshape(data(:, 7),imax,jmax);
        P       = reshape(data(:, 8),imax,jmax);
        theta   = reshape(data(:, 9),imax,jmax);
        rho0    = reshape(data(:,10),imax,jmax);
        P0      = reshape(data(:,11),imax,jmax);
        Pexner  = reshape(data(:,12),imax,jmax);
        theta0  = reshape(data(:,13),imax,jmax);

        % plot the solution
        figure(figSolution);
        handle = contourf(xcoord,ycoord,theta-theta0,'LineColor','none', ...
                          'LevelList',-0.05:0.005:0.5);
        colorbar;
        xlabel('x','FontName','Times','FontSize',20,'FontWeight','normal');
        ylabel('y','FontName','Times','FontSize',20,'FontWeight','normal');
        set(gca,'FontSize',14,'FontName','Times');
        % write plot to files
        print(figSolution,'-depsc',['Contour_',sprintf('%05d',n-1),'.eps']);
    end
else
    filename = 'op.dat';
    if (~exist(filename,'file'))
        break;
    end
    fprintf('Plotting %s.\n',filename);
    % read in the solution
    data    = load(filename);
    rho     = reshape(data(:, 5),imax,jmax);
    uvel    = reshape(data(:, 6),imax,jmax);
    vvel    = reshape(data(:, 7),imax,jmax);
    P       = reshape(data(:, 8),imax,jmax);
    theta   = reshape(data(:, 9),imax,jmax);
    rho0    = reshape(data(:,10),imax,jmax);
    P0      = reshape(data(:,11),imax,jmax);
    Pexner  = reshape(data(:,12),imax,jmax);
    theta0  = reshape(data(:,13),imax,jmax);
    
    % plot the solution
    figure(figSolution);
    handle = contourf(xcoord,ycoord,theta-theta0,'LineColor','none', ...
                      'LevelList',-0.05:0.005:0.5);
    colorbar;
    xlabel('x','FontName','Times','FontSize',20,'FontWeight','normal');
    ylabel('y','FontName','Times','FontSize',20,'FontWeight','normal');
    set(gca,'FontSize',14,'FontName','Times');
    % write plot to files
    print(figSolution,'-depsc',['Contour.eps']);
end
