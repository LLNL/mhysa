clear all;
close all;

for i=0:10000
    index = sprintf('%05d',i);
    fname = strcat('op_',index,'.dat');
    if (~exist(fname,'file'))
        break;
    end
    fprintf('Plotting file %s.\n',fname);
    
    %% Load the data
    data = load(fname); 
    % data is a 2D array with the following columns: i,j,k,x,y,z,u

    %% Get the dimensions
    imax = max(data(:,1)) + 1;
    jmax = max(data(:,2)) + 1;
    kmax = max(data(:,3)) + 1;

    %% Extract spatial coordinates
    x = data(1:imax,4);
    y = data(1:imax:(imax*jmax),5);
    z = data(1:(imax*jmax):(imax*jmax*kmax),6);

    %% Extract the solution and reshape it to a 3D array
    u = reshape(data(:,7),imax,jmax,kmax);
    
    %% open the figure window
    figure(i+1);
    set(i+1, 'Position', [100, 100, 600, 800]);
    

    %% Plot the 3D solution
    subplot(3,1,1);
    p = patch(isosurface(y,x,z,u,0.5));
    isonormals(y,x,z,u,p);
    set(p,'FaceColor','red','EdgeColor','none');
    daspect([1 1 1]);
    view(3);
    camlight;
    lighting gouraud;
    grid on;
    xlabel('y');
    ylabel('x');
    zlabel('z');
    axis([min(y),max(y),min(x),max(x),min(z),max(z)]);

    %% Extract a mid-slice in the x-y plane and plot it
    u_xy = u(:,:,kmax/2)';
    subplot(3,1,2);
    contourf(x,y,u_xy);
    xlabel('x');
    ylabel('y');
    colorbar;
    grid on;
    axis([min(x),max(x),min(y),max(y)]);

    %% Extract a 1D line along x through the middle and plot it
    u_x = u(:,jmax/2,kmax/2);
    subplot(3,1,3);
    plot(x,u_x,'-bo','Linewidth',2)
    xlabel('x');
    ylabel('u');
    legend('u(x,y=0,z=0)');
    grid on;
    axis([min(x),max(x),-0.1,1.1]);
end