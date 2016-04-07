clear all;
close all;

%% input filename
fname = input('Enter solution filename: ','s');
    
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

%% Extract the density perturbation and reshape it to a 3D array
u = reshape(data(:,7),imax,jmax,kmax);

%% open the figure window
figure(1);
set(1, 'Position', [100, 100,900, 800]);
figure(2);
set(2, 'Position', [100, 100,900, 800]);

%% Plot the 3D solution
figure(1);
p = patch(isosurface(y,x,z,u,-0.0025));
isonormals(y,x,z,u,p);
hold on;
set(p,'FaceColor','red','EdgeColor','none');
p = patch(isosurface(y,x,z,u,-0.002));
isonormals(y,x,z,u,p);
hold on;
set(p,'FaceColor','green','EdgeColor','none','FaceAlpha',0.2);
p = patch(isosurface(y,x,z,u,-0.0015));
isonormals(y,x,z,u,p);
set(p,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.1);
daspect([1 1 1]);
view(45,-10);
camlight;
lighting gouraud;
grid on;
xlabel('y');
ylabel('x');
zlabel('z');
axis([min(y),max(y),min(x),max(x),min(z),max(z)]);

%% Extract a mid-slice in the x-z plane and plot it
u_xz = squeeze(u(:,int32(jmax/2),:))';
% contour levels
v = -0.0034:(0.003/10):-0.0004;
figure(2);
contour(x,z,u_xz,v,'linewidth',2);
xlabel('x');
ylabel('z');
colorbar;
grid on;
axis([min(x),max(x),min(z),max(z)]);