% Script to run HyPar to simulate the following:
% Case: Inertia Gravity Wave
% Model: 2D Navier-Stokes
% Copy this file ("Run.m") to an empty folder and 
% execute it in MATLAB.

clear all;
close all;

% remove all useless files
system('rm -rf *.dat *.inp *.log *.bin *.eps INIT PP');

% Ask for path to HyPar source directory
hypar_path = input('Enter path to HyPar source: ','s');

% Compile the code to generate the initial solution
cmd = ['gcc ',hypar_path, ...
    '/Examples/2D/NavierStokes2D/InertiaGravityWave/aux/init.c ', ...
    '-lm -o INIT'];
system(cmd);
% Compile the code to postprocess the output
cmd = ['gcc ',hypar_path, ...
    '/Examples/2D/NavierStokes2D/InertiaGravityWave/aux/PostProcess.c ', ...
    '-lm -o PP'];
system(cmd);
% find the HyPar binary
hypar = [hypar_path,'/bin/HyPar'];

% add the Matlab scripts directory in HyPar to path
path(path,strcat(hypar_path,'/Examples/Matlab/'));

% Get the default
[~,~,~,~,~,~,~,~,~,~,~,~,par_scheme,~, ...
 ~,~,~,~, ~,input_mode, ...
 output_mode,n_io,~,~,~,~,~,~,~, ...
 ~,~,~,~,~,p,rc,xi,wtol,lutype,norm,maxiter,atol,rtol, ...
 verbose] = SetDefaults();
% for NavierStokes2D, par_type must be 2-stage nonconservative
par_type = 'nonconservative-2stage';

% set problem specific input parameters
ndims = 2;
nvars = 4;
iproc = [12 1];
ghost = 3;

% set grid size;
N = [1200 50];

% specify spatial discretization scheme
hyp_scheme      = 'weno5';
hyp_int_type    = 'components';
hyp_flux_split  = 'no';
% parameters controlling the WENO-type schemes
mapped  = 0;
borges  = 0;
yc      = 1;
nl      = 0;
eps     = 1e-6;

% time integration
dt      = 0.5;
t_final = 3000.0;
niter   = int32(t_final/dt);

% set physical model and related parameters
model     = 'navierstokes2d';
gamma     = 1.4;                  % specific heat ratio
upw       = 'rusanov';            % choice of upwinding scheme
Prandtl   = 0.72;                 % Prandtl number
Reynolds  = -1.0;                 % Inviscid flow
Minf      = 1.0;                  % reference Mach number
grav      = [0.0 9.8];            % gravitational force vector
rho_ref   = 1.1612055171196529;   % reference altitude density
p_ref     = 100000.0;             % reference altitude pressure
HB        = 3;                    % type of hydrostatic balance
BV        = 0.01;                 % Brunt-Vaisala frequency
GasConst  = 287.058;              % Universal gas constant

% other options
cons_check      = 'no';
screen_op_iter  = 10;
op_format       = 'binary';
op_overwrite    = 'yes';
file_op_iter    = 240;
ip_type         = 'binary';

% set time-integration scheme
ts      = 'rk';
tstype  = 'ssprk3';

petsc_flags = ' ';
% set PETSc time-integration flags (comment to turn off)
% petsc_flags = [petsc_flags, '-use-petscts '];
% petsc_flags = [petsc_flags, '-ts_type ',ts,' '];
% petsc_flags = [petsc_flags, '-ts_',ts,'_type ',tstype,' '];
% petsc_flags = [petsc_flags,' -ts_dt ',num2str(dt,'%1.16e'),' '];
% petsc_flags = [petsc_flags,' -ts_final_time ',num2str(t_final,'%f'),' '];
% petsc_flags = [petsc_flags, '-ts_adapt_type none '];
% petsc_flags = [petsc_flags, '-hyperbolic_implicit '];
% petsc_flags = [petsc_flags, '-snes_type newtonls '];
% petsc_flags = [petsc_flags, '-snes_rtol 1e-10 '];
% petsc_flags = [petsc_flags, '-snes_atol 1e-10 '];
% petsc_flags = [petsc_flags, '-ksp_type gmres '];
% petsc_flags = [petsc_flags, '-ksp_rtol 1e-10 '];
% petsc_flags = [petsc_flags, '-ksp_atol 1e-10 '];
% petsc_flags = [petsc_flags, '-log_summary'];

% set boundaries
nb = 4;
bctype = ['periodic '; ...
          'periodic '; ...
          'slip-wall'; ...
          'slip-wall'];
bcdim     = [0; 0; 1; 1;];
face      = [1; -1; 1; -1];
limits    = [0 0 0 10000; 0 0 0 10000; 0 300000 0 0; 0 300000 0 0];
WallVel1  = '0.0  0.0';
WallVel2  = '0.0  0.0';

% set the commands to run the executables
nproc = 1;
for i = 1:max(size(iproc))
    nproc = nproc * iproc(i);
end
init_exec = './INIT > init.log 2>&1';
hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ',hypar, ...
               ' ',petsc_flags];

% write the input files for HyPar
WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,ts,tstype, ...
    hyp_scheme,hyp_flux_split,hyp_int_type,par_type,par_scheme, ...
    dt,cons_check,screen_op_iter,file_op_iter,op_format,ip_type, ...
    input_mode,output_mode,n_io,op_overwrite,model);
WriteBoundaryInp(nb,bctype,bcdim,face,limits,WallVel1,WallVel2);
WritePhysicsInp_NavierStokes2D(gamma,upw,Prandtl,Reynolds,Minf,grav, ...
                               rho_ref,p_ref,HB,BV,GasConst);
WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
% Generate initial solution
system(init_exec);
% Run the simulation
system(hypar_exec);

% Postprocess the output files
fprintf('Simulation finished, postprocessing solution output.\n');
fid = fopen('pp.inp','w');
fprintf(fid,'0\n');
fclose(fid);
system('./PP < pp.inp 2>&1 > pp.log');

if (strcmp(op_overwrite,'no'))
    flag = 1;
else
    flag = 0;
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
height = width * max(10*zlen/xlen,0.1);

% open figure window
scrsz = get(0,'ScreenSize');
figSolution = figure('Position',[1 scrsz(4)/2 width height]);
figCrossSec = figure('Position',[1 scrsz(4)/2 width 0.67*width]);

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
                          'LevelList',-0.0015:0.0005:0.003);
        colorbar;
        xlabel('x','FontName','Times','FontSize',20,'FontWeight','normal');
        ylabel('y','FontName','Times','FontSize',20,'FontWeight','normal');
        set(gca,'FontSize',14,'FontName','Times');
        % plot the cross section
        figure(figCrossSec);
        plot(xcoord(:,1),(theta(:,jmax/2)-theta0(:,jmax/2)),'-ko','LineWidth',2, ...
             'MarkerSize',10);
        xlabel('x','FontName','Times','FontSize',20,'FontWeight','normal');
        ylabel('\Delta\theta','FontName','Times','FontSize',20,'FontWeight','normal');
        set(gca,'FontSize',14,'FontName','Times');
        axis([0 300000 0 0.003]);

        % write plot to files
        print(figSolution,'-depsc',['Contour_',sprintf('%05d',n-1),'.eps']);
        print(figCrossSec,'-depsc',['CrossSc_',sprintf('%05d',n-1),'.eps']);
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
                      'LevelList',-0.0015:0.0005:0.003);
    colorbar;
    xlabel('x','FontName','Times','FontSize',20,'FontWeight','normal');
    ylabel('y','FontName','Times','FontSize',20,'FontWeight','normal');
    set(gca,'FontSize',14,'FontName','Times');
    % plot the cross section
    figure(figCrossSec);
    plot(xcoord(:,1),(theta(:,jmax/2)-theta0(:,jmax/2)),'-k.','LineWidth',2, ...
         'MarkerSize',10);
    xlabel('x','FontName','Times','FontSize',20,'FontWeight','normal');
    ylabel('\Delta\theta','FontName','Times','FontSize',20,'FontWeight','normal');
    set(gca,'FontSize',14,'FontName','Times');
    axis([0 300000 -0.002 0.003]);
    grid on;

    % write plot to files
    print(figSolution,'-depsc','Contour.eps');
    print(figCrossSec,'-depsc','CrossSc.eps');
end

% create the directory to save the solution and plots in
timestamp = clock;
dumpname = strcat(sprintf('%04d',timestamp(1)),'_',...
                  sprintf('%02d',timestamp(2)),'_',...
                  sprintf('%02d',timestamp(3)),'_',...
                  sprintf('%02d',timestamp(4)),'_',...
                  sprintf('%02d',timestamp(5)));
if (~exist(dumpname,'file'))
      system(['rm -rf ',dumpname]);
end
mkdir(dumpname);
system(['mv *.inp op* *.eps ',dumpname,'/']);

% clean up
system('rm -rf *.dat *.inp *.log *.bin *.eps INIT PP');
