% Script to run HyPar to simulate the following:
% Case: Inertia Gravity Wave
% Model: 2D Navier-Stokes
% Copy this file ("Run.m") to an empty folder and 
% execute it in MATLAB.

clear all;
close all;

% remove all useless files
system('rm -rf *.dat *.inp *.log INIT');

% Ask for path to HyPar source directory
hypar_path = input('Enter path to HyPar source: ','s');

% Compile the code to generate the initial solution
cmd = ['gcc ',hypar_path, ...
    '/Examples/2D/NavierStokes2D/InertiaGravityWave/aux/init.c ', ...
    '-lm -o INIT'];
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
hyp_scheme = 'weno5';
hyp_int_type = 'components';
hyp_flux_split = 'no';
% parameters controlling the WENO-type schemes
mapped = 1;
borges = 0;
yc = 0;
nl = 0;
eps = 1e-6;

% time integration
dt = 0.5;
t_final = 3000.0;
niter = int32(t_final/dt);

% set physical model and related parameters
model = 'navierstokes2d';
gamma = 1.4;
grav  = [0.0 9.8];
upw   = 'rusanov';
rho_ref = 1.1612055171196529;
p_ref = 100000.0;
HB = 3;
BV = 0.01;
GasConst = 287.058;

% other options
cons_check='no';
screen_op_iter = 10;
op_format = 'text';
op_overwrite = 'yes';
file_op_iter = 240;
ip_type = 'binary';

% set time-integration scheme
ts = 'rk';
tstype = 'ssprk3';

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
bcdim = [0; 0; 1; 1;];
face = [1; -1; 1; -1];
limits = [0 0 0 10000; 0 0 0 10000; 0 300000 0 0; 0 300000 0 0];
WallVel1 = '0.0  0.0';
WallVel2 = '0.0  0.0';

% set the commands to run the executables
nproc = 1;
for i = 1:max(size(iproc))
    nproc = nproc * iproc(i);
end
init_exec = './INIT > init.log 2>&1';
hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ',hypar, ...
               ' ',petsc_flags];
clean_exec = 'rm -rf *.inp *.dat *.log';

% write the input files for HyPar
WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,ts,tstype, ...
    hyp_scheme,hyp_flux_split,hyp_int_type,par_type,par_scheme, ...
    dt,cons_check,screen_op_iter,file_op_iter,op_format,ip_type, ...
    input_mode,output_mode,n_io,op_overwrite,model);
WriteBoundaryInp(nb,bctype,bcdim,face,limits,WallVel2,WallVel2);
WritePhysicsInp_NavierStokes2D(gamma,grav,upw,rho_ref,p_ref,HB,BV,GasConst);
WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
% Generate initial solution
system(init_exec);
% Run the simulation
system(hypar_exec);

% set dimensions for postprocessing
imax = N(1);
jmax = N(2);
npoints = imax * jmax;

% get physical domain dimensions from any one solution file
if (strcmp(op_overwrite,'no'))
    data = load('op_00000.dat');
else
    data = load('op.dat');
end
xcoord  = reshape(data(:,3),imax,jmax);
ycoord  = reshape(data(:,4),imax,jmax);
xlen = max(xcoord(:,1)) - min(xcoord(:,1));
zlen = max(ycoord(1,:)) - min(ycoord(1,:));

% Get screen size
scrsz = get(0,'ScreenSize');
width = min(scrsz(3),1000);
height = width * max(10*zlen/xlen,0.1);

% figure out how many solution files are there
nfiles = niter / file_op_iter + 1;
% open figure window
scrsz = get(0,'ScreenSize');
figSolution = figure('Position',[1 scrsz(4)/2 width height]);
figCrossSec = figure('Position',[1 scrsz(4)/2 width 0.67*width]);

if (strcmp(op_overwrite,'no'))
    for n = 1:nfiles
        filename = ['op_',sprintf('%05d',n-1),'.dat'];
        if (~exist(filename,'file'))
            break;
        end
        fprintf('Plotting %s.\n',filename);
        % read in the solution
        data = load(filename);
        rho     = reshape(data(:,5),imax,jmax);
        rhou    = reshape(data(:,6),imax,jmax);
        rhov    = reshape(data(:,7),imax,jmax);
        energy  = reshape(data(:,8),imax,jmax);

        % reference values
        Cp = gamma/(gamma-1.0) * GasConst;
        T_ref = p_ref / (GasConst * rho_ref);
        Pexner = 1.0 + ...
            ((grav(2)*grav(2))/(Cp*T_ref*BV*BV))*(exp(-BV*BV*ycoord/grav(2))-1.0);
        P0 = p_ref * Pexner.^(gamma/(gamma-1.0));
        rho0 = rho_ref * Pexner.^(1.0/(gamma-10.0));
        theta0 = T_ref * exp(BV*BV*ycoord/grav(2));

        % primitive flow variables
        uvel = rhou ./ rho;
        vvel = rhov ./ rho;
        theta = (energy-0.5*rho.*(uvel.*uvel+vvel.*vvel))./(Pexner.*rho) ...
                * ((gamma-1.0)/GasConst);
        P = (gamma-1.0) * (energy - 0.5*rho.*(uvel.*uvel+vvel.*vvel));

        % set title string
        title_str = ['Space: ',hyp_scheme,'; Time: ',ts,' (',tstype,'); ', ...
            'Grid Size: ',num2str(N,'%d'),'; Time step: ',num2str(dt,'%f'), ...
            ' Time: ',num2str((n-1)*file_op_iter*dt,'%1.4e')];
        % plot the solution
        figure(figSolution);
        handle = contourf(xcoord,ycoord,theta-theta0,'LineColor','none', ...
                          'LevelList',-0.0015:0.0005:0.003);
        colorbar;
        xlabel('x','FontName','Times','FontSize',20,'FontWeight','normal');
        ylabel('y','FontName','Times','FontSize',20,'FontWeight','normal');
        set(gca,'FontSize',14,'FontName','Times');
        title(title_str);
        % plot the cross section
        figure(figCrossSec);
        plot(xcoord(:,1),(theta(:,jmax/2)-theta0(:,jmax/2)),'-ko','LineWidth',2, ...
             'MarkerSize',10);
        xlabel('x','FontName','Times','FontSize',20,'FontWeight','normal');
        ylabel('\Delta\theta','FontName','Times','FontSize',20,'FontWeight','normal');
        set(gca,'FontSize',14,'FontName','Times');
        axis([0 300000 0 0.003]);
        title(title_str);
    end
else
    filename = 'op.dat';
    if (~exist(filename,'file'))
        break;
    end
    fprintf('Plotting %s.\n',filename);
    % read in the solution
    data = load(filename);
    rho     = reshape(data(:,5),imax,jmax);
    rhou    = reshape(data(:,6),imax,jmax);
    rhov    = reshape(data(:,7),imax,jmax);
    energy  = reshape(data(:,8),imax,jmax);

    % reference values
    Cp = gamma/(gamma-1.0) * GasConst;
    T_ref = p_ref / (GasConst * rho_ref);
    Pexner = 1.0 + ...
        ((grav(2)*grav(2))/(Cp*T_ref*BV*BV))*(exp(-BV*BV*ycoord/grav(2))-1.0);
    P0 = p_ref * Pexner.^(gamma/(gamma-1.0));
    rho0 = rho_ref * Pexner.^(1.0/(gamma-10.0));
    theta0 = T_ref * exp(BV*BV*ycoord/grav(2));

    % primitive flow variables
    uvel = rhou ./ rho;
    vvel = rhov ./ rho;
    theta = (energy-0.5*rho.*(uvel.*uvel+vvel.*vvel))./(Pexner.*rho) ...
            * ((gamma-1.0)/GasConst);
    P = (gamma-1.0) * (energy - 0.5*rho.*(uvel.*uvel+vvel.*vvel));

    % set title string
    title_str = ['Space: ',hyp_scheme,'; Time: ',ts,' (',tstype,'); ', ...
        'Grid Size: ',num2str(N,'%d'),'; Time step: ',num2str(dt,'%f'), ...
        ' Time: ',num2str(t_final,'%1.4e')];
    % plot the solution
    figure(figSolution);
    handle = contourf(xcoord,ycoord,theta-theta0,'LineColor','none', ...
                      'LevelList',-0.0015:0.0005:0.003);
    colorbar;
    xlabel('x','FontName','Times','FontSize',20,'FontWeight','normal');
    ylabel('y','FontName','Times','FontSize',20,'FontWeight','normal');
    set(gca,'FontSize',14,'FontName','Times');
    title(title_str);
    % plot the cross section
    figure(figCrossSec);
    plot(xcoord(:,1),(theta(:,jmax/2)-theta0(:,jmax/2)),'-k.','LineWidth',2, ...
         'MarkerSize',10);
    xlabel('x','FontName','Times','FontSize',20,'FontWeight','normal');
    ylabel('\Delta\theta','FontName','Times','FontSize',20,'FontWeight','normal');
    set(gca,'FontSize',14,'FontName','Times');
    axis([0 300000 -0.002 0.003]);
    grid on;
    title(title_str);
end

%clean up
system(clean_exec);

% clean up
system('rm -rf INIT');
