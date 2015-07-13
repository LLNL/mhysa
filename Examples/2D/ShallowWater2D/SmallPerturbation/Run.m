% Script to run HyPar to simulate the following:
% Case: Small Perturbation
% Model: 2D shallow water equations
% Copy this file ("Run.m") to an empty folder and 
% execute it in MATLAB.

clear all;
close all;

% Ask for path to HyPar source directory
hypar_path = input('Enter path to HyPar source: ','s');

% Compile the code to generate the initial solution
cmd = ['gcc ',hypar_path, ...
    '/Examples/2D/ShallowWater2D/SmallPerturbation/aux/init.c ', ...
    '-lm -o INIT'];
system(cmd);
% find the HyPar binary
hypar = [hypar_path,'/bin/HyPar'];

% add the Matlab scripts directory in HyPar to path
path(path,strcat(hypar_path,'/Examples/Matlab/'));

% Get the default
[~,~,~,~,~,~,~,~,~, ...
 hyp_flux_split,~,par_type,par_scheme,~,~, ...
 ~,~,~,~,input_mode, ...
 output_mode,n_io,~,~,~,~,~,~,~, ...
 mapped,borges,yc,nl,eps,p,rc,xi,wtol,lutype,norm,maxiter,atol,rtol, ...
 verbose] = SetDefaults();

% set problem specific input parameters
ndims = 2;
nvars = 3;
iproc = [1 1];
ghost = 3;

% set grid size;
N = [201 101];

% specify spatial discretization scheme
hyp_scheme      = 'weno5';
hyp_int_type    = 'components';

% time integration
dt      = 0.002;
t_final = 0.6;
niter   = int32(t_final/dt);

% set physical model and related parameters
model   = 'shallow-water-2d';
grav    = 9.812;
bt_type = 0;
upw     = 'roe';

% other options
cons_check      = 'no';
screen_op_iter  = 1;
op_format       = 'text';
op_overwrite    = 'no';
file_op_iter    = niter/5;
ip_type         = 'binary';

% set time-integration scheme
ts      = 'rk';
tstype  = '44';
use_petsc = 0;

if (strcmp(ts,'arkimex'))
    hyp_flux_flag = '-hyperbolic_implicit'; % implicit treatment
else
    hyp_flux_flag = ' ';
end

petsc_flags = ' ';
if (use_petsc)
    % set PETSc time-integration
    petsc_flags = [petsc_flags, '-use-petscts '];
    petsc_flags = [petsc_flags, '-ts_type ',ts,' '];
    petsc_flags = [petsc_flags, '-ts_',ts,'_type ',tstype,' '];
    petsc_flags = [petsc_flags,' -ts_dt ',num2str(dt,'%1.16e'),' '];
    petsc_flags = [petsc_flags,' -ts_final_time ',num2str(t_final,'%f'),' '];
    petsc_flags = [petsc_flags, '-ts_adapt_type none '];
    if (strcmp(ts,'arkimex'))
        % set flags for implicit time integration
        petsc_flags = [petsc_flags, hyp_flux_flag,' '];
        petsc_flags = [petsc_flags, '-snes_type newtonls '];
        petsc_flags = [petsc_flags, '-snes_rtol 1e-10 '];
        petsc_flags = [petsc_flags, '-snes_atol 1e-10 '];
        petsc_flags = [petsc_flags, '-ksp_type gmres '];
        petsc_flags = [petsc_flags, '-ksp_rtol 1e-10 '];
        petsc_flags = [petsc_flags, '-ksp_atol 1e-10 '];
    end
    petsc_flags = [petsc_flags, '-log_summary'];
end

% set boundaries
nb = 4;
bctype = ['extrapolate'; ...
          'extrapolate'; ...
          'extrapolate'; ...
          'extrapolate'];
bcdim     = [0; 0; 1; 1;];
face      = [1; -1; 1; -1];
limits    = [0 0 0 1.0; 0 0 0 1.0; 0 2.0 0 0; 0 2.0 0 0];

% set the commands to run the executables
nproc = 1;
for i = 1:max(size(iproc))
    nproc = nproc * iproc(i);
end
init_exec  = './INIT > init.log 2>&1';
hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ',hypar, ...
               ' ',petsc_flags];

% write the input files for HyPar
WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,ts,tstype, ...
    hyp_scheme,hyp_flux_split,hyp_int_type,par_type,par_scheme, ...
    dt,cons_check,screen_op_iter,file_op_iter,op_format,ip_type, ...
    input_mode,output_mode,n_io,op_overwrite,model);
WriteBoundaryInp(nb,bctype,bcdim,face,limits);
WritePhysicsInp_ShallowWater2D(grav,bt_type,upw);
WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
% Generate initial solution
system(init_exec);
% Run the simulation
system(hypar_exec);

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

for i=0:5
    % read solution
    fname = strcat('op_0000',sprintf('%1d',i),'.dat');
    fprintf('Plotting %s.\n',fname);
    data = load(fname);
    h = reshape(data(:,5),imax,jmax);

    % plot
    figure(i+2);
    contourf(xcoord,ycoord,h+b,30);
    grid on;
    colorbar;
end

% clean up
system('rm -rf *.dat *.inp *.log *.bin INIT');
