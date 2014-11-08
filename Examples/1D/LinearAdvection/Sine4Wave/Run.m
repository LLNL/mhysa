% Script to run HyPar to simulate the following:
% Case: Fourth Power of Sine Wave
% Model: Linear Advection
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
    '/Examples/1D/LinearAdvection/Sine4Wave/aux/init.c ', ...
    '-lm -o INIT'];
system(cmd);
% find the HyPar binary
hypar = [hypar_path,'/bin/HyPar'];

% add the Matlab scripts directory in HyPar to path
path(path,strcat(hypar_path,'/Examples/Matlab/'));

% Get the default
[~,~,~,~,~,~,~,~,~,hyp_flux_split,hyp_int_type,par_type,par_scheme,~, ...
 cons_check,screen_op_iter,~,~, ip_type,input_mode, ...
 output_mode,n_io,~,~,nb,bctype,dim,face,limits, ...
 mapped,borges,yc,nl,eps,p,rc,xi,wtol,lutype,norm,maxiter,atol,rtol, ...
 verbose] = SetDefaults();

% set problem specific input parameters
ndims = 1;
nvars = 1;
iproc = 1;
ghost = 3;

% set grid size;
N = 160;

% specify spatial discretization scheme
hyp_scheme = 'crweno5';

% for spatial convergence, use very small time step
dt = 0.005;
t_final = 2.0;
niter = int32(t_final/dt);

% set physical model and related parameters
model = 'linear-advection-diffusion-reaction';
adv = 1.0;

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

% set solution output to text file
op_format = 'text';
op_overwrite = 'no';
file_op_iter = 2*niter;

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
WriteBoundaryInp(nb,bctype,dim,face,limits);
WritePhysicsInp_LinearADR(adv);
WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
% Generate initial solution
system(init_exec);
% Run the simulation
system(hypar_exec);

% read in the initial solution
data = load('op_00000.dat');
uinit = data(:,3);
% read in the final solution
data = load('op_00001.dat');
index = data(:,1);
xcoord = data(:,2);
u = data(:,3);

%clean up
system(clean_exec);
% open figure window
scrsz = get(0,'ScreenSize');
figSolution = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% set title string
title_str = ['Space: ',hyp_scheme,'; Time: ',ts,' (',tstype,'); ', ...
    'Grid Size: ',num2str(N,'%d'),'; Time step: ',num2str(dt,'%f')];

% plot the solution
plot(xcoord,uinit,'-k','linewidth',2,'MarkerSize',10);
hold on;
plot(xcoord,u,'-bo','linewidth',1,'MarkerSize',10);
hold on;
xlabel('x','FontName','Times','FontSize',20,'FontWeight','normal');
ylabel('u(x)','FontName','Times','FontSize',20,'FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');
legend('Initial','Final','Location','northeast');
title(title_str);
grid on;

% clean up
system('rm -rf INIT');

