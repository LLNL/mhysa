% Script to run HyPar to simulate the following:
% Case: Density Sine Wave
% Model: 1D Euler
% Copy this file ("Run.m") to an empty folder and 
% execute it in MATLAB.

clear all;
close all;

% remove all useless files
system('rm -rf *.dat *.inp *.log INIT');

% Ask for path to HyPar source directory
hypar_path = input('Enter path to HyPar source: ','s');

% Compile the code to generate the initial and exact solution
cmd = ['gcc ',hypar_path, ...
    '/Examples/1D/Euler1D/DensitySineWave/aux/exact.c -lm ', ...
    '-o EXACT'];
system(cmd);
% find the HyPar binary
hypar = [hypar_path,'/bin/HyPar'];

% add the Matlab scripts directory in HyPar to path
path(path,strcat(hypar_path,'/Examples/Matlab/'));

% Get the default
[~,~,~,~,~,~,~,~,~,~,~,par_type,par_scheme,~, ...
 cons_check,screen_op_iter,~,~, ~,input_mode, ...
 output_mode,n_io,~,~,~,~,~,~,~, ...
 mapped,borges,yc,nl,eps,p,rc,xi,wtol,lutype,norm,maxiter,atol,rtol, ...
 verbose] = SetDefaults();

% set problem specific input parameters
ndims = 1;
nvars = 3;
iproc = 1;
ghost = 3;

% set grid size;
N = 40;

% specify spatial discretization scheme
hyp_scheme = 'crweno5';
hyp_int_type = 'components';

% specify dt, final time
dt = 0.01;
t_final = 1.0;
niter = int32(t_final/dt);

% set physical model and related parameters
model = 'euler1d';
gamma = 1.4;
grav  = 0.0;
upw   = 'roe';

% set time-integration scheme
ts = 'arkimex';
tstype = '3';
use_petsc = 1; % 1 - specified time-integration method is from PETSc
               % 0 - specified time-integration method is native

if (strcmp(ts,'arkimex'))
%-------------------------------------------------------------------------%
%   if 'arkimex' time-integrator is used, the following options
%   can be chosen from:-

%   use split hyperbolic flux form (acoustic and entropy modes)?
    hyp_flux_split = 'yes';
%   treat acoustic waves implicitly, and entropy waves explicitly
    hyp_flux_flag  = '-hyperbolic_f_explicit -hyperbolic_df_implicit';
%   treat acoustic and entropy waves implicitly
%     hyp_flux_flag  = '-hyperbolic_f_implicit -hyperbolic_df_implicit';
%   treat acoustic and entropy waves explicitly
%     hyp_flux_flag  = '-hyperbolic_f_explicit -hyperbolic_df_explicit';

%   or no splitting?
%     hyp_flux_split = 'no';
%     hyp_flux_flag = '-hyperbolic_implicit'; % implicit treatment
%     hyp_flux_flag = '-hyperbolic_explicit'; % explicit treatment
%------------------------------------------------------------------------%
else
    hyp_flux_split = 'no';
    hyp_flux_flag = ' ';
end

petsc_flags = ' ';
if (use_petsc)
    % set PETSc time-integration flags 
    petsc_flags = [petsc_flags, '-use-petscts '];
    petsc_flags = [petsc_flags, '-ts_type ',ts,' '];
    petsc_flags = [petsc_flags, '-ts_',ts,'_type ',tstype,' '];
    petsc_flags = [petsc_flags,' -ts_dt ',num2str(dt,'%1.16e'),' '];
    petsc_flags = [petsc_flags,' -ts_final_time ',num2str(t_final,'%f'),' '];
    petsc_flags = [petsc_flags, '-ts_adapt_type none '];
    if (strcmp(ts,'arkimex'))
        % if using ARKIMEX time-integrator, some additional options
        petsc_flags = [petsc_flags, hyp_flux_flag, ' '];
        petsc_flags = [petsc_flags, '-snes_type newtonls '];
        petsc_flags = [petsc_flags, '-snes_rtol 1e-10 '];
        petsc_flags = [petsc_flags, '-snes_atol 1e-10 '];
        petsc_flags = [petsc_flags, '-ksp_type gmres '];
        petsc_flags = [petsc_flags, '-ksp_rtol 1e-10 '];
        petsc_flags = [petsc_flags, '-ksp_atol 1e-10 '];
    end
    % petsc_flags = [petsc_flags, '-log_summary'];
end

% set solution output to text file
op_format = 'text';
op_overwrite = 'no';
file_op_iter = 2*niter;
% set initial solution file type to ascii
ip_type = 'ascii';

% set boundaries
nb = 2;
bctype = ['periodic'; 'periodic'];
bcdim = [0; 0];
face = [1; -1];
limits = [0 0; 0 0];

% set the commands to run the executables
nproc = 1;
for i = 1:max(size(iproc))
    nproc = nproc * iproc(i);
end
exact_exec = './EXACT > exact.log 2>&1';
hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ',hypar, ...
               ' ',petsc_flags];
clean_exec = 'rm -rf *.inp *.dat *.log';

% write the input files for HyPar
WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,ts,tstype, ...
    hyp_scheme,hyp_flux_split,hyp_int_type,par_type,par_scheme, ...
    dt,cons_check,screen_op_iter,file_op_iter,op_format,ip_type, ...
    input_mode,output_mode,n_io,op_overwrite,model);
WriteBoundaryInp(nb,bctype,bcdim,face,limits);
WritePhysicsInp_Euler1D(gamma,grav,upw);
WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
% Generate initial solution
system(exact_exec);
% Run the simulation
system(hypar_exec);

% read in the initial solution
data = load('op_00000.dat');
rhoinit = data(:,3);
% read in the final solution
data = load('op_00001.dat');
index = data(:,1);
xcoord = data(:,2);
rho = data(:,3);

%clean up
system(clean_exec);
% open figure window
scrsz = get(0,'ScreenSize');
figSolution = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% set title string
title_str = ['Space: ',hyp_scheme,'; Time: ',ts,' (',tstype,'); ', ...
    'Grid Size: ',num2str(N,'%d'),'; Time step: ',num2str(dt,'%f')];

% plot the solution
plot(xcoord,rhoinit,'-k','linewidth',2);
hold on;
plot(xcoord,rho,'-bo','linewidth',1,'MarkerSize',8);
hold on;
xlabel('x','FontName','Times','FontSize',20,'FontWeight','normal');
ylabel('Density','FontName','Times','FontSize',20,'FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');
title(title_str);
legend('Initial','Final','Location','Best');
grid on;

% clean up
system('rm -rf EXACT');
