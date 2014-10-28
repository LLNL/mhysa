% Script to test the time convergence rate
% for a smooth solution of the linear advection 
% equation

clear all;
close all;

% remove all useless files
system('rm -rf *.dat *.inp *.log INIT');

fprintf('Spatial convergence test on a smooth solution ');
fprintf('to the linear advection equation.\n');

% Ask for path to HyPar source directory
path = input('Enter path to HyPar source: ','s');

% Compile the code to generate the initial solution
cmd = ['g++ ',path,'/Examples/1D/LinearAdvection/SineWave/aux/init.C ', ...
       '-o INIT'];
system(cmd);
% find the HyPar binary
hypar = [path,'/bin/HyPar'];

% Get the default
[~,~,~,~,~,niter,~,~,~, ...
    hyp_flux_split,hyp_int_type,par_type,par_scheme,~,cons_check, ...
    screen_op_iter,file_op_iter,~, ip_type,input_mode, ...
    output_mode,n_io,op_overwrite,~,nb,bctype,dim,face,limits, ...
    mapped,borges,yc,nl,eps,p,rc,xi,wtol,lutype,norm,maxiter,atol,rtol, ...
    verbose] = SetDefaults();

% set problem specific input parameters
ndims = 1;
nvars = 1;
iproc = 1;
ghost = 3;

% specify a nice, high-order spatial discretization scheme
hyp_scheme = 'weno5';

% for spatial convergence, use very small time step
dt = [0.0005 0.001 0.002 0.003 0.005 0.01 0.02];

% set physical model and related parameters
model = 'linear-advection-diffusion-reaction';
adv = 1.0;

% PETSc's time integration method
ts = 'arkimex';
tstype = '4';

petsc_flags = ' ';
% set PETSc time-integration flags (comment to turn off)
petsc_flags = [petsc_flags, '-use-petscts '];
petsc_flags = [petsc_flags, '-ts_type ',ts,' '];
petsc_flags = [petsc_flags, '-ts_',ts,'_type ',tstype,' '];
petsc_flags = [petsc_flags, '-ts_adapt_type none '];
petsc_flags = [petsc_flags, '-hyperbolic_implicit '];
petsc_flags = [petsc_flags, '-snes_type newtonls '];
petsc_flags = [petsc_flags, '-snes_rtol 1e-10 '];
petsc_flags = [petsc_flags, '-snes_atol 1e-10 '];
petsc_flags = [petsc_flags, '-ksp_type gmres '];
petsc_flags = [petsc_flags, '-ksp_rtol 1e-10 '];
petsc_flags = [petsc_flags, '-ksp_atol 1e-10 '];
petsc_flags = [petsc_flags, '-log_summary'];

% turn off solution output to file
op_format = 'none';

% set number of grid refinement levels
ref_levels = max(size(dt));

% set initial grid size;
N = 320;

% set the commands to run the executables
nproc = 1;
for i = 1:max(size(iproc))
    nproc = nproc * iproc(i);
end
init_exec = './INIT > init.log 2>&1';
clean_exec = 'rm -rf *.inp *.dat *.log';

% preallocate arrays for dx, error and wall times
Errors = zeros(ref_levels,3);
Walltimes = zeros(ref_levels,2);

% convergence analysis
for r = 1:ref_levels
    fprintf('\t%2d:  dt=%1.16e\n',r,dt(r));
    niter = 1000000;
    petscdt = [' -ts_dt ',num2str(dt(r),'%1.16e'),' '];
    petscft = ' -ts_final_time 1.0 ';
    % Write out the input files for HyPar
    WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,ts,tstype, ...
        hyp_scheme,hyp_flux_split,hyp_int_type,par_type,par_scheme, ...
        dt(r),cons_check,screen_op_iter,file_op_iter,op_format,ip_type, ...
        input_mode,output_mode,n_io,op_overwrite,model);
    WriteBoundaryInp(nb,bctype,dim,face,limits);
    WritePhysicsInp_LinearADR(adv);
    WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
    WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
    % Generate the initial and exact solution
    system(init_exec);
    system('ln -sf initial.inp exact.inp');
    % Run HyPar
    hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ',hypar, ...
                   ' ',petsc_flags,petscdt,petscft,' > run.log 2>&1'];
    system(hypar_exec);
    % Read in the errors
    [Errors(r,:),Walltimes(r,:)] = ReadErrorDat(ndims);
    % Clean up
    system(clean_exec);    
end

% open figure window
scrsz = get(0,'ScreenSize');
figErrDt   = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% set title string
title_str = ['Time convergence plot for ''',ts,' (',tstype,')'''];

% plot errors
loglog(dt,Errors(:,1),'-ko','linewidth',2,'MarkerSize',10);
hold on;
loglog(dt,Errors(:,2),'-ks','linewidth',2,'MarkerSize',10);
hold on;
loglog(dt,Errors(:,3),'-kd','linewidth',2,'MarkerSize',10);
hold on;
xlabel('dt','FontName','Times','FontSize',20,'FontWeight','normal');
ylabel('Error','FontName','Times','FontSize',20,'FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');
legend(['L1  ';'L2  ';'Linf'],'Location','northwest');
title(title_str);
grid on;

% clean up
system('rm -rf INIT');

