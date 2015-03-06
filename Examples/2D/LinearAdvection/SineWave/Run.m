% Script to run HyPar to simulate the following:
% Case: Advection of a Sine Wave
% Model: 2D Linear Avection
% Copy this file ("Run.m") to an empty folder and 
% execute it in MATLAB.

clear all;
close all;

% remove all useless files
system('rm -rf *.dat *.inp *.log *.bin *.eps EXACT');

% Ask for path to HyPar source directory
hypar_path = input('Enter path to HyPar source: ','s');

% Compile the code to generate the initial and exact solutions
cmd = ['g++ ',hypar_path, ...
    '/Examples/2D/LinearAdvection/SineWave/aux/exact.c ', ...
    '-o EXACT'];
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
par_type = 'nonconservative-2stage';

% set problem specific input parameters
ndims = 2;
nvars = 1;
iproc = [1 1];
ghost = 3;

% set grid size;
N = [20 20];

% specify spatial discretization scheme
hyp_scheme      = 'weno5';
hyp_int_type    = 'components';
% parameters controlling the WENO-type schemes
mapped  = 0;
borges  = 0;
yc      = 1;
nl      = 0;
eps     = 1e-6;

% time integration
dt      = 0.04;
t_final = 1.0;
niter   = int32(t_final/dt);

% set physical model and related parameters
model     = 'linear-advection-diffusion-reaction';
advection = [1.0 1.0];

% other options
cons_check      = 'yes';
screen_op_iter  = 1;
op_format       = 'text';
op_overwrite    = 'no';
file_op_iter    = int32((t_final/10.0)/dt);
ip_type         = 'ascii';

% set time-integration scheme
ts      = 'rk';
tstype  = 'ssprk3';
use_petsc = 0;

hyp_flux_split = 'no'; % this model has no splitting defined
if (strcmp(ts,'arkimex'))
%     hyp_flux_flag = '-hyperbolic_implicit'; % implicit treatment
%     hyp_flux_flag = '-hyperbolic_explicit'; % explicit treatment
else
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
        petsc_flags = [petsc_flags, hyp_flux_flag, ' '];
        petsc_flags = [petsc_flags, '-snes_type newtonls '];
        petsc_flags = [petsc_flags, '-snes_rtol 1e-10 '];
        petsc_flags = [petsc_flags, '-snes_atol 1e-10 '];
        petsc_flags = [petsc_flags, '-ksp_type gmres '];
        petsc_flags = [petsc_flags, '-ksp_rtol 1e-10 '];
        petsc_flags = [petsc_flags, '-ksp_atol 1e-10 '];
    end
%     petsc_flags = [petsc_flags, '-log_summary'];
end

% set boundaries
nb = 4;
bctype = ['periodic'; ...
          'periodic'; ...
          'periodic'; ...
          'periodic'];
bcdim  = [0; 0; 1; 1;];
face   = [1; -1; 1; -1];
limits = [0 0 0.0 1.0; 0 0 0.0 1.0; 0.0 1.0 0 0; 0.0 1.0 0 0];

% set the commands to run the executables
nproc = 1;
for i = 1:max(size(iproc))
    nproc = nproc * iproc(i);
end
exact_exec = './EXACT > exact.log 2>&1';
hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ',hypar, ...
               ' ',petsc_flags];

% write the input files for HyPar
WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,ts,tstype, ...
    hyp_scheme,hyp_flux_split,hyp_int_type,par_type,par_scheme, ...
    dt,cons_check,screen_op_iter,file_op_iter,op_format,ip_type, ...
    input_mode,output_mode,n_io,op_overwrite,model);
WriteBoundaryInp(nb,bctype,bcdim,face,limits);
WritePhysicsInp_LinearADR(advection);
WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
% Generate initial and exact solution
system(exact_exec);
% Run the simulation
system(hypar_exec);

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
height = width * max(zlen/xlen,0.1);

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

        % plot the solution
        figure(figSolution);
        handle = contourf(xcoord,ycoord,rho,50,'LineStyle','none');
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
   
    % plot the solution
    figure(figSolution);
    handle = surf(xcoord,ycoord,rho,'LineStyle','none');
    view(0,90);
    colorbar;
    xlabel('x','FontName','Times','FontSize',20,'FontWeight','normal');
    ylabel('y','FontName','Times','FontSize',20,'FontWeight','normal');
    set(gca,'FontSize',14,'FontName','Times');
    
    % write plot to files
    print(figSolution,'-depsc','Contour.eps');
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
system('rm -rf *.dat *.inp *.log *.bin *.eps EXACT');
