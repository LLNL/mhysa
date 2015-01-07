% Script to test the convergence in time for the
% inertia-gravity wave problem, governed by the 
% 2D Navier-Stokes equations (with gravity)

clear all;
close all;
system('rm -rf *.dat *.inp *.log *.bin INIT PP');

fprintf('Time convergence test for the inertial gravity wave.\n');
fprintf('Governing equations: 2D Navier-Stokes equations.\n');

% Ask for path to HyPar source directory
hypar_path = input('Enter path to HyPar source: ','s');

% Add to MATLAB path
path(path,strcat(hypar_path,'/Examples/Matlab/'));

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

% Get the default set of parameters
[~,~,~,~,~,~,~,~,~, ...
    ~,~,~,par_scheme,~,~, ...
    ~,~,~, ~,input_mode, ...
    output_mode,n_io,~,~,~,~,~,~,~, ...
    ~,~,~,~,~,p,rc,xi,wtol,lutype,norm,maxiter,atol,rtol, ...
    verbose] = SetDefaults();

% time integration methods to test
ts = [ ...
        'arkimex'; ...
        'arkimex'; ...
        'arkimex'; ...
        'rk     '; ...
        'rk     '; ...
        'rk     '  ...
     ];
tstype = [
            'a2    '; ...
            '3     '; ...
            '4     '; ...
            '22    '; ...
            'ssprk3'; ...
            '44    '  ...
         ];
use_petsc = [1,1,1,0,0,0];
implicit  = [1,1,1,0,0,0];
schemes   = [1,2,3,4,5,6];

% final time
t_final = 10.0;

% time step sizes
dt_implicit = [0.01 0.02 0.05 0.1 0.2 0.5 1.0 2.0 5.0];
dt_explicit = [0.01 0.02 0.05 0.1 0.2 0.5];

% plotting styles
style = [ ...
          '-ko'; ...
          '-ks'; ...
          '-k^'; ...
          '-bo'; ...
          '-bs'; ...
          '-b^'; ...
        ];
    
% maximum expected error
maxerr = 10.0;

% set up problem specific parameters
ndims = 2;                  % number of space dimensions
nvars = 4;                  % number of variables in state vector
iproc = [12 1];             % number of processors in each dimension
ghost = 3;                  % number of ghost points
N = [1200 50];              % grid dimensions

% specify spatial discretization scheme details
hyp_scheme      = 'weno5';
hyp_int_type    = 'components';
hyp_flux_split  = 'no';
par_type        = 'nonconservative-2stage';
% parameters controlling the WENO-type schemes
mapped  = 0;
borges  = 0;
yc      = 0;
nl      = 1;
eps     = 1e-6;

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
file_op_iter    = 999999;
ip_type         = 'binary';
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

% Generate or find the reference solution
RefFlag = input('Generate reference solution? ','s');
if (strcmp(RefFlag,'yes'))
    fprintf('Generating reference solution...\n');
    % small time step for reference solution
    dt_ref = 0.05 * min(min(dt_implicit),min(dt_explicit));    
    % use explicit RK 4-stage, 4th order
    ts_ref = 'rk';              
    tstype_ref = '44';          
    niter = int32(t_final/dt_ref);
    % Write out the input files for HyPar
    WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,ts_ref, ...
        tstype_ref,hyp_scheme,hyp_flux_split,hyp_int_type, ...
        par_type,par_scheme, dt_ref,cons_check,screen_op_iter, ...
        file_op_iter,op_format,ip_type,input_mode,output_mode, ...
        n_io,op_overwrite,model);
    WriteBoundaryInp(nb,bctype,bcdim,face,limits,WallVel1,WallVel2);
    WritePhysicsInp_NavierStokes2D(gamma,upw,Prandtl,Reynolds,Minf,grav, ...
                                   rho_ref,p_ref,HB,BV,GasConst);
    WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
    WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
    % Generate the initial solution
    system(init_exec);
    % Run HyPar
    hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ', ...
                  hypar,' > run.log 2>&1'];
    system(hypar_exec);
    % postprocess the binary output file
    fid = fopen('pp.inp','w');
    fprintf(fid,'0\n');
    fclose(fid);
    system('./PP < pp.inp 2>&1 > pp.log');
    % convert the output to the input format
    system(['gcc ',hypar_path,'/Extras/BinaryOPToInitialSolution.c ', ...
            '-o BINOP2INP']);
    fid = fopen('bin.inp','w');
    fprintf(fid,'op.bin');
    fclose(fid);
    system('./BINOP2INP < bin.inp 2>&1 > conv.log');
    system('mv solution.inp reference.bin');
    system('mv run.log reference.log');
    % save the reference solution and log in a separate directory
    dir_name = strcat('refsoln_',sprintf('%04d_',N),hyp_scheme, ...
                      '_',sprintf('%04.1f',t_final));
    mkdir(dir_name);
    system(['mv op.dat reference.bin reference.log ',dir_name]);
    system(['ln -sf ',dir_name,'/reference.bin reference.bin']);
    % clean up
    system('rm -rf *.inp *.log *.dat BINOP2INP op.bin');
else
    refpath = input('Enter path to reference solution: ','s');
    if (~strcmp(refpath,'./')) && (~strcmp(refpath,'.'))
        if (~exist([refpath,'/reference.bin'],'file'))
            fprintf('Error: reference solution not found in %s.\n', ...
                    refpath);
            return;
        end
        system(['ln -sf ',refpath,'/reference.bin reference.bin']);
    end
end

% open figure window
scrsz = get(0,'ScreenSize');
figErrDt   = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
figErrCost = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% initialize legend string
legend_str = char(zeros(size(schemes,1),(2+size(ts,2)+size(tstype,2))));

% create the directory to dump all the log and solutions files in
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

count = 1;
MinErr  = zeros(size(schemes,2),1);
MaxErr  = zeros(size(schemes,2),1);
MinCost = zeros(size(schemes,2),1);
MaxCost = zeros(size(schemes,2),1);
for j = schemes
    % set legend string
    name_str = [ts(j,:),'(',tstype(j,:),')'];
    legend_str(count,:) = name_str;

    % set PETSc time-integration flags if this is a PETSc TI method
    if (use_petsc(j))
        petsc_flags1 = sprintf('%s', ...
            '-use-petscts ', ...
            '-ts_type ',strtrim(ts(j,:)),' ', ...
            '-ts_',strtrim(ts(j,:)),'_type ',strtrim(tstype(j,:)),' ', ...
            '-ts_adapt_type none ', ...
            ' ');
        if (implicit(j))
            petsc_flags2 = sprintf('%s', ...
                '-hyperbolic_implicit ', ...
                '-snes_type newtonls ', ...
                '-snes_rtol 1e-6 ', ...
                '-snes_atol 1e-6 ', ...
                '-snes_stol 1e-16 ', ...
                '-ksp_type gmres ', ...
                '-ksp_rtol 1e-6 ', ...
                '-ksp_atol 1e-6 ', ...
                '-ksp_max_it 1000 ', ...
                '-snes_max_it 1000 ', ...
                '-ksp_monitor ', ...
                '-snes_monitor ', ...
                ' ');
        else
            petsc_flags2 = ' ';
        end
        petsc_flags = sprintf('%s',petsc_flags1,petsc_flags2);
    else
        petsc_flags = ' ';
    end
    
    % set number of grid refinement levels
    if (implicit(j))
        dt = dt_implicit;
    else
        dt = dt_explicit;
    end
    ref_levels = max(size(dt));
    
    % preallocate arrays for dx, error and wall times
    Errors    = zeros(ref_levels,3);
    Walltimes = zeros(ref_levels,2);
    FCounts   = zeros(ref_levels,1);
    
    % convergence analysis
    for r = 1:ref_levels
        fprintf('\t%s %2d:  dt=%1.16e\n',[ts(j,:),' ',tstype(j,:)], ...
                r,dt(r));
        niter = int32(t_final/dt(r));
        if (strcmp(petsc_flags,' '))
            petscdt = ' ';
            petscft = ' ';
            petscms = ' ';
        else
            petscdt = [' -ts_dt ',num2str(dt(r),'%1.16e'),' '];
            petscft = [' -ts_final_time ',num2str(t_final,'%f'),' '];
            petscms = [' -ts_max_steps ',num2str(100*niter,'%d'),' '];
        end
        % Write out the input files for HyPar
        WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,strtrim(ts(j,:)), ...
            strtrim(tstype(j,:)),hyp_scheme,hyp_flux_split,hyp_int_type, ...
            par_type,par_scheme, dt(r),cons_check,screen_op_iter, ...
            file_op_iter,op_format,ip_type,input_mode,output_mode, ...
            n_io,op_overwrite,model);
        WriteBoundaryInp(nb,bctype,bcdim,face,limits,WallVel1,WallVel2);
        WritePhysicsInp_NavierStokes2D(gamma,upw,Prandtl,Reynolds,Minf,grav, ...
                                       rho_ref,p_ref,HB,BV,GasConst);
        WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
        WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
        % Generate the initial
        system(init_exec);
        % Let the reference solution be the exact solution
        system('ln -sf reference.bin exact.inp');
        % Run HyPar
        hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ', ...
                      hypar,' ',petsc_flags,petscdt,petscft,petscms, ...
                      ' > run.log 2>&1'];
        system(hypar_exec);
        % postprocess the binary output file
        fid = fopen('pp.inp','w');
        fprintf(fid,'0\n');
        fclose(fid);
        system('./PP < pp.inp 2>&1 > pp.log');
        % save the text file
        soln_name = strcat('op_',ts(j,:),'_',tstype(j,:),'_dt', ...
                           sprintf('%05.2f',dt(r)),'.dat');
        system(['mv op.dat ',dumpname,'/',soln_name]);
        % save the run log
        log_name = strcat('run_',ts(j,:),'_',tstype(j,:),'_dt', ...
                           sprintf('%05.2f',dt(r)),'.log');
        system(['mv run.log ',dumpname,'/',log_name]);
        % Read in the errors
        [Errors(r,:),Walltimes(r,:)] = ReadErrorDat(ndims);
        % Read in function counts
        [~,FCounts(r),~,~,~,~,~,~] = ReadFunctionCounts();
        % Clean up
        system('rm -rf *.inp *.log *.dat op.*');
    end
    
    % Isolate the L2 Error
    L2Errors = Errors(:,2);

    % To be used in setting axis limits
    MinErr(count)  = max(1e-16,min(L2Errors));
    MaxErr(count)  = min(maxerr,max(L2Errors));
    MinCost(count) = min(FCounts((FCounts>0)&(L2Errors<maxerr)));
    MaxCost(count) = max(FCounts);
    
    % plot errors
    figure(figErrDt);
    loglog(dt,Errors(:,2),style(j,:),'linewidth',1,'MarkerSize',10);
    hold on;
    % plot cost
    figure(figErrCost);
    loglog(FCounts,Errors(:,2),style(j,:),'linewidth',1, ...
           'MarkerSize',10);
    hold on;

    count = count+1;
end

figure(figErrDt);
xlabel('dt','FontName','Times','FontSize',20,'FontWeight','normal');
ylabel('Error','FontName','Times','FontSize',20,'FontWeight','normal');
set(gca,'FontSize',12,'FontName','Times');
legend(legend_str,'Location','northwest');
axis([min(dt_implicit)/2 2*max(dt_implicit) min(MinErr)/2 10*max(MaxErr)]);
grid on;
hold off;

figure(figErrCost);
ylabel('Error (L_2)','FontName','Times','FontSize',20, ...
       'FontWeight','normal');
xlabel('Number of RHS function calls','FontName','Times', ...
       'FontSize',20,'FontWeight','normal');
set(gca,'FontSize',12,'FontName','Times');
legend(legend_str,'Location','northeast');
axis([min(MinCost)/2 2*max(MaxCost) min(MinErr)/2 10*max(MaxErr)]);
grid on;
hold off;

% print plots to file
filename = strcat('Error_Timestep_', ...
                  sprintf('%04d',N(1)),'_', ...
                  sprintf('%04d',N(2)),'_', ...
                  hyp_scheme,'.eps');
print(figErrDt,'-depsc2',filename);
filename = strcat('Error_Cost_', ...
                  sprintf('%04d',N(1)),'_', ...
                  sprintf('%04d',N(2)),'_', ...
                  hyp_scheme,'.eps');
print(figErrCost,'-depsc2',filename);
system(['mv *.eps ',dumpname,'/']);

% clean up
system('rm -rf INIT PP reference.bin');
