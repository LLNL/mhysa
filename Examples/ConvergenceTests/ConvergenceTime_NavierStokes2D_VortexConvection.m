% Script to test the convergence in time for the
% isentropic vortex convection problem, governed 
% by the 2D Navier-Stokes equations

clear all;
close all;
system('rm -rf *.dat *.inp *.log *.bin *.eps INIT');

fprintf('Time convergence test for isentropic vortex convection.\n');
fprintf('Governing equations: 2D Navier-Stokes equations.\n');

% system-specific MPI run commands and arguments
% set these according to the system you are using
% if running in serial, set serial_flag to 1
serial_flag = 0;
mpiexec = 'mpiexec';
% mpiexec = 'srun';
mpi_args = ' '; % no additional arguments
% mpi_args = '-p pdebug';


% Ask for path to Mhysa source directory
mhysa_path = input('Enter path to Mhysa source: ','s');

% Add to MATLAB path
path(path,strcat(mhysa_path,'/Examples/Matlab/'));
% Compile the code to generate the initial solution
cmd = ['gcc ',mhysa_path, ...
       '/Examples/SingleSpecies/2D_IsentropicVortexConvection/aux/exact.c ', ...
       '-lm -o INIT'];
system(cmd);
% find the Mhysa binary
mhysa = [mhysa_path,'/bin/mhysa'];

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
            '2e    '; ...
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
t_final = 1.0;

% time step sizes
dt_implicit = [0.002 0.005 0.01 0.02 0.05 0.1];
dt_explicit = [0.002 0.005 0.01 0.02];

% plotting styles
style = [ ...
          '-ko'; ...
          '-ks'; ...
          '-kd'; ...
          '-bo'; ...
          '-bs'; ...
          '-bd'; ...
        ];

%-------------------------------------------------------------------------%
% for the 'arkimex' time-integrators, the following options
% can be chosen from ('rk' time-integrators will ignore this):-

hyp_flux_split = 'no';
hyp_flux_flag = '-hyperbolic_implicit'; % implicit treatment
% hyp_flux_flag = '-hyperbolic_explicit'; % explicit treatment
%------------------------------------------------------------------------%
    
% maximum expected error
maxerr = 1e-6;

% set up problem specific parameters
ndims = 3;                 % number of space dimensions
nvars = 5;                 % number of variables in state vector
iproc = [4 4 1];           % number of processors in each dimension
ghost = 3;                 % number of ghost points
N = [128 128 3];           % grid dimensions

% specify spatial discretization scheme details
hyp_scheme      = 'weno5';
hyp_int_type    = 'components';
par_type        = 'nonconservative-2stage';
% parameters controlling the WENO-type schemes
mapped  = 0;
borges  = 0;
yc      = 0;
nl      = 1;
eps     = 1e-6;

% set physical model and related parameters
model     = 'navierstokes3d';
gamma     = 1.4;                  % specific heat ratio
upw       = 'rusanov';            % choice of upwinding scheme
Prandtl   = 0.72;                 % Prandtl number
Reynolds  = -1.0;                 % Inviscid flow
Minf      = 1.0;                  % reference Mach number
nspecies  = 1;                    % number of species
nvibeng   = 0;                    % number of vibrational energies
% other options
cons_check      = 'yes';
screen_op_iter  = 10;
op_overwrite    = 'yes';
file_op_iter    = 999999;
ip_type         = 'binary';
% set boundaries
nb = 6;
bctype = ['periodic'; ...
          'periodic'; ...
          'periodic'; ...
          'periodic'; ...
          'periodic'; ...
          'periodic'];
bcdim     = [0; 0; 1; 1; 2; 2];
face      = [1; -1; 1; -1; 1; -1];
limits    = [   0 0 0 10.0 -1000 1000; ...
                0 0 0 10.0 -1000 1000; ...
                0 10.0 0 0 -1000 1000; ...
                0 10.0 0 0 -1000 1000; ...
                0 10.0 0 10.0 0 0;
                0 10.0 0 10.0 0 0];

% set the commands to run the executables
nproc = 1;
for i = 1:max(size(iproc))
    nproc = nproc * iproc(i);
end
init_exec  = './INIT  > init.log  2>&1';

% Generate or find the reference solution
RefFlag = input('Generate reference solution? ','s');
if (strcmp(RefFlag,'yes'))
    fprintf('Generating reference solution (this will take some time)...\n');
    % small time step for reference solution
    dt_ref = 0.05 * min(min(dt_implicit),min(dt_explicit));    
    % use explicit RK 4-stage, 4th order
    ts_ref = 'rk';              
    tstype_ref = '44';          
    niter = int32(t_final/dt_ref);
    % set reference solution output type as binary
    op_format = 'binary';
    % Write out the input files for Mhysa
    WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,ts_ref, ...
        tstype_ref,hyp_scheme,hyp_flux_split,hyp_int_type, ...
        par_type,par_scheme, dt_ref,cons_check,screen_op_iter, ...
        file_op_iter,op_format,ip_type,input_mode,output_mode, ...
        n_io,op_overwrite,model);
    WriteBoundaryInp(nb,bctype,bcdim,face,limits);
    WritePhysicsInp_NavierStokes3D(gamma,upw,Prandtl,Reynolds,Minf,...
                                   nspecies,nvibeng);
    WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
    % Generate the initial solution
    system(init_exec);
    % Run Mhysa
    if (serial_flag == 1)
        mhysa_exec = [mhysa,' > run.log 2>&1 '];
    else
        mhysa_exec = [mpiexec,' -n ',num2str(nproc),' ',mpi_args,' ', ...
                      mhysa,' > run.log 2>&1 '];
    end
    system(mhysa_exec);
    % convert the output to the input format
    system(['gcc ',mhysa_path,'/Extras/BinaryOPToInitialSolution.c ', ...
            '-o BINOP2INP']);
    fid = fopen('bin.inp','w');
    fprintf(fid,'op.bin');
    fclose(fid);
    system('./BINOP2INP < bin.inp 2>&1 > conv.log && rm bin.inp');
    system('mv solution.inp reference.bin');
    system('mv run.log reference.log');
    % save the reference solution and log in a separate directory
    dir_name = strcat('refsoln_',sprintf('%03d_',N),hyp_scheme, ...
                      '_',sprintf('%04.1f',t_final));
    if (exist(dir_name,'file'))
        fprintf('Removing existing directory %s.\n',dir_name);
        system(['rm -rf ',dir_name]);
    end
    mkdir(dir_name);
    system(['mv op.bin *.inp reference.bin reference.log ',dir_name]);
    system(['rm -rf ',dir_name,'/initial.inp']);
    % clean up
    system('rm -rf *.inp *.log *.dat *.bin BINOP2INP');
    % create a shortcut in the current folder for the reference solution
    system(['ln -sf ',dir_name,'/reference.bin reference.bin']);
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
                hyp_flux_flag, ' ', ...
                '-snes_type newtonls ', ...
                '-snes_rtol 1e-10 ', ...
                '-snes_atol 1e-10 ', ...
                '-snes_stol 1e-16 ', ...
                '-ksp_type gmres ', ...
                '-ksp_rtol 1e-10 ', ...
                '-ksp_atol 1e-10 ', ...
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
    
    % set time step sizes
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
        % set output format as text (for plotting)
        op_format = 'text';
        % Write out the input files for Mhysa
        WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,strtrim(ts(j,:)), ...
            strtrim(tstype(j,:)),hyp_scheme,hyp_flux_split,hyp_int_type, ...
            par_type,par_scheme, dt(r),cons_check,screen_op_iter, ...
            file_op_iter,op_format,ip_type,input_mode,output_mode, ...
            n_io,op_overwrite,model);
        WriteBoundaryInp(nb,bctype,bcdim,face,limits);
        WritePhysicsInp_NavierStokes3D(gamma,upw,Prandtl,Reynolds,Minf,...
                                       nspecies,nvibeng);
        WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
        % Generate the initial solution
        system(init_exec);
        % create the sym link to the exact solution
        system('ln -sf reference.bin exact.inp');
        % Run Mhysa
        if (serial_flag == 1)
            mhysa_exec = [mhysa,' ',petsc_flags,petscdt,petscft,petscms, ...
                          ' > run.log 2>&1 '];
        else
            mhysa_exec = [mpiexec,' -n ',num2str(nproc),' ',mpi_args,' ', ...
                          mhysa,' ',petsc_flags,petscdt,petscft,petscms, ...
                          ' > run.log 2>&1 '];
        end
        system(mhysa_exec);
        % save the solution
        soln_name = strcat('op_',ts(j,:),'_',tstype(j,:),'_dt', ...
                           sprintf('%06.3f',dt(r)),'.dat');
        system(['mv op.dat ',dumpname,'/',soln_name]);
        % save the run log
        log_name = strcat('run_',ts(j,:),'_',tstype(j,:),'_dt', ...
                           sprintf('%06.3f',dt(r)),'.log');
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
% axis([min(dt_implicit)/2 2*max(dt_implicit) min(MinErr)/2 4*max(MaxErr)]);
axis tight;
grid on;
hold off;

figure(figErrCost);
ylabel('Error (L_2)','FontName','Times','FontSize',20, ...
       'FontWeight','normal');
xlabel('Number of RHS function calls','FontName','Times', ...
       'FontSize',20,'FontWeight','normal');
set(gca,'FontSize',12,'FontName','Times');
legend(legend_str,'Location','northeast');
% axis([min(MinCost)/2 2*max(MaxCost) min(MinErr)/2 4*max(MaxErr)]);
axis tight;
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
system('rm -rf INIT *.dat *.log *.inp *.bin *.core');
