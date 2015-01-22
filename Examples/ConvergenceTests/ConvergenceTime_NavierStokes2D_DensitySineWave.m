% Script to test the convergence in time for the
% inviscid advection of a density wave, governed 
% by the 2D Navier-Stokes equations

clear all;
close all;

fprintf('Time convergence test for inviscid advection of density wave.\n');
fprintf('Governing equations: 2D Navier-Stokes equations.\n');

% Ask for path to HyPar source directory
hypar_path = input('Enter path to HyPar source: ','s');

% Add to MATLAB path
path(path,strcat(hypar_path,'/Examples/Matlab/'));
% Compile the code to generate the initial and exact solutions
cmd = ['g++ ',hypar_path, ...
    '/Examples/2D/NavierStokes2D/DensitySineWave/aux/exact.C ', ...
    '-o EXACT'];
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
t_final = 0.1;

% time step sizes
dt_implicit = [0.0001 0.0002 0.0005 0.001 0.002 0.005];
dt_explicit = [0.0001 0.0002 0.0005 0.001];

% plotting styles
style = [ ...
          '-ko'; ...
          '-ks'; ...
          '-kd'; ...
          '-bo'; ...
          '-bs'; ...
          '-bd'; ...
        ];
    
% maximum expected error
maxerr = 10.0;

% set up problem specific parameters
ndims = 2;                 % number of space dimensions
nvars = 4;                 % number of variables in state vector
iproc = [2 2];             % number of processors in each dimension
ghost = 3;                 % number of ghost points
N = [160 160];             % grid dimensions

%-------------------------------------------------------------------------%
% for the 'arkimex' time-integrators, the following options
% can be chosen from ('rk' time-integrators will ignore this):-

% use split hyperbolic flux form (acoustic and entropy modes)?
% hyp_flux_split = 'yes';
% treat acoustic waves implicitly, and entropy waves explicitly
% hyp_flux_flag  = '-hyperbolic_f_explicit -hyperbolic_df_implicit';
% treat acoustic and entropy waves implicitly
% hyp_flux_flag  = '-hyperbolic_f_implicit -hyperbolic_df_implicit';
% treat acoustic and entropy waves explicitly
% hyp_flux_flag  = '-hyperbolic_f_explicit -hyperbolic_df_explicit';

% or no splitting?
hyp_flux_split = 'no';
hyp_flux_flag = '-hyperbolic_implicit'; % implicit treatment
% hyp_flux_flag = '-hyperbolic_explicit'; % explicit treatment
%------------------------------------------------------------------------%

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
model     = 'navierstokes2d';
gamma     = 1.4;                  % specific heat ratio
upw       = 'roe';                % choice of upwinding scheme
Prandtl   = 0.72;                 % Prandtl number
Reynolds  = -1.0;                 % Inviscid flow
Minf      = 1.0;                  % reference Mach number
grav      = [0.0 0.0];            % gravitational force vector
rho_ref   = 1.0;                  % reference altitude density
p_ref     = 1.0;                  % reference altitude pressure
HB        = 0;                    % type of hydrostatic balance
BV        = 0.0;                  % Brunt-Vaisala frequency
GasConst  = 1.0;                  % Universal gas constant
% other options
cons_check      = 'no';
screen_op_iter  = 10;
op_overwrite    = 'yes';
file_op_iter    = 999999;
ip_type         = 'ascii';
op_format       = 'none';
% set boundaries
nb = 4;
bctype = ['periodic'; ...
          'periodic'; ...
          'periodic'; ...
          'periodic'];
bcdim     = [0; 0; 1; 1;];
face      = [1; -1; 1; -1];
limits    = [0 0 0 1.0; 0 0 0 1.0; 0 1.0 0 0; 0 1.0 0 0];

% set the commands to run the executables
nproc = 1;
for i = 1:max(size(iproc))
    nproc = nproc * iproc(i);
end
exact_exec  = './EXACT  > exact.log  2>&1';

% open figure windows
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
                '-snes_rtol 1e-8 ', ...
                '-snes_atol 1e-8 ', ...
                '-snes_stol 1e-16 ', ...
                '-ksp_type gmres ', ...
                '-ksp_rtol 1e-8 ', ...
                '-ksp_atol 1e-8 ', ...
                '-ksp_max_it 1000 ', ...
                '-snes_max_it 1000 ', ...
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
        WriteBoundaryInp(nb,bctype,bcdim,face,limits);
        WritePhysicsInp_NavierStokes2D(gamma,upw,Prandtl,Reynolds,Minf,grav, ...
                                       rho_ref,p_ref,HB,BV,GasConst);
        WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
        WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
        % Generate the initial and exact solution
        system(exact_exec);
        % Run HyPar
        hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ', ...
                      hypar,' ',petsc_flags,petscdt,petscft,petscms, ...
                      ' > run.log 2>&1'];
        system(hypar_exec);
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
axis([min(dt_implicit)/2 2*max(dt_implicit) min(MinErr)/2 4*max(MaxErr)]);
grid on;
hold off;

figure(figErrCost);
ylabel('Error (L_2)','FontName','Times','FontSize',20, ...
       'FontWeight','normal');
xlabel('Number of RHS function calls','FontName','Times', ...
       'FontSize',20,'FontWeight','normal');
set(gca,'FontSize',12,'FontName','Times');
legend(legend_str,'Location','northeast');
axis([min(MinCost)/2 2*max(MaxCost) min(MinErr)/2 4*max(MaxErr)]);
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
system('rm -rf EXACT *.dat *.log *.inp *.bin');
