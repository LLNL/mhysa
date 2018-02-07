% Script to test the time convergence rate
% for a smooth solution of the 1D Euler equations

clear all;
close all;

% system-specific MPI run commands and arguments
% set these according to the system you are using
% if running in serial, set serial_flag to 1
serial_flag = 1;
mpiexec = 'mpiexec';
% mpiexec = 'srun';
mpi_args = ' '; % no additional arguments
% mpi_args = '-p pdebug';

% remove all useless files
system('rm -rf *.dat *.inp *.log EXACT');

fprintf('Time convergence test on a smooth solution ');
fprintf('to the 1D Euler equations.\n');

% Ask for path to Mhysa source directory
mhysa_path = input('Enter path to Mhysa source: ','s');

% Add to MATLAB path
path(path,strcat(mhysa_path,'/Examples/Matlab/'));

% Compile the code to generate the initial solution
cmd = ['gcc ',mhysa_path, ...
       '/Examples/SingleSpecies/1D_DensitySineWaveAdvection/aux/exact.c -lm ', ...
       '-o EXACT'];
system(cmd);
% find the Mhysa binary
mhysa = [mhysa_path,'/bin/mhysa'];

% Get the default
[~,~,~,~,~,niter,~,~,~, ...
    ~,hyp_int_type,par_type,par_scheme,~,cons_check, ...
    screen_op_iter,file_op_iter,~, ~,input_mode, ...
    output_mode,n_io,op_overwrite,~,nb,bctype,dim,face,limits, ...
    mapped,borges,yc,nl,eps,p,rc,xi,wtol,lutype,norm,maxiter,atol,rtol, ...
    verbose] = SetDefaults();

% Set initial solution file type to ASCII
ip_type = 'ascii';

% set problem specific input parameters
ndims = 1;
nvars = 3;
iproc = 1;
ghost = 3;

% specify a nice, high-order spatial discretization scheme
hyp_scheme = 'weno5';

% for spatial convergence, use very small time step
dt_ex = [0.0001 0.0002 0.0005 0.001 0.002];
dt_im = [0.0001 0.0002 0.0005 0.001 0.002 0.005];
t_final = 0.5;

% maximum expected error
maxerr = 1.0;

% set physical model and related parameters
model = 'euler1d';
gamma = 1.4;
upw   = 'rusanov';
nspecies = 1;
nvibeng = 0;

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
is_implicit = [1,1,1,0,0,0];
schemes = 1:size(ts,1);

% turn off solution output to file
op_format = 'none';

% set grid size;
N = 320;

hyp_flux_split = 'no';

%-------------------------------------------------------------------------%
% for the 'arkimex' time-integrators, the following options
% can be chosen from ('rk' time-integrators will ignore this):-

hyp_flux_flag = '-hyperbolic_implicit'; % implicit treatment
% hyp_flux_flag = '-hyperbolic_explicit'; % explicit treatment
%------------------------------------------------------------------------%

% set the commands to run the executables
nproc = 1;
for i = 1:max(size(iproc))
    nproc = nproc * iproc(i);
end
exact_exec = './EXACT > exact.log 2>&1';
clean_exec = 'rm -rf *.inp *.dat *.log';

% open figure window
scrsz = get(0,'ScreenSize');
figErrDt   = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
figErrCost = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% initialize legend string
legend_str = char(zeros(size(schemes,1),(2+size(ts,2)+size(tstype,2))));

% plotting styles
style = [ ...
          '-ro'; ...
          '-rs'; ...
          '-rd'; ...
          '-bp'; ...
          '-bv'; ...
          '-b^'; ...
        ];

count = 1;
MinErr  = zeros(size(schemes,2),1);
MaxErr  = zeros(size(schemes,2),1);
MinCost = zeros(size(schemes,2),1);
MaxCost = zeros(size(schemes,2),1);
for j = schemes
    % set PETSc time-integration flags (comment to turn off)
    if (use_petsc(j))
        if (is_implicit(j))
            petsc_flags = sprintf('%s', ...
                '-use-petscts ', ...
                '-ts_type ',strtrim(ts(j,:)),' ', ...
                '-ts_',strtrim(ts(j,:)),'_type ',strtrim(tstype(j,:)),' ', ...
                '-ts_adapt_type none ', ...
                hyp_flux_flag,' ', ...
                '-snes_type newtonls ', ...
                '-snes_rtol 1e-10 ', ...
                '-snes_atol 1e-12 ', ...
                '-ksp_type gmres ', ...
                '-ksp_rtol 1e-12 ', ...
                '-ksp_atol 1e-12 ', ...
                '-ksp_max_it 10000 ', ...
                '-snes_max_it 10000 ', ...
                ' ');
        else
            petsc_flags = sprintf('%s', ...
                '-use-petscts ', ...
                '-ts_type ',strtrim(ts(j,:)),' ', ...
                '-ts_',strtrim(ts(j,:)),'_type ',strtrim(tstype(j,:)),' ', ...
                '-ts_adapt_type none ', ...
                ' ');
        end
    else
        petsc_flags = ' ';
    end
    
    if (is_implicit(j))
        dt = dt_im;
    else
        dt = dt_ex;
    end
    % set number of grid refinement levels
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
        % Write out the input files for Mhysa
        WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,strtrim(ts(j,:)), ...
            strtrim(tstype(j,:)),hyp_scheme,hyp_flux_split,hyp_int_type, ...
            par_type,par_scheme, dt(r),cons_check,screen_op_iter, ...
            file_op_iter,op_format,ip_type,input_mode,output_mode, ...
            n_io,op_overwrite,model);
        WriteBoundaryInp(nb,bctype,dim,face,limits);
        WritePhysicsInp_Euler1D(gamma,upw,nspecies,nvibeng);
        WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
        WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
        % Generate the initial and exact solution
        system(exact_exec);
        % Run MHYSA
        if (serial_flag == 1)
            mhysa_exec = [mhysa,' ',petsc_flags,petscdt,petscft,petscms, ...
                          ' > run.log 2>&1 '];
        else
            mhysa_exec = [mpiexec,' -n ',num2str(nproc),' ',mpi_args,' ', ...
                          mhysa,' ',petsc_flags,petscdt,petscft,petscms, ...
                          ' > run.log 2>&1 '];
        end
        system(mhysa_exec);
        % Read in the errors
        [Errors(r,:),Walltimes(r,:)] = ReadErrorDat(ndims);
        % Read in function counts
        [~,FCounts(r),~,~,~,~,~,~] = ReadFunctionCounts();
        % Clean up
        system(clean_exec);    
    end
    
    % Isolate the L2 Error
    L2Errors = Errors(:,2);
    
    % To be used in setting axis limits
    MinErr(count)  = min(L2Errors(L2Errors>0));
    MaxErr(count)  = max(L2Errors(L2Errors<maxerr));
    MinCost(count) = min(FCounts((FCounts>0)&(L2Errors<maxerr)));
    MaxCost(count) = max(FCounts);

    % plot errors
    figure(figErrDt);
    loglog(dt,Errors(:,2),style(count,:),'linewidth',2,'MarkerSize',10);
    hold on;
    % plot cost
    figure(figErrCost);
    loglog(FCounts,Errors(:,2),style(count,:),'linewidth',2, ...
           'MarkerSize',10);
    hold on;
    
    % set legend string
    name_str = [ts(j,:),'(',tstype(j,:),')'];
    legend_str(count,:) = name_str;
    
    count = count+1;
end

xmin = 0.5*min([dt_im dt_ex]);
xmax = 2.0*max([dt_im dt_ex]);

figure(figErrDt);
xlabel('dt','FontName','Times','FontSize',20,'FontWeight','normal');
ylabel('Error','FontName','Times','FontSize',20,'FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');
legend(legend_str,'Location','northwest');
axis([xmin xmax min(MinErr)/2 2*max(MaxErr)]);
grid on;
hold off;

figure(figErrCost);
ylabel('Error (L_2)','FontName','Times','FontSize',20, ...
       'FontWeight','normal');
xlabel('Number of RHS function calls','FontName','Times', ...
       'FontSize',20,'FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');
legend(legend_str,'Location','northeast');
axis([min(MinCost)/2 2*max(MaxCost) min(MinErr)/2 2*max(MaxErr)]);
grid on;
hold off;

% clean up
system('rm -rf EXACT');

