% Script to test the time convergence rate
% for a smooth solution of the linear advection 
% equation

clear all;
close all;

% remove all useless files
system('rm -rf *.dat *.inp *.log INIT');

fprintf('Time convergence test on a smooth solution ');
fprintf('to the linear advection equation.\n');

% Ask for path to HyPar source directory
hypar_path = input('Enter path to HyPar source: ','s');

% Add to MATLAB path
path(path,strcat(hypar_path,'/Examples/Matlab/'));

% Compile the code to generate the initial solution
cmd = ['g++ ',hypar_path, ...
       '/Examples/1D/LinearAdvection/SineWave/aux/init.C ', ...
       '-o INIT'];
system(cmd);
% find the HyPar binary
hypar = [hypar_path,'/bin/HyPar'];

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
hyp_scheme = 'crweno5';

% for spatial convergence, use very small time step
dt = [0.0001 0.0002 0.0005 0.001 0.002 0.005 ...
      0.01 0.02 0.05 0.1 0.2 0.5 1.0];
t_final = 1.0;

% maximum expected error
maxerr = 10.0;

% set physical model and related parameters
model = 'linear-advection-diffusion-reaction';
adv = 1.0;

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
schemes = 1:size(ts,1);

% turn off solution output to file
op_format = 'none';

% set grid size;
N = 320;

% set the commands to run the executables
nproc = 1;
for i = 1:max(size(iproc))
    nproc = nproc * iproc(i);
end
init_exec = './INIT > init.log 2>&1';
clean_exec = 'rm -rf *.inp *.dat *.log';

% open figure window
scrsz = get(0,'ScreenSize');
figErrDt   = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
figErrCost = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% initialize legend string
legend_str = char(zeros(size(schemes,1),(2+size(ts,2)+size(tstype,2))));

% plotting styles
style = [ ...
          '-ko'; ...
          '-ks'; ...
          '-kd'; ...
          '-kp'; ...
          '-kv'; ...
          '-k^'; ...
        ];

count = 1;
MinErr  = zeros(size(schemes,2),1);
MaxErr  = zeros(size(schemes,2),1);
MinCost = zeros(size(schemes,2),1);
MaxCost = zeros(size(schemes,2),1);
for j = schemes
    % set PETSc time-integration flags (comment to turn off)
    if (use_petsc(j)) 
        petsc_flags = sprintf('%s', ...
            '-use-petscts ', ...
            '-ts_type ',strtrim(ts(j,:)),' ', ...
            '-ts_',strtrim(ts(j,:)),'_type ',strtrim(tstype(j,:)),' ', ...
            '-ts_adapt_type none ', ...
            '-hyperbolic_implicit ', ...
            '-snes_type newtonls ', ...
            '-snes_rtol 1e-6 ', ...
            '-snes_atol 1e-6 ', ...
            '-ksp_type gmres ', ...
            '-ksp_rtol 1e-6 ', ...
            '-ksp_atol 1e-6 ', ...
            '-ksp_max_it 1000 ', ...
            '-snes_max_it 1000 ', ...
                   ' ');
    else
        petsc_flags = ' ';
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
        % Write out the input files for HyPar
        WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,strtrim(ts(j,:)), ...
            strtrim(tstype(j,:)),hyp_scheme,hyp_flux_split,hyp_int_type, ...
            par_type,par_scheme, dt(r),cons_check,screen_op_iter, ...
            file_op_iter,op_format,ip_type,input_mode,output_mode, ...
            n_io,op_overwrite,model);
        WriteBoundaryInp(nb,bctype,dim,face,limits);
        WritePhysicsInp_LinearADR(adv);
        WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
        WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
        % Generate the initial and exact solution
        system(init_exec);
        system('ln -sf initial.inp exact.inp');
        % Run HyPar
        hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ', ...
                      hypar,' ',petsc_flags,petscdt,petscft,petscms, ...
                      ' > run.log 2>&1'];
        system(hypar_exec);
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
    MinErr(j)  = min(L2Errors(L2Errors>0));
    MaxErr(j)  = max(L2Errors(L2Errors<maxerr));
    MinCost(j) = min(FCounts((FCounts>0)&(L2Errors<maxerr)));
    MaxCost(j) = max(FCounts);

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

figure(figErrDt);
xlabel('dt','FontName','Times','FontSize',20,'FontWeight','normal');
ylabel('Error','FontName','Times','FontSize',20,'FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');
legend(legend_str,'Location','northwest');
axis([min(dt)/2 2*max(dt) min(MinErr)/2 2*max(MaxErr)]);
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
system('rm -rf INIT');

