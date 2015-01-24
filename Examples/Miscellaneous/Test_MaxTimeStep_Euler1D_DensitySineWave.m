% Script to test the maximum stable time step size
% for various time-integration methods.
% Case: Density Sine Wave
% Model: 1D Euler

clear all;
close all;

fprintf('Maximum stable time step test on a smooth solution ');
fprintf('to the 1D Euler equations.\n');

% Ask for path to HyPar source directory
hypar_path = input('Enter path to HyPar source: ','s');

% Ask for initial value of dt and tolerance
dt_init = input('Enter initial dt: ');
tolerance = input('Enter step size tolerance: ');

% Add to MATLAB path
path(path,strcat(hypar_path,'/Examples/Matlab/'));

% Compile the code to generate the initial solution
cmd = ['gcc ',hypar_path, ...
       '/Examples/1D/Euler1D/DensitySineWave/aux/exact.c -lm ', ...
       '-o EXACT'];
system(cmd);
% find the HyPar binary
hypar = [hypar_path,'/bin/HyPar'];

% Get the default
[~,~,~,~,~,~,~,~,~, ...
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

% set final time
t_final = 1.0;

% maximum expected error
maxerr = 1.0;

% set physical model and related parameters
model = 'euler1d';
gamma = 1.4;
grav  = 0.0;
upw   = 'roe';

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
order = [2,3,4,2,3,4];
use_petsc = [1,1,1,0,0,0];
is_implicit = [1,1,1,0,0,0];
schemes = 1:size(ts,1);

% plotting styles
style = [ ...
          '-ko'; ...
          '-ks'; ...
          '-kd'; ...
          '-kp'; ...
          '-kv'; ...
          '-k^'; ...
        ];

% turn off solution output to file
op_format = 'none';

% set grid size;
N = 80;

%-------------------------------------------------------------------------%
% for the 'arkimex' time-integrators, the following options
% can be chosen from ('rk' time-integrators will ignore this):-

% use split hyperbolic flux form (acoustic and entropy modes)?
hyp_flux_split = 'yes';
% treat acoustic waves implicitly, and entropy waves explicitly
% hyp_flux_flag  = '-hyperbolic_f_explicit -hyperbolic_df_implicit';
% treat acoustic and entropy waves implicitly
hyp_flux_flag  = '-hyperbolic_f_implicit -hyperbolic_df_implicit';
% treat acoustic and entropy waves explicitly
% hyp_flux_flag  = '-hyperbolic_f_explicit -hyperbolic_df_explicit';

% or no splitting?
% hyp_flux_split = 'no';
% hyp_flux_flag = '-hyperbolic_implicit'; % implicit treatment
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

% set maximum number of data points
ref_levels = 1000;

count = 1;
MinDt   = zeros(size(schemes,2),1);
MaxDt   = zeros(size(schemes,2),1);
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
                '-snes_rtol 1e-8 ', ...
                '-snes_atol 1e-8 ', ...
                '-ksp_type gmres ', ...
                '-ksp_rtol 1e-8 ', ...
                '-ksp_atol 1e-8 ', ...
                '-ksp_max_it 1000 ', ...
                '-snes_max_it 1000 ', ...
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

    % set dt to its initial value
    dt = dt_init;
    dt_factor = 1.0;
    r = 1;

    % set max time step size to final time
    dt_max = t_final;
    
    % preallocate arrays for dt, error, wall times and function counts
    TimeStep  = zeros(ref_levels,1);
    Errors    = zeros(ref_levels,3);
    Walltimes = zeros(ref_levels,2);
    FCounts   = zeros(ref_levels,1);

    % run simulation with initial dt
    fprintf('\t%s, dt=%1.6e, factor=%8.6f: ',[ts(j,:),' ',tstype(j,:)],dt,dt_factor);
    niter = int32(t_final/dt);
    if (strcmp(petsc_flags,' ')) 
        petscdt = ' ';
        petscft = ' ';
        petscms = ' ';
    else
        petscdt = [' -ts_dt ',num2str(dt,'%1.16e'),' '];
        petscft = [' -ts_final_time ',num2str(t_final,'%f'),' '];
        petscms = [' -ts_max_steps ',num2str(100*niter,'%d'),' '];
    end
    WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,strtrim(ts(j,:)), ...
        strtrim(tstype(j,:)),hyp_scheme,hyp_flux_split,hyp_int_type, ...
        par_type,par_scheme, dt,cons_check,screen_op_iter, ...
        file_op_iter,op_format,ip_type,input_mode,output_mode, ...
        n_io,op_overwrite,model);
    WriteBoundaryInp(nb,bctype,dim,face,limits);
    WritePhysicsInp_Euler1D(gamma,grav,upw);
    WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
    WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
    system(exact_exec);
    hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ', ...
                  hypar,' ',petsc_flags,petscdt,petscft,petscms, ...
                  ' > run.log 2>&1' ...
                  ];
    system(hypar_exec);
    [err,wt] = ReadErrorDat(ndims);
    [~,fcounts,~,~,~,~,~,~] = ReadFunctionCounts();
    system(clean_exec);
    if (~min(isfinite(err)))
        fprintf('failed.\n');
        continue;
    end
    fprintf('passed.\n');
    TimeStep(r)     = dt;
    Errors(r,:)     = err;
    Walltimes(r,:)  = wt;
    FCounts(r)      = fcounts;
    r = r+1;

    while ((dt_factor > tolerance) && (r < ref_levels))
        dt_new = dt * (1.0+dt_factor);
        while (dt_new > dt_max)
            dt_factor = 0.5 * dt_factor;
            dt_new = dt * (1.0+dt_factor);
        end

        % estimate error from previous error based on theoretical order
        err_theoretical = Errors(r-1,2) * (dt_new/dt)^order(j);

        fprintf('\t%s, dt=%1.6e, factor=%8.6f: ',[ts(j,:),' ',tstype(j,:)],dt_new,dt_factor);
        niter = floor(t_final/dt_new);
        if (strcmp(petsc_flags,' ')) 
            petscdt = ' ';
            petscft = ' ';
            petscms = ' ';
        else
            petscdt = [' -ts_dt ',num2str(dt_new,'%1.16e'),' '];
            petscft = [' -ts_final_time ',num2str((niter*dt_new),'%f'),' '];
            petscms = [' -ts_max_steps ',num2str(100*niter,'%d'),' '];
        end
        WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,strtrim(ts(j,:)), ...
            strtrim(tstype(j,:)),hyp_scheme,hyp_flux_split,hyp_int_type, ...
            par_type,par_scheme, dt_new,cons_check,screen_op_iter, ...
            file_op_iter,op_format,ip_type,input_mode,output_mode, ...
            n_io,op_overwrite,model);
        WriteBoundaryInp(nb,bctype,dim,face,limits);
        WritePhysicsInp_Euler1D(gamma,grav,upw);
        WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
        WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
        system(exact_exec);
        hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ', ...
                      hypar,' ',petsc_flags,petscdt,petscft,petscms, ...
                      ' > run.log 2>&1' ...
                      ];
        system(hypar_exec);
        [err,wt] = ReadErrorDat(ndims);
        [~,fcounts,~,~,~,~,~,~] = ReadFunctionCounts();
        system(clean_exec);
        if ((~min(isfinite(err))) || (err(2)/err_theoretical > 2))
            fprintf('failed.\n');
            dt_max = dt_new;
            dt_factor = 0.5 * dt_factor;
        else
            fprintf('passed.\n');
            dt = dt_new;
            TimeStep(r) = dt;
            Errors(r,:) = err;
            Walltimes(r,:) = wt;
            FCounts(r) = fcounts;
            r = r+1;
        end
    end
    
    % Isolate the L2 Error
    L2Errors = Errors(:,2);
    
    % To be used in setting axis limits
    MinDt(count)   = min(TimeStep(1:r-1));
    MaxDt(count)   = max(TimeStep(1:r-1));
    MinErr(count)  = min(L2Errors(1:r-1));
    MaxErr(count)  = max(L2Errors(1:r-1));
    MinCost(count) = min(FCounts (1:r-1));
    MaxCost(count) = max(FCounts (1:r-1));

    % plot errors
    figure(figErrDt);
    loglog(TimeStep(1:r-1),L2Errors(1:r-1),style(count,:),'linewidth',1,'MarkerSize',8);
    hold on;
    % plot cost
    figure(figErrCost);
    loglog(FCounts(1:r-1),L2Errors(1:r-1),style(count,:),'linewidth',1, ...
           'MarkerSize',8);
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
axis([min(MinDt)/2 max(MaxDt)*2 min(MinErr)/2 min(2*max(MaxErr),maxerr)]);
grid on;
hold off;

figure(figErrCost);
ylabel('Error (L_2)','FontName','Times','FontSize',20, ...
       'FontWeight','normal');
xlabel('Number of RHS function calls','FontName','Times', ...
       'FontSize',20,'FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');
legend(legend_str,'Location','northeast');
axis([min(MinCost)/2 2*max(MaxCost) min(MinErr)/2 min(2*max(MaxErr),maxerr)]);
grid on;
hold off;

% print figures to file
print(figErrDt,'-depsc2','figErrDt.eps');
print(figErrCost,'-depsc2','figErrCost.eps');

% clean up
system('rm -rf EXACT');

