% Script to test the maximum stable time step size
% for various time-integration methods.
% Case: Isentropic Vortex Convection
% Model: 2D Navier-Stokes

clear all;
close all;

fprintf('Maximum stable time step test on a smooth solution ');
fprintf('to the 2D Navier-Stokes equations.\n');

% Ask for path to HyPar source directory
hypar_path = input('Enter path to HyPar source: ','s');

% Ask for initial value of dt and tolerance
dt_init = input('Enter initial dt: ');
tolerance = input('Enter step size tolerance: ');

% Add to MATLAB path
path(path,strcat(hypar_path,'/Examples/Matlab/'));

% Compile the code to generate the initial solution
cmd = ['gcc ',hypar_path, ...
       '/Examples/2D/NavierStokes2D/InviscidVortexConvection/aux/init.c ', ...
       '-lm -o INIT'];
system(cmd);
% find the HyPar binary
hypar = [hypar_path,'/bin/HyPar'];

% Get the default
[~,~,~,~,~,~,~,~,~, ...
    ~,~,~,par_scheme,~,~, ...
    ~,~,~, ~,input_mode, ...
    output_mode,n_io,~,~,~,~,~,~,~, ...
    ~,~,~,~,~,p,rc,xi,wtol,lutype,norm,maxiter,atol,rtol, ...
    verbose] = SetDefaults();

% time integration methods to test
% do not use native time-integrators, use only PETSc ones
% if the final time is not an integer multiple of the  time step 
% size being tried, native time-integrators will not yield a 
% solution at that exact final time --> the error will not be
% the true error.
ts = [ ...
        'arkimex'; ...
        'arkimex'; ...
        'arkimex'; ...
        'rk     '; ...
        'rk     '; ...
        'rk     '; ...
        'rk     '; ...
        'rk     '  ...
     ];
tstype = [
            '2e '; ...
            '3  '; ...
            '4  '; ...
            '2a '; ...
            '3  '; ...
            '3bs'; ...
            '4  ';  ...
            '5f '   ...
         ];
order = [2,3,4,2,3,3,4,5];
use_petsc = [1,1,1,1,1,1,1,1];
is_implicit = [1,1,1,0,0,0,0,0];
schemes = 1:size(ts,1);

% set final time
t_final = 1.0;

% plotting styles
style = [ ...
          '-bo'; ...
          '-bs'; ...
          '-bd'; ...
          '-ko'; ...
          '-ks'; ...
          '-kp'; ...
          '-kv'; ...
          '-k^'; ...
        ];

% check
if (size(style,1) < size(schemes,2))
    fprintf('Error: number of plotting styles defined less then number of methods to be tested.\n');
    return;
end

% set problem specific input parameters
ndims = 2;
nvars = 4;
iproc = [4 4];
ghost = 3;
% set grid size;
N = [128 128];

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
hyp_scheme      = 'crweno5';
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
ip_type         = 'binary';
op_format       = 'none';

% maximum expected error
maxerr = 1.0;

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
init_exec = './INIT > init.log 2>&1';
clean_exec = 'rm -rf *.inp *.dat *.log';

% Generate or find the reference solution
RefFlag = input('Generate reference solution? ','s');
if (strcmp(RefFlag,'yes'))
    fprintf('Generating reference solution...\n');
    % small time step for reference solution
    dt_ref = 0.05 * dt_init;    
    % use explicit RK 4-stage, 4th order
    ts_ref = 'rk';              
    tstype_ref = '44';          
    niter = int32(t_final/dt_ref);
    % set reference solution output type as binary
    op_format = 'binary';
    % Write out the input files for HyPar
    WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,ts_ref, ...
        tstype_ref,hyp_scheme,hyp_flux_split,hyp_int_type, ...
        par_type,par_scheme, dt_ref,cons_check,screen_op_iter, ...
        file_op_iter,op_format,ip_type,input_mode,output_mode, ...
        n_io,op_overwrite,model);
    WriteBoundaryInp(nb,bctype,bcdim,face,limits);
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
    % convert the output to the input format
    system(['gcc ',hypar_path,'/Extras/BinaryOPToInitialSolution.c ', ...
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

    % preallocate arrays for dt, error, wall times and function counts
    TimeStep  = zeros(ref_levels,1);
    Errors    = zeros(ref_levels,3);
    Walltimes = zeros(ref_levels,2);
    FCounts   = zeros(ref_levels,1);

    % run simulation with initial dt
    fprintf('\t%s, dt=%1.6e, factor=%8.6f: ',[ts(j,:),' ',tstype(j,:)],dt,dt_factor);
    dt_max = t_final;
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
    WriteBoundaryInp(nb,bctype,bcdim,face,limits);
    WritePhysicsInp_NavierStokes2D(gamma,upw,Prandtl,Reynolds,Minf,grav, ...
                                   rho_ref,p_ref,HB,BV,GasConst);
    WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
    WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
    system(init_exec);
    system('ln -sf reference.bin exact.inp');
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
        niter = int32(t_final/dt_new);
        if (strcmp(petsc_flags,' ')) 
            petscdt = ' ';
            petscft = ' ';
            petscms = ' ';
        else
            petscdt = [' -ts_dt ',num2str(dt_new,'%1.16e'),' '];
            petscft = [' -ts_final_time ',num2str(t_final,'%f'),' '];
            petscms = [' -ts_max_steps ',num2str(100*niter,'%d'),' '];
        end
        WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,strtrim(ts(j,:)), ...
            strtrim(tstype(j,:)),hyp_scheme,hyp_flux_split,hyp_int_type, ...
            par_type,par_scheme, dt_new,cons_check,screen_op_iter, ...
            file_op_iter,op_format,ip_type,input_mode,output_mode, ...
            n_io,op_overwrite,model);
        WriteBoundaryInp(nb,bctype,bcdim,face,limits);
        WritePhysicsInp_NavierStokes2D(gamma,upw,Prandtl,Reynolds,Minf,grav, ...
                                       rho_ref,p_ref,HB,BV,GasConst);
        WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
        WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
        system(init_exec);
        system('ln -sf reference.bin exact.inp');
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
system('rm -rf INIT reference.bin');

