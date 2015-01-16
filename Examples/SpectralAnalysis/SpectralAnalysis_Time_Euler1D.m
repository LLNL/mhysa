% Spectral Analysis Script: Compute and compare the solution spectrum
% for various time integration schemes and time step sizes
% Case: Density Sum-of-Sine Waves
% Model: 1D Euler

clear all;
close all;

% remove all useless files
system('rm -rf *.dat *.inp *.log INIT');

% Ask for path to HyPar source directory
hypar_path = input('Enter path to HyPar source: ','s');

% Compile the code to generate the exact solution
cmd = ['g++ ',hypar_path, ...
    '/Examples/1D/Euler1D/DensitySumOfSineWaves/aux/exact.C ', ...
    '-o EXACT'];
system(cmd);
% Compile the code to do a Fourier transform of the solution
fft_include =' '; % specify path to FFTW3 include, if not in default path
fft_lib = ' -lfftw3'; % specify path to FFTW3 lib, if not in default path
cmd = ['g++ ',hypar_path, ...
    '/Examples/1D/Euler1D/DensitySumOfSineWaves/aux/fourier.C ', ...
    fft_include,' ',fft_lib,' -o FOURIER'];
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
N = 256;
max_wavenumber = 128;

% specify spatial discretization scheme
hyp_scheme = 'crweno5';
hyp_int_type = 'components';

% specify dt
dt_ex = [0.001 0.002 0.004 0.005];
dt_im = [0.001 0.004 0.01  0.05 ];
t_final = 1.0;

% set physical model and related parameters
model = 'euler1d';
gamma = 1.4;
grav  = 0.0;
upw   = 'roe';

% time integration methods to test
ts = [ ...
        'arkimex'; ...
        'rk     '  ...
     ];
tstype = [
            '3     '; ...
            'ssprk3'  ...
         ];
use_petsc = [1,0];
is_implicit = [1,0];
schemes = 1:size(ts,1);

% specify plotting styles
style = [ '-bo'; '-bs'; '-b^'; '-bd'; '-ro'; '-rs'; '-r^'; '-rd'];

%-------------------------------------------------------------------------%
% if 'arkimex' time-integrator is used, the following options
% can be chosen from:-

% use split hyperbolic flux form (acoustic and entropy modes)?
hyp_flux_split = 'yes';
% treat acoustic waves implicitly, and entropy waves explicitly
hyp_flux_flag  = '-hyperbolic_f_explicit -hyperbolic_df_implicit';
% treat acoustic and entropy waves implicitly
% hyp_flux_flag  = '-hyperbolic_f_implicit -hyperbolic_df_implicit';
% treat acoustic and entropy waves explicitly
% hyp_flux_flag  = '-hyperbolic_f_explicit -hyperbolic_df_explicit';

% or no splitting?
% hyp_flux_split = 'no';
% hyp_flux_flag = '-hyperbolic_implicit'; % implicit treatment
% hyp_flux_flag = '-hyperbolic_explicit'; % explicit treatment
%------------------------------------------------------------------------%

% set solution output to text file
op_format = 'text';
op_overwrite = 'no';
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
exact_exec = './EXACT < max_wavenumber.inp > exact.log 2>&1';
fourier_exec = './FOURIER > fourier.log 2>&1';
clean_exec = 'rm -rf *.inp *.dat *.log';

% open figure windows
scrsz = get(0,'ScreenSize');
figSpectrum = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

count = 0;
max_methods = size(ts,1) * max(size(dt_im,2),size(dt_ex,2));
legend_str = char(zeros(max_methods,23));
legend_str(1,:) = 'Exact                  ';

for j = schemes
    % set PETSc time-integration flags, if the method is a PETSc one
    if (use_petsc(j))
        if (is_implicit(j))
            petsc_flags = sprintf('%s', ...
                '-use-petscts ', ...
                '-ts_type ',strtrim(ts(j,:)),' ', ...
                '-ts_',strtrim(ts(j,:)),'_type ',strtrim(tstype(j,:)),' ', ...
                '-ts_adapt_type none ', ...
                hyp_flux_flag,' ', ...
                '-snes_type newtonls ', ...
                '-snes_rtol 1e-3 ', ...
                '-snes_atol 1e-3 ', ...
                '-ksp_type gmres ', ...
                '-ksp_rtol 1e-3 ', ...
                '-ksp_atol 1e-3 ', ...
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
    
    if (is_implicit(j))
        dt = dt_im;
    else
        dt = dt_ex;
    end
    
    % find out number of dt levels
    ref_levels = max(size(dt));

    % Spectral analysis
    for r = 1:ref_levels
        fprintf('\t%s %2d:  dt=%1.16e\n',[ts(j,:),' ',tstype(j,:)], ...
                r,dt(r));
        niter = int32(t_final/dt(r));
        file_op_iter = niter;
        if (strcmp(petsc_flags,' ')) 
            petscdt = ' ';
            petscft = ' ';
            petscms = ' ';
        else
            petscdt = [' -ts_dt ',num2str(dt(r),'%1.16e'),' '];
            petscft = [' -ts_final_time ',num2str(t_final,'%f'),' '];
            petscms = [' -ts_max_steps ',num2str(100*niter,'%d'),' '];
        end
        % write a file containing the max wavenumber
        fid = fopen('max_wavenumber.inp','w');
        fprintf(fid,'%d\n',max_wavenumber);
        fclose(fid);
        % Write out the input files for HyPar
        WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,strtrim(ts(j,:)), ...
            strtrim(tstype(j,:)),hyp_scheme,hyp_flux_split,hyp_int_type, ...
            par_type,par_scheme, dt(r),cons_check,screen_op_iter, ...
            file_op_iter,op_format,ip_type,input_mode,output_mode, ...
            n_io,op_overwrite,model);
        WriteBoundaryInp(nb,bctype,bcdim,face,limits);
        WritePhysicsInp_Euler1D(gamma,grav,upw);
        WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
        WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
        % Generate the initial and exact solutions
        system(exact_exec);
        % Run HyPar
        hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ', ...
                      hypar,' ',petsc_flags,petscdt,petscft,petscms, ...
                      ' > run.log 2>&1 '];
        system(hypar_exec);
        % create sym links for fourier analysis code
        system('ln -sf op_00000.dat op_exact.dat');
        system('ln -sf op_00001.dat op.dat');
        % run the Fourier analysis code
        system(fourier_exec);
        
        % read in the exact and final spectra
        data = load('spectrum_exact.dat');
        ampl_exact = data(:,2);
        data = load('spectrum_final.dat');
        ampl_comp  = data(:,2);
        wavenumber = data(:,1);
        
        % normalize
        ampl_exact = ampl_exact / sum(ampl_exact);
        ampl_comp  = ampl_comp / sum(ampl_comp);
        
        % Clean up
        system(clean_exec);
        
        % plot the spectrum
        figure(figSpectrum);
        if (~count)
            loglog(wavenumber,ampl_exact,'-k','linewidth',2);
            hold on;
        end
        loglog(wavenumber,ampl_comp,style(count+1,:),'linewidth',1, ...
               'MarkerSize',5);
        hold on;
        name = [sprintf('%7s',ts(j,:)),'-',sprintf('%6s',tstype(j,:)), ...
                ' dt ',sprintf('%5.3f',dt(r))];
        legend_str(count+2,:) = name;
        count = count+1;
    end
end

% embellish the plot
figure(figSpectrum);
xlabel('Wavenumber','FontName','Times','FontSize',20,'FontWeight','normal');
ylabel('Normalized Energy','FontName','Times','FontSize',20,'FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');
legend(legend_str,'Location','SouthWest');
axis([1.0 max_wavenumber-1 ...
      max(min(ampl_comp(1:max_wavenumber-1)),1e-8) 1.0]); 
grid on;
hold off;

% clean up
system('rm -rf EXACT FOURIER');
