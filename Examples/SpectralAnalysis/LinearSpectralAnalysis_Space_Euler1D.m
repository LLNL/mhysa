% Spectral Analysis Script: Compute and compare the solution spectrum
% for various spatial discretizatio schemes
% Case: Density Sum-of-Sine Waves
% Model: 1D Euler

clear all;
close all;

% remove all useless files
system('rm -rf *.dat *.inp *.log EXACT');

% Ask for path to HyPar source directory
hypar_path = input('Enter path to HyPar source: ','s');

% Compile the code to generate the exact solution
cmd = ['gcc ',hypar_path, ...
    '/Examples/1D/Euler1D/DensitySumOfSineWaves/aux/exact.c -lm ', ...
    '-o EXACT'];
system(cmd);
% Compile the code to do a Fourier transform of the solution
fft_include =' '; % specify path to FFTW3 include, if not in default path
fft_lib = ' -lfftw3'; % specify path to FFTW3 lib, if not in default path
cmd = ['gcc ',hypar_path, ...
    '/Examples/1D/Euler1D/DensitySumOfSineWaves/aux/fourier.c -lm ', ...
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
hyp_schemes = [ ...
                'muscl3 '; ...
                'weno5  '; ...
                'crweno5'  ...
              ];
hyp_int_type = 'components';
schemes = 1:size(hyp_schemes,1);

% specify dt and final time
dt = 0.001;
t_final = 1.0;
niter = int32(t_final/dt);

% set physical model and related parameters
model = 'euler1d';
gamma = 1.4;
grav  = 0.0;
upw   = 'roe';

% time integration methods to test
ts = 'rk';
tstype = '44';

% specify plotting styles
style = [ '-go'; '-rs'; '-b^'];

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
file_op_iter = niter;
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
succ_count = 0;
max_methods = size(hyp_schemes,1);
legend_str = char(zeros(max_methods,7));
legend_str(1,:) = 'Exact  ';

ymin = 1.0;
ymax = 1e-16;
for j = schemes
    % set PETSc time-integration flags, if the method is a PETSc one
    if (strcmp(ts,'arkimex'))
        petsc_flags = sprintf('%s', ...
            '-use-petscts ', ...
            '-ts_type ',strtrim(ts),' ', ...
            '-ts_',strtrim(ts),'_type ',strtrim(tstype),' ', ...
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
        petsc_flags = ' ';
    end
    
    % Spectral analysis
    hyp_scheme = hyp_schemes(j,:);
    fprintf('\t%s: ',hyp_scheme);
    if (strcmp(petsc_flags,' ')) 
        petscdt = ' ';
        petscft = ' ';
        petscms = ' ';
    else
        petscdt = [' -ts_dt ',num2str(dt,'%1.16e'),' '];
        petscft = [' -ts_final_time ',num2str(t_final,'%f'),' '];
        petscms = [' -ts_max_steps ',num2str(100*niter,'%d'),' '];
    end
    % write a file containing the max wavenumber
    fid = fopen('max_wavenumber.inp','w');
    fprintf(fid,'%d\n',max_wavenumber);
    fclose(fid);
    % Write out the input files for HyPar
    WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,ts, ...
        tstype,strtrim(hyp_scheme),hyp_flux_split,hyp_int_type, ...
        par_type,par_scheme, dt,cons_check,screen_op_iter, ...
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
    
    % plot exact spectrum if this is first loop
    if (~count)
        figure(figSpectrum);
        loglog(wavenumber(1:max_wavenumber),ampl_exact(1:max_wavenumber), ...
               '-k','linewidth',2);
        hold on;
    end
    
    if (min(isfinite(ampl_comp)))
        fprintf('done.\n');
        % plot the spectrum
        figure(figSpectrum);
        loglog(wavenumber(1:max_wavenumber),ampl_comp(1:max_wavenumber), ...
            style(count+1,:),'linewidth',1, 'MarkerSize',5);
        hold on;
        legend_str(succ_count+2,:) = hyp_scheme;
        succ_count = succ_count + 1;
        
        ymin = min(ymin,min(ampl_comp(1:max_wavenumber)));
        ymax = max(ymax,max(ampl_comp(1:max_wavenumber)));
    else
        fprintf('solution failed.\n');
    end
    count = count+1;
end

% embellish the plot
figure(figSpectrum);
xlabel('Wavenumber','FontName','Times','FontSize',20,'FontWeight','normal');
ylabel('Normalized Energy','FontName','Times','FontSize',20,'FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');
legend(legend_str(1:succ_count+1,:),'Location','SouthWest');
axis([1.0 max_wavenumber max(ymin,1e-12) ymax]); 
grid on;
hold off;

% clean up
system('rm -rf EXACT FOURIER');
