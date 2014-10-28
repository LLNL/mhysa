function [ndims,nvars,size,iproc,ghost,niter,ts,tstype,hyp_scheme, ...
    hyp_flux_split,hyp_int_type,par_type,par_scheme,dt,cons_check, ...
    screen_op_iter,file_op_iter,op_format, ip_type,input_mode, ...
    output_mode,n_io,op_overwrite,model,nb,bctype,dim,face,limits, ...
    mapped,borges,yc,nl,eps,p,rc,xi,tol,lutype,norm,maxiter,atol,rtol, ...
    verbose] = SetDefaults()

%SETDEFAULTS Defaut values of all the input parameters for HyPar

% solver.inp
ndims   = 1;
nvars   = 1;
size    = 100;
iproc   = 1;
ghost   = 3;
niter   = 200;
ts      = 'rk';
tstype  = 'ssprk3';
hyp_scheme = 'weno5';
hyp_flux_split = 'no';
hyp_int_type = 'components';
par_type = 'nonconservative-1.5stage';
par_scheme = '4';
dt = 0.005;
cons_check = 'yes';
screen_op_iter = 1;
file_op_iter = 99999999;
op_format = 'none';
ip_type = 'binary';
input_mode = 'serial';
output_mode = 'serial';
n_io = 1;
op_overwrite = 'yes';
model = 'linear-advection-diffusion-reaction';

% boundary.inp
nb = 2;
bctype = ['periodic';'periodic'];
dim = [0;0];
face = [1,-1];
limits = [0 0; 0 0];

% weno.inp
mapped = 0;
borges = 0;
yc = 1;
nl = 0;
eps = 0.000001;
p = 2.0;
rc = 0.3;
xi = 0.001;
tol= 1e-16;

% lusolver.inp
lutype = 'jacobi';
norm = 1;
maxiter = 10;
atol = 1e-12;
rtol = 1e-10;
verbose = 0;

end

