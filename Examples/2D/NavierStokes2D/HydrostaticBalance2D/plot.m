clear all;
clf;

rho0 = 1.21;
p0   = 1.0;
g    = 1.0;

data = load('op.dat');
x = data(1:51,3);
y = data(1:51:2601,4);
rho = reshape(data(:,5),51,51);
u   = reshape(data(:,6),51,51) ./ rho;
v   = reshape(data(:,7),51,51) ./ rho;
e   = reshape(data(:,8),51,51);

dp = zeros(51,51);
p  = zeros(51,51);
for i = 1:51
    for j = 1:51
        p(i,j)  = 0.4 * (e(i,j) - 0.5*rho(i,j)*(u(i,j)^2+v(i,j)^2));
        dp(i,j) = p(i,j) - p0 * exp( -(rho0*g/p0) * (x(i)+y(j)) );
    end
end

surf(x,y,dp);
view(45,35);