function WritePhysicsInp_Euler1D(gamma,grav,upw)
%WRITEPHYSICSINP Writes the physics.inp file for HyPar 
%                for the 1D Euler model

fid = fopen('physics.inp','w');
fprintf(fid,'begin\n');
fprintf(fid,'\tgamma           %1.16e\n',gamma);
fprintf(fid,'\tgravity         %1.16e\n',grav);
fprintf(fid,'\tupwinding       %s\n',upw);
fprintf(fid,'end\n');
fclose(fid);

end

