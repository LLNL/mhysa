function WritePhysicsInp_Euler1D(gamma,upw,nspecies,nvibeng)
%WRITEPHYSICSINP Writes the physics.inp file for Mhysa 
%                for the 1D Euler model

fid = fopen('physics.inp','w');
fprintf(fid,'begin\n');
fprintf(fid,'\tgamma           %1.16e\n',gamma);
fprintf(fid,'\tnspecies        %d\n',nspecies);
fprintf(fid,'\tnvibeng         %d\n',nvibeng);
fprintf(fid,'\tupwinding       %s\n',upw);
fprintf(fid,'end\n');
fclose(fid);

end

