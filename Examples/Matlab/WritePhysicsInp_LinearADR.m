function WritePhysicsInp_LinearADR(a)
%WRITEPHYSICSINP Writes the physics.inp file for HyPar 
%                for the linear advection model

fid = fopen('physics.inp','w');
fprintf(fid,'begin\n');
fprintf(fid,'\tadvection         ');
fprintf(fid,'%1.16e ',a);
fprintf(fid,'\n');
fprintf(fid,'end\n');
fclose(fid);

end

