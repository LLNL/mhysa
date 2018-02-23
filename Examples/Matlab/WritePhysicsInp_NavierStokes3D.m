function WritePhysicsInp_NavierStokes3D(gamma,upw,Pr,Re,Minf,nspecies,nvibeng)
%WRITEPHYSICSINP Writes the physics.inp file for HyPar 
%                for the 2D Navier-Stokes model

fid = fopen('physics.inp','w');
fprintf(fid,'begin\n');
fprintf(fid,'\tgamma           %1.16e\n',gamma);
fprintf(fid,'\tupwinding       %s\n',upw);
fprintf(fid,'\tnspecies        %d\n',nspecies);
fprintf(fid,'\tnvibeng         %d\n',nvibeng );
fprintf(fid,'\tPr              %1.16e\n',Pr);
fprintf(fid,'\tRe              %1.16e\n',Re);
fprintf(fid,'\tMinf            %1.16e\n',Minf);
fprintf(fid,'end\n');
fclose(fid);

end

