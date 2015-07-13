function WritePhysicsInp_ShallowWater2D(gravity,topo_type,upw)
%WRITEPHYSICSINP Writes the physics.inp file for HyPar 
%                for the 2D shallow water model

fid = fopen('physics.inp','w');
fprintf(fid,'begin\n');
fprintf(fid,'\tgravity         %1.16e\n',gravity);
fprintf(fid,'\ttopography_type %d\n',topo_type);
fprintf(fid,'\tupwinding       %s\n',upw);
fprintf(fid,'end\n');
fclose(fid);

end

