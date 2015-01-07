function WritePhysicsInp_NavierStokes2D(gamma,upw,Pr,Re,Minf,grav, ...
                                        rho_ref,p_ref,HB,BV,GasConst)
%WRITEPHYSICSINP Writes the physics.inp file for HyPar 
%                for the 2D Navier-Stokes model

fid = fopen('physics.inp','w');
fprintf(fid,'begin\n');
fprintf(fid,'\tgamma           %1.16e\n',gamma);
fprintf(fid,'\tupwinding       %s\n',upw);
fprintf(fid,'\tPr              %1.16e\n',Pr);
fprintf(fid,'\tRe              %1.16e\n',Re);
fprintf(fid,'\tMinf            %1.16e\n',Minf);
fprintf(fid,'\tgravity         ');
fprintf(fid,'%1.16e ',grav);
fprintf(fid,'\n');
fprintf(fid,'\trho_ref         %1.16e\n',rho_ref);
fprintf(fid,'\tp_ref           %1.16e\n',p_ref);
fprintf(fid,'\tHB              %d',HB);
if (HB == 3)
    fprintf(fid,'\t%1.16e\n',BV);
else
    fprintf(fid,'\n');
end
fprintf(fid,'\tR               %1.16e\n',GasConst);
fprintf(fid,'end\n');
fclose(fid);

end

