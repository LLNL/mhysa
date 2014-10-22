function WriteLusolverInp(type,norm,maxiter,atol,rtol,verbose)

%WRITELUSOLVERINP Writes the lusolver.inp file for HyPar

fid = fopen('lusolver.inp','w');
fprintf(fid,'begin\n');
fprintf(fid,'\treducedsolvetype            %s\n',type);
fprintf(fid,'\tevaluate_norm               %d\n',norm);
fprintf(fid,'\tmaxiter                     %d\n',maxiter);
fprintf(fid,'\tatol                        %1.16e\n',atol);
fprintf(fid,'\trtol                        %1.16e\n',rtol);
fprintf(fid,'\tverbose                     %d\n',verbose);
fprintf(fid,'end\n');
fclose(fid);

end

