function WriteGLMGEEInp(ee_mode)

%WRITEGLMGEEINP Writes the glm_gee.inp file for HyPar

fid = fopen('glm_gee.inp','w');
fprintf(fid,'begin\n');
fprintf(fid,'\tee_mode              %s\n',ee_mode);
fprintf(fid,'end\n');
fclose(fid);

end

