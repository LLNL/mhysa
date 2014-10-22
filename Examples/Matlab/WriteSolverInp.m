function WriteSolverInp(ndims,nvars,size,iproc,ghost,niter, ...
    ts,tstype,hyp_scheme,hyp_flux_split,hyp_int_type,par_type, ...
    par_scheme,dt,cons_check,screen_op_iter,file_op_iter,op_format, ...
    ip_type,input_mode,output_mode,n_io,op_overwrite,model)

%WRITEHYPARINPUTFILES Writes out the solver.inp file for HyPar

fileID = fopen('solver.inp','w');
fprintf(fileID,'begin\n');
fprintf(fileID,'\tndims              %d\n',ndims);
fprintf(fileID,'\tnvars              %d\n',nvars);
fprintf(fileID,'\tsize               ');
fprintf(fileID,'%d ',size);
fprintf(fileID,'\n');
fprintf(fileID,'\tiproc              ');
fprintf(fileID,'%d ',iproc);
fprintf(fileID,'\n');
fprintf(fileID,'\tghost              %d\n',ghost);
fprintf(fileID,'\tn_iter             %d\n',niter);
fprintf(fileID,'\trestart_iter       0 \n');
fprintf(fileID,'\ttime_scheme        %s\n',ts);
fprintf(fileID,'\ttime_scheme_type   %s\n',tstype);
fprintf(fileID,'\thyp_space_scheme   %s\n',hyp_scheme);
fprintf(fileID,'\thyp_flux_split     %s\n',hyp_flux_split);
fprintf(fileID,'\thyp_interp_type    %s\n',hyp_int_type);
fprintf(fileID,'\tpar_space_type     %s\n',par_type);
fprintf(fileID,'\tpar_space_scheme   %s\n',par_scheme);
fprintf(fileID,'\tdt                 %1.16e\n',dt);
fprintf(fileID,'\tconservation_check %s\n',cons_check);
fprintf(fileID,'\tscreen_op_iter     %d\n',screen_op_iter);
fprintf(fileID,'\tfile_op_iter       %d\n',file_op_iter);
fprintf(fileID,'\top_file_format     %s\n',op_format);
fprintf(fileID,'\tip_file_type       %s\n',ip_type);
fprintf(fileID,'\tinput_mode         %s'  ,input_mode);
if (~strcmp(input_mode,'serial')) 
    fprintf(fileID,' %d\n',n_io);
else
    fprintf(fileID,'\n');
end
fprintf(fileID,'\toutput_mode        %s'  ,output_mode);
if (~strcmp(output_mode,'serial')) 
    fprintf(fileID,' %d\n',n_io);
else
    fprintf(fileID,'\n');
end
fprintf(fileID,'\top_overwrite       %s\n',op_overwrite);
fprintf(fileID,'\tmodel              %s\n',model);
fprintf(fileID,'end\n');
fclose(fileID);

end

