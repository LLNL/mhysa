function [err, wt] = ReadGLMGEEErrorDat(ndims)

%READGLMGEEERRORDAT Reads the glm_err.dat file, if written by HyPar

fid = fopen('glm_err.dat','r');

if (fid == -1)
    % if file doesn't exist, then this simulation
    % likely blew up.
    err = [inf,inf,inf];
    wt  = [0,0];
else
    fscanf(fid,'%f',1);              % read dt
    % now read data to return
    err = fscanf(fid,'%f',3);        % L1,L2 and Linf error
    % done
    fclose(fid);
end

end

