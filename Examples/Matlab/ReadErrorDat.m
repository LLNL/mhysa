function [err, wt] = ReadErrorDat(ndims)

%READERRORDAT Reads the errors.dat file written by HyPar

fid = fopen('errors.dat','r');

if (fid == -1)
    % if file doesn't exist, then this simulation
    % likely blew up.
    err = [inf,inf,inf];
    wt  = [0,0];
else
    fscanf(fid,'%d',ndims);          % read grid sizes
    fscanf(fid,'%d',ndims);          % read number of processors
    fscanf(fid,'%f',1);              % read dt
    % now read data to return
    err = fscanf(fid,'%f',3);        % L1,L2 and Linf error
    wt  = fscanf(fid,'%f',2);        % solver and total wall times
    % done
    fclose(fid);
end

end

