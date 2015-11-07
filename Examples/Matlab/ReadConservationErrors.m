function conserr = ReadConservationErrors(ndims,nvars)

%READERRORDAT Reads the conservation.dat file written by HyPar

conserr = zeros(nvars,1);
fid = fopen('conservation.dat','r');

if (fid == -1)
    % if file doesn't exist, then this simulation
    % likely blew up.
    conserr =  -ones(nvars,1);
else
    fscanf(fid,'%d',ndims);          % read grid sizes
    fscanf(fid,'%d',ndims);          % read number of processors
    fscanf(fid,'%f',1);              % read dt
    % now read data to return
    conserr = fscanf(fid,'%f',nvars);% conservation errors
    % done
    fclose(fid);
end

end

