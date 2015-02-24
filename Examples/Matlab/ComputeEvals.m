% This script computes the eigenvalues of the linearized 
% matrix representations of the right-hand-side functions.
% 
% The code needs to be compiled with the "-Dcompute_rhs_operators"
% flag so that it writes out files that contain the matrices
% in sparse format during a simulation.
% 
% The eigenvalues are written out to text files and can be plotted
% in MATLAB with the script PlotEvals.m (or any other utility).

clear all;

MaxFiles = 10000;
fname_FFunction_root   = 'Mat_FFunction_';
fname_dFFunction_root  = 'Mat_dFFunction_';
fname_FdFFunction_root = 'Mat_FdFFunction_';
fname_SFunction_root   = 'Mat_SFunction_';
fname_extn = '.dat';
nevals = input('Enter number of eigenvalues to compute: ');
opts.tol = 1e-10;

for i=1:MaxFiles
    index = sprintf('%05d',i-1);
    % compute eigenvalues for FFunction matrix, if available
    filename_FFunction = strcat(fname_FFunction_root,index,fname_extn);
    if (exist(filename_FFunction,'file'))
        fprintf('Reading FFunction matrix from %s: ',filename_FFunction);
        fid = fopen(filename_FFunction,'r');
        ndof = fscanf(fid,'%d',1);
        fprintf('ndof = %d, ',ndof);
        A = zeros(ndof,ndof);
        nnz = 0;
        while (~feof(fid))
            coord = fscanf(fid,'%d',2);
            if (max(size(coord)) > 0)
                nnz = nnz + 1;
                A(coord(1),coord(2)) = fscanf(fid,'%f',1);
            end
        end
        fprintf('nnz = %d.\n',nnz);
        fclose(fid);
        FFunction_Mat = sparse(A);
        
        str_nevals   = strcat(sprintf('%05d',nevals),'_');
        FFunction_eval_fname  = strcat('EVals_',str_nevals,filename_FFunction);

        if (~exist(FFunction_eval_fname,'file'))
            fprintf('  Computing %5d eigenvalues of FFunction matrix... ',nevals);
            tic;
            if (nevals < ndof-1)
                lambdaFFunction = eigs(FFunction_Mat,nevals,'lm',opts);
            else
                lambdaFFunction = eig(full(FFunction_Mat));
            end
            waqt = toc;
            fprintf('%f seconds.\n',waqt);
            fprintf('  Saving eigenvalues to file %s.\n',FFunction_eval_fname);
            fid = fopen(FFunction_eval_fname,'w');
            for n = 1:nevals
                fprintf(fid,'%5d %+1.16e %+1.16e\n', ...
                        n,real(lambdaFFunction(n:n)), ...
                        imag(lambdaFFunction(n:n)));
            end
            fclose(fid);
        else
            fprintf('  Skipping eigenvalues of FFunction matrix, file %s already exists.\n', ...
                    FFunction_eval_fname);
        end
        FlagFFunction = 1;
    else
        FlagFFunction = 0;
    end
    % compute eigenvalues for dFFunction matrix, if available
    filename_dFFunction = strcat(fname_dFFunction_root,index,fname_extn);
    if (exist(filename_dFFunction,'file'))
        fprintf('Reading dFFunction matrix from %s: ',filename_dFFunction);
        fid = fopen(filename_dFFunction,'r');
        ndof = fscanf(fid,'%d',1);
        fprintf('ndof = %d, ',ndof);
        A = zeros(ndof,ndof);
        nnz = 0;
        while (~feof(fid))
            coord = fscanf(fid,'%d',2);
            if (max(size(coord)) > 0)
                nnz = nnz + 1;
                A(coord(1),coord(2)) = fscanf(fid,'%f',1);
            end
        end
        fprintf('nnz = %d.\n',nnz);
        fclose(fid);
        dFFunction_Mat = sparse(A);
        
        str_nevals   = strcat(sprintf('%05d',nevals),'_');
        dFFunction_eval_fname  = strcat('EVals_',str_nevals,filename_dFFunction);

        if (~exist(dFFunction_eval_fname,'file'))
            fprintf('  Computing %5d eigenvalues of dFFunction matrix... ',nevals);
            tic;
            if (nevals < ndof-1)
                lambdadFFunction = eigs(dFFunction_Mat,nevals,'lm',opts);
            else
                lambdadFFunction = eig(full(dFFunction_Mat));
            end
            waqt = toc;
            fprintf('%f seconds.\n',waqt);
            fprintf('  Saving eigenvalues to file %s.\n',dFFunction_eval_fname);
            fid = fopen(dFFunction_eval_fname,'w');
            for n = 1:nevals
                fprintf(fid,'%5d %+1.16e %+1.16e\n', ...
                        n,real(lambdadFFunction(n:n)), ...
                        imag(lambdadFFunction(n:n)));
            end
            fclose(fid);
        else
            fprintf('  Skipping eigenvalues of dFFunction matrix, file %s already exists.\n', ...
                    dFFunction_eval_fname);
        end
        FlagdFFunction = 1;
    else
        FlagdFFunction = 0;
    end
    % compute eigenvalues for (FFunction-dFFunction matrix), if available
    filename_FdFFunction = strcat(fname_FdFFunction_root,index,fname_extn);
    if (exist(filename_FdFFunction,'file'))
        fprintf('Reading FdFFunction matrix from %s: ',filename_FdFFunction);
        fid = fopen(filename_FdFFunction,'r');
        ndof = fscanf(fid,'%d',1);
        fprintf('ndof = %d, ',ndof);
        A = zeros(ndof,ndof);
        nnz = 0;
        while (~feof(fid))
            coord = fscanf(fid,'%d',2);
            if (max(size(coord)) > 0)
                nnz = nnz + 1;
                A(coord(1),coord(2)) = fscanf(fid,'%f',1);
            end
        end
        fprintf('nnz = %d.\n',nnz);
        fclose(fid);
        FdFFunction_Mat = sparse(A);
        
        str_nevals   = strcat(sprintf('%05d',nevals),'_');
        FdFFunction_eval_fname  = strcat('EVals_',str_nevals,filename_FdFFunction);

        if (~exist(FdFFunction_eval_fname,'file'))
            fprintf('  Computing %5d eigenvalues of FdFFunction matrix... ',nevals);
            tic;
            if (nevals < ndof-1)
                lambdaFdFFunction = eigs(FdFFunction_Mat,nevals,'lm',opts);
            else
                lambdaFdFFunction = eig(full(FdFFunction_Mat));
            end
            waqt = toc;
            fprintf('%f seconds.\n',waqt);
            fprintf('  Saving eigenvalues to file %s.\n',FdFFunction_eval_fname);
            fid = fopen(FdFFunction_eval_fname,'w');
            for n = 1:nevals
                fprintf(fid,'%5d %+1.16e %+1.16e\n', ...
                        n,real(lambdaFdFFunction(n:n)), ...
                        imag(lambdaFdFFunction(n:n)));
            end
            fclose(fid);
        else
            fprintf('  Skipping eigenvalues of dFFunction matrix, file %s already exists.\n', ...
                    FdFFunction_eval_fname);
        end
        FlagFdFFunction = 1;
    else
        FlagFdFFunction = 0;
    end
    % compute eigenvalues for SFunction matrix, if available
    filename_SFunction = strcat(fname_SFunction_root,index,fname_extn);
    if (exist(filename_SFunction,'file'))
        fprintf('Reading SFunction matrix from %s: ',filename_SFunction);
        fid = fopen(filename_SFunction,'r');
        ndof = fscanf(fid,'%d',1);
        fprintf('ndof = %d, ',ndof);
        A = zeros(ndof,ndof);
        nnz = 0;
        while (~feof(fid))
            coord = fscanf(fid,'%d',2);
            if (max(size(coord)) > 0)
                nnz = nnz + 1;
                A(coord(1),coord(2)) = fscanf(fid,'%f',1);
            end
        end
        fprintf('nnz = %d.\n',nnz);
        fclose(fid);
        SFunction_Mat = sparse(A);
        
        str_nevals   = strcat(sprintf('%05d',nevals),'_');
        SFunction_eval_fname  = strcat('EVals_',str_nevals,filename_SFunction);

        if (~exist(SFunction_eval_fname,'file'))
            fprintf('  Computing %5d eigenvalues of SFunction matrix... ',nevals);
            tic;
            if (nevals < ndof-1)
                lambdaSFunction = eigs(SFunction_Mat,nevals,'lm',opts);
            else
                lambdaSFunction = eig(full(SFunction_Mat));
            end
            waqt = toc;
            fprintf('%f seconds.\n',waqt);
            fprintf('  Saving eigenvalues to file %s.\n',SFunction_eval_fname);
            fid = fopen(SFunction_eval_fname,'w');
            for n = 1:nevals
                fprintf(fid,'%5d %+1.16e %+1.16e\n', ...
                        n,real(lambdaSFunction(n:n)), ...
                        imag(lambdaSFunction(n:n)));
            end
            fclose(fid);
        else
            fprintf('  Skipping eigenvalues of SFunction matrix, file %s already exists.\n', ...
                    SFunction_eval_fname);
        end
        FlagSFunction = 1;
    else
        FlagSFunction = 0;
    end
    if (    (~FlagFFunction) ...
         && (~FlagdFFunction) ...
         && (~FlagFdFFunction) ...
         && (~FlagSFunction) ...
       )
        break;
    end
end
