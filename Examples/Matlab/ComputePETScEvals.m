% This script computes the eigenvalues of the linearized 
% matrix representations of the right-hand-side functions,
% when PETSc time integration is used.
% 
% The code needs to be compiled with the "-Dcompute_rhs_operators"
% flag so that it writes out files that contain the matrices
% in sparse format during a simulation.
% 
% The eigenvalues are written out to text files and can be plotted
% in MATLAB with the script PlotPETScEvals.m (or any other utility).

clear all;

MaxFiles = 10000;
fname_RHSFunctionExpl_root   = 'Mat_RHSFunctionExpl_';
fname_RHSFunctionIMEX_root   = 'Mat_RHSFunctionIMEX_';
fname_IFunctionIMEX_root     = 'Mat_IFunctionIMEX_';
fname_extn = '.dat';
nevals = input('Enter number of eigenvalues to compute: ');
opts.tol = 1e-10;

for i=1:MaxFiles
    index = sprintf('%05d',i-1);
    % compute eigenvalues for RHSFunctionExpl matrix, if available
    filename_RHSFunctionExpl = strcat(fname_RHSFunctionExpl_root,index,fname_extn);
    if (exist(filename_RHSFunctionExpl,'file'))
        fprintf('Reading RHSFunctionExpl matrix from %s: ',filename_RHSFunctionExpl);
        fid = fopen(filename_RHSFunctionExpl,'r');
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
        RHSFunctionExpl_Mat = sparse(A);
        
        str_nevals   = strcat(sprintf('%05d',nevals),'_');
        RHSFunctionExpl_eval_fname  = strcat('EVals_',str_nevals,filename_RHSFunctionExpl);

        if (~exist(RHSFunctionExpl_eval_fname,'file'))
            fprintf('  Computing %5d eigenvalues of RHSFunctionExpl matrix... ',nevals);
            tic;
            if (nevals < ndof-1)
                lambdaRHSFunctionExpl = eigs(RHSFunctionExpl_Mat,nevals,'lm',opts);
            else
                lambdaRHSFunctionExpl = eig(full(RHSFunctionExpl_Mat));
            end
            waqt = toc;
            fprintf('%f seconds.\n',waqt);
            fprintf('  Saving eigenvalues to file %s.\n',RHSFunctionExpl_eval_fname);
            fid = fopen(RHSFunctionExpl_eval_fname,'w');
            for n = 1:nevals
                fprintf(fid,'%5d %+1.16e %+1.16e\n', ...
                        n,real(lambdaRHSFunctionExpl(n:n)), ...
                        imag(lambdaRHSFunctionExpl(n:n)));
            end
            fclose(fid);
        else
            fprintf('  Skipping eigenvalues of RHSFunctionExpl matrix, file %s already exists.\n', ...
                    RHSFunctionExpl_eval_fname);
        end
        FlagRHSFunctionExpl = 1;
    else
        FlagRHSFunctionExpl = 0;
    end
    % compute eigenvalues for RHSFunctionIMEX matrix, if available
    filename_RHSFunctionIMEX = strcat(fname_RHSFunctionIMEX_root,index,fname_extn);
    if (exist(filename_RHSFunctionIMEX,'file'))
        fprintf('Reading RHSFunctionIMEX matrix from %s: ',filename_RHSFunctionIMEX);
        fid = fopen(filename_RHSFunctionIMEX,'r');
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
        RHSFunctionIMEX_Mat = sparse(A);
        
        str_nevals   = strcat(sprintf('%05d',nevals),'_');
        RHSFunctionIMEX_eval_fname  = strcat('EVals_',str_nevals,filename_RHSFunctionIMEX);

        if (~exist(RHSFunctionIMEX_eval_fname,'file'))
            fprintf('  Computing %5d eigenvalues of RHSFunctionIMEX matrix... ',nevals);
            tic;
            if (nevals < ndof-1)
                lambdaRHSFunctionIMEX = eigs(RHSFunctionIMEX_Mat,nevals,'lm',opts);
            else
                lambdaRHSFunctionIMEX = eig(full(RHSFunctionIMEX_Mat));
            end
            waqt = toc;
            fprintf('%f seconds.\n',waqt);
            fprintf('  Saving eigenvalues to file %s.\n',RHSFunctionIMEX_eval_fname);
            fid = fopen(RHSFunctionIMEX_eval_fname,'w');
            for n = 1:nevals
                fprintf(fid,'%5d %+1.16e %+1.16e\n', ...
                        n,real(lambdaRHSFunctionIMEX(n:n)), ...
                        imag(lambdaRHSFunctionIMEX(n:n)));
            end
            fclose(fid);
        else
            fprintf('  Skipping eigenvalues of RHSFunctionIMEX matrix, file %s already exists.\n', ...
                    RHSFunctionIMEX_eval_fname);
        end
        FlagRHSFunctionIMEX = 1;
    else
        FlagRHSFunctionIMEX = 0;
    end
    % compute eigenvalues for IFunctionIMEX matrix, if available
    filename_IFunctionIMEX = strcat(fname_IFunctionIMEX_root,index,fname_extn);
    if (exist(filename_IFunctionIMEX,'file'))
        fprintf('Reading IFunctionIMEX matrix from %s: ',filename_IFunctionIMEX);
        fid = fopen(filename_IFunctionIMEX,'r');
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
        IFunctionIMEX_Mat = sparse(A);
        
        str_nevals   = strcat(sprintf('%05d',nevals),'_');
        IFunctionIMEX_eval_fname  = strcat('EVals_',str_nevals,filename_IFunctionIMEX);

        if (~exist(IFunctionIMEX_eval_fname,'file'))
            fprintf('  Computing %5d eigenvalues of IFunctionIMEX matrix... ',nevals);
            tic;
            if (nevals < ndof-1)
                lambdaIFunctionIMEX = eigs(IFunctionIMEX_Mat,nevals,'lm',opts);
            else
                lambdaIFunctionIMEX = eig(full(IFunctionIMEX_Mat));
            end
            waqt = toc;
            fprintf('%f seconds.\n',waqt);
            fprintf('  Saving eigenvalues to file %s.\n',IFunctionIMEX_eval_fname);
            fid = fopen(IFunctionIMEX_eval_fname,'w');
            for n = 1:nevals
                fprintf(fid,'%5d %+1.16e %+1.16e\n', ...
                        n,real(lambdaIFunctionIMEX(n:n)), ...
                        imag(lambdaIFunctionIMEX(n:n)));
            end
            fclose(fid);
        else
            fprintf('  Skipping eigenvalues of IFunctionIMEX matrix, file %s already exists.\n', ...
                    IFunctionIMEX_eval_fname);
        end
        FlagIFunctionIMEX = 1;
    else
        FlagIFunctionIMEX = 0;
    end
    if (    (~FlagRHSFunctionExpl) ...
         && (~FlagRHSFunctionIMEX) ...
         && (~FlagIFunctionIMEX) ...
       )
        break;
    end
end
