function WriteBoundaryInp(nb,bctype,dim,face,limits,varargin)

%WRITEBOUNDARYINP Writes out the boundary.inp file for HyPar

fileID = fopen('boundary.inp','w');
fprintf(fileID,'%d\n',nb);
count = 1;
for i=1:nb
    fprintf(fileID,'%-25s  ',strtrim(bctype(i,:)));
    fprintf(fileID,'%4d    ',dim(i));
    fprintf(fileID,'%4d    ',face(i));
    fprintf(fileID,'%1.16e ',limits(i,:));
    fprintf(fileID,'\n');
    if (     strcmp(strtrim(bctype(i,:)),'dirichlet'                   ) ...
          || strcmp(strtrim(bctype(i,:)),'sponge'                      ) ...
          || strcmp(strtrim(bctype(i,:)),'slip-wall'                   ) ...
          || strcmp(strtrim(bctype(i,:)),'noslip-wall'                 ) ...
          || strcmp(strtrim(bctype(i,:)),'subsonic-inflow'             ) ...
          || strcmp(strtrim(bctype(i,:)),'subsonic-outflow'            ) ...
          || strcmp(strtrim(bctype(i,:)),'supersonic-inflow'           ) ...
          || strcmp(strtrim(bctype(i,:)),'turbulent-supersonic-inflow' ) ...
       )
        fprintf(fileID,'%s\n',varargin{count});
        count = count+1;
    end
end
fclose(fileID);

end

