function WriteBoundaryInp(nb,bctype,dim,face,limits)

%WRITEBOUNDARYINP Writes out the boundary.inp file for HyPar

fileID = fopen('boundary.inp','w');
fprintf(fileID,'%d\n',nb);
for i=1:nb
    fprintf(fileID,'%-25s  ',strtrim(bctype(i,:)));
    fprintf(fileID,'%4d    ',dim(i));
    fprintf(fileID,'%4d    ',face(i));
    fprintf(fileID,'%1.16e ',limits(i,:));
    fprintf(fileID,'\n');
end
fclose(fileID);

end

