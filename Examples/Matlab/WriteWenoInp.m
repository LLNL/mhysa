function WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,tol)

%WRITEWENOINP Writes the weno.inp file for HyPar

fid = fopen('weno.inp','w');
fprintf(fid,'begin\n');
fprintf(fid,'\tmapped               %d\n',mapped);
fprintf(fid,'\tborges               %d\n',borges);
fprintf(fid,'\tyc                   %d\n',yc);
fprintf(fid,'\tno_limiting          %d\n',nl);
fprintf(fid,'\tepsilon              %1.16e\n',eps);
fprintf(fid,'\tp                    %1.16e\n',p);
fprintf(fid,'\trc                   %1.16e\n',rc);
fprintf(fid,'\txi                   %1.16e\n',xi);
fprintf(fid,'\ttol                  %1.16e\n',tol);
fprintf(fid,'end\n');
fclose(fid);

end

