% This script plots of the eigenvalues computed by the 
% ComputeEvals.m script. It saves the plots to a multi-page
% Postscript file ("Eigenvalues.ps").

clear all;
close all;

MaxFiles = 10000;
fname_FFunction_root   = 'Mat_FFunction_';
fname_dFFunction_root  = 'Mat_dFFunction_';
fname_FdFFunction_root = 'Mat_FdFFunction_';
fname_SFunction_root   = 'Mat_SFunction_';
fname_extn = '.dat';
nevals = input('Enter number of eigenvalues computed: ');

scrsz = get(0,'ScreenSize');
width = scrsz(4)/2;
height = 0.9*width;
EigPlot = figure('Position',[1 scrsz(4)/2 width height]);

nplot = 0;
for i=1:MaxFiles
    index = sprintf('%05d',i-1);
    filename_FFunction   = strcat(fname_FFunction_root  ,index,fname_extn);
    filename_dFFunction  = strcat(fname_dFFunction_root ,index,fname_extn);
    filename_FdFFunction = strcat(fname_FdFFunction_root,index,fname_extn);
    filename_SFunction   = strcat(fname_SFunction_root  ,index,fname_extn);
    str_nevals    = strcat(sprintf('%05d',nevals),'_');
    Feval_fname   = strcat('EVals_',str_nevals,filename_FFunction);
    dFeval_fname  = strcat('EVals_',str_nevals,filename_dFFunction);
    FdFeval_fname = strcat('EVals_',str_nevals,filename_FdFFunction);
    Seval_fname   = strcat('EVals_',str_nevals,filename_SFunction);
    
    xmin = -1.0;
    xmax = 0.0;
    ymin = -1.0;
    ymax = 1.0;

    if (exist(Feval_fname,'file'))
        fprintf('Plotting FFunction eigenvalues from  %s\n',Feval_fname);
        dataF  = load(Feval_fname);
        figure(EigPlot);
        plot(dataF(:,2),dataF(:,3),'bo');
        hold on;
        flagFFunction = 1;
        legend_str = '     F(u) ';
        xmin = min(xmin,min(dataF(:,2)));
        xmax = max(xmax,max(dataF(:,2)));
        ymin = min(ymin,min(dataF(:,3)));
        ymax = max(ymax,max(dataF(:,3)));
    else
        flagFFunction = 0;
    end
    if (exist(dFeval_fname,'file'))
        fprintf('Plotting dFFunction eigenvalues from  %s\n',dFeval_fname);
        datadF  = load(dFeval_fname);
        figure(EigPlot);
        plot(datadF(:,2),datadF(:,3),'rx');
        hold on;
        flagdFFunction = 1;
        legend_str = [legend_str;'     dF(u)'];
        xmin = min(xmin,min(datadF(:,2)));
        xmax = max(xmax,max(datadF(:,2)));
        ymin = min(ymin,min(datadF(:,3)));
        ymax = max(ymax,max(datadF(:,3)));
    else
        flagdFFunction = 0;
    end
    if (exist(FdFeval_fname,'file'))
        fprintf('Plotting FdFFunction eigenvalues from  %s\n',FdFeval_fname);
        dataFdF  = load(FdFeval_fname);
        figure(EigPlot);
        plot(dataFdF(:,2),dataFdF(:,3),'ks');
        hold on;
        flagFdFFunction = 1;
        legend_str = [legend_str;'F(u)-dF(u)'];
        xmin = min(xmin,min(dataFdF(:,2)));
        xmax = max(xmax,max(dataFdF(:,2)));
        ymin = min(ymin,min(dataFdF(:,3)));
        ymax = max(ymax,max(dataFdF(:,3)));
    else
        flagFdFFunction = 0;
    end
    if (exist(Seval_fname,'file'))
        fprintf('Plotting FFunction eigenvalues from  %s\n',Seval_fname);
        dataS  = load(Seval_fname);
        figure(EigPlot);
        plot(dataS(:,2),dataS(:,3),'gx');
        hold on;
        flagSFunction = 1;
        legend_str = [legend_str;'      S(u)'];
        xmin = min(xmin,min(dataS(:,2)));
        xmax = max(xmax,max(dataS(:,2)));
        ymin = min(ymin,min(dataS(:,3)));
        ymax = max(ymax,max(dataS(:,3)));
    else
        flagSFunction = 0;
    end
    if (flagFFunction || flagdFFunction || flagFdFFunction || flagSFunction)
        figure(EigPlot);
        set(gca,'FontSize',10,'FontName','Times');
        xlabel('Real(\lambda)','FontName','Times','FontSize',14, ...
            'FontWeight','normal');
        ylabel('Imag(\lambda)','FontName','Times','FontSize',14, ...
            'FontWeight','normal');
        legend(legend_str,'Location','northwest');
        legend('boxoff');
        title(index);
        axis([xmin, xmax, ymin, ymax]);
        grid on;
        hold off;
        figname = 'Eigenvalues.ps';
        if (nplot)
            print('-dpsc2',figname,'-append',EigPlot);
        else
            print('-dpsc2',figname,EigPlot);
        end
        nplot = nplot + 1;
    end
    if (    (~flagFFunction) ...
         && (~flagdFFunction) ...
         && (~flagFdFFunction) ...
         && (~flagSFunction) ...
       )
        break;
    end
end
