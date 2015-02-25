% This script plots of the eigenvalues computed by the 
% ComputePETScEvals.m script. It saves the plots to a multi-page
% Postscript file ("Eigenvalues.ps").

clear all;
close all;

MaxFiles = 10000;
fname_RHSFunctionExpl_root = 'Mat_RHSFunctionExpl_';
fname_RHSFunctionIMEX_root = 'Mat_RHSFunctionIMEX_';
fname_IFunctionIMEX_root   = 'Mat_IFunctionIMEX_';
fname_extn = '.dat';
nevals = input('Enter number of eigenvalues computed: ');

scrsz = get(0,'ScreenSize');
width = scrsz(4)/2;
height = 0.9*width;
EigPlot = figure('Position',[1 scrsz(4)/2 width height]);

nplot = 0;
for i=1:MaxFiles
    index = sprintf('%05d',i-1);
    filename_RHSFunctionExpl = strcat(fname_RHSFunctionExpl_root,index,fname_extn);
    filename_RHSFunctionIMEX = strcat(fname_RHSFunctionIMEX_root,index,fname_extn);
    filename_IFunctionIMEX   = strcat(fname_IFunctionIMEX_root  ,index,fname_extn);
    str_nevals    = strcat(sprintf('%05d',nevals),'_');
    RHSFunctionExpl_eval_fname = strcat('EVals_',str_nevals,filename_RHSFunctionExpl);
    RHSFunctionIMEX_eval_fname = strcat('EVals_',str_nevals,filename_RHSFunctionIMEX);
    IFunctionIMEX_eval_fname   = strcat('EVals_',str_nevals,filename_IFunctionIMEX  );
    
    xmin = -1.0;
    xmax = 0.0;
    ymin = -1.0;
    ymax = 1.0;

    if (exist(RHSFunctionExpl_eval_fname,'file'))
        fprintf('Plotting RHSFunctionExpl eigenvalues from  %s\n',RHSFunctionExpl_eval_fname);
        dataRHSFunctionExpl  = load(RHSFunctionExpl_eval_fname);
        figure(EigPlot);
        plot(dataRHSFunctionExpl(:,2),dataRHSFunctionExpl(:,3),'bs');
        hold on;
        flagRHSFunctionExpl = 1;
        legend_str = 'F(u) Explicit';
        xmin = min(xmin,min(dataRHSFunctionExpl(:,2)));
        xmax = max(xmax,max(dataRHSFunctionExpl(:,2)));
        ymin = min(ymin,min(dataRHSFunctionExpl(:,3)));
        ymax = max(ymax,max(dataRHSFunctionExpl(:,3)));
    else
        flagRHSFunctionExpl = 0;
    end

    if (exist(IFunctionIMEX_eval_fname,'file'))
        fprintf('Plotting IFunctionIMEX eigenvalues from  %s\n',IFunctionIMEX_eval_fname);
        dataIFunctionIMEX  = load(IFunctionIMEX_eval_fname);
        figure(EigPlot);
        plot(dataIFunctionIMEX(:,2),dataIFunctionIMEX(:,3),'rx');
        hold on;
        flagIFunctionIMEX = 1;
        legend_str = [legend_str;'G(u) IMEX    '];
        xmin = min(xmin,min(dataIFunctionIMEX(:,2)));
        xmax = max(xmax,max(dataIFunctionIMEX(:,2)));
        ymin = min(ymin,min(dataIFunctionIMEX(:,3)));
        ymax = max(ymax,max(dataIFunctionIMEX(:,3)));
    else
        flagIFunctionIMEX = 0;
    end

    if (exist(RHSFunctionIMEX_eval_fname,'file'))
        fprintf('Plotting RHSFunctionIMEX eigenvalues from  %s\n',RHSFunctionIMEX_eval_fname);
        dataRHSFunctionIMEX  = load(RHSFunctionIMEX_eval_fname);
        figure(EigPlot);
        plot(dataRHSFunctionIMEX(:,2),dataRHSFunctionIMEX(:,3),'kx');
        hold on;
        flagRHSFunctionIMEX = 1;
        legend_str = [legend_str;'F(u) IMEX    '];
        xmin = min(xmin,min(dataRHSFunctionIMEX(:,2)));
        xmax = max(xmax,max(dataRHSFunctionIMEX(:,2)));
        ymin = min(ymin,min(dataRHSFunctionIMEX(:,3)));
        ymax = max(ymax,max(dataRHSFunctionIMEX(:,3)));
    else
        flagRHSFunctionIMEX = 0;
    end

    if (flagRHSFunctionExpl || flagRHSFunctionIMEX || flagIFunctionIMEX)
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
    if (    (~flagRHSFunctionExpl) ...
         && (~flagRHSFunctionIMEX) ...
         && (~flagIFunctionIMEX  ) ...
       )
        break;
    end
end
