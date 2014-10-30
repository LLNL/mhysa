function [TimeSteps,Hyp,Par,Sou,RHSFunc,IFunc,IJac,IJFunc] ...
    = ReadFunctionCounts()

% READFUNCTIONCOUNTS Reads and returns the function counts written
%                    out by HyPar

if (exist('function_counts.dat','file'))
    data = load('function_counts.dat');
    TimeSteps   = data(1);
    Hyp         = data(2); 
    Par         = data(3);
    Sou         = data(4);
    RHSFunc     = data(5);
    IFunc       = data(6);
    IJac        = data(7);
    IJFunc      = data(8);
else
    TimeSteps   = 0;
    Hyp         = 0;
    Par         = 0;
    Sou         = 0;
    RHSFunc     = 0;
    IFunc       = 0;
    IJac        = 0;
    IJFunc      = 0;
end

