function [x,fval,info,output]= optimize()
load('guess_polished.mat','guess_polished');
option.derivative_option = 0;
guess=guess_polished(1,:);
xlow=[-10,-10,-10,0,0,0,zeros(1,sum(N)),0.5,0.5,0.5,0.5,0.5,0.5,ones(1,26)*(-3),0,0,0];
xupp=[1,1,1,1,1,1,ones(1,sum(N)),5.0,5.0,5.0,5.0,5.0,5.0,ones(1,26)*(3),5,5,5];
xmul=zeros(1,51);
xstate=zeros(1,51);
Flow=-inf;
Fupp=inf;
Fmul=0;
Fstate=0;
[x,fval,info,output]=snopt(guess,xlow,xupp,xmul,xstate, ...
		     Flow,Fupp,Fmul,Fstate,@compute_error,options);
