function [S] = construct_S(n,sizes)

% n is the number of discretized grids
% size is the total number of inputs and outputs
% this function constructs the error matrix from which we draw the random
% error psi

[x11_sup,x11_prob]=ini_normal(11.32,41.90,n);
[x12_sup,x12_prob]=ini_normal(35.24,66.03,n);
[x13_sup,x13_prob]=ini_normal(788.69,703.06,n);
[x14_sup,x14_prob]=ini_normal(141.39,139.76,n);
[x15_sup,x15_prob]=ini_normal(10.32,3.76,n);
[x21_sup,x21_prob]=ini_normal(147.15,142.52,n);
[x22_sup,x22_prob]=ini_normal(34.82,68.96,n);
[x23_sup,x23_prob]=ini_normal(3.82,1.96,n);
[x31_sup,x31_prob]=ini_normal(220.65,157.49,n);
[x32_sup,x32_prob]=ini_normal(12.65,5.49,n);
[y3_sup,y3_prob]=ini_normal(1777.23,803.24,n);

lnx11_sup=log(x11_sup);
lnx12_sup=log(x12_sup);
lnx13_sup=log(x13_sup);
lnx14_sup=log(x14_sup);
lnx15_sup=log(x15_sup);
lnx21_sup=log(x21_sup);
lnx22_sup=log(x22_sup);
lnx23_sup=log(x23_sup);
lnx31_sup=log(x31_sup);
lnx32_sup=log(x32_sup);
lny3_sup=log(y3_sup);

std_lnx11=calculate_std(lnx11_sup,x11_prob,n);
std_lnx12=calculate_std(lnx12_sup,x12_prob,n);
std_lnx13=calculate_std(lnx13_sup,x13_prob,n);
std_lnx14=calculate_std(lnx14_sup,x14_prob,n);
std_lnx15=calculate_std(lnx15_sup,x15_prob,n);
std_lnx21=calculate_std(lnx21_sup,x21_prob,n);
std_lnx22=calculate_std(lnx22_sup,x22_prob,n);
std_lnx23=calculate_std(lnx23_sup,x23_prob,n);
std_lnx31=calculate_std(lnx31_sup,x31_prob,n);
std_lnx32=calculate_std(lnx32_sup,x32_prob,n);
std_lny3=calculate_std(lny3_sup,y3_prob,n);

a=[std_lnx11,std_lnx12,std_lnx13,std_lnx14,std_lnx15,std_lnx21,std_lnx22,std_lnx23,std_lnx31,std_lnx32,std_lny3,std_lny3,std_lny3];
rho=0.15;
for i=1:sizes
    for j=1:sizes
        if (i==j)
            S(i,j)=a(i)*a(j);
        else
            S(i,j)=a(i)*a(j)*rho;
        end
    end
end


end

