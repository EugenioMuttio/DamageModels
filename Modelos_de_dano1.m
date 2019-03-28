function [rtrial] = Modelos_de_dano1 (MDtype,ce,eps_n1,n)
%**************************************************************************************
%*          Defining damage criterion surface                                        %*
%*                                                                                   %*
%*
%*                          MDtype=  1      : SYMMETRIC                              %*
%*                          MDtype=  2      : ONLY TENSION                           %*
%*                          MDtype=  3      : NON-SYMMETRIC                          %*
%*                                                                                   %*
%*                                                                                   %*
%* OUTPUT:                                                                           %*
%*                          rtrial                                                   %*               
%**************************************************************************************



%**************************************************************************************
if (MDtype==1)      %* Symmetric
rtrial= sqrt(eps_n1*ce*eps_n1')                       ;

elseif (MDtype==2)  %* Only tension 

    sigma_n1=ce*eps_n1';
    for j=1:4
        if sigma_n1(j)<0
            sigma_n1(j)=0;
        end
    end
    
    cei=inv(ce);
    eps_n1m=(cei*sigma_n1)';
    
    rtrial= sqrt(eps_n1m*ce*eps_n1m')                       ;
    
elseif (MDtype==3)  %*Non-symmetric
   
    sigma_n1=eps_n1*ce;
    sum_stress_abs=sum(abs(sigma_n1));
    sum_stress=0;
    for j=1:4
        if sigma_n1(j)<0
            sigma_n1(j)=0;
        end
        sum_stress=sum_stress+sigma_n1(j);
    end
    
    theta_stress=sum_stress/sum_stress_abs;
    coeft=theta_stress+(1-theta_stress)/n;

    rtrial= coeft*sqrt(eps_n1*ce*eps_n1');

end
%**************************************************************************************
return