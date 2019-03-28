function [sigma_n1,hvar_n1,aux_var] = rmap_dano1 (eps_n,eps_n1,hvar_n,Eprop,ce,MDtype,n,delta_t)

%**************************************************************************************
%*                                         *
%*           Integration Algorithm for a isotropic damage model
%*
%*                                                                                    *
%*            [sigma_n1,hvar_n1,aux_var] = rmap_dano1 (eps_n1,hvar_n,Eprop,ce)        *
%*                                                                                    *
%* INPUTS              eps_n1(4)   strain (almansi)    step n+1                       *
%*                                 vector R4    (exx eyy exy ezz)                     *
%*                     hvar_n(6)   internal variables , step n                        *
%*                                 hvar_n(1:4) (empty)                          *
%*                                 hvar_n(5) = r  ; hvar_n(6)=q                       *
%*                     Eprop(:)    Material parameters                                *
%*
%*                     ce(4,4)     Constitutive elastic tensor                        *
%*                                                                                    *
%* OUTPUTS:            sigma_n1(4) Cauchy stress  , step n+1                          *
%*                     hvar_n(6)   Internal variables , step n+1                           *
%*                     aux_var(3)  Auxiliar variables for computing const. tangent tensor  *
%***************************************************************************************


hvar_n1 = hvar_n;
r_n     = hvar_n(5);
q_n     = hvar_n(6);
E       = Eprop(1);
nu      = Eprop(2);
H       = Eprop(3);
sigma_u = Eprop(4);
hard_type = Eprop(5) ;

eta=Eprop(7);
alpha=Eprop(8);
%*************************************************************************************


%*************************************************************************************
%*       initializing                                                %*
 r0 = sigma_u/sqrt(E);
 zero_q=1.d-6*r0;
 A=abs(H);
 q_inf=r0+sign(H)*0.99*r0;
% if(r_n<=0.d0)
%     r_n=r0;
%     q_n=r0;
% end
%*************************************************************************************


%*************************************************************************************
%*       Damage surface                                                              %*
%[rtrial] = Modelos_de_dano1 (MDtype,ce,eps_n1,n);

[tau_eps_n]=Modelos_de_dano1 (MDtype,ce,eps_n,n);%sqrt(eps_n*ce*eps_n');
[tau_eps_n1]=Modelos_de_dano1 (MDtype,ce,eps_n1,n);%sqrt(eps_n1*ce*eps_n1');
[rtrial]=(1-alpha)*tau_eps_n+alpha*tau_eps_n1;
%*************************************************************************************


%*************************************************************************************
%*   Ver el Estado de Carga                                                           %*
%*   --------->    fload=0 : elastic unload                                           %*
%*   --------->    fload=1 : damage (compute algorithmic constitutive tensor)         %*
fload=0;


if(rtrial > r_n)
    %*   Loading

    fload=1;
    delta_r=rtrial-r_n;
    %r_n1= rtrial  ;
    r_n1= ((eta-delta_t*(1-alpha))/(eta+alpha*delta_t))*r_n+(delta_t/(eta+alpha*delta_t))*rtrial;
    if hard_type == 0
        %  Linear
        q_n1= q_n+ H*delta_r;
    else
        % Comment/delete lines below once you have implemented this case
        % *******************************************************
        %menu({'Hardening/Softening exponential law has not been implemented yet. '; ...
        %'Modify file "rmap_dano1" ' ; ...
        %'to include this option'},  ...
        %'STOP');
        %error('OPTION NOT AVAILABLE')
        
        %Hexp=A*(zero_q-r0)/r0*exp(A*(1-r_n/r0));
        
           
        q_n1= q_inf-(q_inf-r0)*exp(A*(1-r_n1/r0));

        
    end

    if(q_n1<zero_q)
        q_n1=zero_q;
    end


else

    %*     Elastic load/unload
    fload=0;
    r_n1= r_n  ;
    q_n1= q_n  ;


end
% Damage variable
% ---------------
dano_n1   = 1.d0-(q_n1/r_n1);
%  Computing stress
%  ****************
sigma_n1  =(1.d0-dano_n1)*ce*eps_n1';

%hold on 
%plot(sigma_n1(1),sigma_n1(2),'bx')

%*************************************************************************************


%*************************************************************************************
%* Updating historic variables                                            %*
%  hvar_n1(1:4)  = eps_n1p;
hvar_n1(5)= r_n1 ;
hvar_n1(6)= q_n1 ;
%*************************************************************************************




%*************************************************************************************
%* Auxiliar variables                                                               %*
aux_var(1) = fload;
aux_var(2) = q_n1/r_n1;
aux_var(3) = (q_n1-H*r_n1)/r_n1^3;
Ce_tan=(1-dano_n1)*ce-aux_var(1)*aux_var(3)*((ce*eps_n1')*(ce*eps_n1')');
Ce_alg=Ce_tan-aux_var(1)*((alpha*delta_t)/(eta+alpha*delta_t)*aux_var(3)*r_n1/(sqrt(sigma_n1'*(inv(ce))*sigma_n1)/(1-dano_n1))*((ce*eps_n1')*(ce*eps_n1')'));
aux_var(4)=Ce_tan(1,1); %Ce_tan
aux_var(5)=Ce_alg(1,1); %Ce_alg
%*************************************************************************************
 










