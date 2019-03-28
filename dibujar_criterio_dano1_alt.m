function hplot = dibujar_criterio_dano1(ce,nu,q,tipo_linea,MDtype,n)
%*************************************************************************************
%*                 PLOT DAMAGE SURFACE CRITERIUM: ISOTROPIC MODEL                             %*
%*                                                                                  %*
%*      function [ce] = tensor_elastico (Eprop, ntype)                    %*
%*                                                                                  %*
%*      INPUTS                                                       %*
%*                                                                                  %*
%*                    Eprop(4)    vector de propiedades de material                 %*
%*                                      Eprop(1)=  E------>modulo de Young          %*
%*                                      Eprop(2)=  nu----->modulo de Poisson        %*
%*                                      Eprop(3)=  H----->modulo de Softening/hard. %*
%*                                      Eprop(4)=sigma_u----->tensi�n �ltima        %*
%*                     ntype                                 %*
%*                                 ntype=1  plane stress                            %*
%*                                 ntype=2  plane strain                            %*
%*                                 ntype=3  3D                                      %*
%*                     ce(4,4)     Constitutive elastic tensor  (PLANE S.       )    %*
%*                     ce(6,6)                                  ( 3D)                %*
%*************************************************************************************


%*************************************************************************************
%*        Inverse ce                                                                %*
ce_inv=inv(ce);
c11=ce_inv(1,1);
c22=ce_inv(2,2);
c12=ce_inv(1,2);
c21=c12;
c14=ce_inv(1,4);
c24=ce_inv(2,4);
%**************************************************************************************







%**************************************************************************************
% POLAR COORDINATES
if MDtype==1
    tetha=[0:0.01:2*pi];
    %**************************************************************************************
    %* RADIUS
    D=size(tetha);                       %*  Range
    m1=cos(tetha);                       %*
    m2=sin(tetha);                       %*
    Contador=D(1,2);                     %*
    
    
    radio = zeros(1,Contador) ;
    s1    = zeros(1,Contador) ;
    s2    = zeros(1,Contador) ;
               
    
    for i=1:Contador
        
        radio(i)= q/sqrt([m1(i) m2(i) 0 nu*(m1(i)+m2(i))]*ce_inv*[m1(i) m2(i) 0 ...
            nu*(m1(i)+m2(i))]');
        
        s1(i)=radio(i)*m1(i);
        s2(i)=radio(i)*m2(i);  
        
    end
    hplot =plot(s1,s2,tipo_linea);

    
elseif MDtype==2
    % Comment/delete lines below once you have implemented this case
    % *******************************************************
    %menu({'Damage surface "ONLY-TENSION" has not been implemented yet. '; ...
    %    'Modify files "Modelos_de_dano1" and "dibujar_criterio_dano1"' ; ...
    %    'to include this option'},  ...
    %    'STOP');
    %error('OPTION NOT AVAILABLE')
    
    tetha=[0:0.01:2*pi];
    %**************************************************************************************
    %* RADIUS
    D=size(tetha);                       %*  Range
    m1=cos(tetha);                       %*
    m2=sin(tetha);                       %*
    Contador=D(1,2);                     %*
    
    
    radio = zeros(1,Contador) ;
    s1    = zeros(1,Contador) ;
    s2    = zeros(1,Contador) ;
    m1_p    = zeros(1,Contador) ;
    m2_p    = zeros(1,Contador) ;

    
    for i=1:Contador
        m1_p(i)=m1(i);
        m2_p(i)=m2(i);
        if m1(i)<0
            m1_p(i)=0;
        end
        if m2(i)<0
            m2_p(i)=0;
        end
        radio(i)= q/sqrt([m1_p(i) m2_p(i) 0 nu*(m1_p(i)+m2_p(i))]*ce_inv*[m1(i) m2(i) 0 ...
            nu*(m1(i)+m2(i))]');
        
        s1(i)=radio(i)*m1(i);
        s2(i)=radio(i)*m2(i);  
        
    end
    
    
    hplot =plot(s1,s2,tipo_linea);
    
     
  
  
    
elseif MDtype==3
    % Comment/delete lines below once you have implemented this case
    % *******************************************************
%     menu({'Damage surface "NON-SYMMETRIC" has not been implemented yet. '; ...
%         'Modify files "Modelos_de_dano1" and "dibujar_criterio_dano1"' ; ...
%         'to include this option'},  ...
%         'STOP');
%     error('OPTION NOT AVAILABLE')
    
    tetha=[0:0.01:2*pi];
    %**************************************************************************************
    %* RADIUS
    D=size(tetha);                       %*  Range
    m1=cos(tetha);                       %*
    m2=sin(tetha);                       %*
    Contador=D(1,2);                     %*
    
    
    radio = zeros(1,Contador) ;
    s1    = zeros(1,Contador) ;
    s2    = zeros(1,Contador) ;
    
    for i=1:Contador
        sum_stress=0;
        vec_stress=[m1(i) m2(i) nu*(m1(i)+m2(i))];
        sum_stress_abs=sum(abs(vec_stress));
        for j=1:3
            if vec_stress(j)<0
                vec_stress(j)=0;
            end
            sum_stress=sum_stress+vec_stress(j);
        end
        
        theta_stress=sum_stress/sum_stress_abs;
        coeft=theta_stress+(1-theta_stress)/n;
        
        radio(i)= q/(coeft*(sqrt([m1(i) m2(i) 0 nu*(m1(i)+m2(i))]*ce_inv*[m1(i) m2(i) 0 ...
            nu*(m1(i)+m2(i))]')));
        
        s1(i)=radio(i)*m1(i);
        s2(i)=radio(i)*m2(i);  
        
    end
    hplot =plot(s1,s2,tipo_linea);
end
%**************************************************************************************



%**************************************************************************************
return



