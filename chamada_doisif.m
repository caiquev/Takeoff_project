%FUNÇÃO MODEFRONTIER
X_position = 0
% while X_position<55
%    
    angulo_incidencia_asa = 3;
    angulo_incidencia_emp = -4 ;
    beta1 = 0 ;
    beta2 = -19 ;
    empenagem = 1 ;
    mass_acft = 18.15 ;
    tdp_x = 0.14;
    
%         mass_acft = mass_acft +0.05 ;
    [M_tdp,Vstall,X_speed,MTOW,liftoff,alfa,X_position,alfaestol,R_c] = doisif(mass_acft,angulo_incidencia_asa,angulo_incidencia_emp,tdp_x,empenagem,beta1,beta2)
    
    
% end


