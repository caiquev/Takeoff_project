function [Cl,Cd,Cm] = coeficientes_asa(alfa,X_speed)

%     Cl =  -0.0005*alfa^2 + 0.0931*alfa + 0.9406 ; %asa Bruxo perfil goe na ponta
%     Cm =  -0.0047*alfa - 0.3495    ; %asa Bruxo perfil goe na ponta
%     Cd = -2*10^(-6)*alfa^4 - 2*10^(-5)*alfa^3 + 0.0012*alfa^2 + 0.001*alfa + 0.0527  ; %asa Bruxo perfil goe na ponta
%         
%         Cl = 0.0911*alfa + 1.0995 ; % asa VANT 12.5
%         Cm = 5*10^(-5)*alfa^2 - 0.0048*alfa - 0.4308 ;  % asa VANT 12.5
%         Cd = 0.0006*alfa^2 + 0.0036*alfa + 0.0727 ; % asa VANT 12.5

% Cl = 0.0886*alfa + 1.0357 ; %VANT GOE NA PONTA Estol 12.5
% Cm = -0.0044*alfa - 0.4037 ; %VANT GOE NA PONTA 12.5 
% Cd = -2*10^(-5)*alfa^3 + 0.0007*alfa^2 + 0.0034*alfa + 0.0602 ; %VANT GOE NA PONTA 12.5

%  Cl = 0.0947*alfa + 0.9589 ; %BRUXO
%  Cm = -0.0049*alfa - 0.3676 ;  %BRUXO
%  Cd = 0.0011*alfa^2 + 0.0016*alfa + 0.0705 ;  %BRUXO

    Cm =  -0.0046*alfa - 0.4166 ; %Madrugada 15
%      Cd = -3E-5*alfa^3 + 0.0009*alfa^2 + 0.0018*alfa + 0.0712  ; %Madrugada 15
   Cl =  -0.0004*alfa^2 + 0.0879*alfa + 1.0808 ; %Madrugada 15
   
   if X_speed <= 5
   
    Cd = -0.0104*X_speed^2 + 0.0043*X_speed + 0.3448 ; % CD para 3 graus de incidencia da asa Variando a velocidade.

   else 
       
     Cd =  0.0003*X_speed^2 - 0.0073*X_speed + 0.1331 ; % CD para 3 graus de incidencia da asa Variando a velocidade.
     
   end
end