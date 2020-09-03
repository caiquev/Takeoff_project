function [Cl,Cd,Cm] = coeficientes_empenagem(alfa,beta,empenagem,X_speed)

switch empenagem
    
    case 2 %Perfil themyto
        
        Cl = 0.0632*alfa + 0.2118 + 0.065*(beta+12) ; %Rô 1.15
        Cm =  -0.0002*alfa^2 - 5E-05*alfa - 0.002 + 0.0117*(beta+12) ; %Rô 1.15
        %
        %             Cl = -0.0315 ; %Rô 1.07
        %             Cm = -0.0069 ; %Rô 1.07
        
        if X_speed <= 1
            
            Cd =  0.00184484939789*alfa^2 + 0.01301310473961*alfa + 0.06607719794145 ; %rô 1.15
            %                 Cd = 0.0378409 ;
            
        elseif   (1 <X_speed)   && (X_speed <= 3 )
            
            Cd_menor =  0.00184484939789*alfa^2 + 0.01301310473961*alfa + 0.06607719794145 ; %rô 1.15
            %                 Cd_menor = 0.0378409 ;
            
            
            Cd_maior =  0.00000010755257*alfa^6 + 0.00000198379714*alfa^5 - 0.00000103685269*alfa^4 - 0.00013379869815*alfa^3 + 0.00021481140413*alfa^2 + 0.00409231065635*alfa + 0.02953319060325; %rô 1.15
            %                 Cd_maior = 0.02577871 ;
            
            
            
            dif = abs((3-X_speed)/2) ;
            
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif (3 < X_speed ) && (X_speed <= 5 )
            
            Cd_menor = 0.00000010755257*alfa^6 + 0.00000198379714*alfa^5 - 0.00000103685269*alfa^4 - 0.00013379869815*alfa^3 + 0.00021481140413*alfa^2 + 0.00409231065635*alfa + 0.02953319060325; %rô 1.15
            %                 Cd_menor = 0.02577871 ;
            
            Cd_maior =  0.00038760352577*alfa^2 + 0.00253495557974*alfa + 0.02164873823553 ; %rô 1.15
            %                 Cd_maior = 0.01917881 ;
            
            dif = abs((5-X_speed)/2) ;
            
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif (5 < X_speed ) && (X_speed <= 7)
            
            Cd_menor = 0.00038760352577*alfa^2 + 0.00253495557974*alfa + 0.02164873823553 ; %rô 1.15
            %                 Cd_menor = 0.01917881 ;
            
            Cd_maior = 0.00037376230329*alfa^2 + 0.00235632275236*alfa + 0.01884708416207; %rô 1.15
            %                 Cd_maior = 0.01917881;
            
            
            dif = abs((7-X_speed)/2) ;
            
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif (9 > X_speed ) && (X_speed > 7)
            
            Cd_menor =  0.00037376230329*alfa^2 + 0.00235632275236*alfa + 0.01884708416207; %rô 1.15
            %                 Cd  = 0.01917881 ;
            Cd_maior =  0.00000036500558*alfa^4 + 0.00000515524382*alfa^3 + 0.00033476280330*alfa^2 + 0.00194044994842*alfa + 0.01793728505472 ; %rô 1.15
            
            dif = abs((5-X_speed)/2) ;
            
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif (9 <= X_speed ) && (X_speed <= 11)
            
            Cd_menor = 0.00000036500558*alfa^4 + 0.00000515524382*alfa^3 + 0.00033476280330*alfa^2 + 0.00194044994842*alfa + 0.01793728505472 ;  %rô 1.15
            Cd_maior =  0.00000146881395*alfa^3 + 0.00038724079691*alfa^2 + 0.00224939029752*alfa + 0.01598144896917 ; %rô 1.15
            dif = abs((5-X_speed)/2) ;
            
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif X_speed > 11
            
            Cd =  0.00036770536124*alfa^2 + 0.00231010773896*alfa + 0.01627466160848 ;
            
        end
        
        
end
end

% Cl
% Cd
% Cm