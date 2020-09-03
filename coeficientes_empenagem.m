function [Cl,Cd,Cm] = coeficientes_empenagem(alfa,beta,empenagem)

% alfa = -4
% empenagem = 1 
% beta = 10

switch empenagem
    
    case 1 %Perfil HS2
        
        
        if beta == -19
        
            
            Cl = 0.0002*alfa^2 + 0.0641*alfa - 0.9618 ;
            Cm = -0.0003*alfa^2 + 0.0033*alfa + 0.2979 ;
            Cd = 1*10^(-6)*alfa^4 - 3e-05*alfa^3 + 0.0004*alfa^2 - 0.009*alfa + 0.1167 ;
            
        
        elseif beta == -15
            
            
            Cl = 0.065*alfa - 0.8003 ; %flap -15
            Cm =  -0.0005*alfa^2 + 0.0053*alfa + 0.2598 ; %flap -15
            Cd = 0.0004*alfa^2 - 0.0087*alfa + 0.0844 ; %flap -15
            
        elseif beta>-15 && beta<-10
            
            
            Cl_maior = 0.0638*alfa - 0.5877 ; %flap -10
            Cm_maior = -0.0004*alfa^2 + 0.0042*alfa + 0.202; %flap -10
            Cd_maior = 0.003*alfa^2 - 0.0064*alfa + 0.0524; %flap -10
            
            Cl_menor = 0.065*alfa - 0.8003 ; %flap -15
            Cm_menor =  -0.0005*alfa^2 + 0.0053*alfa + 0.2598 ; %flap -15
            Cd_menor = 0.0004*alfa^2 - 0.0087*alfa + 0.0844 ; %flap -15
            
            dif = abs((15 - (-beta))/5) ;
            
            Cl = Cl_menor*dif + Cl_maior*(1-dif) ;
            Cm = Cm_menor*dif + Cm_maior*(1-dif) ;
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif beta == -10
            
            
            Cl = 0.0638*alfa - 0.5877 ; %flap -10
            Cm = -0.0004*alfa^2 + 0.0042*alfa + 0.202; %flap -10
            Cd = 0.003*alfa^2 - 0.0064*alfa + 0.0524; %flap -10
            
        elseif beta>-10 && beta<-5
            
            Cl_menor = 0.0638*alfa - 0.5877 ; %flap -10
            Cm_menor = -0.0004*alfa^2 + 0.0042*alfa + 0.202; %flap -10
            Cd_menor = 0.003*alfa^2 - 0.0064*alfa + 0.0524; %flap -10
            
            Cl_maior =  0.063*alfa - 0.3854 ; %flap -5
            Cm_maior = -0.0004*alfa^2 + 0.003*alfa + 0.1459 ; %flap -5
            Cd_maior = 0.0003*alfa^2 - 0.0042*alfa + 0.0298 ; %flap -5
            
            dif = abs((10 - (-beta))/5) ;
            
            Cl = Cl_menor*dif + Cl_maior*(1-dif) ;
            Cm = Cm_menor*dif + Cm_maior*(1-dif) ;
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif beta == -5
            
            Cl =  0.063*alfa - 0.3854 ; %flap -5
            Cm = -0.0004*alfa^2 + 0.003*alfa + 0.1459 ; %flap -5
            Cd = 0.0003*alfa^2 - 0.0042*alfa + 0.0298 ; %flap -5
            
        elseif beta>-5 && beta<0
            
            Cl_maior = 0.0629*alfa - 0.1867 ; %Flap 0
            Cm_maior =  -0.0004*alfa^2 + 0.0019*alfa + 0.0911 ; %Flap 0
            Cd_maior =  0.0004*alfa^2 - 0.002*alfa + 0.0139 ; %Flap 0
            
            Cl_menor =  0.063*alfa - 0.3854 ; %flap -5
            Cm_menor = -0.0004*alfa^2 + 0.003*alfa + 0.1459 ; %flap -5
            Cd_menor = 0.0003*alfa^2 - 0.0042*alfa + 0.0298 ; %flap -5
            
            dif = abs((5 - (-beta))/5) ;
            
            Cl = Cl_menor*dif + Cl_maior*(1-dif) ;
            Cm = Cm_menor*dif + Cm_maior*(1-dif) ;
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif beta == 0
            
            
            Cl = 0.0629*alfa - 0.1867 ; %Flap 0
            Cm =  -0.0004*alfa^2 + 0.0019*alfa + 0.0911 ; %Flap 0
            Cd =  0.0004*alfa^2 - 0.002*alfa + 0.0139 ; %Flap 0
            
            
        elseif beta>0 && beta<5
            
            
            Cl_maior = 0.0602*alfa + 0.0087 ; %flap 5
            Cm_maior =  -0.0003*alfa^2 + 0.001*alfa + 0.0381 ; %flap 5
            Cd_maior = 0.0004*alfa^2 + 0.0001*alfa + 0.008 ; %flap 5
            
            Cl_menor = 0.0629*alfa - 0.1867 ; %Flap 0
            Cm_menor =  -0.0004*alfa^2 + 0.0019*alfa + 0.0911 ; %Flap 0
            Cd_menor =  0.0004*alfa^2 - 0.002*alfa + 0.0139 ; %Flap 0
            
            dif = abs((-beta)/5) ;
            
            Cl = Cl_menor*dif + Cl_maior*(1-dif) ;
            Cm = Cm_menor*dif + Cm_maior*(1-dif) ;
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif beta == 5
            
            
            Cl = 0.0602*alfa + 0.0087 ; %flap 5
            Cm =  -0.0003*alfa^2 + 0.001*alfa + 0.0381 ; %flap 5
            Cd = 0.0004*alfa^2 + 0.0001*alfa + 0.008 ; %flap 5
            
        elseif beta>5 && beta<10
            
            Cl_menor = 0.0602*alfa + 0.0087 ; %flap 5
            Cm_menor =  -0.0003*alfa^2 + 0.001*alfa + 0.0381 ; %flap 5
            Cd_menor = 0.0004*alfa^2 + 0.0001*alfa + 0.008 ; %flap 5
            
            Cl_maior =  0.0632*alfa + 0.2013 ; %Flap 10
            Cm_maior =  -0.0004*alfa^2 - 0.0004*alfa - 0.0134 ; %Flap 10
            Cd_maior = 0.0004*alfa^2 + 0.0022*alfa + 0.0131 ; %Flap 10
            
            dif = abs((5 -beta)/5) ;
            
            Cl = Cl_menor*dif + Cl_maior*(1-dif) ;
            Cm = Cm_menor*dif + Cm_maior*(1-dif) ;
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif beta == 10
            
            Cl =  0.0632*alfa + 0.2013 ; %Flap 10
            Cm =  -0.0004*alfa^2 - 0.0004*alfa - 0.0134 ; %Flap 10
            Cd = 0.0004*alfa^2 + 0.0022*alfa + 0.0131 ; %Flap 10
            
        elseif beta>10 && beta<15
            
            Cl_menor =  0.0632*alfa + 0.2013 ; %Flap 10
            Cm_menor =  -0.0004*alfa^2 - 0.0004*alfa - 0.0134 ; %Flap 10
            Cd_menor = 0.0004*alfa^2 + 0.0022*alfa + 0.0131 ; %Flap 10
            
            Cl_maior = 0.0645*alfa + 0.4075 ; %flap 15
            Cm_maior = -0.0004*alfa^2 - 0.0016*alfa - 0.0645 ;%flap 15
            Cd_maior =  0.0004*alfa^2 + 0.0043*alfa + 0.0274 ; %flap 15
            
            
            dif = abs((10 -beta)/5) ;
            
            Cl = Cl_menor*dif + Cl_maior*(1-dif) ;
            Cm = Cm_menor*dif + Cm_maior*(1-dif) ;
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif beta == 15
            
            Cl = 0.0645*alfa + 0.4075 ; %flap 15
            Cm = -0.0004*alfa^2 - 0.0016*alfa - 0.0645 ;%flap 15
            Cd =  0.0004*alfa^2 + 0.0043*alfa + 0.0274 ; %flap 15
            
            
            
        end
        
    case 2 %Perfil novo
        
        if beta == -15
            
            Cl = 0.0002*alfa^2 + 0.0632*alfa - 0.8336 ; %flap -15
            Cm = -0.0005*alfa^2 + 0.005*alfa + 0.3187 ; %flap -15
            Cd = 2*10^(-05)*alfa^3 - 0.0002*alfa^2 - 0.007*alfa + 0.1547 ; %flap -15
        elseif beta>-15 && beta<-10
            
            
            Cl_maior = 0.0002*alfa^2 + 0.0632*alfa - 0.8336 ; %flap -10
            Cm_maior = -0.0005*alfa^2 + 0.005*alfa + 0.3187; %flap -10
            Cd_maior = 2*10^(-6)*alfa^4 - 6*10^(-05)*alfa^3 + 0.0006*alfa^2 - 0.0075*alfa + 0.1068; %flap -10
            
            Cl_menor = 0.0002*alfa^2 + 0.0632*alfa - 0.8336 ; %flap -15
            Cm_menor =  -0.0005*alfa^2 + 0.005*alfa + 0.3187 ; %flap -15
            Cd_menor = 2*10^(-05)*alfa^3 - 0.0002*alfa^2 - 0.007*alfa + 0.1547 ; %flap -15
            
            dif = abs((15 - (-beta))/5) ;
            
            Cl = Cl_menor*dif + Cl_maior*(1-dif) ;
            Cm = Cm_menor*dif + Cm_maior*(1-dif) ;
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif beta == -10
            
            Cl = 0.0002*alfa^2 + 0.0632*alfa - 0.8336 ; %flap -10
            Cm = -0.0005*alfa^2 + 0.005*alfa + 0.3187; %flap -10
            Cd = 2*10^(-6)*alfa^4 - 6*10^(-05)*alfa^3 + 0.0006*alfa^2 - 0.0075*alfa + 0.1068; %flap -10
            
        elseif beta>-10 && beta<-5
            
            Cl_menor = 0.0002*alfa^2 + 0.0632*alfa - 0.8336 ; %flap -10
            Cm_menor = -0.0005*alfa^2 + 0.005*alfa + 0.3187; %flap -10
            Cd_menor = 2*10^(-6)*alfa^4 - 6*10^(-5)*alfa^3 + 0.0006*alfa^2 - 0.0075*alfa + 0.1068; %flap -10
            
            Cl_maior =  0.063*alfa - 0.3854 ; %flap -5
            Cm_maior = -0.0004*alfa^2 + 0.004*alfa + 0.2522 ; %flap -5
            Cd_maior = 3*10^(-7)*alfa^3 + 0.0004*alfa^2 - 0.0066*alfa + 0.0744 ; %flap -5
            
            dif = abs((10 - (-beta))/5) ;
            
            Cl = Cl_menor*dif + Cl_maior*(1-dif) ;
            Cm = Cm_menor*dif + Cm_maior*(1-dif) ;
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif beta == -5
            
            Cl =  0.063*alfa - 0.3854 ; %flap -5
            Cm = -0.0004*alfa^2 + 0.004*alfa + 0.2522 ; %flap -5
            Cd = 3*10^(-7)*alfa^3 + 0.0004*alfa^2 - 0.0066*alfa + 0.0744 ; %flap -5
            
        elseif beta>-5 && beta<0
            
            Cl_maior = 0.0625*alfa - 0.4038 ; %Flap 0
            Cm_maior =  -0.0004*alfa^2 + 0.0032*alfa + 0.2011 ; %Flap 0
            Cd_maior = -7*10^(-6)*alfa^3 + 0.0004*alfa^2 - 0.0039*alfa + 0.047
            
            Cl_menor =  0.063*alfa - 0.3854 ; %flap -5
            Cm_menor = -0.0004*alfa^2 + 0.004*alfa + 0.2522 ; %flap -5
            Cd_menor = 3*10^(-7)*alfa^3 + 0.0004*alfa^2 - 0.0066*alfa + 0.0744 ; %flap -5
            
            dif = abs((5 - (-beta))/5) ;
            
            Cl = Cl_menor*dif + Cl_maior*(1-dif) ;
            Cm = Cm_menor*dif + Cm_maior*(1-dif) ;
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif beta == 0
            
            
            Cl = 0.0629*alfa - 0.1867 ; %Flap 0
            Cm =   -0.0004*alfa^2 + 0.0032*alfa + 0.2011 ; %Flap 0
            Cd = -7*10^(-6)*alfa^3 + 0.0004*alfa^2 - 0.0039*alfa + 0.047; %Flap 0
        elseif beta>0 && beta<5
            
            
            Cl_maior = 0.0621*alfa - 0.199 ; %flap 5
            Cm_maior =  -0.0004*alfa^2 + 0.0019*alfa + 0.1384 ; %flap 5
            Cd_maior = 10^(-7)*alfa^6 - 2*10^(-7)*alfa^5 - 2*10^(-5)*alfa^4 + 3*10^(-5)*alfa^3 + 0.0011*alfa^2 - 0.0023*alfa + 0.0268 ; %flap 5
            
            Cl_menor = 0.0629*alfa - 0.1867 ; %Flap 0
            Cm_menor =  -0.0004*alfa^2 + 0.0032*alfa + 0.2011 ; %Flap 0
            Cd_menor = -7*10^(-6)*alfa^3 + 0.0004*alfa^2 - 0.0039*alfa + 0.047 %Flap 0
            
            dif = abs((-beta)/5) ;
            
            Cl = Cl_menor*dif + Cl_maior*(1-dif) ;
            Cm = Cm_menor*dif + Cm_maior*(1-dif) ;
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif beta == 5
            
            
            Cl =  0.0621*alfa - 0.199 ; %flap 5
            Cm =   -0.0004*alfa^2 + 0.0019*alfa + 0.1384 ; %flap 5
            Cd = 10^(-7)*alfa^6 - 2*10^(-7)*alfa^5 - 2*10^(-5)*alfa^4 + 3*10^(-5)*alfa^3 + 0.0011*alfa^2 - 0.0023*alfa + 0.0268 ; %flap 5
            
        elseif beta>5 && beta<10
            
            Cl_menor =  0.0621*alfa - 0.199 ; %flap 5
            Cm_menor = -0.0004*alfa^2 + 0.0019*alfa + 0.1384 ; %flap 5
            Cd_menor = 10^(-7)*alfa^6 - 2*10^(-7)*alfa^5 - 2*10^(-5)*alfa^4 + 3*10^(-5)*alfa^3 + 0.0011*alfa^2 - 0.0023*alfa + 0.0268 ; %flap 5
            
            Cl_maior =   0.063*alfa - 0.0027 ; %Flap 10
            Cm_maior =   -0.0004*alfa^2 + 0.0007*alfa + 0.0851; %Flap 10
            Cd_maior = 2*10^(-7)*alfa^6 + 9*10^(-7)*alfa^5 - 2*10^(-5)*alfa^4 - 9*10^(-5)*alfa^3 + 0.0012*alfa^2 + 0.0021*alfa + 0.0187 ; %Flap 10
            
            dif = abs((5 -beta)/5) ;
            
            Cl = Cl_menor*dif + Cl_maior*(1-dif) ;
            Cm = Cm_menor*dif + Cm_maior*(1-dif) ;
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif beta == 10
            
            Cl = 0.063*alfa - 0.0027 ; %Flap 10
            Cm =  -0.0004*alfa^2 + 0.0007*alfa + 0.0851; %Flap 10
            Cd = 2*10^(-7)*alfa^6 + 9*10^(-7)*alfa^5 - 2*10^(-5)*alfa^4 - 9*10^(-5)*alfa^3 + 0.0012*alfa^2 + 0.0021*alfa + 0.0187 ; %Flap 10
            
        elseif beta>10 && beta<15
            
            Cl_menor =  0.063*alfa - 0.0027 ; %Flap 10
            Cm_menor =   -0.0004*alfa^2 + 0.0007*alfa + 0.0851; %Flap 10
            Cd_menor = 2*10^(-7)*alfa^6 + 9*10^(-7)*alfa^5 - 2*10^(-5)*alfa^4 - 9*10^(-5)*alfa^3 + 0.0012*alfa^2 + 0.0021*alfa + 0.0187 ; %Flap 10
            
            Cl_maior = 0.0685*alfa - 1.0779 ; %flap 15
            Cm_maior = -0.0005*alfa^2 + 0.006*alfa + 0.3841 ;%flap 15
            Cd_maior = 2*10^(-7)*alfa^6 + 2*10^(-6)*alfa^5 - 7*10^(-6)*alfa^4 - 0.0002*alfa^3 + 0.0003*alfa^2 + 0.005*alfa + 0.0253 ; %flap 15
            
            
            
            dif = abs((10 -beta)/5) ;
            
            Cl = Cl_menor*dif + Cl_maior*(1-dif) ;
            Cm = Cm_menor*dif + Cm_maior*(1-dif) ;
            Cd = Cd_menor*dif + Cd_maior*(1-dif) ;
            
        elseif beta == 15
            
            Cl = 0.0685*alfa - 1.0779 ; %flap 15
            Cm = -0.0005*alfa^2 + 0.006*alfa + 0.3841 ;%flap 15
            Cd = 2*10^(-7)*alfa^6 + 2*10^(-6)*alfa^5 - 7*10^(-6)*alfa^4 - 0.0002*alfa^3 + 0.0003*alfa^2 + 0.005*alfa + 0.0253 ; %flap 15
            
        end
        
        
        
        
end

end

% Cl
% Cd
% Cm