function [cd] = CD (superficie, perfilsuperficie, v, cl, ~)

%Essa fun��o calcula o CD0 baseado nos c�lculos do cap�tulo 3 do Sadraey
%Como entrada, ele toma o necess�rio da superficie em forma de uma
%estrutura de dados, e tamb�m do perfil da superf�cie sustentadora a ser
%analisada.
%O �ltimo input necess�rio � a velocidade a ser rodada a an�lise, a fim de
%calcular o Re para as superf�cies
%Um input opcional � um CL arbitr�rio, para calcular o arrasto completo da
%superf�cie
%Outro input opcional � h. Se houver esse input, o arrasto induzido ser�
%calculado computando h como sendo a altura da asa relativa ao solo, a fim
%de incluir o efeito solo

M = 0; %Para as velocidades em nosso envelope, M � aprox 0.

% -----------------------------Arrasto parasita da asa
rho = 1.15;
muu = 1.8608e-5;

tcmax = perfilsuperficie.tcmax;
cdminw = perfilsuperficie.cdminw;
s = superficie.s;
mac = superficie.mac;

Swet =  2*(1+(0.5*tcmax))*superficie.s;                              % Area molhada da superf�cie
ftcw = 1+(2.7*tcmax)+(100*(tcmax^4));                                % Fator em funcao da espessura m?xima do perfil na asa (EQ 3.25 do Sadraey)
Re = (rho*v*mac)/muu;                                                % Numero de Reynolds
cfw = 0.455/log10(Re)^2.58;                                          % Fator em funcao do numero de Reynolds
fM = 1-(0.08*(M^1.45));                                              % Fator em funcao do numero de Mach
CD0 = (cfw*ftcw*fM)*(Swet/s)*((cdminw/0.004)^0.4);                   % Arrasto parasita da superf�cie

switch nargin
    
    case 3
        
        cd = CD0;                                                    % Caso seja necess�rio calcular apenas o arrasto parasita
        
    case 4
        e = superficie.e;                                            % Fator de efici�ncia de Oswald
        ar = superficie.ar;                                          % Alongamento da asa
        k = 1/(pi*e*ar);                                             % Constante da polar
            
        cd = CD0+k*cl^2;                                             %Coeficiente de arrasto completo da superf�cie, para o CL especificado
    case 5  
        
        e = superficie.e;                                            % Fator de efici�ncia de Oswald
        ar = superficie.ar;                                          % Alongamento da asa
        k = 1/(pi*e*ar);                                             % Constante da polar
        
        G = ((16*superficie.h/superficie.b))^2/(1+ ((16*superficie.h/superficie.b))^2);    %Efeito solo, baseado em Anderson, EQ 6.77, pg 357.
        k = k * G;
            cd = CD0+k*cl^2;     
end
end