function [M_tdp,Vstall,X_speed,MTOW,liftoff,alfa,X_position,alfaestol,R_c] = doisif(mass_acft,angulo_incidencia_asa,angulo_incidencia_emp,tdp_x,empenagem,beta1,beta2)
%Algoritmo de corrida de decolagem, baseada em uma integração das forças da
%corrida em um timestep fixo. O algoritmo precisa, para funcionar, das
%funções que determinam o coeficiente de sustentação e arrasto, suas
%variações com respeito a alfa, e também da tração dinâmica do motor. A
%saída do algoritmo é valor X_position que determina o comprimento de pista
%percorrido
%
close all
% clear
% clc

% Constantes

dt = 0.001;        % Timestep
ro = 1.07;        % Densidade do ar em jua
mi = 0.06;      % Coef. de atrito
n = 0.9 ;        %Eficiência de cauda
alfaestol = 15 ;
alfa_cl0 = -11 ; %MADRUGADA
%-------------------------------------------------------------------------------


asa = struct ('b', 2.45, 'perc_retan', 0.733, 'cr', 0.59, 'cp', 0.3, 'ar', 4.47, 'lambdac2', 6.8, 's', 1.336, 'mac', 0.559, 'e', 0.9);
emp = struct ('b', 1, 'perc_retan', 0, 'cr', 0.36 , 'cp', 0.235, 'ar', 3.3, 'lambdac2', -6.5, 's', 0.2952, 'mac', 0.301, 'e', 0.9);

W = mass_acft*9.81;

asa.h = 0.13;     % Altura da asa com relação ao solo

emp.h = 0.23;    % Altura da empenagem com relação ao solo

lht = 0.986; %Distância CA da asa até o CA de empenagem

x = 0.13;         %largura compartimento de carga
y = 0.2;          %comprimento compartimento de carga
z = 0.07;         %altura compartimento de carga

%--------------------------------------------------------------------------

% Distâncias relativas para uso dos cálculos dos braços de momento.

vetor_posicao_ca_tdp = [-tdp_x asa.h-0.055];
vetor_posicao_cg_tdp = [0.0535-tdp_x 0.15]; %0.078
vetor_posicao_ca_emp_tdp = [lht+(0.0535-tdp_x) emp.h-0.055];
vetor_posicao_motor_tdp = [0.507 0.18-0.055];


d = sqrt(asa.h^2 + vetor_posicao_cg_tdp(1)^2);      %distancia entre eixo central do compartimento e eixo do tdp
Iy = (1/12)*(mass_acft)*(x^2 + z^2);  %Momento de inércia em torno do eixo y, considerando 80% da massa total da aeronave para a carga paga
Iy_tdp = Iy + (mass_acft)*d^2;   %Eixos paralelos - Teorema Steiner

%--------------------------------------------------------------------------

% Valores iniciais para a iteração

liftoff = 0;      % Setando flag de decolagem como FALSO
time = 0;
X_speed = 0;
X_position = 0;
X_acceleration = 0;
alfa = 0;
pitch_rate = 0;
Pitch_acceleration = 0;
M_tdp = 0;
F_y = 0;
F_thrust = 0;
F_drag = 0;

%--------------------------------------------------------------------------
CL_max = coeficientes_asa(alfaestol,X_speed) ;

Vstall = sqrt(2*W/(ro*asa.s*CL_max)); %Formula da velocidade de estol

i = 0; %Contador
while liftoff == 0 % Início da corrida de decolagem.
    
    
    %------------------------- Cálculo do downwash --------------------------
    
    
    lambda = asa.cp/asa.cr ;
    lambdac4 = atan(((3*asa.cr/4)-(3*asa.cp/4))/(asa.b/2));
    AR_eff_asa = asa.ar/(-(1/0.35^2)*((asa.h/asa.b)-0.35)^2 + 1) ;
    k_a = 1/AR_eff_asa - 1/(1+AR_eff_asa^1.7);
    k_lambda = (10-3*lambda)/7;
    k_h = (1-abs(0.16/asa.b))/((2*lht/asa.b)^(1/3));
    de_dalfa = 4.44*(((k_a*k_lambda*k_h)*(cos(lambdac4)^0.5))^1.19);
    
    e0 = de_dalfa*(alfa_cl0) ;
    
    downwash = e0 + de_dalfa*(alfa+angulo_incidencia_asa) ; % Sabendo que downwash = 0 para o angulo de Cl = 0 da asa.
    
    %-------------------------------------------------------------------------
    
    matriz_rotacao = [cosd(alfa) sind(alfa) ; -sind(alfa) cosd(alfa)];
    vetor_posicao_ca_tdp_rotacionado = vetor_posicao_ca_tdp*matriz_rotacao;
    vetor_posicao_cg_tdp_rotacionado = vetor_posicao_cg_tdp*matriz_rotacao;
    vetor_posicao_ca_emp_tdp_rotacionado = vetor_posicao_ca_emp_tdp*matriz_rotacao;
    %
    asa.h = 0.045+ vetor_posicao_ca_tdp_rotacionado(2);
    emp.h = 0.045+ vetor_posicao_ca_emp_tdp_rotacionado(2);
    
    if X_speed >= 1.05*Vstall %Velocidade mínima para iniciar o comando do profundor
        
        %perfilemp = struct ('a0', 7.037, 'angulo0', 7.05 , 'tcmax', 0.0914, 'cdminw', 0.02782, 'clmax', -1.5);
        
        
        [CL_emp,CD_emp,Cm_emp] = coeficientes_empenagem(angulo_incidencia_emp+alfa-downwash,beta2,empenagem) ; %esta função utiliza as curvas obtidas no xflr5 e interpola para encontrar os dados para qualquer deflexão beta
        F_lift_emp = (0.5 * ro * X_speed^2)* n * emp.s * CL_emp ;
        F_drag_emp = (0.5 * ro * X_speed^2)* n * emp.s * CD_emp ;
        
        [CL_asa,~,Cm_asa] = coeficientes_asa(angulo_incidencia_asa+alfa,X_speed) ; %curvas no xflr5
        CD_asa = -3E-5*alfa^3 + 0.0009*alfa^2 + 0.0018*alfa + 0.0712  ; %Madrugada 15
        
    else
        
        [CL_emp,~,Cm_emp] = coeficientes_empenagem(angulo_incidencia_emp+alfa-downwash,beta1,empenagem) ;
        
        CD_emp =  -3E-05*X_speed^3 + 0.0008*X_speed^2 - 0.0073*X_speed + 0.0341 ;
        
        F_lift_emp = (0.5 * ro * X_speed^2)* n * emp.s * CL_emp ;
        F_drag_emp = (0.5 * ro * X_speed^2)* n * emp.s * CD_emp ;
        
            [CL_asa,CD_asa,Cm_asa] = coeficientes_asa(angulo_incidencia_asa+alfa,X_speed) ; %curvas no xflr5
    end
    

    
    
    time = time+dt ;
    
    X_speed = X_speed + X_acceleration*dt;   %Integração da velocidade
    X_position = X_position + X_speed*dt;   %Integração da posição
    X_acceleration = (F_thrust*cosd(alfa) -(F_drag)  - (-F_y*mi) )/mass_acft; %Resultantes das forças no eixo horizontal, divido pela massa da aeronave.
    
    
    %------------------------- Pitch Rate -------------------------------------
    
    
    alfa = alfa + pitch_rate*dt; %Integração do ângulo
    
    if (alfa+angulo_incidencia_asa) > alfaestol
        alfa =  alfaestol-angulo_incidencia_asa;
        
    end
    
    Pitch_acceleration = M_tdp/Iy_tdp;              %Aceleração angular
    pitch_rate = pitch_rate + rad2deg(Pitch_acceleration)*dt; %Velocidade angular
    
    
    
    %     F_thrust = -4e-6*X_speed^4 + 0.0004*X_speed^3 - 0.0445*X_speed^2 - 0.0631*X_speed + 41.453; %Tração dinâmica. Resultados obtidos experimentalmente, e formula obtida através de um fit-data no excel.
    F_thrust = -4e-6*X_speed^4 + 0.0004*X_speed^3 - 0.0391*X_speed^2 - 0.0554*X_speed + 36.43 ; % Tração rô 1.07
%     F_thrust =  -4E-06*X_speed^4 + 0.0004*X_speed^3 - 0.0425*X_speed^2 - 0.0602*X_speed + 39.573 ; %Tração 1.15
    F_lift_asa = 0.5 * ro * X_speed^2 * asa.s * CL_asa;
    
    F_drag_asa = 0.5 * ro * X_speed^2 * asa.s * CD_asa ;
    
    
    %----------------------%Somatório dos momentos-----------------------------
    
    M_tdp = -F_lift_asa*vetor_posicao_ca_tdp_rotacionado(1)...
        + F_drag_asa*vetor_posicao_ca_tdp_rotacionado(2)...
        + W*vetor_posicao_cg_tdp_rotacionado(1)...
        + (-F_lift_emp)*vetor_posicao_ca_emp_tdp_rotacionado(1)...
        + F_drag_emp*vetor_posicao_ca_emp_tdp_rotacionado(2)...
        - F_thrust*(vetor_posicao_motor_tdp(2))...
        + mass_acft*X_acceleration*(vetor_posicao_cg_tdp_rotacionado(2))...
        + Cm_asa*ro*asa.s*asa.mac*X_speed^2/2 + n*Cm_emp*ro*emp.s*emp.mac*X_speed^2/2;   %n é a eficiência de cauda
    
    
    
    
    i = i+1;
    Momento(i) =M_tdp;
    
    
    
    
    %--------------------------------------------------------------------------
    
    F_drag = F_drag_asa+F_drag_emp;
    F_y = F_lift_asa + F_thrust*sind(alfa) - W - (-F_lift_emp)  ; %Resultantes das forças no eixo vertical.
    
    % Plotagens
    Velocidade(i) = X_speed;
    Posicao(i) = X_position;
    Aceleracao(i) = X_acceleration;
    Forca_vertical(i) = F_y;
    Arrasto(i) = F_drag_asa;
    Forca_empenagem(i) = F_lift_emp;
    Angulo(i) = alfa ;
    Time(i) = time;
    Arrasto(i) = F_drag ;
    Atrito(i) = -F_y*mi ;
    
    if F_y > 0 && X_speed >= 1.05*Vstall %&& M_tdp >= 0  %Condição de decolagem
        flag = 'ok';
        liftoff = 1;
        X_position;
        
    end
    if X_position > 60
        flag = 'err';
        break
    end
    
    if M_tdp <= 0 && alfa<=0 %Se a aeronave não tiver momento para iniciar a rotação
        M_tdp = 0; %Ela permanece sem momento. A normal da bequilha neutraliza o momento picador (em direção ao solo)
        alfa = 0 ;
    elseif alfa < 0
        alfa = 0 ;
    end
    
end

MTOW  = mass_acft;


R_c = X_speed*(F_thrust*cosd(alfa) -(F_drag))/(mass_acft*9.81); % Razão de subida




% % Gráficos -----------------------------
plot (Posicao, Forca_vertical)
title 'Forca Y'
figure()
plot (Posicao, Velocidade)
title ' Velocidade'
figure()
plot (Time, Momento)
title ' momento'
figure()
plot (Posicao, Angulo)
title ' Alfa'
figure()
plot (Posicao, Arrasto)
title ' Arrasto'
figure()
plot (Posicao, Atrito)
title ' Atrito'
figure()
plot (Velocidade, Arrasto)
title 'Velocidade X Arrasto'


end

