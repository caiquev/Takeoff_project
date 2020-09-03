function [M_tdp,Vstall,X_speed,MTOW,liftoff,alfa,X_position,alfaestol,R_c] = doisif(mass_acft,angulo_incidencia_asa,angulo_incidencia_emp,tdp_x,empenagem,beta1,beta2)
%Algoritmo de corrida de decolagem, baseada em uma integra��o das for�as da
%corrida em um timestep fixo. O algoritmo precisa, para funcionar, das
%fun��es que determinam o coeficiente de sustenta��o e arrasto, suas
%varia��es com respeito a alfa, e tamb�m da tra��o din�mica do motor. A
%sa�da do algoritmo � valor X_position que determina o comprimento de pista
%percorrido
%
close all
% clear
% clc

% Constantes
vento = 0 ; % Velocidade do vento [m/s]
dt = 0.001;        % Timestep [s]
ro = 1.15;        % Densidade do ar [kg/m3]
mi = 0.07;      % Coef. de atrito
n = sqrt(0.9) ;        %Efici�ncia de cauda
inclinacao_pista = 3 ; %Inclina��o da pista de corrida de decolagem em SJC (PRECISA PEGAR ESSE DADO DA PISTA DE JUA D� PRA PEGAR PELO GOOGLEMAPS?!?!) [graus]
alfaestol = 13.5 ; %ASA AR 5 [graus]
% alfaestol = 14 ; %Asa Oficina

% alfaestol = 12.5 ; %VANT
% alfaestol = 14.5 ; %ASA MADRUGADA
% alfaestol  = 13.5 ; %asa BRUXO


%--------------------------------------------------
%Valores iniciais para a itera��o



i = 1; %Contador
liftoff = 0;      %  Flag de decolagem como FALSO
time = 0;         % tempo de simula��o [s]
X_speed = 0;      %Come�ar com um valor muito baixo pra evitar alguns erros por causa da aerod n�o estacion�ria.
X_position = 0;   %Posi��o Inicial do Avi�o [m]
X_acceleration = 0; %Acelera��o inicial do avi�o [m/s2]
alfa = 0;         %  Nivel do aviao eixo lateral  [graus]
pitch_rate = 0;   % Velocidade de arfagem [graus/s]
Pitch_acceleration = 0; % Acelera��o angular de arfagem [graus/s2]
M_tdp = 0;        % Momento Resultante em rela��o ao TDP [N*m]
F_y = -mass_acft*9.81 ; %Valor inicial da For�a resultante no eixo y (for�a peso) .
F_thrust = 0;     % Tra��o Inicial [N]
F_drag = 0;       %For�a de Arrasto da asa/emp, posteriormente atualizado com valores para aeronave completa. [N]
R_c = 0 ;         % Raz�o de subida logo ap�s a decolagem 
%--------------------------------------------------

% Dados Geom�tricos da Asa e Horizontal


%b=envergadura[m] / perc_retan = percentual retangular [0 a 1] / cr = corda
%na raiz [m] / cp = corda na ponta [m] / ar = raz�o de aspecto / lambdac2 =
%afilamento do meio da corda / s = �rea da superf�cie [m2] / mac = corda
%m�dia aerodin�mica [m] / e = efici�ncia de cauda.


asa = struct ('b', 2.45, 'perc_retan', 0.733, 'cr', 0.50, 'cp', 0.3,'ar', 5, 'lambdac2', 9.75 , 's', 1.2, 'mac', 0.503, 'e', sqrt(0.9)); 

emp = struct ('b', 1, 'perc_retan', 0, 'cr', 0.36 , 'cp', 0.235, 'ar', 3.3, 'lambdac2', -6.5, 's', 0.2952, 'mac', 0.301, 'e', sqrt(0.9));



W = mass_acft*9.81; % For�a peso [N]

asa.h = 0.13;       % Altura da asa com rela��o ao solo [m]

emp.h = 0.23;       % Altura da empenagem com rela��o ao solo [m]

lht = 0.988; %Dist�ncia CA da asa at� o CA de empenagem [m]
vht = 0.52 ; % Volume de cauda 

%--------------------------------------------------------------------------
%Estima as dimens�es do compartimento de carga para estimar um valor de
%momento de in�rcia.

% x = 0.13;         %largura compartimento de carga [m]
% y = 0.2;          %comprimento compartimento de carga [m]
% z = 0.07;         %altura compartimento de carga [m]
%-----------------------------------------------------------------

raio_roda = 0.055; % Raio da roda do avi�o [m]

%--- Dist�ncias dos pontos de aplica��o de for�a ao tdp nos eixos x e y.

vetor_posicao_ca_tdp = [-tdp_x asa.h-raio_roda]; % C.A.- Trem de Pouso [m].
vetor_posicao_cg_tdp = [0.064-tdp_x 0.13]; %0.077 � o CA-CG AR 5 0.06 AR asa toni-TDP [m].
vetor_posicao_ca_emp_tdp = [lht+(0.06-tdp_x) emp.h-raio_roda]; % C.A. da SH-TDP [m].
vetor_posicao_motor_tdp = [0.507 0.18-0.055]; % Motor-TDP [m].


d = sqrt(asa.h^2 + vetor_posicao_cg_tdp(1)^2);      %distancia resultante do c.g. ao tdp no plano x-y.
Iy = 0.409434721;  %Momento de in�rcia do compartimento de carga em rela��o ao c.g. [kg*m2]
Iy_tdp = Iy + (mass_acft)*d^2;   %Momento de in�rcia do compartimento de carga em rela��o ao centro da roda do tdp principal. [kg*m2] Eixos paralelos - Teorema Steiner


%------------------------- C�lculo do downwash --------------------------

CL_asa_inc = coeficientes_asa(angulo_incidencia_asa,X_speed) ; % Coeficiente de sustenta��o da asa no �ngulo de inci�ncia da asa durante a corrida.

lambda = asa.cp/asa.cr ; % Taper ratio 
lambdac4 = atan(((3*asa.cr/4)-(3*asa.cp/4))/(asa.b/2));  % One Quarter chord Taper ratio 
AR_eff_asa = asa.ar/(-(1/0.35^2)*((asa.h/asa.b)-0.35)^2 + 1) ; % Raz�o de aspecto virtual da asa devido o efeito solo.

k_a = 1/AR_eff_asa - 1/(1+AR_eff_asa^1.7); % Par�metros da equa��o que calcula o downwash causado pela asa
k_lambda = (10-3*lambda)/7; % Par�metros da equa��o que calcula o downwash causado pela asa
k_h = (1-abs(0.16/asa.b))/((2*lht/asa.b)^(1/3)); % Par�metros da equa��o que calcula o downwash causado pela asa
de_dalfa = 4.44*(((k_a*k_lambda*k_h)*(cos(lambdac4)^0.5))^1.19) ; % Derivada do downwash em rela��o a alfa [graus?]

e0 = 2*CL_asa_inc*57.3/(pi*AR_eff_asa) ; %Downwash para alfa = 0. [graus?]

downwash = e0 + de_dalfa*(alfa) ; % Downwash resultante. [graus?]

%-------------------------------------------------------------------------


[CL_emp_decolagem,~,~] = coeficientes_empenagem(angulo_incidencia_emp+(alfaestol-angulo_incidencia_asa)+downwash,beta2,empenagem,X_speed); %CL da empenagem quando o avi�o arfou completamente


CL_max = coeficientes_asa(alfaestol,X_speed) + (CL_emp_decolagem*vht); % CL m�x do avi�o.

Vstall = sqrt(2*W/(ro* asa.s * CL_max)); % Velocidade de stall [m/s].


while liftoff == 0 % In�cio da corrida de decolagem.
    
    
    matriz_rotacao = [cosd(alfa) sind(alfa) ; -sind(alfa) cosd(alfa)]; %Matriz de rota��o utilizada para girar a aeronave.
    
    vetor_posicao_ca_tdp_rotacionado = vetor_posicao_ca_tdp*matriz_rotacao; % Vetor posi��o CA-tdp rotacionado.
    vetor_posicao_cg_tdp_rotacionado = vetor_posicao_cg_tdp*matriz_rotacao; % Vetor posi��o cg-tdp rotacionado.
    vetor_posicao_ca_emp_tdp_rotacionado = vetor_posicao_ca_emp_tdp*matriz_rotacao; % Vetor posi��o CA da SH-tdp rotacionado.
    %
    asa.h = 0.045+ vetor_posicao_ca_tdp_rotacionado(2); % Altura do C.A da asa em rela��o ao solo. [m]
    emp.h = 0.045+ vetor_posicao_ca_emp_tdp_rotacionado(2); % Altura do C.A da SH em rela��o ao solo. [m]
    
    
    time = time+dt ; % Contando o tempo de simula��o [s]
    
    F_x = F_thrust*cosd(alfa) -(F_drag)  - (-F_y*mi) - W*sind(inclinacao_pista) ; %For�a resultante eixo x [N] 
    
    X_acceleration = F_x/mass_acft; %Acelera��o instant�nea [m/s2]
    
        X_speed = X_speed + X_acceleration*dt ;  %Velocidade instant�nea [m/s2]
        X_position = X_position + X_speed*dt;   %Posi��o instant�nea [m]

    
    
    if  X_speed >= 1.05*Vstall %Velocidade m�nima para iniciar o comando do profundor
        
        %Caso - profundor defletido
        
        [CL_emp,CD_emp,Cm_emp] = coeficientes_empenagem(angulo_incidencia_emp+alfa-downwash,beta2,empenagem,X_speed+vento) ; %Fun��o que retorna coeficientes aerod. da SH (graus,inteiro, m/s) o valor inteiro � utilizado para testar diferentes perfis.
        F_lift_emp = (0.5 * ro * (X_speed+vento)^2)* n * emp.s * CL_emp ; %Sustenta��o instant�nea da empenagem  [N]
        F_drag_emp = (0.5 * ro * (X_speed+vento)^2)* n * emp.s * CD_emp ; %Arrasto instant�neo da empenagem  [N]
        
        [CL_asa,CD_asa,Cm_asa] = coeficientes_asa(angulo_incidencia_asa+alfa,X_speed+vento) ; %Fun��o que retorna coeficientes aerod. da SH (graus, m/s) sem a fun��o para armazenar diferentes perfis.
        
    else
        
        %Caso - Profundor incidencia de corrida
        [CL_emp,CD_emp,Cm_emp] = coeficientes_empenagem(angulo_incidencia_emp+alfa-downwash,beta1,empenagem,X_speed+vento) ;
        
        F_lift_emp = (0.5 * ro * (X_speed)^2)* n * emp.s * CL_emp ;
        F_drag_emp = (0.5 * ro * (X_speed)^2)* n * emp.s * CD_emp ;
        
        [CL_asa,CD_asa,Cm_asa] = coeficientes_asa(angulo_incidencia_asa+alfa,X_speed+vento) ; %curvas no xflr5
        
    end
    
    
%Tra��o din�mica. Resultados obtidos experimentalmente, e formula obtida atrav�s de um fit-data no excel.
    
% F_thrust = (-4e-6*(X_speed(i)+vento)^4 + 0.0004*(X_speed(i)+vento)^3 - 0.0391*(X_speed(i)+vento)^2 - 0.0554*(X_speed(i)+vento) + 36.43)*0.95 ; % Tra��o r� 1.07

%     F_thrust =  (-4E-06*(X_speed+vento)^4 + 0.0004*(X_speed+vento)^3 - 0.0425*(X_speed)^2 - 0.0602*X_speed + 39.573)*0.95 ; %Tra��o 1.15 
    F_thrust = -0.00000132046044*(X_speed)^6 + 0.00012517370160*(X_speed)^5 - 0.00420610352997*(X_speed)^4 + 0.05838308056186*(X_speed)^3 - 0.29837734119295*(X_speed)^2 - 0.28406427868039*(X_speed) + 39.94010752332670 ; %Tra��o 1.15 

    
% F_thrust = (-1E-05*(X_speed(i)+vento)^4 + 0.0007*(X_speed(i)+vento)^3 - 0.046*(X_speed(i)+vento)^2 - 0.0405*(X_speed(i)+vento) + 40.126)*0.95 ; %Tra��o 1.15 

    F_lift_asa = 0.5 * ro * (X_speed+vento)^2 * asa.s * CL_asa;  %Sustenta��o instant�nea da asa  [N]
    
    F_drag_asa = 0.5 * ro * (X_speed+vento)^2 * asa.s * CD_asa ; %Arrasto instant�nea da asa  [N]
    

    
    
    %------------------------- Pitch Rate -------------------------------------
   %Limitador m�gico do �ngulo do avi�o para travar o avi�o no angulo de CL m�x que zera o momento resultante e a velocidade angular da aeronave. 
  
    if (alfa+angulo_incidencia_asa) > alfaestol %Verifica se est� no angulo m�ximo possivel fisicamente.
        
        
        alfa =  alfaestol-angulo_incidencia_asa; %Faz travar no angulo m�ximo pois encosta no ch�o.
        pitch_rate = 0 ; %Para de girar.
        M_tdp =  0;
        
        
    else %Est� livre para girar.
        
        
        %----------------------%Momento resultante em rela��o ao tdp-----------------------------
        
        M_tdp = -F_lift_asa*vetor_posicao_ca_tdp_rotacionado(1)...
            + F_drag_asa*vetor_posicao_ca_tdp_rotacionado(2)...
            + W*vetor_posicao_cg_tdp_rotacionado(1)...
            + (-F_lift_emp)*vetor_posicao_ca_emp_tdp_rotacionado(1)...
            + F_drag_emp*vetor_posicao_ca_emp_tdp_rotacionado(2)...
            - F_thrust*cosd(alfa)*(vetor_posicao_motor_tdp(2))...
            + F_thrust*sind(alfa)*(vetor_posicao_motor_tdp(1))...
            + mass_acft*X_acceleration*(vetor_posicao_cg_tdp_rotacionado(2))...
            + Cm_asa*ro*asa.s*asa.mac*X_speed^2/2 + n*Cm_emp*ro*emp.s*emp.mac*X_speed^2/2  ;   %n � a efici�ncia de cauda
        
        %-------------------------------------------------------------------------
        
        
        if M_tdp <= 0  %Se a aeronave n�o tiver momento para iniciar a rota��o
            
            if alfa <= 0
            
                M_tdp = 0; %Ela permanece sem momento. A normal da bequilha neutraliza o momento picador (em dire��o ao solo)
                alfa = 0 ;
                pitch_rate = 0 ; 
                
            end 
            
        end 
           
       
        
        
        %---------------- Depois de verificar se era possivel a
        %aeronave girar devido aos limites f�sicos (bater no ch�o).
        %Incrementa o angulo. Caso n�o fosse poss�vel girar a aeronave M_tdp recebeu 0
        %
        
        alfa = alfa + pitch_rate*dt; % Angulo instant�neo do avi�o [graus]
        
        Pitch_acceleration = M_tdp/Iy_tdp;              %Acelera��o angular instantanea [rad/s2]
        pitch_rate = pitch_rate + rad2deg(Pitch_acceleration)*dt; %Velocidade angular instantanea [graus/s]
        
        
        
    end
    
    
    F_drag = (F_drag_asa+F_drag_emp)*1.07 ;%Aumento de 7% do arrasto total devido aos Ber�o,Fuselagem,TDP,Boom. Valor baseado em historico das aeronaves F-carranca.
    F_y = F_lift_asa + F_thrust*sind(alfa) - W - (-F_lift_emp)  ; %Resultantes das for�as no eixo vertical.
    
%     Plotagens
    Velocidade(i) = X_speed;
    Posicao(i) = X_position;
    Aceleracao(i) = X_acceleration;
    Forca_vertical(i) = F_y;
    Arrasto(i) = F_drag;
    Arrasto_asa(i) = F_drag_asa ;
    Forca_empenagem(i) = F_lift_emp;
    Angulo(i) = alfa ;
    Time(i) = time;
    Atrito(i) = -F_y*mi ;
    Arrasto_emp(i) = F_drag_emp ;
    Momento(i) = M_tdp ;
    forca_emp(i) = F_lift_emp;
%     
    if  F_y >=0  
        if X_speed > 1.05*Vstall
           
        flag = 'ok';
        liftoff = 1
        X_position
        X_speed
    
        end
    end
    
    if X_position > 57
        flag = 'err'
        break
    end
    
    i = i+1;
    
end

MTOW  = mass_acft;

% Ay = F_y/mass_acft
% Ax = (F_thrust*cosd(alfa) -(F_drag))/mass_acft;
% R_c = X_speed(i)*(F_thrust*cosd(alfa) -(F_drag))/(mass_acft*9.81) - Ax*X_speed(i)/9.81  % Raz�o de subida segundo Anderson Eq 6.61
% X_acceleration
% Vstall
%
% (F_thrust*cosd(alfa) -(F_drag))/(mass_acft*9.81)
% Ax*X_speed(i)/9.81

% % Gr�ficos -----------------------------
% plot (Posicao, Forca_vertical)
% title 'Forca Y'
% figure()
% plot (Posicao, Velocidade)
% title ' Velocidade'
% figure()
% plot (Time, Momento)
% title ' momento'
% figure()
% plot (Posicao, Angulo)
% title ' Alfa'
% figure()
% plot (Posicao, Arrasto)
% title ' Arrasto'
% figure()
% plot (Posicao, Atrito)
% title ' Atrito'
% figure()
% plot (Velocidade, Arrasto)
% title 'Velocidade X Arrasto'
% figure()
% plot (Velocidade, Posicao)
% title 'Velocidade x Posicao'
% figure()
% plot (Time, Angulo)
% title ' Alfa x time'
% figure()
% plot (Velocidade, Arrasto_asa)
% title 'Arrasto asa'
% figure()
% plot (Velocidade, Arrasto_emp)
% title 'Arrasto emp'
% figure()
% plot (Velocidade, Time)
% title 'Velocidade (m/s)'
% ylabel 'tempo (s)'
% figure()
% plot (Time, Aceleracao)
% title 'Acelera��o x time'
% figure()
% plot (Time, forca_emp)
% title 'time x forca_emp'
end

