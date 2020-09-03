function [cl, a, alfaestol, clmax] = CL (superficie, perfil, alfa, ~)

% Essa função toma como inputs o necessário para fazer a inclinação CLxAlfa
% da superficie, baseado na eq. baseada no DATCOM, como visto em Airplane Design
% and Aerodynamics - Roskam (pg. 104), que é dito como válida para qualquer
% AR e enflechamento
% 
% Como input, ele necessariamente toma uma estrutura "superficie" que deve
% conter pelo menos a razão de aspecto e o enflechamento de meio corda
% (half_chord_angle) sob "superficie.arw" e "superficie.lambdac2"
% Também tem como input obrigatório alguns dados do perfil, em uma estrutura
% "perfil", no caso, o a0 sob "perfil.a0" e alfa do CL = 0 sob 
% "perfil.angulo0"
% alfa seria um argumento opcional que serviria para calcular o CL de certo
% ângulo a partir de 0. Se não houver este argumento, então alfa = 0

% Como saídas, a função pode dar o CL para o alfa dado, a inclinação
% CLxAlfa calculada, e o CLmax pelo método da seção crítica. NÃO USAR a
% terceira saída simultaneamenete com o efeito solo, a não ser que o
% usuário saiba o que está fazendo

%--------------------------------------------------------------------------
% Ultima alteração realizada (14/02, Caio Guimarães) - Adicionado a
% possibilidade de calcular o CL da asa no solo, em uma determinada altura
% h. Para realizar esse calculo, basta entrar com o "h" como argumento.
% Ler efeito_solo.m para mais detalhes
%--------------------------------------------------------------------------


if nargin < 3 %Caso não seja explicitado um alfa
    alfa = 0;
    deltaa0g = 0;
end
        
ar = superficie.ar;
half_chord_angle = superficie.lambdac2;
a0 = perfil.a0;
angulo0 = perfil.angulo0;

k = a0/(2*pi);

a = (2*pi*ar)/(2 + (((ar^2/k^2)*(1 + ((tand(half_chord_angle))^2))+4))^0.5);
if angulo0 < 0 % É necessario consertar isso. Problemas quando o perfil é invertido.
    cl = (abs(angulo0)+alfa)*deg2rad(a); %Contabilizando o CL a partir do angulo em que CL é igual a zero. Assumindo linearidade do CLxAlfa
else           % Quando o perfil é invertido. Consertar essa linha.
    cl = (-angulo0+alfa)*deg2rad(a);
end

placeholder = a;

if nargin == 4 % Caso a asa esteja no efeito solo
    
    deltaa0g = perfil.tcmax*(-0.1177*(1/(superficie.h/superficie.mac)^2)+3.5655*(1/(superficie.h/superficie.mac))); %Variação no ângulo de cl = 0, baseado na equação 10.10, página 449 do ROSKAM.
    
    %Esses 'fatores foram obtidos a partir de uma "linearização" exponencial,
    %obtida pelo MATLAB, da curva descrita em ROSKAM (Fig 10.11, pg 448)
    
    termo1_efeito_solo = 0.8375;
    termo2_efeito_solo = 0.09214;
    termo3_efeito_solo = -0.7859;
    termo4_efeito_solo = -3.695;
    
    Razao_alongamento_efetivo = termo1_efeito_solo*exp...
        (termo2_efeito_solo*2*superficie.h/superficie.b) + termo3_efeito_solo*exp(termo4_efeito_solo*2*superficie.h/superficie.b);
    
    ar_efetivo = ar/Razao_alongamento_efetivo; %Alongamento efetivo da superfície
    a = (2*pi*ar_efetivo)/(2 + (((ar_efetivo^2/k^2)*(1 + ((tand(half_chord_angle))^2))+4))^0.5);
    
    cl = cl*a/placeholder - a*deg2rad(deltaa0g);

end
    
    
    
if nargout >= 3 %Calculo do CLmax pelo metodo da seção crítica
    clmax = distribuicao_multhopp(superficie, perfil, deg2rad(placeholder));
    alfaestol = clmax/(a/57.3) - abs(angulo0) - deltaa0g;
end

a = placeholder;
a = a/57.3;

end


