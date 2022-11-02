
Nx = 213;

Ske_t = zeros(6);
Ske_v = zeros(6);

Sme_t = zeros(6);
Sme_v = zeros(6);

SK = zeros (3 * Nx);
SM = zeros (3 * Nx);

rho = 7600;
E = 210*10^9;
shift = 0;


for index = 1:Nx-1-6
    
    i1 = index + shift;
    
    if ismember(i1,[71, 116, 131, 162, 179, 196])
        shift = shift + 1;
    end
    
    i1 = index + shift;
    i2 = index + shift + 1;
   
    
    %identifica o caso em que cada limite se encontra
    if i1 < 71
       caso = 1;
    
    elseif i1< 116
        caso =2;
    
    elseif i1 < 162
        caso = 3;
    
    else
        caso = 4;
        
    end
    
    
    if caso == 1
        %definicao das matrizes dos elementos
        Ske_t = zeros(6);
        Ske_v = zeros(6);
        Sme_t = zeros(6);
        Sme_v = zeros(6);
        
        L = 0.5;
        A = pi*(0.07); %pi*(0.4^2) - pi*(0.3^2)
        I = 0.0137444679;
        
        var_t = E*A/L;
        Ske_t(1,1) = var_t;
        Ske_t(1,4) = -var_t;
        Ske_t(4,1) = -var_t;
        Ske_t(4,4) = var_t;
        
        var_viga = E*I/L^3;
        Ske_v(2,2) = var_viga*(12);
        Ske_v(2,3) = var_viga*(6*L);
        Ske_v(2,5) = var_viga*(-12);
        Ske_v(2,6) = var_viga*(6*L);
        Ske_v(3,2) = var_viga*(6*L);
        Ske_v(3,3) = var_viga*(4*L^2);
        Ske_v(3,5) = var_viga*(-6*L);
        Ske_v(3,6) = var_viga*(2*L^2);
        
        Ske_v(5,2) = var_viga*(-12);
        Ske_v(5,3) = var_viga*(-6*L);
        Ske_v(5,5) = var_viga*(12);
        Ske_v(5,6) = var_viga*(-6*L);
        Ske_v(6,2) = var_viga*(6*L);
        Ske_v(6,3) = var_viga*(2*L^2);
        Ske_v(6,5) = var_viga*(-6*L);
        Ske_v(6,6) = var_viga*(4*L^2);
        
        Ske = Ske_t + Ske_v;
        
        var_t = rho*A*L/6;
        Sme_t(1,1) = 2*var_t;
        Sme_t(1,4) = 1*var_t;
        Sme_t(4,1) = 1*var_t;
        Sme_t(4,4) = 2*var_t;
        
        var_viga = rho*A*L;
        rg = sqrt(I/A);
        Sme_v(2,2) = var_viga*(13/35 + (6*rg^2)/(5*L^2) );
        Sme_v(2,3) = var_viga*(11*L/210 + (rg^2)/(10*L));
        Sme_v(2,5) = var_viga*(9/70 + (-6*rg^2)/(5*L^2));
        Sme_v(2,6) = var_viga*(-13*L/420 + (rg^2)/(10*L));
        Sme_v(3,2) = var_viga*(11*L/210 + (rg^2)/(10*L));
        Sme_v(3,3) = var_viga*(L^2/105 + (2*rg^2)/(15));
        Sme_v(3,5) = var_viga*(13*L/420 + (-rg^2)/(10*L));
        Sme_v(3,6) = var_viga*(-L^2/140 + (-rg^2)/(30));
        
        Sme_v(5,2) = var_viga*(9/70 + (-6*rg^2)/(5*L^2));
        Sme_v(5,3) = var_viga*(13*L/420 + (-rg^2)/(10*L));
        Sme_v(5,5) = var_viga*(13/35 + (6*rg^2)/(5*L^2) );
        Sme_v(5,6) = var_viga*(-11*L/210 - (rg^2)/(10*L));
        Sme_v(6,2) = var_viga*(-13*L/420 + (rg^2)/(10*L));
        Sme_v(6,3) = var_viga*(-L^2/140 + (-rg^2)/(30));
        Sme_v(6,5) = var_viga*(-11*L/210 - (rg^2)/(10*L));
        Sme_v(6,6) = var_viga*(L^2/105 + (2*rg^2)/(15));
        
        Sme = Sme_t + Sme_v;
        
        T = gera_rot(90);
        Sme = T' * Sme * T;
        Ske = T' * Ske * T;
        
        SK(3*i1-2:3*i2,3*i1-2:3*i2)=SK(3*i1-2:3*i2,3*i1-2:3*i2)+Ske;
        SM(3*i1-2:3*i2,3*i1-2:3*i2)=SM(3*i1-2:3*i2,3*i1-2:3*i2)+Sme;
        
    elseif caso == 2
        %definicao das matrizes dos elementos
        Ske_t = zeros(6);
        Ske_v = zeros(6);
        Sme_t = zeros(6);
        Sme_v = zeros(6);
        
        
        L = 0.5*sqrt(2);
        A = 0.15; %0.3*0.65 - 0.1*0.45 
        I = 0.00610625;
        
        %trelica
        var_t = E*A/L;
        Ske_t(1,1) = var_t;
        Ske_t(1,4) = -var_t;
        Ske_t(4,1) = -var_t;
        Ske_t(4,4) = var_t;  
        
        %monatgem da parte de rigidez como viga
        var_viga = E*I/L^3;
        Ske_v(2,2) = var_viga*(12);
        Ske_v(2,3) = var_viga*(6*L);
        Ske_v(2,5) = var_viga*(-12);
        Ske_v(2,6) = var_viga*(6*L);
        Ske_v(3,2) = var_viga*(6*L);
        Ske_v(3,3) = var_viga*(4*L^2);
        Ske_v(3,5) = var_viga*(-6*L);
        Ske_v(3,6) = var_viga*(2*L^2);
        
        Ske_v(5,2) = var_viga*(-12);
        Ske_v(5,3) = var_viga*(-6*L);
        Ske_v(5,5) = var_viga*(12);
        Ske_v(5,6) = var_viga*(-6*L);
        Ske_v(6,2) = var_viga*(6*L);
        Ske_v(6,3) = var_viga*(2*L^2);
        Ske_v(6,5) = var_viga*(-6*L);
        Ske_v(6,6) = var_viga*(4*L^2);
        
        Ske = Ske_t + Ske_v; %soma matriz de rididez
        
        
        var_t = rho*A*L/6;
        Sme_t(1,1) = 2*var_t;
        Sme_t(1,4) = 1*var_t;
        Sme_t(4,1) = 1*var_t;
        Sme_t(4,4) = 2*var_t;
        
        %matriz de massa como viga
        var_viga = rho*A*L;
        rg = sqrt(I/A);
        Sme_v(2,2) = var_viga*(13/35 + (6*rg^2)/(5*L^2) );
        Sme_v(2,3) = var_viga*(11*L/210 + (rg^2)/(10*L));
        Sme_v(2,5) = var_viga*(9/70 + (-6*rg^2)/(5*L^2));
        Sme_v(2,6) = var_viga*(-13*L/420 + (rg^2)/(10*L));
        Sme_v(3,2) = var_viga*(11*L/210 + (rg^2)/(10*L));
        Sme_v(3,3) = var_viga*(L^2/105 + (2*rg^2)/(15));
        Sme_v(3,5) = var_viga*(13*L/420 + (-rg^2)/(10*L));
        Sme_v(3,6) = var_viga*(-L^2/140 + (-rg^2)/(30));
        
        Sme_v(5,2) = var_viga*(9/70 + (-6*rg^2)/(5*L^2));
        Sme_v(5,3) = var_viga*(13*L/420 + (-rg^2)/(10*L));
        Sme_v(5,5) = var_viga*(13/35 + (6*rg^2)/(5*L^2) );
        Sme_v(5,6) = var_viga*(-11*L/210 - (rg^2)/(10*L));
        Sme_v(6,2) = var_viga*(-13*L/420 + (rg^2)/(10*L));
        Sme_v(6,3) = var_viga*(-L^2/140 + (-rg^2)/(30));
        Sme_v(6,5) = var_viga*(-11*L/210 - (rg^2)/(10*L));
        Sme_v(6,6) = var_viga*(L^2/105 + (2*rg^2)/(15));
        
       %montagem matriz total
       Sme = Sme_t + Sme_v;
       
       
       T = gera_rot(45);
       Sme = T' * Sme * T;
       Ske = T' * Ske * T;
       
       SK(3*i1-2:3*i2,3*i1-2:3*i2)=SK(3*i1-2:3*i2,3*i1-2:3*i2)+Ske;
       SM(3*i1-2:3*i2,3*i1-2:3*i2)=SM(3*i1-2:3*i2,3*i1-2:3*i2)+Sme;

        
    elseif caso == 3
        %definicao das matrizes dos elementos
        Ske_t = zeros(6);
        Ske_v = zeros(6);
        Sme_t = zeros(6);
        Sme_v = zeros(6);        
        
        L = 0.5;
        A = 0.2; %0.85*0.35 - 0.65*0.15 
        I = 0.0144791667;
        
        %matriz rigidez treliça
        var_t = E*A/L;
        Ske_t(1,1) = var_t;
        Ske_t(1,4) = -var_t;
        Ske_t(4,1) = -var_t;
        Ske_t(4,4) = var_t;
        
        %matriz rigidez viga
        var_viga = E*I/L^3;
        Ske_v(2,2) = var_viga*(12);
        Ske_v(2,3) = var_viga*(6*L);
        Ske_v(2,5) = var_viga*(-12);
        Ske_v(2,6) = var_viga*(6*L);
        Ske_v(3,2) = var_viga*(6*L);
        Ske_v(3,3) = var_viga*(4*L^2);
        Ske_v(3,5) = var_viga*(-6*L);
        Ske_v(3,6) = var_viga*(2*L^2);
        
        Ske_v(5,2) = var_viga*(-12);
        Ske_v(5,3) = var_viga*(-6*L);
        Ske_v(5,5) = var_viga*(12);
        Ske_v(5,6) = var_viga*(-6*L);
        Ske_v(6,2) = var_viga*(6*L);
        Ske_v(6,3) = var_viga*(2*L^2);
        Ske_v(6,5) = var_viga*(-6*L);
        Ske_v(6,6) = var_viga*(4*L^2);
        
        %soma das duas partes
        Ske = Ske_t + Ske_v;
        
        %matriz de massa treliça
        var_t = rho*A*L/6;
        
        Sme_t(1,1) = 2*var_t;
        Sme_t(1,4) = var_t;
        Sme_t(4,1) = var_t;
        Sme_t(4,4) = 2*var_t;
        
        %matriz de massa como viga
        var_viga = rho*A*L;
        rg = sqrt(I/A);
        Sme_v(2,2) = var_viga*(13/35 + (6*rg^2)/(5*L^2) );
        Sme_v(2,3) = var_viga*(11*L/210 + (rg^2)/(10*L));
        Sme_v(2,5) = var_viga*(9/70 + (-6*rg^2)/(5*L^2));
        Sme_v(2,6) = var_viga*(-13*L/420 + (rg^2)/(10*L));
        Sme_v(3,2) = var_viga*(11*L/210 + (rg^2)/(10*L));
        Sme_v(3,3) = var_viga*(L^2/105 + (2*rg^2)/(15));
        Sme_v(3,5) = var_viga*(13*L/420 + (-rg^2)/(10*L));
        Sme_v(3,6) = var_viga*(-L^2/140 + (-rg^2)/(30));
        
        Sme_v(5,2) = var_viga*(9/70 + (-6*rg^2)/(5*L^2));
        Sme_v(5,3) = var_viga*(13*L/420 + (-rg^2)/(10*L));
        Sme_v(5,5) = var_viga*(13/35 + (6*rg^2)/(5*L^2) );
        Sme_v(5,6) = var_viga*(-11*L/210 - (rg^2)/(10*L));
        Sme_v(6,2) = var_viga*(-13*L/420 + (rg^2)/(10*L));
        Sme_v(6,3) = var_viga*(-L^2/140 + (-rg^2)/(30));
        Sme_v(6,5) = var_viga*(-11*L/210 - (rg^2)/(10*L));
        Sme_v(6,6) = var_viga*(L^2/105 + (2*rg^2)/(15));
        
        %soma matriz de massa
        Sme = Sme_t + Sme_v;
        
        SK(3*i1-2:3*i2,3*i1-2:3*i2)=SK(3*i1-2:3*i2,3*i1-2:3*i2)+Ske;
        SM(3*i1-2:3*i2,3*i1-2:3*i2)=SM(3*i1-2:3*i2,3*i1-2:3*i2)+Sme;
       
    elseif caso == 4
        %definicao das matrizes dos elementos
        Ske_t = zeros(6);
        Ske_v = zeros(6);
        Sme_t = zeros(6);
        Sme_v = zeros(6);        
        
        L = 0.5;
        A = 0.07; %0.2*0.35
        I = 0.0002333333;
        
        %matriz rigidez treliça
        var_t = E*A/L;
        Ske_t(1,1) = var_t;
        Ske_t(1,4) = -var_t;
        Ske_t(4,1) = -var_t;
        Ske_t(4,4) = var_t;
        
        %matriz rigidez viga
        var_viga = E*I/L^3;
        Ske_v(2,2) = var_viga*(12);
        Ske_v(2,3) = var_viga*(6*L);
        Ske_v(2,5) = var_viga*(-12);
        Ske_v(2,6) = var_viga*(6*L);
        Ske_v(3,2) = var_viga*(6*L);
        Ske_v(3,3) = var_viga*(4*L^2);
        Ske_v(3,5) = var_viga*(-6*L);
        Ske_v(3,6) = var_viga*(2*L^2);
        
        Ske_v(5,2) = var_viga*(-12);
        Ske_v(5,3) = var_viga*(-6*L);
        Ske_v(5,5) = var_viga*(12);
        Ske_v(5,6) = var_viga*(-6*L);
        Ske_v(6,2) = var_viga*(6*L);
        Ske_v(6,3) = var_viga*(2*L^2);
        Ske_v(6,5) = var_viga*(-6*L);
        Ske_v(6,6) = var_viga*(4*L^2);
        
        %soma das duas partes
        Ske = Ske_t + Ske_v;
        
        %matriz de massa treliça
        var_t = rho*A*L/6;
        
        Sme_t(1,1) = 2*var_t;
        Sme_t(1,4) = var_t;
        Sme_t(4,1) = var_t;
        Sme_t(4,4) = 2*var_t;
        
        %matriz de massa como viga
        var_viga = rho*A*L;
        rg = sqrt(I/A);
        Sme_v(2,2) = var_viga*(13/35 + (6*rg^2)/(5*L^2) );
        Sme_v(2,3) = var_viga*(11*L/210 + (rg^2)/(10*L));
        Sme_v(2,5) = var_viga*(9/70 + (-6*rg^2)/(5*L^2));
        Sme_v(2,6) = var_viga*(-13*L/420 + (rg^2)/(10*L));
        Sme_v(3,2) = var_viga*(11*L/210 + (rg^2)/(10*L));
        Sme_v(3,3) = var_viga*(L^2/105 + (2*rg^2)/(15));
        Sme_v(3,5) = var_viga*(13*L/420 + (-rg^2)/(10*L));
        Sme_v(3,6) = var_viga*(-L^2/140 + (-rg^2)/(30));
        
        Sme_v(5,2) = var_viga*(9/70 + (-6*rg^2)/(5*L^2));
        Sme_v(5,3) = var_viga*(13*L/420 + (-rg^2)/(10*L));
        Sme_v(5,5) = var_viga*(13/35 + (6*rg^2)/(5*L^2) );
        Sme_v(5,6) = var_viga*(-11*L/210 - (rg^2)/(10*L));
        Sme_v(6,2) = var_viga*(-13*L/420 + (rg^2)/(10*L));
        Sme_v(6,3) = var_viga*(-L^2/140 + (-rg^2)/(30));
        Sme_v(6,5) = var_viga*(-11*L/210 - (rg^2)/(10*L));
        Sme_v(6,6) = var_viga*(L^2/105 + (2*rg^2)/(15));
        
        %soma matriz de massa
        Sme = Sme_t + Sme_v;
        
        SK(3*i1-2:3*i2,3*i1-2:3*i2)=SK(3*i1-2:3*i2,3*i1-2:3*i2)+Ske;
        SM(3*i1-2:3*i2,3*i1-2:3*i2)=SM(3*i1-2:3*i2,3*i1-2:3*i2)+Sme;
    
    end 
end

%ajustando os valores dos nós em comum que sao repetidos 
%41 - 72 // CASO 2
%55 - 117 //CASO 3
%71 - 132 //CASO 3
%86 - 131 - 163 //CASO 3
%102 - 162 - 180 
%116 - 197
matriz_corr = [41 72; 55 117; 71 132; 86 131; 86 163; 102 162; 102 180; 116 197];

for tup = 1:length(matriz_corr(:,1))

    i1 = matriz_corr(tup,1);
    i1_new = matriz_corr(tup,2);

    %arruma diagonal principal
    SK(3*i1-2:3*i1,3*i1-2:3*i1)=SK(3*i1-2:3*i1,3*i1-2:3*i1)+ SK(3*i1_new-2:3*i1_new,3*i1_new-2:3*i1_new);
    SM(3*i1-2:3*i1,3*i1-2:3*i1)=SM(3*i1-2:3*i1,3*i1-2:3*i1)+ SM(3*i1_new-2:3*i1_new,3*i1_new-2:3*i1_new);
    
    %arruma interacoes de nos distantes
    SK(3*i1-2:3*i1,3*(i1_new+1)-2:3*(i1_new+1))=SK(3*i1-2:3*i1,3*(i1_new+1)-2:3*(i1_new+1)) + SK(3*i1_new-2:3*i1_new,3*(i1_new+1)-2:3*(i1_new+1));
    SM(3*i1-2:3*i1,3*(i1_new+1)-2:3*(i1_new+1))=SM(3*i1-2:3*i1,3*(i1_new+1)-2:3*(i1_new+1)) + SM(3*i1_new-2:3*i1_new,3*(i1_new+1)-2:3*(i1_new+1));

    %arruma interacoes de nos distantes 2
    SK(3*(i1_new+1)-2:3*(i1_new+1),3*i1-2:3*i1)=SK(3*(i1_new+1)-2:3*(i1_new+1), 3*i1-2:3*i1) + SK(3*(i1_new+1)-2:3*(i1_new+1),3*i1_new-2:3*i1_new);
    SM(3*(i1_new+1)-2:3*(i1_new+1),3*i1-2:3*i1)=SM(3*(i1_new+1)-2:3*(i1_new+1), 3*i1-2:3*i1) + SM(3*(i1_new+1)-2:3*(i1_new+1),3*i1_new-2:3*i1_new);

end

vector_delete = [];

for value=matriz_corr(:,2)

    vector_delete = [vector_delete, value*3-2];
    vector_delete = [vector_delete, value*3-1];
    vector_delete = [vector_delete, value*3];
    
end

vector_delete = reshape(vector_delete, [24,1]);


%deletando as colunas dos nos repetidos
SK(vector_delete,:) = [];
SM(vector_delete,:) = [];

SK(:,vector_delete) = [];
SM(:,vector_delete) = [];


%ajuste para as condições de contorno (engaste no ponto 1)
SK(:,1) = zeros(length(SK),1);
SK(1,:) = zeros(1,length(SK));
SK(1,1) = 1;

SK(2,2:length(SK)) = zeros(1,length(SK)-1);
SK(2:length(SK),2) = zeros(length(SK)-1,1);
SK(2,2) = 1;

SK(3,3:length(SK)) = zeros(1,length(SK)-2);
SK(3:length(SK),3) = zeros(length(SK)-2,1);
SK(3,3) = 1;

%matriz de massa
SM(:,1) = zeros(length(SM),1);
SM(1,:) = zeros(1,length(SM));
SM(1,1) = 1;

SM(2,2:length(SM)) = zeros(1,length(SM)-1);
SM(2:length(SM),2) = zeros(length(SM)-1,1);
SM(2,2) = 1;

SM(3,3:length(SM)) = zeros(1,length(SM)-2);
SM(3:length(SM),3) = zeros(length(SM)-2,1);
SM(3,3) = 1;


%ANALISE MODAL:
%obtendo os autovalores
[D,R]=eig(SK,SM);

%rescalando os autovalores
for i = 1:length(SK)
 valuei=D(:,i)'*SM*D(:,i);
 D(:,i)=D(:,i)/sqrt(valuei);
end

freq = sqrt(sort(unique(diag(R)))) / (2 * pi);
freq = sqrt(diag(R))/(2*pi);

function [T] = gera_rot(theta)

    ang = theta*pi/180;

    T = [cos(ang) sin(ang) 0 0 0 0;
         -sin(ang) cos(ang) 0 0 0 0;
         0 0 1 0 0 0;
         0 0 0 cos(ang) sin(ang) 0;
         0 0 0 -sin(ang) cos(ang) 0;
         0 0 0 0 0 1];

end





