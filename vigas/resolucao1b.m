Nx = input('número de nós = ');
leng = input('comprimento da barra = ');
EI = 1;
A = 1;
Row = 1;
Row = Row * A;
Nelx = Nx - 1;
X= [ ];
for i = 1:Nx
 X(i) = leng * (i - 1) / Nelx;
end
%construindo a matriz de massa e rigidez
SK = zeros (2 * Nx);
SM = zeros (2 * Nx);
Ske = zeros (4);
Sme = zeros (4);
for Nel = 1:Nelx
    L = X(Nel + 1) - X(Nel);
    i1 = Nel;
    i2 = Nel + 1;
    Ske(1,1) = 12*EI/L^3;
    Sme(1,1) = 13*L*Row/35;
    Ske(1,2) = 6 * EI / L ^ 2;
    Sme(1,2) = 11 * L ^2 * Row / 210;
    Ske(1,3) = -12 * EI / L ^ 3;
    Sme(1,3) = 9 * L * Row / 70;
    Ske(1,4) = 6 * EI / L ^ 2;
    Sme(1,4) = -13 * L ^2 * Row / 420;
    Ske(2,1) = 6 * EI / L ^ 2;
    Sme(2,1) = 11 * L ^2 * Row / 210;
    Ske(2,2) = 4 * EI / L;
    Sme(2,2) = L ^3 * Row / 105;
    Ske(2,3) = -6 * EI / L ^ 2;
    Sme(2,3) = 13 * L ^2 * Row / 420;
    Ske(2,4) = 2 * EI / L;
    Sme(2,4) = - (L^3*Row)/140;
    Ske(3,1) = -12 * EI / L ^ 3;
    Sme(3,1) = 9 * L * Row / 70;
    Ske(3,2) = -6 * EI / L ^ 2;
    Sme(3,2) = 13 * L ^2 * Row / 420;
    Ske(3,3) = 12 * EI / L ^ 3;
    Sme(3,3) = 13 * L * Row / 35;
    Ske(3,4) = -6 * EI / L ^ 2;
    Sme(3,4) = -11 * L ^2 * Row / 210;
    Ske(4,1) = 6 * EI / L ^ 2;
    Sme(4,1) = -13 * L ^2 * Row / 420;
    Ske(4,2) = 2 * EI / L;
    Sme(4,2) = - L ^3 * Row / 140;
    Ske(4,3) = -6 * EI / L ^ 2;
    Sme(4,3) = -11 * L ^2 * Row / 210;
    Ske(4,4) = 4 * EI / L;
    Sme(4,4) = L ^3 * Row / 105;
    SK(2*i1-1:2*i2,2*i1-1:2*i2)=SK(2*i1-1:2*i2,2*i1-1:2*i2)+Ske;
    SM(2*i1-1:2*i2,2*i1-1:2*i2)=SM(2*i1-1:2*i2,2*i1-1:2*i2)+Sme;
end
%condições de contorno (suporte simples)
SK(:,1)=zeros(2*Nx,1);
SK(1,:)=zeros(1,2*Nx);
SK (1,1)=1;
SK(:,2*Nx-1)=zeros(2*Nx,1);
SK(2*Nx-1,:)=zeros(1,2*Nx);
SK(2*Nx-1,2*Nx-1)=1;
SM(:,1)=zeros(2*Nx,1);
SM(1,:)=zeros(1,2*Nx);
SM(1,1)=1;
SM(:,2*Nx-1)=zeros(2*Nx,1);
SM(2*Nx-1,:)=zeros(1,2*Nx);
SM(2*Nx-1,2*Nx-1)=1;
% Problema de autovalor. Os autovetores são normalizados de forma que a norma 
[D,R]=eig(SK,SM)
%reescalando os autovetores tal que X'*SM*X=I
for i = 1:2*Nx
    valuei=D(:,i)'*SM*D(:,i);
    D(:,i)=D(:,i)/sqrt(valuei);
end
%plota os modos de vibração
for i = 1:2*Nx
    i
    xf = [ ];
    Df = [ ];
    ipx = 20;
    for Nel = 1:Nelx
        de(1) = D(2*Nel-1,i);
        de(2) = D(2*Nel,i);
        de(3) = D(2*Nel+1,i);
        de(4) = D(2*Nel+2,i);
        for ip = 1:ipx+1
            L = X(Nel+1)-X(Nel);
            S= L*(ip-1)/ipx;
            N(1)=1-3*S^2/L^2+2*S^3/L^3;
            N(2)=S-2*S^2/L+S^3/L^2;
            N(3)=3*S^2/L^2-2*S^3/L^3;
            N(4)=-S^2/L+S^3/L^2;
            Df(ipx*(Nel-1)+ip)=de*N';
            xf(ipx*(Nel-1)+ip)=X(Nel)+S;
        end
    end
    plot(xf,Df)
    title('autovetor')
    xlabel('D')
    ylabel('deflexao')
    pause
    end