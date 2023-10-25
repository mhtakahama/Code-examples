clear all; echo off; close all force; clc; format long; %clear another variables
%% Doutorado em Engenharia Mec�nica PPGEM
%Universidade Tecnol�gica Federal do Paran�
%Campus Curitiba
%Aluno: Marcos Takahama
%Professor Paulo Henrique Dias dos Santos
%Disciplina:CFD
%03/10/2020
% Ex_4_Aula_4_CFD.f90
%% 1 - Inputs

%input parameters
prompt = {'Valor da Temperatura 1','Valor da Temperatura 2','n�mero de VC'};
dlg_title = 'Condi��es de Contorno';
num_lines = 1;
defaultans = {'150','50','50'};
N=50; %this will control the width of the inputdlg
answer = inputdlg(prompt,dlg_title,[1, length(dlg_title)+N],defaultans);
%input

tp1=str2num(answer{1});
tp2=str2num(answer{2});
nv=str2num(answer{3});

%% 2-M�todo VF
r1=0.5;
r2=1;
dr = (r2-r1)/nv;
k = 25;
q = 1000;
T_est = 25;
erro = 1e-4;
kmax = 8*10^5;

%Estimativa do campo inicial da vari�vel
for i=1:nv
    T(i) = T_est;
end

%C�lculo da malha:
r(1) = r1+dr/2;

for i=2:nv
    r(i) = r(i-1) + dr;
end

fprintf('Informa��es da Malha')
r

%C�lculo temperatura atrav�s do MVF
for i=1:nv
    if (i==1)
        Ae(i)=(r(i)+dr/2);
        Aw(i)=0;
        Sp(i)=-(2*(r(i)-dr/2));
        Ap(i)=Ae(i)+Aw(i)-Sp(i);
        Su(i)= 2*(r(i)-dr/2)*tp1+q/k*(dr/2)^2;
        
        E(i,i)=Ap(i);
        E(i,i+1)=-Ae(i);
    else if((i>1)&&(i<nv))
        Ae(i)=r(i)+dr/2;
        Aw(i)=r(i)-dr/2;
        Sp(i)=0;
        Ap(i)=Ae(i)+Aw(i)-Sp(i);
        Su(i)= q/k*(dr/2)^2;
            
        E(i,i-1)=-Aw(i);
        E(i,i)=Ap(i);
        E(i,i+1)=-Ae(i);
        else
        Ae(i)=0;
        Aw(i)=(r(i)-dr/2);
        Sp(i)=-(2*(r(i)+dr/2));
        Ap(i)=Ae(i)+Aw(i)-Sp(i);
        Su(i)= 2*(r(i)+dr/2)*tp2+q/k*(dr/2)^2;

        E(i,i-1)=-Aw(i);
        E(i,i)=Ap(i);
        end
    end
end
fprintf('Informa��es da Matriz e vetor')
E
Su

% Subrotina para o c�lculo dos coeficientes:
[residue,it_final,T]=jacobi(Ae,Aw,Su,Ap,nv,T,erro,kmax)

%% plot de resultados

pontosr=linspace(r1,r2,10000);
pontosy=50+100/log(0.5)*log(pontosr);
gcf=figure;
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);%Maximize window

plot(pontosr,pontosy,'b')
hold on
plot(r,T,'r-o')

title(['Exerc�cio 4'],'FontSize',20) %Legend options
ylabel(['Temperatura (�C)'])
xlabel('Comprimento (m)')
legend('Solu��o Anal�tica','Solu��o Calculada')
%% Functions

function [residue,it_final,T]=jacobi(Ae,Aw,Su,Ap,N,T,erro,kmax)
RMS=0;
%1) Estimativa inicial da vari�vel
for i=1:N
    M(i,1) = T(i); %A(i) � oriunda do programa principal (T)
end

%2) In�cio do processo iterativo
k=1; %Inicializa��o do contador de itera��es
while (k<kmax)
    %C�lculo da vari�vel com os valores da itera��o anterior
    for i=1:N
        if(i==1) %Primeiro Volume de Controle
            M(i,2)=(Ae(i)*M(i+1,1)+Su(i))/Ap(i);
        elseif((i>1)&&(i<N)) %Volumes de Controle Internos
            M(i,2)=(Aw(i)*M(i-1,1)+Ae(i)*M(i+1,1)+Su(i))/Ap(i);
        else %�ltimo Volume de Controle
            M(i,2)=(Aw(i)*M(i-1,1)+Su(i))/Ap(i);
        end
    end
    
    %3) C�lculo do Res�duo
    for i=1:N
        if((i>1)&&(i<N))
            %Res�duo Local
            Ri(i)=abs(Su(i)+Ae(i)*M(i+1,2)+Aw(i)*M(i-1,2)-Ap(i)*M(i,2) );
            % Res�duo da Itera��o
            RMS = RMS + Ri(i)^2.0;
        end
    end
    RMS = RMS^0.5;
    
    %4) Avan�o da itera��o
    for i=1:N
        M(i,1)=M(i,2);
    end
    
    %5) Verifica��o da converg�ncia
    if(RMS<=erro+1) %CONVERG�NCIA:
        it_final = k; %armazena a �ltima posi��o iterativa
        k=kmax; %for�a o fim da itera��es
        residue = RMS-1;
        
    else %SEM CONVERG�NCIA:
        k=k+1; %incremento de k
        it_final = k; %armazena a �ltima itera��o ->kmax
        residue = RMS-1;
    end
end %Finaliza��o do contador de itera��es

%6) Resultado que retorna para o programa principal
for i=1:N
    T(i)= M(i,2);
end
end