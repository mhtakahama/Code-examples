clear all; close all; clc

M=load('treinamento.txt'); vet_entrada= M(:,1:13)'; vet_desejado=M(:,14:16)'; [i1,j1]=size(M);

x=minmax(vet_entrada);net=newff(x,[5 3],{'logsig' 'logsig'},'trainlm');
net.trainParam.epochs=500; net.trainParam.goal=1e-6; net.trainParam.lr=0.5; net.trainParam.show=5;

net=train(net,vet_entrada,vet_desejado); %3o Passo, Treinar a rede

%Testar Rede
T=load('teste.txt'); vet_teste_entrada= T(:,1:13)'; vet_teste_desejado=T(:,14:16)';
[i2,j2]=size(T);

vet_saida=sim(net, vet_teste_entrada); %4o Passo Testar a Rede
[i3,j3]=size(vet_saida);

for k=1:j3 %substituir numeros reais por inteiros
    for i=1:i3
        if vet_saida(i,k)>0.5
            vet_saida(i,k)=1;
            if i==1;
                vet_resp(k)={'A'};
            end
            if i==2;
                vet_resp(k)={'B'};
            end
            if i==3;
                vet_resp(k)={'C'};
            end
        else
            vet_saida(i,k)=0;
        end
        
        if vet_saida(i,k)-vet_teste_desejado(i,k)~=0 %comparação entre vetor de saída e vetor desejado
            fprintf('Houve erro na amostra %d, com o parâmetro %d \n',k,i);
        end
    end
end

if vet_saida(i,k)==vet_teste_desejado(i,k)
            fprintf('O resultado obtido é compatível com o desejado\n');
            
        else
           fprintf('O resultado obtido não é compatível com o desejado\n');
end

vet_resp'
