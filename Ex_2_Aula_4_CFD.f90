    !****************************************************************************
    ! Doutorado em Engenharia Mecânica PPGEM
    ! Universidade Tecnológica Federal do Paraná
    ! Campus Curitiba
    ! Aluno: Marcos Takahama
    ! Professor: Paulo Henrique Dias dos Santos
    ! Disciplina: CFD
    ! Setembro de 2020
    !****************************************************************************
    ! Encontre uma nova solução para o problema anterior mas agora considere que a
    ! barra possui uma geração de energia volumétrica constante de 1000 W/m3
    !
    !Ex_2_Aula_4_CFD.f90
    !
    ! SUBROUTINES:
    ! jacobi
    !****************************************************************************
    !1 - Nome do codigo e entrada de parâmetros
    program Ex_2_Aula_4_CFD

    !Declaração de variáveis; integer=integer*4 (4bytes); real=real*4
    !parameter=valor constante
    !dimension=declaração de vetor ou matriz

    implicit none !Variáveis implicitas
    integer, parameter:: nv=5 !Número de Vomules de Controles
    integer:: i, it_final
    real:: L, dx, T_est, erro, kmax, k, residue, q
    real, dimension(nv):: Aw, Ap, Ae, Su, Sp
    real, dimension(nv):: T, x
    real, dimension(nv,nv):: E
    real:: Tp1, Tp2
    !Inputs de usuário
    write(*,*) 'Qual o valor da temperatura da parede em x=0'
    read(*,*) Tp1
    write(*,*) 'Qual o valor da temperatura da parede x=L'
    read(*,*) Tp2
    !****************************************************************************
    !2 - Método VF
    L = 1.0
    dx = L/nv
    k = 25
    q = 1000
    T_est = 25
    erro = 0.0001
    kmax = 8*10**5

    !Estimativa do campo inicial da variável
    do i=1,nv !for
        T(i) = T_est
    end do

    !Cálculo da malha:
    x(1) = 0.5D0*dx !D0 = double precision
    do i=2,nv
        x(i) = x(i-1) + dx
    end do

    write(*,*) 'Posicoes da malha'
    write(*,*) x !Escreve a posição da malha
    pause

    !Cálculo temperatura através do MVF
    do i=1,nv
        if (i==1) then
            Ae(i)=1
            Aw(i)=0
            Sp(i)=-2
            Ap(i)=Ae(i)+Aw(i)-Sp(i)
            Su(i)= 2*Tp1

            E(i,i)=Ap(i)
            E(i,i+1)=-Ae(i)
        else if((i>1).and.(i<nv)) then
            Ae(i)=1
            Aw(i)=1
            Sp(i)=0
            Ap(i)=Ae(i)+Aw(i)-Sp(i)
            Su(i)= q*(dx**2)/k

            E(i,i-1)=-Aw(i)
            E(i,i)=Ap(i)
            E(i,i+1)=-Ae(i)
        else
            Ae(i)=0
            Aw(i)=1
            Sp(i)=-2
            Ap(i)=Ae(i)+Aw(i)-Sp(i)
            Su(i)= 2*Tp2

            E(i,i-1)=-Aw(i)
            E(i,i)=Ap(i)
        end if
    end do

    write(*,*) 'Matriz A'
    write(*,*) E !Escreve a matriz de A
    write(*,*) 'Vetor B'
    write(*,*) Su !Escreve o vetor B
    pause
    !Subrotina para o cálculo dos coeficientes:
    call jacobi(Ae,Aw,Su,Ap,nv,T,erro,kmax,it_final,residue)

    write(*,*) 'erro calculado'
    write(*,*) residue !Escreve o numero de iterações pelo método de jacobi

    write(*,*) 'numero de iteracoes para convergencia'
    write(*,*) it_final !Escreve o numero de iterações pelo método de jacobi

    write(*,*) 'Temperatura calculada nos pontos nodais'
    write(*,*) T !Escreve temperaturas calculados

    pause
    !****************************************************************************
    !3 - "Plot resultados"

    !Resultados plotados no TecPlot:
    !Temperatura através de MVF:
150 format(2x,30(F20.12,1X))
    open(15,file="Tnum.plt")
    write (15,*) 'Title ="Tnum"'
    write (15,*) "Variables = 'x', 'T'"
    write (15,*) 'ZONE T = "CAMPO" i=',nv
    do i=1,nv
        write(15,150) x(i), T(i)
    end do
    close(15)

    !****************************************************************************
    !4 - Subrotinas
    Contains !As subrotinas são separadas do programa principal por esse comando
    !Subrotina do Método Jacobi
    !   call jacobi(Ae,Aw,Su,Ap,nv,T,erro,kmax)
    !****************************************************************************
    !Método de Jacobi
    SUBROUTINE jacobi(Ae,Aw,Su,Ap,N,A,erro,kmax,it_final,residue)
    implicit none
    integer:: i,k,it_final
    integer, intent(in):: N
    real, dimension(N):: Ae, Aw, Su, Ap
    real, dimension(N):: T
    real, dimension(1:N,1:2):: M
    real, dimension(N):: Ri, A
    real :: RMS, erro, kmax, residue
    !it_final: armazena a última iteração convergido ou não
    !M: matriz problema (posição x iteração)
    !Ri: resíduo do ponto i
    !RMS: resíduo médio quadrático da iteração k+1
    !1) Estimativa inicial da variável
    do i=1,N
        M(i,1) = A(i) !A(i) é oriunda do programa principal
    enddo
    RMS=0
    !2) Início do processo iterativo
    k=1 !Inicialização do contador de iterações
    do while (k<kmax)
        !Cálculo da variável com os valores da iteração anterior
        do i=1,N
            if(i==1) then !Primeiro Volume de Controle
                M(i,2)=(Ae(i)*M(i+1,1)+Su(i))/Ap(i)
            elseif((i>1).and.(i<N)) then !Volumes de Controle Internos
                M(i,2)=(Aw(i)*M(i-1,1)+Ae(i)*M(i+1,1)+Su(i))/Ap(i)
            else !Último Volume de Controle
                M(i,2)=(Aw(i)*M(i-1,1)+Su(i))/Ap(i)
            endif
        enddo

        !3) Cálculo do Resíduo
        do i=1,N
            if((i>1).and.(i<N)) then
                !Resíduo Local
                Ri(i)=abs(Su(i)+Ae(i)*M(i+1,2)+Aw(i)*M(i-1,2)-Ap(i)*M(i,2))
                !write(*,*) i,  Ri(i)
                !pause
                ! Resíduo da Iteração
                RMS=RMS+Ri(i)**2.0D0  !usa o residuo maximo calculado em cada nó
            endif
        enddo
        RMS = RMS**0.5D0

        !4) Avanço da iteração
        do i=1,N
            M(i,1)=M(i,2)
        enddo

        !5) Verificação da convergência
        if(RMS<erro+1) then !CONVERGÊNCIA: RMS tende a 1
            it_final = k !armazena a última posição iterativa
            k=kmax !força o fim da iterações
            residue = RMS-1
        else !SEM CONVERGÊNCIA:
            k=k+1 !incremento de k
            it_final = k !armazena a última iteração ->kmax
            residue = RMS-1
            !write(*,*) k, RMS
            !pause
        endif
    enddo !Finalização do contador de iterações

    !6) Resultado que retorna para o programa principal
    do i=1,N
        A(i)= M(i,2)
    enddo
    END SUBROUTINE jacobi

    end program Ex_2_Aula_4_CFD

