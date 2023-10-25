    !****************************************************************************
    ! Doutorado em Engenharia Mec�nica PPGEM
    ! Universidade Tecnol�gica Federal do Paran�
    ! Campus Curitiba
    ! Aluno: Marcos Takahama
    ! Professor: Paulo Henrique Dias dos Santos
    ! Disciplina: CFD
    ! Setembro de 2020
    !****************************************************************************
    ! Resolva o exerc�cio, agora em coordenadas cil�ndricas com um raio interno de 
    ! 0.5m e raio externo de 1m
    !
    !Ex_3_Aula_4_CFD.f90
    !
    ! SUBROUTINES:
    ! jacobi
    !****************************************************************************
    !1 - Nome do codigo e entrada de par�metros
    program Ex_3_Aula_4_CFD

    !Declara��o de vari�veis; integer=integer*4 (4bytes); real=real*4
    !parameter=valor constante
    !dimension=declara��o de vetor ou matriz

    implicit none !Vari�veis implicitas
    integer, parameter:: nv=50 !N�mero de Vomules de Controles
    integer:: i, it_final
    real:: dr, T_est, erro, kmax, k, residue, q,Tp1, Tp2
    real, parameter:: r2=1, r1=0.5
    real, dimension(nv):: Aw, Ap, Ae, Su, Sp, T, r
    real, dimension(nv,nv):: E
    !Inputs de usu�rio
    write(*,*) 'Qual o valor da temperatura da parede em x=0'
    read(*,*) Tp1
    write(*,*) 'Qual o valor da temperatura da parede x=L'
    read(*,*) Tp2
    !****************************************************************************
    !2 - M�todo VF
        write(*,*) r2, r1

    dr = (r2-r1)/nv
    k = 25
    q = 1000
    T_est = 25
    erro = 0.0001
    kmax = 8*10**5

    !Estimativa do campo inicial da vari�vel
    do i=1,nv !for
        T(i) = T_est
    end do

    !C�lculo da malha:
    r(1) = r1 + dr/2
    do i=2,nv
        r(i) = r(i-1) + dr
    end do

    write(*,*) 'Posicoes da malha'
    write(*,*) r !Escreve a posi��o da malha
    pause

    !C�lculo temperatura atrav�s do MVF
    do i=1,nv
        if (i==1) then
            Ae(i)=(r(i)+dr/2)
            Aw(i)=0
            Sp(i)=-(2*(r(i)-dr/2))
            Ap(i)=Ae(i)+Aw(i)-Sp(i)
            Su(i)= 2*(r(i)-dr/2)*tp1

            E(i,i)=Ap(i)
            E(i,i+1)=-Ae(i)
        else if((i>1).and.(i<nv)) then
            Ae(i)=r(i)+dr/2
            Aw(i)=r(i)-dr/2
            Sp(i)=0
            Ap(i)=Ae(i)+Aw(i)-Sp(i)
            Su(i)= 0

            E(i,i-1)=-Aw(i)
            E(i,i)=Ap(i)
            E(i,i+1)=-Ae(i)
        else
            Ae(i)=0
            Aw(i)=(r(i)-dr/2)
            Sp(i)=-(2*(r(i)+dr/2))
            Ap(i)=Ae(i)+Aw(i)-Sp(i)
            Su(i)= 2*(r(i)+dr/2)*tp2

            E(i,i-1)=-Aw(i)
            E(i,i)=Ap(i)
        end if
    end do

    write(*,*) 'Matriz A'
    write(*,*) E !Escreve a matriz de A
    write(*,*) 'Vetor B'
    write(*,*) Su !Escreve o vetor B
    pause
    !Subrotina para o c�lculo dos coeficientes:
    call jacobi(Ae,Aw,Su,Ap,nv,T,erro,kmax,it_final,residue)

    write(*,*) 'erro calculado'
    write(*,*) residue !Escreve o numero de itera��es pelo m�todo de jacobi

    write(*,*) 'numero de iteracoes para convergencia'
    write(*,*) it_final !Escreve o numero de itera��es pelo m�todo de jacobi

    write(*,*) 'Temperatura calculada nos pontos nodais'
    write(*,*) T !Escreve temperaturas calculados

    pause
    !****************************************************************************
    !3 - "Plot resultados"

    !Resultados plotados no TecPlot:
    !Temperatura atrav�s de MVF:
150 format(2x,30(F20.12,1X))
    open(15,file="Tnum.plt")
    write (15,*) 'Title ="Tnum"'
    write (15,*) "Variables = 'x', 'T'"
    write (15,*) 'ZONE T = "CAMPO" i=',nv
    do i=1,nv
        write(15,150) r(i), T(i)
    end do
    close(15)

    !****************************************************************************
    !4 - Subrotinas
    Contains !As subrotinas s�o separadas do programa principal por esse comando
    !Subrotina do M�todo Jacobi
    !   call jacobi(Ae,Aw,Su,Ap,nv,T,erro,kmax)
    !****************************************************************************
    !M�todo de Jacobi
    SUBROUTINE jacobi(Ae,Aw,Su,Ap,N,A,erro,kmax,it_final,residue)
    implicit none
    integer:: i,k,it_final
    integer, intent(in):: N
    real, dimension(N):: Ae, Aw, Su, Ap
    real, dimension(N):: T
    real, dimension(1:N,1:2):: M
    real, dimension(N):: Ri, A
    real :: RMS, erro, kmax, residue
    !it_final: armazena a �ltima itera��o convergido ou n�o
    !M: matriz problema (posi��o x itera��o)
    !Ri: res�duo do ponto i
    !RMS: res�duo m�dio quadr�tico da itera��o k+1
    !1) Estimativa inicial da vari�vel
    do i=1,N
        M(i,1) = A(i) !A(i) � oriunda do programa principal
    enddo
    RMS=0
    !2) In�cio do processo iterativo
    k=1 !Inicializa��o do contador de itera��es
    do while (k<kmax)
        !C�lculo da vari�vel com os valores da itera��o anterior
        do i=1,N
            if(i==1) then !Primeiro Volume de Controle
                M(i,2)=(Ae(i)*M(i+1,1)+Su(i))/Ap(i)
            elseif((i>1).and.(i<N)) then !Volumes de Controle Internos
                M(i,2)=(Aw(i)*M(i-1,1)+Ae(i)*M(i+1,1)+Su(i))/Ap(i)
            else !�ltimo Volume de Controle
                M(i,2)=(Aw(i)*M(i-1,1)+Su(i))/Ap(i)
            endif
        enddo

        !3) C�lculo do Res�duo
        do i=1,N
            if((i>1).and.(i<N)) then
                !Res�duo Local
                Ri(i)=abs(Su(i)+Ae(i)*M(i+1,2)+Aw(i)*M(i-1,2)-Ap(i)*M(i,2))
                !write(*,*) i,  Ri(i)
                !pause
                ! Res�duo da Itera��o
                RMS=RMS+Ri(i)**2.0D0  !usa o residuo maximo calculado em cada n�
            endif
        enddo
        RMS = RMS**0.5D0

        !4) Avan�o da itera��o
        do i=1,N
            M(i,1)=M(i,2)
        enddo

        !5) Verifica��o da converg�ncia
        if(RMS<erro+1) then !CONVERG�NCIA: RMS tende a 1
            it_final = k !armazena a �ltima posi��o iterativa
            k=kmax !for�a o fim da itera��es
            residue = RMS-1
        else !SEM CONVERG�NCIA:
            k=k+1 !incremento de k
            it_final = k !armazena a �ltima itera��o ->kmax
            residue = RMS-1
            !write(*,*) k, RMS
            !pause
        endif
    enddo !Finaliza��o do contador de itera��es

    !6) Resultado que retorna para o programa principal
    do i=1,N
        A(i)= M(i,2)
    enddo
    END SUBROUTINE jacobi

    end program Ex_3_Aula_4_CFD

