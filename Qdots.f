        program Atomico U finito
        
C      Emplea las Funciones de Green para el modelo atomico de la impureza de Anderson, considerando
C      como canal de conducci�n, un nanotubo de carbono, tipo zig-zag met�lico (n=3) y una banda
C      correspondiente a un canal bal�stico.
C      Calcula propiedades termoel�ctricas de transporte (G,S,K,ZT y Wideman-Franz), adem�s
C      de la densidad de estados y los n�meros de ocupaci�n.
C      Se impone la regla de suma de Fridel en el proceso de dispersi�n, para determinar
C      energ�a del cumulante at�mico, asociado a la banda de conducci�n.

      implicit complex*16 (x-z)
      implicit real*8 (a-h,o-w)
        
      INTEGER DT1,DT2,DT3,DT38,DT11,OUT
      PARAMETER (DT1=1,DT2=2,DT38=38,DT11=11,OUT=3)
      DIMENSION U(17),RF(17),EN(16),RC(13),HVEC(1001),XKJ(501,501),
     *RFD(16),RFS(16),UD(16),US(16)
      CHARACTER(100) :: Tconsole, Archivo1, Archivo2
        EXTERNAL FMU
        EXTERNAL FMU2
c        EXTERNAL PIMGCOEF1
        EXTERNAL Fmin
        
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
        COMMON/RETAN/E00,EFF,E1,E2,DEL,DEL0,ZUI,EPSIL
        COMMON/ET/ETA,OI,OF,NT
        COMMON/DELTA/DELTA
        COMMON/RESI/U,EN,RF,RC
        COMMON/OCC/TXFA,TDS,TXFC,TXCA,TOT,TSOMA,FRIED,R0F,R0C
        COMMON/BRETAN/GNT,GNF,GNQ,GDS,GNR,COMPLETEZA
        COMMON/COMPLETEZA/ANXF,ANDS,ANSOMA
        COMMON/OCC2/TXF,TXC,TXUP,TXDOWN,TXDD,TXVAC,TSOMAF,TSOMA2
        COMMON/OCC3/TX11,TX22,TX33,TX44,TX13,TX31,FSUM,TXCAD
        COMMON/OCUP_EXAT/ANDD,ANSS,ANDDS,AVNSS,SOMA,ANDM13,ANDM24
        COMMON/OCUP_ATOMICO/ANVAC,ANSSS,ANDDD
        COMMON/RF22/DELTA22
        COMMON/QD1/EF1,VQ1,EL1,TXF1,U1
        COMMON/QD2/EF2,VQ2,EL2,TXF2,U2
       
      IF(COMMAND_ARGUMENT_COUNT().NE.1)THEN
        WRITE(*,*)'ERROR, ONE COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
        STOP
      ENDIF
      
      CALL GET_COMMAND_ARGUMENT(1,Tconsole)
      
      OPEN (UNIT=DT1,FILE="DATA-"//trim(Tconsole)//".dat",
     *     STATUS='UNKNOWN')
      OPEN (UNIT=DT2,FILE="DATA+"//trim(Tconsole)//".dat",
     *     STATUS='UNKNOWN')
c     Chequea si se recibio un argumento por consola
     
     
      
        pi=3.14159265358979323D0        
        CC=1.0d0/pi
C       ESTOS SON LOS DATOS PARA UN CANAL BALISTICO

        DELTA=1.0d-2
        D=1.0D0
        V2=2.0D0*D*DELTA*CC
        V=DSQRT(V2)
        VQ2=V
        EF2=-10.0d0*DELTA
        EF=EF2

        ETA=1.0d-4
        
        COU=20.0D0*DELTA
        U2=COU

        AMU=0.D0*DELTA

        EQ=AMU
        EL2=EQ
        EL1=AMU
C     EL1=3.0D-1*DELTA
       Read(Tconsole,*)T  
        T=T*DELTA
        BETA=1.0D0/T

        EPS=1.0D-8

C       HIBRIZACIÓN PARA SEGUNDO QD

        VQ1=VQ2

C REPULSION COULOMBIANA PARA SEGUNDO QD

        U1=1.0d0*U2
        EF2=(-U2/2.0D0)+AMU

        AEKI2=-1.0D-4+AMU
        AEKF2=1.0D-4+AMU
        
        AEKI1=-4.5D-3+AMU
        AEKF1=4.5D-3+AMU

        TXF1=1.0d0
        TXF2=1.0d0

        
C     ESTOS  SON PARA EL NANOTUB0

C       ZWE=DCMPLX(AMU,1.0D-3)
C        Dx = 0.d0
C       Dy = 0.d0

C      CALL zig_mauro(ZWE,Dx,Dy,ZG0,ZG0N)
C      RHOC=DIMAG(ZGO)*CC

C      DELTA=0.01D

C       V2=DELTA*CC/RHOC
C       V=DSQRT(V2)CLAVE	

        
c      LA TEMPERATURA

        LT=1

        TI=4.0d0*DELTA
        TF=7.0d0*DELTA
C        WRITE(6,33)
C33      FORMAT(' LT=1',' TI=DELTA',' TF=1.')
C        READ(*,*)LT,TI,TF
        AR=TF-TI
        IF(LT.EQ.1)GO TO 59
        H1=AR/(LT-1)
        DL=(DLOG10(TF)-DLOG10(TI))/(LT-1)
59      CONTINUE
        DO 869 IIQ=1,LT
        DLT=DLOG10(TI)+(IIQ-1)*DL
        T=10**DLT
        BETA=1.D0/T

C CALCULOS segundo QD

c        EF=EF2
c        CALL EXATA_UFINITO
c        CALL OCUP3
      
        IF (T.LT.(1.0D0*DELTA)) THEN
       DEK=Fmin(AEKI2,AEKF2,FMU2,EPS)
      EQ=DEK
C        EQ=AMU+1.0D-12
        ELSE
        EQ=AMU+1.0D-12
        ENDIF
        CONTINUE
        EL2=EQ
        EF=EF2
        V=VQ2
        V2=VQ2*VQ2
        COU=U2

        CALL EXATA_UFINITO
        CALL OCUP3

        TXF2=2.0d0*TXF

        IF (T.GE.(1.0D0*DELTA)) THEN
        DIF=0.0D0
        DIF2=0.0D0
        ELSE
        DIF=DABS(FRIED-R0F)
        DIF2=(DIF/R0F)*1.0D2
        ENDIF
        CONTINUE

        if ((DIF2.LT.(1.0d0)).OR.(DIF.LT.(1.0D-8))) then
         AEKI2=EL2-5.0D-5
         AEKF2=EL2+5.0D-5
         AQD2=1.0D0
         EE2=EL2
         ELSE
         AEKI2=-5.0D-3+EE2
         AEKF2=5.0D-3+EE2
         AQD2=0.0D0
         ENDIF
         CONTINUE

ccc        DELTA22=PIMGCOEF1(AMU)*VQ2*VQ2
ccc        DELTA22=DABS(DELTA22)
c        VQ2=1.0D0/DABS(PIMGCOEF1(AMU))
c        VQ2=DSQRT(VQ2)

c      Energ�a del nivel localizado

C        WRITE(6,4001)
c4001    FORMAT('AEFI=-7.0DELTA',5X,'AEFF=DELTA',5X,' NEF=1')
C        READ(*,*)AEFI,AEFF,NEF

       AEFI=(-U1/2.0D0)+AMU
c       AEFI=-22.0D0*DELTA
        AEFF=AMU-45.0D0*DELTA
C       AEFF=-22*DELTA
        NEF=1
        EQ=AMU
C       EQ=2.0D-1*DELTA
        EL1=EQ

        IF(NEF.EQ.1)GO TO 4000
        DEF=(AEFF-AEFI)/(NEF-1)
4000    CONTINUE

        DO 5002 IB=1,NEF
        AEF=AEFI+(IB-1)*DEF
c        AEF=-40.0D0*DELTA
        EF=AEF
        EF1=EF
           
c      CALCULOS PARA primer QD
        
c        CALL EXATA_UFINITO
c        CALL OCUP3 
 

        IF (T.LT.(1.0D0*DELTA)) THEN
        DEK=Fmin(AEKI1,AEKF1,FMU,EPS)
        EQ=DEK
        ELSE
        EQ=AMU+1.0D-12
        ENDIF
        CONTINUE
        EL1=EQ 
        V=VQ1
        V2=V*V 
        COU=U1
       
        CALL EXATA_UFINITO
        CALL OCUP3
        TXF1=2.0d0*TXF

        IF (T.GE.(1.0D0*DELTA)) THEN
        DIFB=0.0D0
        DIF2B=0.0D0
        ELSE
        DIFB=DABS(FRIED-R0F)
        DIF2BB=(DIFB/R0F)*1.0D2
        DIF2B=DIF2BB+DABS((EL2-EL1)/EL2)*1.0D2
        DIF2B=DIF2B/2.0D0
        ENDIF
        CONTINUE
        
        CALL FASE(FASELO,FASECON,FASEQ1,FASEQ2)
         FASELO=FASELO+5.0D-1
         FASECON=FASECON+5.0D-1
         FASEQ1=FASEQ1+5.0D-1
         FASEQ2=FASEQ2+5.0D-1
c        FRIED2=FRIED*PI*1.0d-2
c        FRIED=FRIED2/(PI*DELTA22)
        
         if ((DIF2BB.LT.(1.0d0)).OR.(DIFB.LT.(1.0D-8))) then
         AEKI1=EL1-5.0D-3
         AEKF1=EL1+5.0D-3
         AQD1=1.0D0
         EE=EL1
         ELSE
         AEKI1=EE-5.0D-2
         AEKF1=EE+5.0D-2
         AQD1=0.0D0
         ENDIF
         CONTINUE
              

870     FORMAT(17F16.8)       
        CALL THERMO(G20,S0,AK0,TEM0,WFR0,G2,S,AK,TEM,WFR,E)
          
C       DIF3=2.0D0*(DABS((FASEQ1-TXF1*0.5D0)/TXF1))*1.0D2
C       DIF3B=2.0D0*(DABS((FASEQ2-TXF2*0.5D0)/TXF2))*1.0D2

       G2=G2/2.0d0

       WRITE(6,870)EF1/DELTA,EL1/DELTA,EL2/DELTA,DIF2,DIF2BB,
     *       G2,S,AK,TEM,WFR,TXF1,TXF2,FASELO,FASEQ2,FASECON,T/DELTA,E


C      IF ((AQD2.EQ.(1.0d0)).AND.(AQD1.EQ.(1.0D0))) THEN
      IF ((DIF2.LT.(1.0d0)).AND.(DIF2BB.LT.(5.0D0))) THEN

180   FORMAT(19F16.8)
c      WRITE(6,171)T,EF,FRIED,R0F,DIF
c171   FORMAT('T=',F10.8,1X,'EF=',F10.8,1X,'FRIED=',F10.8,1X,
c     *       'ROF=',F10.8,1X,'DIF=',F10.8) 
      WRITE(DT1,180)EF1/DELTA,EL1/DELTA,EL2/DELTA,DIF2,DIF2BB,G2,S,AK,
     *TEM,WFR,TXF1,TXF2,FASELO,FASEQ2,FASEQ1,FASECON,E,T/DELTA,EF2/DELTA
    
       ELSE
       CONTINUE
       ENDIF


      EFF=AEFI-(IB-1)*DEF
C      EFF=-2.0D0*DELTA-(IB-1)*DEF
c       EFF=2.0D0*DELTA
       EF=EFF
       EF1=EFF
       
C       IF (T.LT.(1.0d0*DELTA)) THEN
C        DEK=Fmin(AEKI1,AEKF1,FMU,EPS)
C        EQ=DEK
C        ELSE
C        EQ=AMU
C        ENDIF
C        CONTINUE
        EL2=-EL2
        EL1=-EL1
        EQ=EL1
        V=VQ1
        V2=V*V
        COU=U1

        CALL EXATA_UFINITO
        CALL OCUP3
        TXF1=2.0d0*TXFF

        IF (T.GE.(1.0D0*DELTA)) THEN
        DIFB=0.0D0
        DIF2B=0.0D0
        ELSE
        DIFB=DABS(FRIED-R0F)
        DIF2BB=(DIFB/R0F)*1.0D2
        DIF2B=DIF2BB+DABS((EL2-EL1)/EL2)*1.0D2
        DIF2B=DIF2B/2.0D0
        ENDIF
        CONTINUE

        CALL FASE(FASELO,FASECON,FASEQ1,FASEQ2)
         FASELO=FASELO+5.0D-1
         FASECON=FASECON+5.0D-1
         FASEQ1=FASEQ1+5.0D-1
         FASEQ2=FASEQ2+5.0D-1


         if ((DIF2BB.LT.(1.0d0)).OR.(DIFB.LT.(1.0D-8))) then
         AEKI1=EL1-5.0D-3
         AEKF1=EL1+5.0D-3
         AQD1=1.0D0
         EE=EL1
         ELSE
         AEKI1=EE-5.0D-2
         AEKF1=EE+5.0D-2
         AQD1=0.0D0
         ENDIF
         CONTINUE

C870     FORMAT(17F16.8)
        CALL THERMO(G20,S0,AK0,TEM0,WFR0,G2,S,AK,TEM,WFR,E)

C       DIF3=2.0D0*(DABS((FASEQ1-TXF1*0.5D0)/TXF1))*1.0D2
C       DIF3B=2.0D0*(DABS((FASEQ2-TXF2*0.5D0)/TXF2))*1.0D2

       G2=G2/2.0d0

       WRITE(6,870)EF1/DELTA,EL1/DELTA,EL2/DELTA,DIF2,DIF2BB,
     *       G2,S,AK,TEM,WFR,TXF1,TXF2,FASELO,FASEQ2,FASECON,T/DELTA,E


C       IF ((AQD2.EQ.(1.0d0)).AND.(AQD1.EQ.(1.0D0))) THEN
      IF ((DIF2.LT.(1.0d0)).AND.(DIF2BB.LT.(5.0D0))) THEN

C180   FORMAT(19F16.8)
 1      WRITE(DT2,180)EF1/DELTA,EL1/DELTA,EL2/DELTA,DIF2,DIF2BB,G2,S,AK,
     *TEM,WFR,TXF1,TXF2,FASELO,FASEQ2,FASEQ1,FASECON,E,T/DELTA,EF2/DELTA

       ELSE
       CONTINUE
       ENDIF
       
      EQ=-EQ
      EL1=-EL1
      EL2=-EL2

c4002    CONTINUE
5002    CONTINUE
     
869   CONTINUE

      A=1.0d-2
      B=4.0d0
      NTT=1000
       
C Esto es para ekl ajust


c      CALL GRAF(A,B,NTT)
c      CALL Ajus(A,B,NTT)
      
      CLOSE(DT1)
      
      CLOSE(DT2)
      
      end
         
C      Fin del programa principal
C*********************************************************************
C*********************************************************************

        Function ALateral2(EW)
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
        External Ferm

        FER=FERM(EW)
	DFER=FER*(1.D0-FER)*BETA

        ALateral2=2.0D0*DFER*(EW**(2.0D0))

        RETURN
        END

        Function ALateral0(EW)
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
        External Ferm

        FER=FERM(EW)
	DFER=FER*(1.D0-FER)*BETA

        ALateral0=2.0D0*DFER

        RETURN
        END
      

      SUBROUTINE Ajus(A,B,NTT)
*******************************************************************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      
****************************************************************************
      INTEGER DT31
      PARAMETER (DT31=31)

      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
      COMMON/OCC2/TXF,TXC,TXUP,TXDOWN,TXDD,TXVAC,TSOMAF,TSOMA2

      EXTERNAL ALateral0, ALateral2

      OPEN (UNIT=DT31,FILE='2QDS-Universal.dat',STATUS='UNKNOWN')


      PI=1.0D0/CC
      A0=-0.00935763d0
      A2=-9.90825d-5
      A1=62.1537d0
      A3=0.853788d0
      A4=-62.1795d0
      A5=0.0400315d0
      A6=0.0247632d0
      A7=0.1081824d0
      A8=1.45984d0
      A9=-0.958835d0

      A02=6.11746D28
      A12=1.85047D29
      A22=-2.01215D0
      A32=-0.352106D0
      A42=1.72775D28
      A52=-0.243646D0


C para la Conductancia electrica entre 1.0d-5*Tk y 1000Tk
c      T=1.0d0/BETA
c      TR=T/Tk   Para Pi/3, Tk=850 mK

       TK=0.87037d0

c       DELTA=PI/2.0d0
c       DELTA=1.0002d0*PI/2.0D0
c        DELTA=0.848209d0*PI/2.0d0
        DELTA=0.d0
        AR=B-A

        EPSIL=1.0D-8
        AMU=0.0D0
        GL0=2.0d0

        IF(NTT.EQ.1)GO TO 59
        H1=AR/(NTT-1)
        DL=(DLOG10(B)-DLOG10(A))/(NTT-1)
59      CONTINUE
        DO 80 IIQ=1,NTT
        DLT=DLOG10(A)+(IIQ-1)*DL

       TR=10**DLT

c        TR=(8.0d-1/(3.3d0*1.16D0))
c        TR=1.05d0
        BETA=1.0d0/TR

       AI=-80.0D0*TR
       AF=0.0D0*TR

C      TR=A+(K-1)*DM1

c     AQUI HAY QUE DESCRIBIR EL VALOR DE TK!!!

c iNTEGRAL PARA EL CASO LATERAL
       
C      GL0=QSUBA(AI,AF,EPSIL,NPTS,ICHECK,RELERR,ALateral0)
C      GL0=2.0d0*GL0
     
      AL2L0=QSUBA(AI,AF,EPSIL,NPTS,ICHECK,RELERR,ALateral2)
      AL2L0=2.0d0*AL2L0*BETA*BETA

      Gs=(A0/(A1*DEXP(A2*(TR**(A3)))+A4))+A5
     *   +(A6/(DEXP(A7*TR**(A8))+A9))

      Gs=4.0d0*Gs/2.00856d0
      Gs=(2.0d0/1.99994891d0)*Gs
   
      AL2S=A52+(A02/(A12*DEXP(A22*((TR)**(A32)))+A42))
      AL2S=2.0D0*AL2S
  
      A01L=TR*TR*((10.64D0**(1.0D0/0.728D0))-1.0D0)+1.0D0
      A01L=(1.0D0/A01L)**(0.728D0)

      AL10T=(2.0D0*PI*PI*PI*TR*0.923d0/3.0D0)*A01L


      AA=4.0D0*(Gs-1.0D0)/(AL10T*AL10T)
      AA=AA*(AL2S-PI*PI/6.0D0)
      BB=((Gs-2.0d0)*PI*PI/6.0d0)+AL2S
      BB=BB*(-4.0d0/(AL10T*AL10T))
      CCC=2.0d0*PI*PI/(3.0d0*AL10T*AL10T)
      DD=BB*BB-4.0d0*(AA+1.0d0)*(CCC-1.0d0)
      DD=DSQRT(DD)

      Raiz1=(-BB-DD)/(2.0d0*(AA+1.0d0))
      Raiz2=(-BB+DD)/(2.0d0*(AA+1.0d0))

c      DO 80 K=1,NTT
c      Delta=A+(K-1)*H1
           
c      Delta=0.59364d0*Deltaa+0.2895682d0
c      Delta=Delta*Pi/2.0d0
   
      GD=-Gs*DCOS(2.0d0*Delta)+(DCOS(2.0d0*Delta)+1.0d0)
    
C    Gs para el caso lateral

      GDL=(GL0-GD)
      GDL=GDL*2.0d0/1.91991994d0
      
C    L2/(T*T) 
      AL2=-(AL2S-(PI*PI/6.0D0))*DCOS(2.0D0*Delta)+(PI*PI/6.0D0)

c   AL2SIMETRICO -CASO LATERAL

      AL2L=AL2L0-AL2
      
C    L1/T
      AL1T=AL10T*DSIN(Delta)*DCOS(Delta)

c   AL1sIMETRICO CASO LATERAL

      AL1TL=-AL1T

      SD=-AL1T/GD
      SDL=-AL1TL/GDL

      AKD=AL2*TR-(AL1T*AL1T*TR/GD)
      AKDL=AL2L*TR-(AL1TL*AL1TL*TR/GDL)
     

      EPSI=(AL1T*AL1T)/(GD*AL2)
      EPSIL=(AL1TL*AL1TL)/(GDL*AL2L)
      AZT=EPSI/(1.0D0-EPSI)
      AZTL=EPSIL/(1.0D0-EPSIL)
      TRR=TR*TK
      
*************************************
c      WRITE(DT31,5900)TR,Raiz1,Raiz2
      WRITE(DT31,5900)TR,TRR,GD,GDL,SD,SDL
     * ,(3.0d0*AKD/(TR*Pi*Pi))
     * ,(3.0d0*AKDL/(Pi*Pi*TR))
     * ,AZT,AZTL,EPSI,EPSIL,GD*SD/TRR
80    CONTINUE
5900  FORMAT(13F15.8)
      CLOSE(DT31)
**********************************************************************************

      RETURN
      END

        FUNCTION ZGCOEF1(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/QD1/EF1,VQ1,EL1,TXF1,U1
        COMMON/QD2/EF2,VQ2,EL2,TXF2,U2
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU

        EF=EF1
        V=VQ1
        V2=VQ1*VQ1
        V21=V2       
        COU=U1
	ETTA=1.D-4
        EQ=EL1

        ZW=DCMPLX(EW,ETTA)

c   ESTO ES PARA los 2 qds inmersos.
        CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
     *ZGFXATSS,ZG00,ZM)
************************************************
        ZGQD1=ZGF
        ZGC1=ZM
        ZGCOEF1=ZGC1*ZGC1*ZGC1*V21*ZGQD1

C Esto es para el primer QD lateral

CCC      ZGCOEF1=ZGC1*ZGC1+ZGC1*ZGC1*ZGC1*V21*ZGQD1

C esto es para el segundo QD lateral

C        EF=EF2
C        V=VQ2
C        V2=VQ2*VQ2
C        V22=V2
C        COU=U2
C	ETTA=1.D-4
C        EQ=EL2

C        ZW=DCMPLX(EW,ETTA)

C       CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
C     *ZGFXATSS,ZG00,ZM)
C       ZGQD2=ZGF

c     ZGC00=ZGC1*ZGC1V21*ZGQD1
c     ZGCOEF1=(ZGC00/(V22*ZGQD2))+ZGC00*ZGC00

C esto es para los dos QDs acoplados uno arriba del otro

c       ZGCOEF1=ZGC1*ZGC1*v21*ZGQD1*ZGQD1+
c     *   (ZGC1*ZGC1*V22*ZGQD2/(V22*ZGQD2))

C esto es para dos QDs laterales
c         ZGCOEF1=(ZGC1*ZGC1*V21*ZGQD1+ZGC1)
c         ZGCOEF1=(ZGCOEF1/(V22*ZGQD2))+ZGCOEF1*ZGCOEF1
        RETURN
        END
  

        FUNCTION PIMGCOEF1(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/QD1/EF1,VQ1,EL1,TXF1,U1
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU

        EXTERNAL ZGCOEF1

c        EF=EF1
c        V=VQ1
c        V2=VQ1*VQ1
c        COU=U1
c        EQ=EL1

       FAS=DIMAG(ZGCOEF1(EW))/DREAL(ZGCOEF1(EW))
       FAS=DATAN(FAS)
        FAS=FAS/2.0D0

CC Esto es una prueba

      pi=1.0d0/CC
c      FAS=FAS+0.5D0*pi

      PIMG=DREAL(ZGCOEF1(EW))*DREAL(ZGCOEF1(EW))+
     *    DIMAG(ZGCOEF1(EW))*DIMAG(ZGCOEF1(EW))
        PIMG=DSQRT(PIMG)
        PIMGCOEF1=(DSQRT(PIMG))*DSIN(FAS)

c        B2=DREAL(ZGCOEF1(EW))*DREAL(ZGCOEF1(EW))+	
c     *    DIMAG(ZGCOEF1(EW))*DIMAG(ZGCOEF1(EW))
c        B2=DSQRT(B2)+DREAL(ZGCOEF1(EW))
c        B2=DSQRT(B2)
c        PIM=2.0d0*DIMAG(ZGCOEF1(EW))/B2


c        A=-DREAL(ZGCOEF1(EW))*0.5D0
c        B=DREAL(ZGCOEF1(EW))*DREAL(ZGCOEF1(EW))+
c     *    DIMAG(ZGCOEF1(EW))*DIMAG(ZGCOEF1(EW))
c        B=DSQRT(B)*0.5D0
c        PIMGCOEF1=DSQRT(A+B)
c        PIMGCOEF1=PIMGCOEF1+PIM
c        PIMGCOEF1=PIMGCOEF1/2.0D0
        
        RETURN
        END

         FUNCTION ZGCON(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/QD1/EF1,VQ1,EL1,TXF1,U1

        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU

        EXTERNAL ZGCOEF1
C        EXTERNAL PIMGCOEF1

c        EF=EF1
c        V=VQ1
c        V2=VQ1*VQ1
c        COU=U1
c	ETTA=1.D-4
c        EQ=EL1

       FAS=DIMAG(ZGCOEF1(EW))/DREAL(ZGCOEF1(EW))
       FAS=DATAN(FAS)
        FAS=FAS/2.0D0
cc PRUEBA
        pi=1.0d0/CC
c        FAS=FAS+0.5D0*pi

      PIM2=DREAL(ZGCOEF1(EW))*DREAL(ZGCOEF1(EW))+
     *    DIMAG(ZGCOEF1(EW))*DIMAG(ZGCOEF1(EW))
        PIM2=DSQRT(PIM2)
        PIM2=DSQRT(PIM2)
        PIM=PIM2*DSIN(FAS)

        PRE=PIM2*DCOS(FAS)

C        PRE=DREAL(ZGCOEF1(EW))*DREAL(ZGCOEF1(EW))+
C     *    DIMAG(ZGCOEF1(EW))*DIMAG(ZGCOEF1(EW))
C        PRE=(DSQRT(PRE)+DREAL(ZGCOEF1(EW)))*0.5D0
C        PRE=DSQRT(PRE)
        
C        A=-DREAL(ZGCOEF1(EW))*0.5D0
C        B=DREAL(ZGCOEF1(EW))*DREAL(ZGCOEF1(EW))+
C     *    DIMAG(ZGCOEF1(EW))*DIMAG(ZGCOEF1(EW))
C        B=DSQRT(B)*0.5D0
C        B=A+B
         
C        PREAL=DREAL(ZGCOEF1(EW))+B
C        PREAL=DSQRT(PREAL)
C        PREAL=(PREAL+PRE)/2.0D0
C        PIMAG=PIMGCOEF1(EW)
        ZGCON=DCMPLX(PRE,PIM)
        
        RETURN
        END

        FUNCTION FASELOCAL(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU

        EXTERNAL ZGLOCAL
        EXTERNAL FERM
        
        A=DIMAG(ZGLOCAL(EW))
        B=DREAL(ZGLOCAL(EW))
        FER=FERM(EW)
	DFER=FER*(1.D0-FER)*BETA 
        
        FASELOCAL=DFER*(DATAN(A/B))
        
        RETURN
        END

        FUNCTION FASEQD1(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
        COMMON/QD1/EF1,VQ1,EL1,TXF1,U1
        
        EXTERNAL FERM

        EF=EF1
        V=VQ1
        V2=VQ1*VQ1
        COU=U1
	ETTA=1.D-4
        EQ=EL1
 
        ZW=DCMPLX(EW,ETTA)

        CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
     *ZGFXATSS,ZG00,ZM) 
        
        A=DIMAG(ZGF)
        B=DREAL(ZGF)
        FER=FERM(EW)
	DFER=FER*(1.D0-FER)*BETA 
        
        FASEQD1=DFER*(DATAN(A/B))
        
        RETURN
        END

         FUNCTION FASEQD2(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
C        COMMON/QD1/EF1,VQ1,EL1,TXF1,U1
        COMMON/QD2/EF2,VQ2,EL2,TXF2,U2
        
        EXTERNAL FERM

        EF=EF2
        V=VQ2
        V2=VQ2*VQ2
        COU=U2
	ETTA=1.D-4
        EQ=EL2        

        ZW=DCMPLX(EW,ETTA)
        CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
     *ZGFXATSS,ZG00,ZM) 
        
        A=DIMAG(ZGF)
        B=DREAL(ZGF)
        FER=FERM(EW)
	DFER=FER*(1.D0-FER)*BETA 
        
        FASEQD2=DFER*(DATAN(A/B))
        
        RETURN
        END

        FUNCTION FASEC(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU

        EXTERNAL ZGCOEF1
        EXTERNAL FERM
        EXTERNAL ZGCON
        

        A=DIMAG(ZGCON(EW))
        B=DREAL(ZGCON(EW))
        FER=FERM(EW)
	DFER=FER*(1.D0-FER)*BETA 
        
        FASEC=(DFER*(DATAN(A/B)))
        
        RETURN
        END

         
        SUBROUTINE FASE(FASELO,FASECON,FASEQ1,FASEQ2)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU

        EXTERNAL FASELOCAL
        EXTERNAL FASEC
        EXTERNAL FASEQD1
        EXTERNAL FASEQD2

        EPSIL=1.0D-10
        T=1.0d0/BETA

        A=(-15.0D0*T)+AMU
        B=(15.0D0*T)+AMU
        
        FASELO=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,FASELOCAL)*CC
        FASECON=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,FASEC)*CC
        FASEQ1=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,FASEQD1)*CC
        FASEQ2=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,FASEQD2)*CC
        
        RETURN
        END

        

        FUNCTION ZGLOCAL(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/QD2/EF2,VQ2,EL2,TXF2,U2
        COMMON/QD1/EF1,VQ1,EL1,TXF1,U1
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU

        EXTERNAL ZGCOEF1

c        EF=EF1
c        V=VQ1
c        V2=VQ1*VQ1
c        COU=U1
	ETTA=1.0D-4
c        EQ=EL1
        
        Z1=ZGCOEF1(EW)

        EF=EF2
        V=VQ2
        V2=V*V
        COU=U2
        EQ=EL2
        
        ZW=DCMPLX(EW,ETTA)
        CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
     *ZGFXATSS,ZG00,ZM)
************************************************
        ZGLOCAL=Z1*VQ2*VQ2*ZGF
        RETURN
        END
  

       SUBROUTINE EXATA_UFINITO
C******************************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      INTEGER DT61
      PARAMETER (DT61=61)
      DIMENSION U(17),RF(17),EN(16),RC(13),DXP(16),AMINU(17),FER(17),
     *UD(16),RFD(16),RFM13(16),UM13(16),RFM24(16),UM24(16),US(16),
     *RFS(16),FERD(16),FERS(16),FER13(16),FER24(16),EN1(16),US3(16),
     *RFS3(16),Z3(13),ZZ3(3),A(7),FED(16),FED2(16),RFS2(16),US2(16),
     *RC11(16),UC11(16),RC22(16),UC22(16),RC33(16),UC33(16),
     *RC44(16),UC44(16)
      COMMON /RESI/U,EN,RF,RC
      COMMON/RESI2/UD,RFD,US,RFS
      COMMON/RESI3/UM13,RFM13,UM24,RFM24,US3,RFS3,US2,RFS2
      COMMON/RESIC/UC11,RC11,UC22,RC22,UC33,RC33,UC44,RC44
      COMMON/DADOS/AEF,AEQ,AMU,VV,VV2,BETA,DS,CC,D,COU
      COMMON/ET/ETA,OI,OF,NT
      COMMON/AMI/EMIN,KK
      COMMON/OCUP_EXAT/ANDD,ANSS,ANDDS,ANVAC,SOMA,ANDM13,ANDM24
      COMMON/OCUP_ATOMICO/ANVACATO,ANSSATO,ANDDATO
      COMMON/DISCRI/DISC
      COMMON/CORREC/DEK
      OPEN (UNIT=DT61,FILE='aresiduo.dat',STATUS='UNKNOWN')
      DATA C1 /2.8284271247D0/
      DATA PI /3.14159265358979323D0/
      T=1.D0/BETA
      AVG4=PI*VV/(2.D0*D)
      V=VV*AVG4
      V2=V*V
      EF=AEF-AMU
      EQ=AEQ-AMU
      DELTA=DSQRT((EF-EQ)*(EF-EQ)+4.0D0*V2)
      DELTALI=DSQRT((EF+COU-EQ)*(EF+COU-EQ)+4.0D0*V2)
      TANGFIDEN=EQ-EF+DELTA
      TANGFI=2.D0*V/TANGFIDEN
      DEN1=DSQRT(1.D0+TANGFI*TANGFI)
      SENFI=TANGFI/DEN1
      SENFI2=SENFI*SENFI
      COSFI=1.D0/DEN1
      COSFI2=COSFI*COSFI
      TANGTETADEN=EF+COU-EQ-DELTALI
      TANGTETA=2.D0*V/TANGTETADEN
      DEN2=DSQRT(1.D0+TANGTETA*TANGTETA)
      SENTETA=TANGTETA/DEN2
      SENTETA2=SENTETA*SENTETA
      COSTETA=1.D0/DEN2
      COSTETA2=COSTETA*COSTETA
      AA=EQ+EF
      AA2=AA*AA
      AA3=AA2*AA
      AB=2.D0*EF+COU
	AB2=AB*AB
	AB3=AB2*AB
	AC=2.D0*EQ
	AC2=AC*AC
	AC3=AC2*AC
	ARAUX=-3.D0*AA2*AB -3.D0*AA2*AC -3.D0*AB2*AA -3.D0*AB2*AC
     *-3D0*AC2*AA -3D0*AC2*AB +12.D0*AA*AB*AC +2.D0*AA3+2.D0*AB3+
     *2.D0*AC3 +36.D0*V2*AA -18.D0*V2*AB -18.D0*V2*AC
	CZ=1/54.D0
	AR=CZ*ARAUX


        AQAUX=12.D0*V2+AA2+AB2+AC2-AA*AB-AA*AC-AB*AC
	CZQ=-1/9.D0
	AQ=CZQ*AQAUX 
        AMENOSAQ=-AQ
	AMENOSAQ3=AMENOSAQ*AMENOSAQ*AMENOSAQ
	RAIZQ=DSQRT(AMENOSAQ)
	RAIZQ3=DSQRT(AMENOSAQ3)
        ATETA1=ACOS(AR/RAIZQ3)

	ARGCOS1=ATETA1/3.D0
	ARGCOSAUX2=2.D0*PI/3
	ARGCOS2=ARGCOSAUX2+ARGCOS1
	ARGCOSAUX3=4.D0*PI/3
	ARGCOS3=ARGCOSAUX3+ARGCOS1
        DISC=AQ*AQ*AQ+AR*AR
C................ AUTOVALORES DE ENERGIA..........
       DO 21 I=1,16
	EN(I)=0
21    CONTINUE
       EN(1)=0.D0
       EN(2)=0.5D0*(EF+EQ-DELTA)
       EN(3)=EN(2)
       EN(4)=0.5D0*(EF+EQ+DELTA)
       EN(5)=EN(4)
       EN(6)=EF+EQ
       EN(7)=EN(6)
       EN(8)=EN(6)


C CALCULO NUMERICO DA EQ. DO 3 GRAU
      NDEG=3
      A1=EQ+EF
      B1=2.D0*EF+COU
      C1=2.D0*EQ
      A(1)=1.D0
      A(2)=-(A1+B1+C1)
      A(3)=A1*B1+A1*C1+B1*C1-4.D0*V2
      A(4)=2.D0*V2*(B1+C1)-A1*B1*C1
        CALL ROOTS(A,Z3,NDEG,IER)
*********************************
        ZZ3(1)=Z3(1)
        ZZ3(2)=Z3(2)
        ZZ3(3)=Z3(3)
C        WRITE(6,10)'XXXXXXXXXXXXXXXXXXXX=',ZZ3(1),ZZ3(2),ZZ3(3)
C10      FORMAT(20F15.8,//)
        EN(9)=DREAL(ZZ3(1))
        EN(10)=DREAL(ZZ3(2))
        EN(11)=DREAL(ZZ3(3))
C FIM DO CALCULO NUMERICO
C	EF+COU/3.D0+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       EN(9)=2.0D0*RAIZQ*COS(ARGCOS1)
C       EN(10)=2.0D0*RAIZQ*COS(ARGCOS2)
C       EN(11)=2.0D0*RAIZQ*COS(ARGCOS3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       EN(12)=0.5D0*(3.D0*EQ+3.D0*EF+COU+DELTALI)
       EN(13)=EN(12)
       EN(14)=0.5D0*(3.D0*EQ+3.D0*EF+COU-DELTALI)
       EN(15)=EN(14)
       EN(16)=2.D0*EQ+2.D0*EF+COU

       AMINU(1)=DMIN1(EN(2),EN(1),EN(7),EN(4),EN(8),EN(5))
       AMINU(2)=DMIN1(EN(5),EN(1),EN(8),EN(3))
       AMINU(3)=DMIN1(EN(10),EN(2))
       AMINU(4)=MIN(EN(11),EN(2))
       AMINU(5)=DMIN1(EN(9),EN(2))
       AMINU(6)=DMIN1(EN(10),EN(4))
       AMINU(7)=DMIN1(EN(12),EN(6),EN(16),EN(14))
       AMINU(8)=DMIN1(EN(12),EN(9))
       AMINU(9)=DMIN1(EN(12),EN(10))
       AMINU(10)=DMIN1(EN(12),EN(11))
       AMINU(11)=DMIN1(EN(14),EN(10))
       AMINU(12)=DMIN1(EN(16),EN(12))
       AMINU(13)=DMIN1(EN(9),EN(4))
       AMINU(14)=DMIN1(EN(11),EN(4))
       AMINU(15)=DMIN1(EN(14),EN(6))
       AMINU(16)=DMIN1(EN(14),EN(9))
       AMINU(17)=DMIN1(EN(14),EN(11))

	ENAUX9=(EN(9)-2.D0*EF-COU)*(EN(9)-2.D0*EF-COU)
	EIAUX9=1.D0/ENAUX9
	EMAUX9=(EN(9)-2.D0*EQ)*(EN(9)-2.D0*EQ)
	EJAUX9=1.D0/EMAUX9
	AAAUX9=DSQRT(2.D0+4.D0*V2*(EIAUX9+EJAUX9))
        AA9=1.D0/AAAUX9

        CADA9=EN(9)-2.D0*EQ
	CAUX9=2.D0*V/CADA9
	CC9=CAUX9*AA9

	BADA9=EN(9)-2.D0*EF-COU
	BAUX9=2.D0*V/BADA9
	BB9=BAUX9*AA9


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	ENAUX10=(EN(10)-2.D0*EF-COU)*(EN(10)-2.D0*EF-COU)
	EIAUX10=1.D0/ENAUX10
	EMAUX10=(EN(10)-2.D0*EQ)*(EN(10)-2.D0*EQ)
	EJAUX10=1.D0/EMAUX10
	AAAUX10=SQRT(2.D0+4.D0*V2*(EIAUX10+EJAUX10))
        AA10=1.D0/AAAUX10

        CADA10=EN(10)-2.D0*EQ
	CAUX10=2.D0*V/CADA10
	CC10=CAUX10*AA10

	BADA10=EN(10)-2.D0*EF-COU
	BAUX10=2.D0*V/BADA10
	BB10=BAUX10*AA10

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	ENAUX11=(EN(11)-2.D0*EF-COU)*(EN(11)-2.D0*EF-COU)
	EIAUX11=1.D0/ENAUX11
	EMAUX11=(EN(11)-2.D0*EQ)*(EN(11)-2.D0*EQ)
	EJAUX11=1.D0/EMAUX11
	AAAUX11=DSQRT(2.D0+4.D0*V2*(EIAUX11+EJAUX11))
        AA11=1.D0/AAAUX11

        CADA11=EN(11)-2.D0*EQ
	CAUX11=2.D0*V/CADA11
	CC11=CAUX11*AA11

	BADA11=EN(11)-2.D0*EF-COU
	BAUX11=2.D0*V/BADA11
	BB11=BAUX11*AA11

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        EMIN=EN(1)
        DO 50 I=2,16
        IF(EMIN.GT.EN(I))EMIN=EN(I)
50      CONTINUE
        KK=0
        DO 51 I=1,16
        KK=KK+1
        IF(EMIN.EQ.EN(I))GO TO 52
51      CONTINUE
52      CONTINUE

C.......... FUNCAO DE PARTICAO............
C
      AZ=0.D0
      DO 100 I=1,16
      EN1(I)=EN(I)
      EN(I)=EN(I)-EMIN
      AZ=AZ+DEXP(-BETA*EN(I))
C      WRITE(DT91,70)I,EN(I),EN1(I)
100   CONTINUE
C........ DIFERENCAS DE ENERGIAS............
       U(1)=EN(2)-EN(1)
       U(2)=EN(5)-EN(1)
       U(3)=EN(10)-EN(2)
       U(4)=EN(11)-EN(2)
       U(5)=EN(9)-EN(2)
       U(6)=EN(10)-EN(4)
       U(7)=EN(12)-EN(6)
       U(8)=EN(12)-EN(9)
       U(9)=EN(12)-EN(10)
       U(10)=EN(12)-EN(11)
       U(11)=EN(14)-EN(10)
       U(12)=EN(14)-EN(11)
       U(13)=EN(9)-EN(4)
       U(14)=EN(11)-EN(4)
       U(15)=EN(14)-EN(6)
       U(16)=EN(14)-EN(9)
C       U(17)=EN(14)-EN(11)

C............ EXPONENCIAIS .........................
        DO 61 I=1,16
        DXP(I) = DEXP(-BETA*EN(I))
61     CONTINUE

CC............... RESIDUOS .........................

C.................SEM TERMOS DE U....................

      RF(1)=COSFI2*(DXP(1)+DXP(2)+1.5D0*DXP(6)+1.5D0*DXP(4))/AZ
      RF(2)=SENFI2*(DXP(1)+DXP(4)+1.5D0*DXP(2)+1.5D0*DXP(6))/AZ

C.................COM TERMOS DE U....................

      DESP31=DXP(3)+DXP(10)
	ACEM3=(AA10*SENFI)*(AA10*SENFI)
	BCEM3=(BB10*COSFI)*(BB10*COSFI)
      RF(3)=DESP31*(ACEM3+BCEM3)/AZ

      DESP41=DXP(3)+DXP(11)
	ACEM4=(AA11*SENFI)*(AA11*SENFI)
	BCEM4=(BB11*COSFI)*(BB11*COSFI)
      RF(4)=DESP41*(ACEM4+BCEM4)/AZ

	DESP51=DXP(3)+DXP(9)
	ACEM5=(AA9*SENFI)*(AA9*SENFI)
	BCEM5=(BB9*COSFI)*(BB9*COSFI)
      RF(5)=DESP51*(ACEM5+BCEM5)/AZ

	DESP61=DXP(4)+DXP(10)
	ACEM6=(AA10*COSFI)*(AA10*COSFI)
	BCEM6=(BB10*SENFI)*(BB10*SENFI)
      RF(6)=DESP61*(ACEM6+BCEM6)/AZ

	DESP71=DXP(8)+DXP(12)
	DESP72=DXP(15)+DXP(16)
	ACEM7=1.5d0*SENTETA2
	BCEM7=SENTETA2
      RF(7)=(DESP71*ACEM7+DESP72*BCEM7)/AZ

	DESP81=DXP(9)+DXP(12)
	ACEM8=(CC9*COSTETA)*(CC9*COSTETA)
	BCEM8=(AA9*SENTETA)*(AA9*SENTETA)
      RF(8)=DESP81*(ACEM8+BCEM8)/AZ

	DESP91=DXP(10)+DXP(12)
	ACEM9=(CC10*COSTETA)*(CC10*COSTETA)
	BCEM9=(AA10*SENTETA)*(AA10*SENTETA)
      RF(9)=DESP91*(ACEM9+BCEM9)/AZ

	DESP101=DXP(11)+DXP(12)
	ACEM10=(CC11*COSTETA)*(CC11*COSTETA)
	BCEM10=(AA11*SENTETA)*(AA11*SENTETA)
      RF(10)=DESP101*(ACEM10+BCEM10)/AZ

	DESP111=DXP(10)+DXP(14)
	ACEM11=(CC10*SENTETA)*(CC10*SENTETA)
	BCEM11=(AA10*COSTETA)*(AA10*COSTETA)
      RF(11)=DESP111*(ACEM11+BCEM11)/AZ

      	DESP171=DXP(11)+DXP(15)
	ACEM17=(AA11*COSTETA)*(AA11*COSTETA)
	BCEM17=(CC11*SENTETA)*(CC11*SENTETA)
      RF(12)=DESP171*(ACEM17+BCEM17)/AZ

	DESP131=DXP(5)+DXP(9)
	ACEM13=(AA9*COSFI)*(AA9*COSFI)
	BCEM13=(BB9*SENFI)*(BB9*SENFI)
      RF(13)=DESP131*(ACEM13+BCEM13)/AZ

	DESP141=DXP(5)+DXP(11)
	ACEM14=(AA11*COSFI)*(AA11*COSFI)
	BCEM14=(BB11*SENFI)*(BB11*SENFI)
      RF(14)=DESP141*(ACEM14+BCEM14)/AZ
        DESP121=DXP(13)+DXP(16)
	ACEM12=COSTETA2
	DESP151=DXP(8)+DXP(14)
	ACEM15=1.5D0*COSTETA2
	RF(15)=(DESP151*ACEM15)/AZ+(DESP121*ACEM12)/AZ

	DESP161=DXP(9)+DXP(15)
	ACEM16=(AA9*COSTETA)*(AA9*COSTETA)
	BCEM16=(CC9*SENTETA)*(CC9*SENTETA)
      RF(16)=DESP161*(ACEM16+BCEM16)/AZ

          SRF=0.D0
	  DO 56 I=1,17
	  SRF=SRF+RF(I)
C        WRITE(91,454)I,RF(I),U(I),EN(I),AMINU(I),EMIN
C454     FORMAT('I=',I5,3X,'RESF=',F20.10,3X,'POLOS (Ui)=',F20.10,F15.8/,
C     *'EN=',F20.10,3X,'AMINU=',F20.10,3X,'EMIN=',F20.10)
56      CONTINUE


C	WRITE(DT91,70)T,RF(1),RF(2),RF(3),RF(4),RF(5),RF(6),
C     *RF(7),RF(8),RF(9),RF(10),RF(11),RF(12),RF(13),
C     *RF(14),RF(15),RF(16),RF(17)
70    FORMAT(2F15.8,I10)

CCCCCCCCCCCCCCCCCCCC ANIQUILA SPIN UP m13CCCCCCCCCCCCCCCCCCCCCCC
      UM13(1)=EN(9)-EN(3)
      UM13(2)=EN(10)-EN(3)
      UM13(3)=EN(11)-EN(3)
      UM13(4)=EN(9)-EN(5)
      UM13(5)=EN(10)-EN(5)
      UM13(6)=EN(11)-EN(5)
      UM13(7)=EN(12)-EN(9)
      UM13(8)=EN(12)-EN(10)
      UM13(9)=EN(12)-EN(11)
      UM13(10)=EN(14)-EN(9)
      UM13(11)=EN(14)-EN(10)
      UM13(12)=EN(14)-EN(11)


      RFM13(1)=(BB9*COSFI)*(-AA9*SENFI)*(DXP(3)+DXP(9))/AZ
      RFM13(2)=(BB10*COSFI)*(-AA10*SENFI)*(DXP(3)+DXP(10))/AZ
      RFM13(3)=(BB11*COSFI)*(-AA11*SENFI)*(DXP(3)+DXP(11))/AZ
      RFM13(4)=(BB9*SENFI)*(AA9*COSFI)*(DXP(5)+DXP(9))/AZ
      RFM13(5)=(BB10*SENFI)*(AA10*COSFI)*(DXP(5)+DXP(10))/AZ
      RFM13(6)=(BB11*SENFI)*(AA11*COSFI)*(DXP(5)+DXP(11))/AZ
      RFM13(7)=(-AA9*SENTETA)*(CC9*COSTETA)*(DXP(12)+DXP(9))/AZ
      RFM13(8)=(-AA10*SENTETA)*(CC10*COSTETA)*(DXP(12)+DXP(10))/AZ
      RFM13(9)=(-AA11*SENTETA)*(CC11*COSTETA)*(DXP(12)+DXP(11))/AZ
      RFM13(10)=(AA9*COSTETA)*(CC9*SENTETA)*(DXP(14)+DXP(9))/AZ
      RFM13(11)=(AA10*COSTETA)*(CC10*SENTETA)*(DXP(14)+DXP(10))/AZ
      RFM13(12)=(AA11*COSTETA)*(CC11*SENTETA)*(DXP(14)+DXP(11))/AZ


CCCCCCCCCCCCCCCCCCCC ANIQUILA SPIN DOWN m24 CCCCCCCCCCCCCCCCCCCCCCC
      UM24(1)=EN(9)-EN(2)
      UM24(2)=EN(10)-EN(2)
      UM24(3)=EN(11)-EN(2)
      UM24(4)=EN(9)-EN(4)
      UM24(5)=EN(10)-EN(4)
      UM24(6)=EN(11)-EN(4)
      UM24(7)=EN(13)-EN(9)
      UM24(8)=EN(13)-EN(10)
      UM24(9)=EN(13)-EN(11)
      UM24(10)=EN(15)-EN(9)
      UM24(11)=EN(15)-EN(10)
      UM24(12)=EN(15)-EN(11)


      RFM24(1)=(BB9*COSFI)*(AA9*SENFI)*(DXP(2)+DXP(9))/AZ
      RFM24(2)=(BB10*COSFI)*(AA10*SENFI)*(DXP(2)+DXP(10))/AZ
      RFM24(3)=(BB11*COSFI)*(AA11*SENFI)*(DXP(2)+DXP(11))/AZ
      RFM24(4)=(BB9*SENFI)*(-AA9*COSFI)*(DXP(4)+DXP(9))/AZ
      RFM24(5)=(BB10*SENFI)*(-AA10*COSFI)*(DXP(4)+DXP(10))/AZ
      RFM24(6)=(BB11*SENFI)*(-AA11*COSFI)*(DXP(4)+DXP(11))/AZ
      RFM24(7)=(AA9*SENTETA)*(CC9*COSTETA)*(DXP(13)+DXP(9))/AZ
      RFM24(8)=(AA10*SENTETA)*(CC10*COSTETA)*(DXP(13)+DXP(10))/AZ
      RFM24(9)=(AA11*SENTETA)*(CC11*COSTETA)*(DXP(13)+DXP(11))/AZ
      RFM24(10)=(-AA9*COSTETA)*(CC9*SENTETA)*(DXP(15)+DXP(9))/AZ
      RFM24(11)=(-AA10*COSTETA)*(CC10*SENTETA)*(DXP(15)+DXP(10))/AZ
      RFM24(12)=(-AA11*COSTETA)*(CC11*SENTETA)*(DXP(15)+DXP(11))/AZ


C         nao mecher daqui p baixo esta correto!!!!

CCCCCCCCCCCCCCCCCCCC OCUPACAO DUPLA m44CCCCCCCCCCCCCCCCCCCCCCC
      UD(1)=EN(9)-EN(2)
      UD(2)=EN(10)-EN(2)
      UD(3)=EN(11)-EN(2)
      UD(4)=EN(9)-EN(4)
      UD(5)=EN(10)-EN(4)
      UD(6)=EN(11)-EN(4)
      UD(7)=EN(13)-EN(9)
      UD(8)=EN(13)-EN(10)
      UD(9)=EN(13)-EN(11)
      UD(10)=EN(15)-EN(9)
      UD(11)=EN(15)-EN(10)
      UD(12)=EN(14)-EN(11)
      UD(13)=EN(12)-EN(6)
      UD(14)=EN(14)-EN(6)
      UD(15)=EN(16)-EN(12)
      UD(16)=EN(16)-EN(14)
      RFD(1)=(BB9*COSFI)*(BB9*COSFI)*(DXP(2)+DXP(9))/AZ
      RFD(2)=(BB10*COSFI)*(BB10*COSFI)*(DXP(2)+DXP(10))/AZ
      RFD(3)=(BB11*COSFI)*(BB11*COSFI)*(DXP(2)+DXP(11))/AZ
      RFD(4)=(BB9*SENFI)*(BB9*SENFI)*(DXP(4)+DXP(9))/AZ
      RFD(5)=(BB10*SENFI)*(BB10*SENFI)*(DXP(4)+DXP(10))/AZ
      RFD(6)=(BB11*SENFI)*(BB11*SENFI)*(DXP(4)+DXP(11))/AZ
      RFD(7)=(AA9*SENTETA)*(AA9*SENTETA)*(DXP(13)+DXP(9))/AZ
      RFD(8)=(AA10*SENTETA)*(AA10*SENTETA)*(DXP(13)+DXP(10))/AZ
      RFD(9)=(AA11*SENTETA)*(AA11*SENTETA)*(DXP(13)+DXP(11))/AZ
      RFD(10)=(AA9*COSTETA)*(AA9*COSTETA)*(DXP(15)+DXP(9))/AZ
      RFD(11)=(AA10*COSTETA)*(AA10*COSTETA)*(DXP(15)+DXP(10))/AZ
      RFD(12)=(AA11*COSTETA)*(AA11*COSTETA)*(DXP(15)+DXP(11))/AZ
      RFD(13)=1.5D0*SENTETA*SENTETA*(DXP(6)+DXP(12))/AZ
      RFD(14)=1.5D0*COSTETA*COSTETA*(DXP(6)+DXP(14))/AZ
      RFD(15)=COSTETA*COSTETA*(DXP(12)+DXP(16))/AZ
      RFD(16)=SENTETA*SENTETA*(DXP(14)+DXP(16))/AZ
CCCCCCCCCCCCCCCCCCCC OCUPACAO SIMPLES m11CCCCCCCCCCCCCCCCCCCCC
      US(1)=EN(9)-EN(5)
      US(2)=EN(10)-EN(5)
      US(3)=EN(11)-EN(5)
      US(4)=EN(9)-EN(3)
      US(5)=EN(10)-EN(3)
      US(6)=EN(11)-EN(3)
      US(7)=EN(14)-EN(9)
      US(8)=EN(14)-EN(10)
      US(9)=EN(14)-EN(11)
      US(10)=EN(12)-EN(9)
      US(11)=EN(12)-EN(10)
      US(12)=EN(12)-EN(11)
      US(13)=EN(6)-EN(2)
      US(14)=EN(6)-EN(4)
      US(15)=EN(2)-EN(1)
      US(16)=EN(4)-EN(1)
      RFS(1)=(AA9*COSFI)*(AA9*COSFI)*(DXP(5)+DXP(9))/AZ
      RFS(2)=(AA10*COSFI)*(AA10*COSFI)*(DXP(5)+DXP(10))/AZ
      RFS(3)=(AA11*COSFI)*(AA11*COSFI)*(DXP(5)+DXP(11))/AZ

      RFS(4)=(AA9*SENFI)*(AA9*SENFI)*(DXP(3)+DXP(9))/AZ
      RFS(5)=(AA10*SENFI)*(AA10*SENFI)*(DXP(3)+DXP(10))/AZ
      RFS(6)=(AA11*SENFI)*(AA11*SENFI)*(DXP(3)+DXP(11))/AZ

      RFS(7)=(CC9*SENTETA)*(CC9*SENTETA)*(DXP(14)+DXP(9))/AZ
      RFS(8)=(CC10*SENTETA)*(CC10*SENTETA)*(DXP(14)+DXP(10))/AZ
      RFS(9)=(CC11*SENTETA)*(CC11*SENTETA)*(DXP(14)+DXP(11))/AZ

      RFS(10)=(CC9*COSTETA)*(CC9*COSTETA)*(DXP(12)+DXP(9))/AZ
      RFS(11)=(CC10*COSTETA)*(CC10*COSTETA)*(DXP(12)+DXP(10))/AZ
      RFS(12)=(CC11*COSTETA)*(CC11*COSTETA)*(DXP(12)+DXP(11))/AZ

      RFS(13)=1.5D0*SENFI*SENFI*(DXP(2)+DXP(6))/AZ
      RFS(14)=1.5D0*COSFI*COSFI*(DXP(6)+DXP(4))/AZ
      RFS(15)=COSFI*COSFI*(DXP(1)+DXP(2))/AZ
      RFS(16)=SENFI*SENFI*(DXP(1)+DXP(4))/AZ
C      DO 105 IJ=1,16
C      WRITE(DT91,70)IJ,US(IJ),RFS(IJ),UD(IJ),RFD(IJ)
C105   CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    m33    CCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      US3(1)=EN(9)-EN(3)
      US3(2)=EN(10)-EN(3)
      US3(3)=EN(11)-EN(3)

      US3(4)=EN(9)-EN(4)
      US3(5)=EN(10)-EN(4)
      US3(6)=EN(11)-EN(4)

      US3(7)=EN(12)-EN(8)

      US3(8)=EN(12)-EN(9)
      US3(9)=EN(12)-EN(10)
      US3(10)=EN(12)-EN(11)

      US3(11)=EN(14)-EN(8)

      US3(12)=EN(14)-EN(9)
      US3(13)=EN(14)-EN(10)
      US3(14)=EN(14)-EN(11)

      US3(15)=EN(16)-EN(13)

      US3(16)=EN(16)-EN(15)

      RFS3(1)=(BB9*COSFI)*(BB9*COSFI)*(DXP(3)+DXP(9))/AZ
      RFS3(2)=(BB10*COSFI)*(BB10*COSFI)*(DXP(3)+DXP(10))/AZ
      RFS3(3)=(BB11*COSFI)*(BB11*COSFI)*(DXP(3)+DXP(11))/AZ

      RFS3(4)=(BB9*SENFI)*(BB9*SENFI)*(DXP(4)+DXP(9))/AZ
      RFS3(5)=(BB10*SENFI)*(BB10*SENFI)*(DXP(4)+DXP(10))/AZ
      RFS3(6)=(BB11*SENFI)*(BB11*SENFI)*(DXP(4)+DXP(11))/AZ

      RFS3(7)=1.5D0*SENTETA*SENTETA*(DXP(8)+DXP(12))/AZ

      RFS3(8)=(AA9*SENTETA)*(AA9*SENTETA)*(DXP(12)+DXP(9))/AZ
      RFS3(9)=(AA10*SENTETA)*(AA10*SENTETA)*(DXP(12)+DXP(10))/AZ
      RFS3(10)=(AA11*SENTETA)*(AA11*SENTETA)*(DXP(12)+DXP(11))/AZ

      RFS3(11)=1.5D0*COSTETA*COSTETA*(DXP(14)+DXP(8))/AZ

      RFS3(12)=(AA9*COSTETA)*(AA9*COSTETA)*(DXP(14)+DXP(9))/AZ
      RFS3(13)=(AA10*COSTETA)*(AA10*COSTETA)*(DXP(14)+DXP(10))/AZ
      RFS3(14)=(AA11*COSTETA)*(AA11*COSTETA)*(DXP(14)+DXP(11))/AZ

      RFS3(15)=COSTETA*COSTETA*(DXP(16)+DXP(13))/AZ
      RFS3(16)=SENTETA*SENTETA*(DXP(16)+DXP(15))/AZ

C      DO 105 IJ=1,17
C      WRITE(DT61,70)U(IJ),RF(IJ),IJ
C105   CONTINUE
C      WRITE(DT61,75)EF,RF(1),RF(2),RF(3),RF(4),RF(5),RF(6),RF(7),RF(8),
C     *RF(9),RF(10),RF(11),RF(12),RF(13),RF(14),RF(15),RF(16),DEK
C75    FORMAT(30F15.8)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    m22    CCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      US2(1)=EN(3)-EN(1)
      US2(2)=EN(5)-EN(1)
      US2(3)=EN(7)-EN(3)
      US2(4)=EN(7)-EN(5)

      US2(5)=EN(9)-EN(2)
      US2(6)=EN(10)-EN(2)
      US2(7)=EN(11)-EN(2)

      US2(8)=EN(9)-EN(4)
      US2(9)=EN(10)-EN(4)
      US2(10)=EN(11)-EN(4)

      US2(11)=EN(13)-EN(9)
      US2(12)=EN(13)-EN(10)
      US2(13)=EN(13)-EN(11)

      US2(14)=EN(15)-EN(9)
      US2(15)=EN(15)-EN(10)
      US2(16)=EN(15)-EN(11)

      RFS2(1)=COSFI*COSFI*(DXP(3)+DXP(1))/AZ
      RFS2(2)=SENFI*SENFI*(DXP(5)+DXP(1))/AZ
      RFS2(3)=1.5D0*SENFI*SENFI*(DXP(7)+DXP(3))/AZ
      RFS2(4)=1.5D0*COSFI*COSFI*(DXP(7)+DXP(5))/AZ

      RFS2(5)=(AA9*SENFI)*(AA9*SENFI)*(DXP(9)+DXP(2))/AZ
      RFS2(6)=(AA10*SENFI)*(AA10*SENFI)*(DXP(10)+DXP(2))/AZ
      RFS2(7)=(AA11*SENFI)*(AA11*SENFI)*(DXP(11)+DXP(2))/AZ

      RFS2(8)=(AA9*COSFI)*(AA9*COSFI)*(DXP(9)+DXP(4))/AZ
      RFS2(9)=(AA10*COSFI)*(AA10*COSFI)*(DXP(10)+DXP(4))/AZ
      RFS2(10)=(AA11*COSFI)*(AA11*COSFI)*(DXP(11)+DXP(4))/AZ

      RFS2(11)=(CC9*COSTETA)*(CC9*COSTETA)*(DXP(13)+DXP(9))/AZ
      RFS2(12)=(CC10*COSTETA)*(CC10*COSTETA)*(DXP(13)+DXP(10))/AZ
      RFS2(13)=(CC11*COSTETA)*(CC11*COSTETA)*(DXP(13)+DXP(11))/AZ

      RFS2(14)=(CC9*SENTETA)*(CC9*SENTETA)*(DXP(15)+DXP(9))/AZ
      RFS2(15)=(CC10*SENTETA)*(CC10*SENTETA)*(DXP(15)+DXP(10))/AZ
      RFS2(16)=(CC11*SENTETA)*(CC11*SENTETA)*(DXP(15)+DXP(11))/AZ
Cccccccccccccccccccccccc ELETRONS DE CONDUCAO cccccccccccccccccccccc
      UC11(1)=EN(2)-EN(1)
      UC11(2)=EN(4)-EN(1)
      UC11(3)=EN(6)-EN(2)
      UC11(4)=EN(6)-EN(4)
      UC11(5)=EN(9)-EN(3)
      UC11(6)=EN(10)-EN(3)
      UC11(7)=EN(11)-EN(3)
      UC11(8)=EN(9)-EN(5)
      UC11(9)=EN(10)-EN(5)
      UC11(10)=EN(11)-EN(5)
      UC11(11)=EN(14)-EN(9)
      UC11(12)=EN(14)-EN(10)
      UC11(13)=EN(14)-EN(11)
      UC11(14)=EN(12)-EN(9)
      UC11(15)=EN(12)-EN(10)
      UC11(16)=EN(12)-EN(11)

      RC11(1)=SENFI*SENFI*(DXP(1)+DXP(2))/AZ
      RC11(2)=COSFI*COSFI*(DXP(1)+DXP(4))/AZ
      RC11(3)=1.5D0*COSFI*COSFI*(DXP(2)+DXP(6))/AZ
      RC11(4)=1.5D0*SENFI*SENFI*(DXP(4)+DXP(6))/AZ

      RC11(5)=(AA9*COSFI)*(AA9*COSFI)*(DXP(9)+DXP(3))/AZ
      RC11(6)=(AA10*COSFI)*(AA10*COSFI)*(DXP(10)+DXP(3))/AZ
      RC11(7)=(AA11*COSFI)*(AA11*COSFI)*(DXP(11)+DXP(3))/AZ

      RC11(8)=(AA9*SENFI)*(AA9*SENFI)*(DXP(9)+DXP(5))/AZ
      RC11(9)=(AA10*SENFI)*(AA10*SENFI)*(DXP(10)+DXP(5))/AZ
      RC11(10)=(AA11*SENFI)*(AA11*SENFI)*(DXP(11)+DXP(5))/AZ

      RC11(11)=(BB9*COSTETA)*(BB9*COSTETA)*(DXP(14)+DXP(9))/AZ
      RC11(12)=(BB10*COSTETA)*(BB10*COSTETA)*(DXP(14)+DXP(10))/AZ
      RC11(13)=(BB11*COSTETA)*(BB11*COSTETA)*(DXP(14)+DXP(11))/AZ

      RC11(14)=(BB9*SENTETA)*(BB9*SENTETA)*(DXP(12)+DXP(9))/AZ
      RC11(15)=(BB10*SENTETA)*(BB10*SENTETA)*(DXP(12)+DXP(10))/AZ
      RC11(16)=(BB11*SENTETA)*(BB11*SENTETA)*(DXP(12)+DXP(11))/AZ
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC RC22 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      UC22(1)=EN(3)-EN(1)
      UC22(2)=EN(5)-EN(1)
      UC22(3)=EN(7)-EN(3)
      UC22(4)=EN(7)-EN(5)
      UC22(5)=EN(9)-EN(2)
      UC22(6)=EN(10)-EN(2)
      UC22(7)=EN(11)-EN(2)
      UC22(8)=EN(9)-EN(4)
      UC22(9)=EN(10)-EN(4)
      UC22(10)=EN(11)-EN(4)
      UC22(11)=EN(13)-EN(9)
      UC22(12)=EN(13)-EN(10)
      UC22(13)=EN(13)-EN(11)
      UC22(14)=EN(15)-EN(9)
      UC22(15)=EN(15)-EN(10)
      UC22(16)=EN(15)-EN(11)

      RC22(1)=SENFI*SENFI*(DXP(1)+DXP(3))/AZ
      RC22(2)=COSFI*COSFI*(DXP(1)+DXP(5))/AZ
      RC22(3)=1.5D0*COSFI*COSFI*(DXP(3)+DXP(7))/AZ
      RC22(4)=1.5D0*SENFI*SENFI*(DXP(5)+DXP(7))/AZ
      RC22(5)=(AA9*COSFI)*(AA9*COSFI)*(DXP(9)+DXP(2))/AZ
      RC22(6)=(AA10*COSFI)*(AA10*COSFI)*(DXP(10)+DXP(2))/AZ
      RC22(7)=(AA11*COSFI)*(AA11*COSFI)*(DXP(11)+DXP(2))/AZ
      RC22(8)=(AA9*SENFI)*(AA9*SENFI)*(DXP(9)+DXP(4))/AZ
      RC22(9)=(AA10*SENFI)*(AA10*SENFI)*(DXP(10)+DXP(4))/AZ
      RC22(10)=(AA11*SENFI)*(AA11*SENFI)*(DXP(11)+DXP(4))/AZ
      RC22(11)=(BB9*SENTETA)*(BB9*SENTETA)*(DXP(13)+DXP(9))/AZ
      RC22(12)=(BB10*SENTETA)*(BB10*SENTETA)*(DXP(13)+DXP(10))/AZ
      RC22(13)=(BB11*SENTETA)*(BB11*SENTETA)*(DXP(13)+DXP(11))/AZ
      RC22(14)=(BB9*COSTETA)*(BB9*COSTETA)*(DXP(15)+DXP(9))/AZ
      RC22(15)=(BB10*COSTETA)*(BB10*COSTETA)*(DXP(15)+DXP(10))/AZ
      RC22(16)=(BB11*COSTETA)*(BB11*COSTETA)*(DXP(15)+DXP(11))/AZ
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC RC33 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      UC33(1)=EN(9)-EN(3)
      UC33(2)=EN(10)-EN(3)
      UC33(3)=EN(11)-EN(3)
      UC33(4)=EN(9)-EN(5)
      UC33(5)=EN(10)-EN(5)
      UC33(6)=EN(11)-EN(5)
      UC33(7)=EN(12)-EN(9)
      UC33(8)=EN(12)-EN(10)
      UC33(9)=EN(12)-EN(11)
      UC33(10)=EN(14)-EN(9)
      UC33(11)=EN(14)-EN(10)
      UC33(12)=EN(14)-EN(11)
      UC33(13)=EN(12)-EN(8)
      UC33(14)=EN(14)-EN(8)
      UC33(15)=EN(16)-EN(13)
      UC33(16)=EN(16)-EN(15)

      RC33(1)=(CC9*SENFI)*(CC9*SENFI)*(DXP(9)+DXP(3))/AZ
      RC33(2)=(CC10*SENFI)*(CC10*SENFI)*(DXP(10)+DXP(3))/AZ
      RC33(3)=(CC11*SENFI)*(CC11*SENFI)*(DXP(11)+DXP(3))/AZ
      RC33(4)=(CC9*COSFI)*(CC9*COSFI)*(DXP(9)+DXP(5))/AZ
      RC33(5)=(CC10*COSFI)*(CC10*COSFI)*(DXP(10)+DXP(5))/AZ
      RC33(6)=(CC11*COSFI)*(CC11*COSFI)*(DXP(11)+DXP(5))/AZ
      RC33(7)=(AA9*COSTETA)*(AA9*COSTETA)*(DXP(12)+DXP(9))/AZ
      RC33(8)=(AA10*COSTETA)*(AA10*COSTETA)*(DXP(12)+DXP(10))/AZ
      RC33(9)=(AA11*COSTETA)*(AA11*COSTETA)*(DXP(12)+DXP(11))/AZ
      RC33(10)=(AA9*SENTETA)*(AA9*SENTETA)*(DXP(14)+DXP(9))/AZ
      RC33(11)=(AA10*SENTETA)*(AA10*SENTETA)*(DXP(14)+DXP(10))/AZ
      RC33(12)=(AA11*SENTETA)*(AA11*SENTETA)*(DXP(14)+DXP(11))/AZ
      RC33(13)=1.5D0*COSTETA*COSTETA*(DXP(12)+DXP(8))/AZ
      RC33(14)=1.5D0*SENTETA*SENTETA*(DXP(14)+DXP(8))/AZ
      RC33(15)=SENTETA*SENTETA*(DXP(13)+DXP(16))/AZ
      RC33(16)=COSTETA*COSTETA*(DXP(15)+DXP(16))/AZ

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC RC44 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      UC44(1)=EN(9)-EN(2)
      UC44(2)=EN(10)-EN(2)
      UC44(3)=EN(11)-EN(2)
      UC44(4)=EN(9)-EN(4)
      UC44(5)=EN(10)-EN(4)
      UC44(6)=EN(11)-EN(4)
      UC44(7)=EN(13)-EN(9)
      UC44(8)=EN(13)-EN(10)
      UC44(9)=EN(13)-EN(11)
      UC44(10)=EN(15)-EN(9)
      UC44(11)=EN(15)-EN(10)
      UC44(12)=EN(15)-EN(11)
      UC44(13)=EN(12)-EN(6)
      UC44(14)=EN(14)-EN(6)
      UC44(15)=EN(16)-EN(12)
      UC44(16)=EN(16)-EN(14)

      RC44(1)=(CC9*SENFI)*(CC9*SENFI)*(DXP(9)+DXP(2))/AZ
      RC44(2)=(CC10*SENFI)*(CC10*SENFI)*(DXP(10)+DXP(2))/AZ
      RC44(3)=(CC11*SENFI)*(CC11*SENFI)*(DXP(11)+DXP(2))/AZ
      RC44(4)=(CC9*COSFI)*(CC9*COSFI)*(DXP(9)+DXP(4))/AZ
      RC44(5)=(CC10*COSFI)*(CC10*COSFI)*(DXP(10)+DXP(4))/AZ
      RC44(6)=(CC11*COSFI)*(CC11*COSFI)*(DXP(11)+DXP(4))/AZ
      RC44(7)=(AA9*COSTETA)*(AA9*COSTETA)*(DXP(13)+DXP(9))/AZ
      RC44(8)=(AA10*COSTETA)*(AA10*COSTETA)*(DXP(13)+DXP(10))/AZ
      RC44(9)=(AA11*COSTETA)*(AA11*COSTETA)*(DXP(13)+DXP(11))/AZ
      RC44(10)=(AA9*SENTETA)*(AA9*SENTETA)*(DXP(15)+DXP(9))/AZ
      RC44(11)=(AA10*SENTETA)*(AA10*SENTETA)*(DXP(15)+DXP(10))/AZ
      RC44(12)=(AA11*SENTETA)*(AA11*SENTETA)*(DXP(15)+DXP(11))/AZ
      RC44(13)=1.5D0*COSTETA*COSTETA*(DXP(12)+DXP(6))/AZ
      RC44(14)=1.5D0*SENTETA*SENTETA*(DXP(14)+DXP(6))/AZ
      RC44(15)=SENTETA*SENTETA*(DXP(16)+DXP(12))/AZ
      RC44(16)=COSTETA*COSTETA*(DXP(16)+DXP(14))/AZ


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC FIM DAS FG CCCCCCCCCCCCCCCCCCCC
        AN22=0.D0
        AN33=0.D0
        ANDD=0.D0
	  ANSS=0.D0
	  ANVAC=0.D0
        DO 22 I=1,16
	  FERS(I)=FERM(US(I))
	  ANS=FERS(I)*RFS(I)
        ANSS=ANSS+ANS
        AN2=FERM(US2(I))*RFS2(I)
	  ANVACC=(1.D0-FERS(I))*RFS(I)
	  ANVAC=ANVAC+ANVACC

	  AN22=AN22+AN2

	  FED(I)=FERM(US3(I))
        AN3=FED(I)*RFS3(I)
        AN33=AN33+AN3

 	  FERD(I)=FERM(UD(I))
        AND=FERD(I)*RFD(I)
        ANDD=ANDD+AND

22      CONTINUE
C          IF(ANSS.LT.0.5D0)GO TO 18
C           ANDD=0.5D0*ANDD
C           ANSS=0.5D0*ANSS
C           AVNSS=0.5D0*AVNSS
C           ANDDS=0.5D0*ANDDS
C18        CONTINUE

        DO 23 K=1,12
	  FER13(K)=FERM(UM13(K))
	  AM13=FER13(K)*RFM13(K)
	  FER24(K)=FERM(UM24(K))
	  AM24=FER24(K)*RFM24(K)
	ANDM13=ANDM13+AM13
	ANDM24=ANDM24+AM24
23      CONTINUE

      GATUP=ANSS+AN33
C      GATUP=ANSS
C      GATDOWN=AN22+ANDD
      GATDOWN=AN22
      ATOT=GATUP+GATDOWN+ANVAC
C	WRITE(DT91,*)ANSS,AN33,ANDD,ANVAC,ATOT,ANDM13
C        WRITE(6,*)'TTTTTTTTTT=',ANSS,AN33,ANDD,ANVAC,ATOT,ANDM13
C         write(6,*)'SSSSSSSSSSSSSSS=',ANSS,AN33,ANDD,ANVAC,ATOT
C	WRITE(6,24)ANSS,ANDM13,AN33,AN22,ANDM24,ANDD,ANVAC,GATDOWN,
C    *GATUP,ATOT
C24    FORMAT('GAT11=',F13.10,2X,'GAT13=',F13.10,2X,'GAT33=',F13.10,2X,/,
C     *'GAT22=',F13.10,2X,'GAT24=',F13.10,2X,'GAT44=',F13.10,2X,/,
C     *'ANVAC=',F13.10,2X,'GATDOWN=',F13.10,2X,'GATUP=',F13.10,2X,
C     *'ATOT=',F13.10)
        SS=ANVAC+2.D0*ANSS+ANDD
        ANVACATO=ANVAC
        ANSSATO=ANSS
        ANDDATO=ANDD
C	WRITE(DT91,70)ANVAC,ANSS,ANDD,SS
C	WRITE(6,*)'XXXXXXXXXXXXXXXXX=',ANVAC,ANSS,ANDD,SS
       RETURN
       END
C********************  FIM DA SUB EXATA_UFINITO ********************

        SUBROUTINE ROOTS(A,Z,NDEG,IER)
C       *****************************
        REAL*8 A(11),RZ(20)
        COMPLEX*16 Z(10),ZA(10)
        EQUIVALENCE (ZA(1),RZ(1))
C
C       Calcula las raices de un polinomio de  grado NDEG < 10
C       con coeficientes reales
C
        CALL ZRPOLY(A,NDEG,RZ,IER)
        DO 20 I=1,NDEG
        Z(I) = ZA(I)
   20   CONTINUE
        RETURN
        END
***********************************************************************
      SUBROUTINE GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
     *ZGFXATSS,ZG00,ZM)
*******************************************************************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      REAL*8 MU,MUC
      DIMENSION U(17),RF(17),EN(16),RC(13)
      DIMENSION UU(17),MU(17),UN(16),MUC(13),
     *RFD(16),RFS(16),UD(16),US(16),UM13(16),
     *UM24(16),RFM13(16),RFM24(16),US3(16),RFS3(16),US2(16),RFS2(16),
     *UC11(16),RC11(16),UC22(16),RC22(16),UC33(16),RC33(16),UC44(16),
     *RC44(16)
      COMMON /RESI/U,EN,RF,RC
      COMMON/RESI2/UD,RFD,US,RFS
      COMMON/RESI3/UM13,RFM13,UM24,RFM24,US3,RFS3,US2,RFS2
      COMMON/RESIC/UC11,RC11,UC22,RC22,UC33,RC33,UC44,RC44
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
      COMMON/OCUP_EXAT/ANDD,ANSS,ANDDS,AVNSS,SOMA,ANDM13,ANDM24
      COMMON/OCUP_ATOMICO/ANVACATO,ANSSATO,ANDDATO
c      COMMON/DELTA/DELTA

        DATA PI /3.14159265358979323D0/
        GAMMA=2.D0*D/PI
        DELTA=0.5D0*PI*V2/D
        GAMMA2=GAMMA*GAMMA
        D2=D*D
        ZW2=ZW*ZW
C        ETTA=0.0031D0
        ETTA=1.D-6
        AVG4=PI*V/(2.D0*D)
        VG2=V*V*AVG4*AVG4
        BW=DREAL(ZW)
        ZWD=DCMPLX(BW,ETTA)
C      ZGFT=RF(1)/(ZWD-U(1))+RF(2)/(ZWD-U(2))+RF(3)/(ZWD-U(3))+
C     *RF(4)/(ZWD-U(4))+RF(5)/(ZWD-U(5))+RF(6)/(ZWD-U(6))+
C     *RF(7)/(ZWD-U(7))+RF(8)/(ZWD-U(8))+RF(9)/(ZWD-U(9))+
C     *RF(10)/(ZWD-U(10))+RF(11)/(ZWD-U(11))+RF(12)/(ZWD-U(12))+
C     *RF(13)/(ZWD-U(13))+RF(14)/(ZWD-U(14))+RF(15)/(ZWD-U(15))+
C     *RF(16)/(ZWD-U(16))+RF(17)/(ZWD-U(17))

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC DUPLA E SIMPLES CCCCCCCCCCCC
      ZGFXATDD=RFD(1)/(ZWD-UD(1))+RFD(2)/(ZWD-UD(2))+RFD(3)/(ZWD-UD(3))+
     *RFD(4)/(ZWD-UD(4))+RFD(5)/(ZWD-UD(5))+RFD(6)/(ZWD-UD(6))+
     *RFD(7)/(ZWD-UD(7))+RFD(8)/(ZWD-UD(8))+RFD(9)/(ZWD-UD(9))+
     *RFD(10)/(ZWD-UD(10))+RFD(11)/(ZWD-UD(11))+RFD(12)/(ZWD-UD(12))+
     *RFD(13)/(ZWD-UD(13))+RFD(14)/(ZWD-UD(14))+RFD(15)/(ZWD-UD(15))+
     *RFD(16)/(ZWD-UD(16))
        ZGFXATDD=-ZGFXATDD
CC
      ZGFXATSS=RFS(1)/(ZWD-US(1))+RFS(2)/(ZWD-US(2))+RFS(3)/(ZWD-US(3))+
     *RFS(4)/(ZWD-US(4))+RFS(5)/(ZWD-US(5))+RFS(6)/(ZWD-US(6))+
     *RFS(7)/(ZWD-US(7))+RFS(8)/(ZWD-US(8))+RFS(9)/(ZWD-US(9))+
     *RFS(10)/(ZWD-US(10))+RFS(11)/(ZWD-US(11))+RFS(12)/(ZWD-US(12))+
     *RFS(13)/(ZWD-US(13))+RFS(14)/(ZWD-US(14))+RFS(15)/(ZWD-US(15))+
     *RFS(16)/(ZWD-US(16))
           ZGFXATSS=-ZGFXATSS
CCCCCCCCCCCCCCCCCCCCCCCCCCCC M13, M24 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ZGFXM13=RFM13(1)/(ZWD-UM13(1))+RFM13(2)/(ZWD-UM13(2))+
     *RFM13(3)/(ZWD-UM13(3))+RFM13(4)/(ZWD-UM13(4))+
     *RFM13(5)/(ZWD-UM13(5))+RFM13(6)/(ZWD-UM13(6))+
     *RFM13(7)/(ZWD-UM13(7))+RFM13(8)/(ZWD-UM13(8))+
     *RFM13(9)/(ZWD-UM13(9))+RFM13(10)/(ZWD-UM13(10))+
     *RFM13(11)/(ZWD-UM13(11))+RFM13(12)/(ZWD-UM13(12))+
     *RFM13(13)/(ZWD-UM13(13))+RFM13(14)/(ZWD-UM13(14))+
     *RFM13(15)/(ZWD-UM13(15))+RFM13(16)/(ZWD-UM13(16))
	 ZGFXM13=-ZGFXM13
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ZGFXM24=RFM24(1)/(ZWD-UM24(1))+RFM24(2)/(ZWD-UM24(2))+
     *RFM24(3)/(ZWD-UM24(3))+RFM24(4)/(ZWD-UM24(4))+
     *RFM24(5)/(ZWD-UM24(5))+RFM24(6)/(ZWD-UM24(6))+
     *RFM24(7)/(ZWD-UM24(7))+RFM24(8)/(ZWD-UM24(8))+
     *RFM24(9)/(ZWD-UM24(9))+RFM24(10)/(ZWD-UM24(10))+
     *RFM24(11)/(ZWD-UM24(11))+RFM24(12)/(ZWD-UM24(12))+
     *RFM24(13)/(ZWD-UM24(13))+RFM24(14)/(ZWD-UM24(14))+ 
     *RFM24(15)/(ZWD-UM24(15))+RFM24(16)/(ZWD-UM24(16))
      ZGFXM24=-ZGFXM24
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC M33 CCCCCCCCCCCCCCCCCCCCCC
      ZGF33=RFS3(1)/(ZWD-US3(1))+RFS3(2)/(ZWD-US3(2))+
     *RFS3(3)/(ZWD-US3(3))+RFS3(4)/(ZWD-US3(4))+RFS3(5)/(ZWD-US3(5))+
     *RFS3(6)/(ZWD-US3(6))+RFS3(7)/(ZWD-US3(7))+RFS3(8)/(ZWD-US3(8))+
     *RFS3(9)/(ZWD-US3(9))+RFS3(10)/(ZWD-US3(10))+
     *RFS3(11)/(ZWD-US3(11))+RFS3(12)/(ZWD-US3(12))+
     *RFS3(13)/(ZWD-US3(13))+RFS3(14)/(ZWD-US3(14))+
     *RFS3(15)/(ZWD-US3(15))+RFS3(16)/(ZWD-US3(16))
           ZGF33=-ZGF33
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC M22 CCCCCCCCCCCCCCCCCCCCCCCCCC
      ZGF22=RFS2(1)/(ZWD-US2(1))+RFS2(2)/(ZWD-US2(2))+
     *RFS2(3)/(ZWD-US2(3))+RFS2(4)/(ZWD-US2(4))+RFS2(5)/(ZWD-US2(5))+
     *RFS2(6)/(ZWD-US2(6))+RFS2(7)/(ZWD-US2(7))+RFS2(8)/(ZWD-US2(8))+
     *RFS2(9)/(ZWD-US2(9))+RFS2(10)/(ZWD-US2(10))+
     *RFS2(11)/(ZWD-US2(11))+RFS2(12)/(ZWD-US2(12))+
     *RFS2(13)/(ZWD-US2(13))+RFS2(14)/(ZWD-US2(14))+
     *RFS2(15)/(ZWD-US2(15))+RFS2(16)/(ZWD-US2(16))
           ZGF22=-ZGF22

      ZGF_AT=ZGFXATSS
C      +2.D0*ZGFXM13+ZGF33

CCCCCCCCCCCCCCCCCCCCCCCCCCCCC ELETRONS DE CONDU��O CCCCCCCCCC
      ZGC11=RC11(1)/(ZWD-UC11(1))+RC11(2)/(ZWD-UC11(2))+
     *RC11(3)/(ZWD-UC11(3))+RC11(4)/(ZWD-UC11(4))+RC11(5)/(ZWD-UC11(5))+
     *RC11(6)/(ZWD-UC11(6))+RC11(7)/(ZWD-UC11(7))+RC11(8)/(ZWD-UC11(8))+
     *RC11(9)/(ZWD-UC11(9))+RC11(10)/(ZWD-UC11(10))+
     *RC11(11)/(ZWD-UC11(11))+RC11(12)/(ZWD-UC11(12))+
     *RC11(13)/(ZWD-UC11(13))+RC11(14)/(ZWD-UC11(14))+
     *RC11(15)/(ZWD-UC11(15))+RC11(16)/(ZWD-UC11(16))

      ZGC22=RC22(1)/(ZWD-UC22(1))+RC22(2)/(ZWD-UC22(2))+
     *RC22(3)/(ZWD-UC22(3))+RC22(4)/(ZWD-UC22(4))+RC22(5)/(ZWD-UC22(5))+
     *RC22(6)/(ZWD-UC22(6))+RC22(7)/(ZWD-UC22(7))+RC22(8)/(ZWD-UC22(8))+
     *RC22(9)/(ZWD-UC22(9))+RC22(10)/(ZWD-UC22(10))+
     *RC22(11)/(ZWD-UC22(11))+RC22(12)/(ZWD-UC22(12))+
     *RC22(13)/(ZWD-UC22(13))+RC22(14)/(ZWD-UC22(14))+
     *RC22(15)/(ZWD-UC22(15))+RC22(16)/(ZWD-UC22(16))

      ZGC33=RC33(1)/(ZWD-UC33(1))+RC33(2)/(ZWD-UC33(2))+
     *RC33(3)/(ZWD-UC33(3))+RC33(4)/(ZWD-UC33(4))+RC33(5)/(ZWD-UC33(5))+
     *RC33(6)/(ZWD-UC33(6))+RC33(7)/(ZWD-UC33(7))+RC33(8)/(ZWD-UC33(8))+
     *RC33(9)/(ZWD-UC33(9))+RC33(10)/(ZWD-UC33(10))+
     *RC33(11)/(ZWD-UC33(11))+RC33(12)/(ZWD-UC33(12))+
     *RC33(13)/(ZWD-UC33(13))+RC33(14)/(ZWD-UC33(14))+
     *RC33(15)/(ZWD-UC33(15))+RC33(16)/(ZWD-UC33(16))

      ZGC44=RC44(1)/(ZWD-UC44(1))+RC44(2)/(ZWD-UC44(2))+
     *RC44(3)/(ZWD-UC44(3))+RC44(4)/(ZWD-UC44(4))+RC44(5)/(ZWD-UC44(5))+
     *RC44(6)/(ZWD-UC44(6))+RC44(7)/(ZWD-UC44(7))+RC44(8)/(ZWD-UC44(8))+
     *RC44(9)/(ZWD-UC44(9))+RC44(10)/(ZWD-UC44(10))+
     *RC44(11)/(ZWD-UC44(11))+RC44(12)/(ZWD-UC44(12))+
     *RC44(13)/(ZWD-UC44(13))+RC44(14)/(ZWD-UC44(14))+
     *RC44(15)/(ZWD-UC44(15))+RC44(16)/(ZWD-UC44(16))

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCC BANDA QUADRADA   CCCCCCCCCCCCCCCCCCCCC

      ZARG=(ZWD-D+AMU)/(ZWD+D+AMU)
      ZM=0.5D0*CDLOG(ZARG)/D
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C FG ATOMICAS UP
      ZGA11=ZGFXATSS
      ZGA13=ZGFXM13
      ZGA31=ZGFXM13
      ZGA33=ZGF33

C      ZWW=VG2*ZM
      ZM0=-1.D0/(ZWD-EQ+AMU)
      ZWW=VG2*ZM0
      ZA=1.D0+(ZGA11+ZGA13)*ZWW
      ZB=(ZGA11+ZGA13)*ZWW
      ZC=(ZGA31+ZGA33)*ZWW
      ZD=1.D0+(ZGA31+ZGA33)*ZWW
      ZMDELTA=ZA*ZD-ZB*ZC
      ZM11=((ZD*ZGA11)-(ZB*ZGA31))/ZMDELTA
      ZM13=((ZD*ZGA13)-(ZB*ZGA33))/ZMDELTA
      ZM31=((ZA*ZGA31)-(ZC*ZGA11))/ZMDELTA
      ZM33=((ZA*ZGA33)-(ZC*ZGA13))/ZMDELTA
      ZWWV=V2*ZM
      ZDELTA1A=1.D0-ZWWV*(ZM11+ZM33+ZM31+ZM13)
      ZWWA=ZWWV*((ZM11*ZM33)-(ZM13*ZM31))
      ZG11=(ZM11-ZWWA)/ZDELTA1A
      ZG13=(ZM13+ZWWA)/ZDELTA1A
      ZG31=(ZM31+ZWWA)/ZDELTA1A
      ZG33=(ZM33-ZWWA)/ZDELTA1A
C O CALCULO DA AREA � FEITO COM O ZGF ABAIXO
      ZGF=-(ZG11+ZG13+ZG31+ZG33)
      ZGC=ZM+V2*ZGF/(D2-ZW2)
      ZGC=ZGC
      ZGFC=V*ZM*ZGF

      ZM=-ZM 
C     PARA ACOPLAMIENTO LATERAL
C      ZG00=ZM*(1.0D0+ZM*V2*ZGF)

C Para acoplamiento inmerso
           
      ZG00=ZM*ZM*V2*ZGF
C      FIM

c     CONSIDERO SISTEMA 2 QDS, 1 INMERSO Y OTRO LATERAL
c      ZG00=ZGC*ZGC*V2*ZGF+((ZGC*ZGC*V2*ZGF)**(2))*V2*ZGF


CCCCCCCCCCCCCCCCCCCCCCCC rede mod. atomico CCCCCCCCCCCCCCCCCCCCCCCC
      C1=0.5D0*V2/D
      ZCADAUXA=(ZM11*ZM33)-(ZM13*ZM31)
      ZCADAUXB=(ZM11+ZM33+ZM31+ZM13)
      ZARGCAD=(ZW-D+V2*ZCADAUXB)/(ZW+D+V2*ZCADAUXB)
      ZGFCAD11=C1*CDLOG(ZARGCAD)*(ZM11*ZCADAUXB-ZCADAUXA)+ZM11
      ZGFCAD13=C1*CDLOG(ZARGCAD)*(ZM13*ZCADAUXB+ZCADAUXA)+ZM13
      ZGFCAD31=C1*CDLOG(ZARGCAD)*(ZM31*ZCADAUXB+ZCADAUXA)+ZM31
      ZGFCAD33=C1*CDLOG(ZARGCAD)*(ZM33*ZCADAUXB-ZCADAUXA)+ZM33
      ZGFREDE=ZGFCAD11+ZGFCAD13+ZGFCAD31+ZGFCAD33

      ZG11=-ZG11
      ZG13=-ZG13
      ZG31=-ZG31
      ZG33=-ZG33

C FG ATOMICAS DOWN CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ZGAD22=ZGF22
      ZGAD24=ZGFXM24
      ZGAD42=ZGFXM24
      ZGAD44=ZGFXATDD

C      ZWW=VG2*ZM
      ZM0=-1.D0/(ZWD-EQ+AMU)
      ZWW=VG2*ZM0
      ZAD=1.D0+(ZGAD22-ZGAD24)*ZWW
      ZBD=(ZGAD24-ZGAD22)*ZWW
      ZCD=(ZGAD42-ZGAD44)*ZWW
      ZDD=1.D0+(ZGAD44-ZGAD42)*ZWW
      ZMDDELTA=ZAD*ZDD-ZBD*ZCD
      ZM22=(ZDD*ZGAD22-ZBD*ZGAD42)/ZMDDELTA
      ZM24=(ZDD*ZGAD24-ZBD*ZGAD44)/ZMDDELTA
      ZM42=(ZAD*ZGAD42-ZCD*ZGAD22)/ZMDDELTA
      ZM44=(ZAD*ZGAD44-ZCD*ZGAD24)/ZMDDELTA
      ZWWV=V2*ZM
      ZDELTA1AD=1.D0-ZWWV*(ZM22+ZM44-ZM24-ZM42)
      ZWWAD=ZWWV*(ZM22*ZM44-ZM24*ZM42)
      ZG22=(ZM22-ZWWAD)/ZDELTA1AD
      ZG24=(ZM24-ZWWAD)/ZDELTA1AD
      ZG42=(ZM42-ZWWAD)/ZDELTA1AD
      ZG44=(ZM44-ZWWAD)/ZDELTA1AD
      ZGFB=-(ZG22+ZG44)
      ZG22=-ZG22
      ZG44=-ZG44
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C APROXIMACAO DE CADEIAS PARA U FINITO
      G0F=-(ANVACATO+ANSSATO)/(ZW-EF)
      G0FU=-(ANSSATO+ANDDATO)/(ZW-EF-COU)
      ZDEN=1.D0-V2*ZM*(G0F+G0FU)
      ZGFF_CAD=-G0F*(1.D0-V2*ZM*G0FU)/ZDEN
      ZGFD_CAD=-V2*G0F*G0FU*ZM/ZDEN
      ZGDF_CAD=ZGFC_CAD
      ZGDD_CAD=-G0FU*(1.D0-V2*ZM*G0F)/ZDEN
      ZGFCAD=ZGFF_CAD+ZGFD_CAD+ZGDF_CAD+ZGDD_CAD
      ZGFSS=ZGFF_CAD
      ZGFDD=ZGFU_CAD
CCCCCCCCCCCCCCCCCCCCCCCC rede aproximacao da cadeia CCCCCCCCCCCCCCC
      DSS=ANSSATO+ANVACATO
      DSD=ANSSATO+ANDDATO
      AQ1=DSS*V2/(ZW-EF)
      AQ2=DSD*V2/(ZW-EF-COU)
      AQ=AQ1+AQ2
      ZARGREDECAD=(ZW-D-AQ)/(ZW+D-AQ)
C
      ZREDECNUM1A=-DSS*DSD*V2*CDLOG(ZARGREDECAD)
      ZREDECDEN1A=2.D0*D*(ZW-EF)*(ZW-EF-COU)
      ZREDE1A=ZREDECNUM1A/ZREDECDEN1A
      ZREDECNUM1B=DSS*(AQ*CDLOG(ZARGREDECAD)-2.D0*D)
      ZREDECDEN1B=2.D0*D*(ZW-EF)
      ZREDE1B=ZREDECNUM1B/ZREDECDEN1B
      ZREDE1=ZREDE1B+ZREDE1A
C
      ZREDECNUM2B=DSS*DSD*V2*CDLOG(ZARGREDECAD)
      ZREDECDEN2B=2.D0*D*(ZW-EF)*(ZW-EF-COU)
      ZREDE2B=ZREDECNUM2B/ZREDECDEN2B
      ZREDE2=ZREDE2B
      ZREDE3=ZREDE2
C
      ZREDECNUM4A=-DSS*DSD*V2*CDLOG(ZARGREDECAD)
      ZREDECDEN4A=2.D0*D*(ZW-EF)*(ZW-EF-COU)
      ZREDE4A=ZREDECNUM4A/ZREDECDEN4A
      ZREDECNUM4B=DSD*AQ*CDLOG(ZARGREDECAD)-DSD*2.D0*D
      ZREDECDEN4B=2.D0*D*(ZW-EF-COU)
      ZREDE4B=ZREDECNUM4B/ZREDECDEN4B
      ZREDE4=ZREDE4B+ZREDE4A
CC
      ZTOS=ZREDE1+ZREDE2+ZREDE3+ZREDE4


C CALCULO PELA FG DOS ELETRONS F DE ACORDO COM O TRABALHO DO BWSP13-SP
C          ZATO=ZM11+ZM33+ZM13+ZM31
C          ZATO=(ZWD-EQ+AMU)*(ZGFT)/((ZWD-EQ+AMU)-VG2*(ZGFT))
C          Z1NUM=ZATO
C          Z1DEN=1.D0-V2*ZATO*ZM
C          ZGFT=-Z1NUM/Z1DEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          ZGC=ZM/V2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        RETURN
        END
        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        SUBROUTINE PARAMETRO
****************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        REAL*8 MU
        COMMON/DADOS/EF,EQ,MU,V,V2,BETA,DS,CC,D,COU
        COMMON/RETAN/E00,EFF,E1,E2,DEL,DEL0,ZUI,EPSIL
        COMMON/ET/ETA,OI,OF,NT
        COMMON/DELTA/DELTA
        DATA PI /3.14159265358979323D0/
        CC=1.D0/PI
        ZUI=(0.D0,1.D0)
C       OI=-0.01D0
C       OF=0.01D0

        OI=-0.15D0
        OF=0.25D0
        D=1.0D0
        DELTA=0.01D0
        NT=5001
        WRITE(6,432)
432     FORMAT(' OI=-0.002D0',5X,' OF=0.002D0',5X,' D=1.D0',
     *5X, ' NT=1001')
        READ(*,*)OI,OF,D,NT
        EQ=0.D0
        MU=0.D0
        V=DSQRT(2.D0*D*DELTA/PI)
        V=0.08D0
        EF=-0.07D0
cccccccccccccccccccccccccccccccccccccccccccccccccc

        COU=0.20D0

cccccccccccccccccccccccccccccccccccccccccccccccccc
        T=1.D-4
        ETA=1.D-8
        WRITE(6,433)
433     FORMAT('EF=-0.08D0',5X,'EQ=0.D0',5X,'COU=0.210D0'/,
     *' MU=0.D0',5X,' T=1.D-4',5X,'ETA=1.D-4')
        READ(*,*)EF,EQ,COU,MU,T,ETA
        BETA=1.D0/T
        V2=V*V
        E00=-10.D0
        EFF=10.D0
        DEL=5.D0
        E1=-0.1
        E2=0.1
C        DEL0=PI/(1.04D0*BETA)
        DEL0=1.D-5
        EPSIL=1.D-8
        WRITE(6,434)
434     FORMAT(' E00=-10.D0',5X,' EFF=10.D0',5X,' E1=-0.1D0',5X,' E2
     *0.1D0',/,' DEL=5.D0',5X,' EPSIL=1.D-8')
        READ(*,*)E00,EFF,E1,E2,DEL,EPSIL
        WRITE(3,14)
        WRITE(6,14)
        WRITE(3,16)ETA,E00,EFF,EPSIL,E1,E2,DEL,DEL0
        WRITE(6,16)ETA,E00,EFF,EPSIL,E1,E2,DEL,DEL0
        WRITE(3,15)OI,OF,NT,D,EQ,ENA
        WRITE(6,15)OI,OF,NT,D,EQ,ENA
        WRITE(3,25)MU,EF,V,T,COU
        WRITE(6,25)MU,EF,V,T,COU
16      FORMAT(' ETA=',F15.8,5X,' E00=',F15.8,/,
     *' EFF=',F15.8,5X,' EPSIL=',F15.8,5X,' E1=',F15.8,/,
     *' E2=',F15.8,5X,' DEL=',F15.8,5X,' DEL0=',F15.8,/)
14      FORMAT(' OS PARAMETROS DO PROBLEMA SAO:')
15      FORMAT(' OI=',F10.5,5X,' OF=',F10.5,' NT=',I5,5X,/,
     *' D=',F10.5,5X,' EQ=',F10.5,5X,' ENA=',F10.5)
25      FORMAT(' MU=',F10.5,5X,5X,' EF=',F10.5,/,
     *' V=',F10.5,5X,' T=',E10.5,5X,'COU=',F10.5)
         RETURN
         END
************++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
********************************* FIM DA SUB PARAMETRO****************************
          
      SUBROUTINE OCUP3
*******************************
      IMPLICIT REAL*8 (A-H,O-W)
      IMPLICIT COMPLEX*16 (X-Z)
      REAL*8 MU
      INTEGER DT30,DT31
      PARAMETER (DT30=30,DT31=31)
      EXTERNAL GFF
      EXTERNAL GCC
      EXTERNAL GCAD
      EXTERNAL G11
      EXTERNAL G33
      EXTERNAL G22
      EXTERNAL G44
      EXTERNAL G13
      EXTERNAL G31
      EXTERNAL GVA
      DIMENSION VA(7),VB(7)
      COMMON/DADOS/EF,EQ,MU,V,V2,BETA,DS,CC,D,COU
      COMMON/RETAN/E00,EFF,E1,E2,DEL,DEL0,ZUI,EPSIL
      COMMON/ET/ETA,OI,OF,NT
      COMMON/OCC/TXFA,TDS,TXFC,TXC,TOT,TSOMA,FRIED,RF,RC
      COMMON/OCC2/TXFF,TXCC,TXUP,TXDOWN,TXDD,TXVAC,TSOMAF,TSOMA2
      COMMON/OCC3/TX11,TX22,TX33,TX44,TX13,TX31,FSUM,TXCAD
      COMMON/OCUP_EXAT/ANDD,ANSS,ANDDS,AVNSS,SOMA,ANDM13,ANDM24
      COMMON/CORREC/DEK
      OPEN (UNIT=DT30,FILE='ocupacao_qsuba.dat',STATUS='UNKNOWN')
      OPEN (UNIT=DT31,FILE='integrando.dat',STATUS='UNKNOWN')
      PI=1.0D0/CC
      AMU=MU
      T=1.D0/BETA
      A=-3.2D0
      B=3.2D0
      SFF=0.D0
      SCC=0.D0
      SCAD=0.D0
      S11=0.D0
      S33=0.D0
      S22=0.D0
      S44=0.D0
      S13=0.D0
      S31=0.D0
      SVA=0.D0
      SDS=0.D0

       DD1=0.01D0
       DD2=0.001D0
       DD3=0.01D0

      VA(1)=-3.2D0
      VB(1)=-0.15D0
      VA(2)=-0.15D0
      VB(2)=-0.02d0
      VA(3)=-0.02D0
      VB(3)=0.02D0
      VA(4)=0.02D0
      VB(4)=0.20D0

      VA(5)=0.20D0
      VB(5)=3.0D0
      EPSIL=1.0d-10

          DO 10 I=1,5
          A=VA(I)
          B=VB(I)
        RES1=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,GFF)
        RES10=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,GCC)
        RES2=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,G11)
        RES3=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,G33)
	RES4=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,G22)
	RES5=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,G13)
	RES6=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,GVA)
        RES7=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,G44)
        RES8=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,G31)
        RES9=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,GCAD)
C        SDS=SDS+RES1
        SFF=SFF+RES1
        SCC=SCC+RES10
        S11=S11+RES2
        S33=S33+RES3
	S22=S22+RES4
	S44=S44+RES7
	S13=S13+RES5
	S31=S31+RES8
        SVA=SVA+RES6
        SCAD=SCAD+RES9
C	  WRITE(6,*)'A,B,SSS=',A,B,RES1
10        CONTINUE
c	  ZZ=DCMPLX(0.0D0,ETA)

        TXFF=CC*SFF
        TXCC=-CC*SCC
	TX11=CC*S11
	TX33=CC*S33
	TX22=CC*S22
	TX13=CC*S13
	TX31=CC*S31
	TX44=CC*S44
	TXVAC=CC*SVA
	TXCAD=CC*SCAD
	DELL=0.5D0*PI*V2/D
	DELLI=1.D0/(PI*DELL)
	FRIED=DELLI*SIN(PI*TXFF)*SIN(PI*TXFF)
C	 IF(COU.GT.1.1D0)TXSD=ANDDS
C      IF(TXSS.LT.0.5D0)GO TO 18
C      TXDD=0.5D0*TXDD
C      TXSS=0.5D0*TXSS
C      TXVAC=0.5D0*TXVAC
C      TXSD=0.5D0*TXSD
C18    CONTINUE
C      TSOMA2=TXF+TX22+TXVAC
c      TSOMAF=TXDS
c      TXDD=TX44
c      TXUP=TX11
c      TXDOWN=TX22
c      TSOMA2=TXUP+TXDOWN+TXDD+TXVAC
c      DELL=0.5D0*PI*V2/D
c      DELLI=1.D0/(PI*DELL)
C CALCULAMOS TX13 E TX31 E SUA SOMA SE ANULA - 060508
      ZY=DCMPLX(AMU,ETA)
      CALL GKONDO(ZY,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
     *ZGFXATSS,ZG00,ZM)
      RF=-CC*DIMAG(ZGF)
c      RF=DABS(RF)
      RC=-CC*DIMAG(ZGC)
c      RC=DABS(RC)
      WRITE(DT30,31)T,TXFF,TXCC,TXUP,TXDOWN,TXDD,TXVAC,TX33,TSOMA2,
     *FRIED,RF,TXCAD
C          WRITE(DT30,31)T,MU,DEK,TSOMA,TR2,TXC,TXF,TDS,TXFC,TOT,
C     *R0F,R0C,FRIED
C          WRITE(6,30)DEK,TXF,TX22,TXVAC,TSOMA2
C		TXF,TDS,TXFC,TSOMA,TOT,DEK,R0F,R0C,FRIED
C30        FORMAT('DEK=',F10.8,1X,'TXF=',F10.8,1X,'TX11B=',F10.8,
C     *1X,'TXVAC=',F10.8,1X,'TSOMA2=',F10.8)
C     ,5X,*'TXC=',F15.8,/,'TXF=',F15.8,
C     *'TDS=',F15.8,5X,' NFC=',F15.8,/,'TSOMA=',F15.8,5X,'NTOT=',
C     *F15.8,5X,'DEK=',F15.8,
C     */,'R0F=',F15.8,5X,'R0C=',F15.8,5X,'FRIEDEL=',F15.8)
31        FORMAT(20F15.8)
          RETURN
          END

        FUNCTION GFF(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
	TT=1.D0/BETA
C        ETTA=1.5d0*TT
	ETTA=1.D-4
        ZW=DCMPLX(EW,ETTA)
        CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
     *ZGFXATSS,ZG00,ZM)
************************************************
        ZFE=XFER(ZW)
        GFF=-DIMAG(ZFE*ZGF)
        RETURN
        END
        
        
         FUNCTION GCC(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
	TT=1.D0/BETA
C        ETTA=1.5d0*TT
	ETTA=1.D-4
        ZW=DCMPLX(EW,ETTA)
        CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
     *ZGFXATSS,ZG00,ZM)
************************************************
        ZFE=XFER(ZW)
        GCC=-DIMAG(ZFE*ZGC)
        RETURN
        END
        
        FUNCTION GCAD(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
	TT=1.D0/BETA
C        ETTA=1.5d0*TT
	ETTA=1.D-4
        ZW=DCMPLX(EW,ETTA)
        CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
     *ZGFXATSS,ZG00,ZM)
************************************************
        ZFE=XFER(ZW)
        GCAD=-DIMAG(ZTOS)
        RETURN
        END

        FUNCTION G11(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
	TT=1.D0/BETA
C        ETTA=1.5d0*TT
	ETTA=1.D-4
        ZW=DCMPLX(EW,ETTA)
        CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
     *ZGFXATSS,ZG00,ZM)
************************************************
        ZFE=XFER(ZW)
        G11=-DIMAG(ZFE*ZG11)
C        G11=-DIMAG((1.D0-ZFE)*ZG44)
        RETURN
        END


        FUNCTION G33(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
	TT=1.D0/BETA
C        ETTA=1.5d0*TT
	ETTA=1.D-4
        ZW=DCMPLX(EW,ETTA)
        CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
     *ZGFXATSS,ZG00,ZM)
************************************************
        ZFE=XFER(ZW)
        G33=-DIMAG(ZFE*ZG33)
C        G33=-DIMAG((1.D0-ZFE)*ZG44)
        RETURN
        END


        FUNCTION G22(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
	TT=1.D0/BETA
C        ETTA=1.5d0*TT
	ETTA=1.D-4
        ZW=DCMPLX(EW,ETTA)
        CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
     *ZGFXATSS,ZG00,ZM)
************************************************
        ZFE=1.D0-XFER(ZW)
C        G22=-DIMAG(ZFE*ZG22)
       G22=-DIMAG(ZFE*ZG33)
        RETURN
        END
        
        FUNCTION G44(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
	TT=1.D0/BETA
C        ETTA=1.5d0*TT
	ETTA=1.D-4
        ZW=DCMPLX(EW,ETTA)
        CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
     *ZGFXATSS,ZG00,ZM)
************************************************
        ZFE=XFER(ZW)
        G44=-DIMAG(ZFE*ZG33)
        RETURN
        END

        FUNCTION G13(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
	TT=1.D0/BETA
C        ETTA=1.5d0*TT
	ETTA=1.D-4
        ZW=DCMPLX(EW,ETTA)
        CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
     *ZGFXATSS,ZG00,ZM)
************************************************
        ZFE=XFER(ZW)
        G13=-DIMAG(ZFE*ZG13)
C        G13=-DIMAG(ZG13)
        RETURN
        END
        
        FUNCTION G31(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
	TT=1.D0/BETA
C        ETTA=1.5d0*TT
	ETTA=1.D-4
        ZW=DCMPLX(EW,ETTA)
        CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
     *ZGFXATSS,ZG00,ZM)
************************************************
        ZFE=XFER(ZW)
        G31=-DIMAG((1.D0-ZFE)*ZG13)
C        G31=-DIMAG(ZG13)
        RETURN
        END


        FUNCTION GVA(EW)
*****************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
	TT=1.D0/BETA
C        ETTA=1.5d0*TT
	ETTA=1.D-4
        ZW=DCMPLX(EW,ETTA)
        CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
     *ZGFXATSS,ZG00,ZM)
************************************************
        ZFE=1.D0-XFER(ZW)
        GVA=-DIMAG(ZFE*ZG11)
        RETURN
        END

***********************************************************************
        SUBROUTINE GRAF(OI,OF,NTT)
*************************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)

        EXTERNAL FERM
        EXTERNAL ZGLOCAL
        EXTERNAL PIMGCOEF1
        EXTERNAL ZGCON

        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
        COMMON/OCC/TXFA,TDS,TXFC,TXC,TOT,TSOMA,FRIED,RF,RC
        COMMON/QD2/EF2,VQ2,EL2,TXF2,U2
        COMMON/QD1/EF1,VQ1,EL1,TXF1,U1

        INTEGER DT81
        PARAMETER (DT81=81)
      OPEN (UNIT=DT81,FILE='2QDS-U,T=5.5(21.3).dat',STATUS='UNKNOWN')
C ESTA SUBROTINA CALCULA O GRAFICO DA DENSIDADE DE ESTADOS
C NO EIXO REAL.
        ETA=1.0d-4
        DM1=(OF-OI)/(NTT-1)
        GAM=4.0D0*CC*CC*D*D*VQ2*VQ2/(VQ1*VQ1)

        DO 80 K=1,NTT
        BW=OI+(K-1)*DM1
        ZW=DCMPLX(BW,ETA)

        TRA=DIMAG(ZGLOCAL(BW))*DIMAG(ZGLOCAL(BW))+DREAL(ZGLOCAL(BW))*DREAL(ZGLOCAL(BW))
        TRA=TRA*GAM

        R0F1=-CC*DIMAG(ZGLOCAL(BW))
        R0F1=DABS(R0F1)
        R0F2=CC*CC*0.25D0*R0F1*FERM(BW)
        R0F3=CC*CC*0.25D0*(1.0D0-FERM(BW))*R0F1
        R0F4=-CC*DIMAG(ZGCON(BW))
        R0F4=DABS(R0F4)
        R0F5=CC*CC*0.25D0*R0F4*FERM(BW)
        R0F6=CC*CC*0.25D0*(1.0D0-FERM(BW))*R0F4
        ROF7=TRA


c        CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
c     *ZGFXATSS,ZG00,ZM)
c	R0F1=-CC*DIMAG(ZGF)
c	R0F2=-CC*DIMAG(ZGC)
c	R0F3=-CC*DIMAG(ZGFC)
c	R0F4=-CC*DIMAG(ZG00)

*************************************
        WRITE(DT81,5800)BW*1.0d2,R0F1,R0F2,R0F3,R0F4,R0F5,R0F6,ROF7
80      CONTINUE
5800    FORMAT(8F15.8)
        CLOSE(DT81)
        RETURN

        END

C**********************************************************************
      SUBROUTINE SIMP(N,H,F,RES)
C  ********************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(N)
      M=N-3
      SP=F(2)
      SI=F(3)
      DO 1 I=4,M,2
      SP=SP+F(I)
    1 SI=SI+F(I+1)
      RES=H*(F(1)+2.D0*(SI+2.D0*(SP+F(N-1)))+F(N))/3.D0
      RETURN
      END


******************************************************************************************
         SUBROUTINE THERMO(G20,S0,AK0,TEM0,WFR0,G2,S,AK,TEM,WFR,E)
**********************************************************
         IMPLICIT COMPLEX*16 (X-Z)
         IMPLICIT REAL*8 (A-H,O-W)
         COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
             
              CALL COEF11(AL11,AL11QSUBA)
              G2=AL11QSUBA
              G20=AL11
           
              CALL COEF12(AL12,AL12QSUBA)
******************************
              S=-AL12QSUBA*BETA/AL11QSUBA
              S0=-AL12*BETA/AL11
           
              CALL COEF22(AL22,AL22QSUBA)
******************************

              AK=BETA*(AL22QSUBA-((AL12QSUBA*AL12QSUBA)/AL11QSUBA))

              AK0=BETA*(AL22-((AL12*AL12)/AL11))

              TEM=S*S*G2/(BETA*AK)
              TEM0=S0*S0*G20/(BETA*AK0)

              WFR=3.0D0*BETA*AK*CC*CC/G2
              WFR0=3.0d0*BETA*AK0*CC*CC/G20
              E=AL12QSUBA*AL12QSUBA/(AL11QSUBA*AL22QSUBA)

             
        RETURN
        END

***************************************************************

        FUNCTION AELE11(EW)
*************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)

        EXTERNAL FERM
        EXTERNAL ZGLOCAL
        EXTERNAL PIMGCOEF1
        EXTERNAL ZGCON

        COMMON/ET/ETA,OI,OF,NT
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
        COMMON/RF22/DELTA22
        COMMON/QD2/EF2,VQ2,EL2,TXF2,U2
        COMMON/QD1/EF1,VQ1,EL1,TXF1,U1
 
        DATA PI /3.14159265358979323D0/
        ETA=1.0D-4   
c         DX=0.0D0
c         DY=0.0D0

c         ZWE=DCMPLX(EW,ETA)
C         CALL GKONDO(ZWE,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
C     *   ZGFXATSS,ZG00,ZM)
c         CALL zig_mauro(ZWE,DX,DY,ZG0,ZG0N)
c              AT=DIMAG(ZG0)*DIMAG(ZG0)+DREAL(ZG0)*DREAL(ZG0)
c        AT=(1.0d0/DIMAG(ZG0))
c        AT=CC*CC*2.0d0*D*2.0d0*D
  
c        ZAT=ZM*ZM*V2*ZGF
c         ZN=ZGCON(AMU)
c         AT2=DIMAG(ZN)*DIMAG(ZN)+DREAL(ZN)*DREAL(ZN)
c         GAM=1.0d0/AT2
c        GAM=1.0D0/DABS(DIMAG(ZN))
c        GAM=VQ2*VQ2/DELTA22
c        GAM=GAM*GAM
         GAM=4.0D0*CC*CC*D*D*VQ2*VQ2/(VQ1*VQ1)

        ZAT=ZGLOCAL(EW)
        AT=DIMAG(ZAT)*DIMAG(ZAT)+DREAL(ZAT)*DREAL(ZAT)
        AT=AT*GAM

c        IF (AT.LE.(1.0D0)) THEN
c        AT=AT
c        ELSE
c        DAT=AT-1.0D0
c        AT=AT-1.0D0*DAT
c        ENDIF
c        CONTINUE 

        FER=FERM(EW)
	DFER=FER*(1.D0-FER)
	
          AELE11=2.0d0*DFER*AT*BETA
	  RETURN
	  END

        FUNCTION AELE12(EW)
*************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)

        EXTERNAL FERM
        EXTERNAL ZGLOCAL
        EXTERNAL PIMGCOEF1
        EXTERNAL ZGCON        

        COMMON/ET/ETA,OI,OF,NT
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
        COMMON/RF22/DELTA22
        COMMON/QD2/EF2,VQ2,EL2,TXF2,U2
        COMMON/QD1/EF1,VQ1,EL1,TXF1,U1
       
        DATA PI /3.14159265358979323D0/
        ETA=1.0D-4
c        DX=0.0D0
c         DY=0.0D0
c         ZWE=DCMPLX(EW,ETA)
C      CALL GKONDO(ZWE,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
C     *   ZGFXATSS,ZG00,ZM)
c        ZN=ZGCON(AMU)
c        GAM=DREAL(ZN)*DREAL(ZN)+DIMAG(ZN)*DIMAG(ZN)
c        GAM=1.0d0/GAM
         GAM=4.0D0*CC*CC*D*D*VQ2*VQ2/(VQ1*VQ1)
       
c        GAM=1.0D0/DABS(DIMAG(ZN))
           
c         GAM=VQ2*VQ2/DELTA22
c         GAM=GAM*GAM

         ZAT=ZGLOCAL(EW)

        AT=DIMAG(ZAT)*DIMAG(ZAT)+DREAL(ZAT)*DREAL(ZAT)
        AT=AT*GAM

c        IF (AT.LE.(1.0D0)) THEN
c        AT=AT
c        ELSE
c        DAT=AT-1.0D0
c        AT=AT-1.0D0*DAT
c        ENDIF
c        CONTINUE 

c         CALL zig_mauro(ZWE,DX,DY,ZG0,ZG0N)
c              AT=DIMAG(ZG0)*DIMAG(ZG0)+DREAL(ZG0)*DREAL(ZG0)
c        AT=(1.0d0/DIMAG(ZG0))
c        AT=AT*AT
c        AT=CC*CC*2.0d0*D*2.0d0*D     
        FER=FERM(EW)
	DFER=FER*(1.D0-FER)
	
********************************************************
           AELE12=2.0d0*DFER*(EW-AMU)*AT*BETA

	  RETURN
	  END


        FUNCTION AELE22(EW)
*************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)

        EXTERNAL FERM
        EXTERNAL ZGLOCAL
        EXTERNAL PIMGCOEF1
        EXTERNAL ZGCON

        COMMON/ET/ETA,OI,OF,NT
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
        COMMON/RF22/DELTA22
        COMMON/QD2/EF2,VQ2,EL2,TXF2,U2
        COMMON/QD1/EF1,VQ1,EL1,TXF1,U1
        
        DATA PI /3.14159265358979323D0/
        ETA=1.0D-4
c         DX=0.0D0
c         DY=0.0D0
c        ZWE=DCMPLX(EW,ETA)
C        CALL GKONDO(ZWE,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
C     *   ZGFXATSS,ZG00,ZM)
c        ZN=ZGCON(AMU)
c        GAM=DIMAG(ZN)*DIMAG(ZN)+DREAL(ZN)*DREAL(ZN)
c        GAM=1.0D-1/GAM
         GAM=4.0D0*CC*CC*D*D*VQ2*VQ2/(VQ1*VQ1)

c        GAM=1.0D0/DABS(DIMAG(ZN))
c         GAM=VQ2*VQ2/DELTA22
c         GAM=GAM*GAM
         
         ZAT=ZGLOCAL(EW)

        AT=DIMAG(ZAT)*DIMAG(ZAT)+DREAL(ZAT)*DREAL(ZAT)  
        AT=AT*GAM
        
c        IF (AT.LE.(1.0D0)) THEN
c        AT=AT
c        ELSE
c        DAT=AT-1.0D0
c        AT=AT-1.0D0*DAT
c        ENDIF
c        CONTINUE 

c         CALL zig_mauro(ZWE,DX,DY,ZG0,ZG0N)
c              AT=DIMAG(ZG0)*DIMAG(ZG0)+DREAL(ZG0)*DREAL(ZG0)
c        AT=(1.0d0/DIMAG(ZG0))
c        AT=AT*AT
c        AT=CC*CC*2.0d0*D*2.0d0*D      
        FER=FERM(EW)
	DFER=FER*(1.D0-FER)

          AELE22=2.0d0*DFER*(EW-AMU)*(EW-AMU)*AT*BETA

	  RETURN
	  END


       SUBROUTINE COEF11(AL11,AL11QSUBA)
*********************************************
         IMPLICIT COMPLEX*16 (X-Z)
         IMPLICIT REAL*8 (A-H,O-W)
         INTEGER DT75
         PARAMETER (DT75=75)
         DIMENSION FF11(3001)
         
        EXTERNAL FERM
        EXTERNAL ZGLOCAL
        EXTERNAL PIMGCOEF1
        EXTERNAL ZGCON
        EXTERNAL AELE11

         COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
         COMMON/OCC/TXFA,TDS,TXFC,TXC,TOT,TSOMA,FRIED,RF,RC
         COMMON/RF22/DELTA22
         COMMON/QD2/EF2,VQ2,EL2,TXF2,U2
         COMMON/QD1/EF1,VQ1,EL1,TXF1,U1
         
      OPEN (UNIT=DT75,FILE='AL11.dat',STATUS='UNKNOWN')
         DATA PI /3.14159265358979323D0/

               DELTA=0.01D0
               D=1.0D2*DELTA
               ETA=1.0D-4
               EPSIL=1.0D-10
               CC=1.0d0/PI
               TT=1.D0/BETA
c               ZN=ZGCON(AMU)
               GAM=4.0D0*CC*CC*D*D*VQ2*VQ2/(VQ1*VQ1)
c               GAM=DREAL(ZN)*DREAL(ZN)+DIMAG(ZN)*DIMAG(ZN)
c               GAM=1.0D0/GAM
c               GAM=1.0D0/DABS(DIMAG(ZN))
c               GAM=VQ2*VQ2/DELTA
c               GAM=GAM*GAM

C Esto para Coeficiente de transmisión
c                AI=-D
c                AF=D
               AI=-15.0D0*TT+AMU
               AF=15.0D0*TT+AMU
              NM=3001
C              BETA2=BETA*BETA
              AL120=0.D0
              RES12Q=0.0d0
              RES12=0.0d0
              DO 70 KK=1,NM
              FF11(KK)=0.D0
70             CONTINUE
              CM1=(AF-AI)/(NM-1)
              DO 50 I=1,NM
              EW=AI+(I-1)*CM1
 
C              ZW=DCMPLX(EW,ETA) 
c              FER=1.0d0/(1.0d0+DEXP((EW-AMU)*BETA))

	      FER=FERM(EW)
              DFER=FER*(1.D0-FER)
                           
C              CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
C     *ZGFXATSS,ZG00,ZM)
              ZAT=ZGLOCAL(EW)
*************************************************
              ATRANS=DREAL(ZAT)*DREAL(ZAT)+DIMAG(ZAT)*DIMAG(ZAT)
              ATRANS=ATRANS*GAM 

c             IF (ATRANS.LE.(1.0D0)) THEN
c             ATRANS=ATRANS
c             ELSE
c             DAT=ATRANS-1.0D0
c             ATRANS=ATRANS-1.0D0*DAT
c             ENDIF
c             CONTINUE  
              
              FF1=2.0d0*DFER*ATRANS*BETA
            
              FF11(I)=FF1
              CALL SIMP(NM,CM1,FF11,RES12)
*****************************************
              AL120=AL120+RES12
18            FORMAT(10F15.8)
              WRITE(DT75,18)EW,ATRANS,FF1
50            CONTINUE
C MULTIPLIQUEI POR BETA O AL12
              AL11=AL120
              RES12Q=QSUBA(AI,AF,EPSIL,NPTS,ICHECK,RELERR,AELE11)
              AL11QSUBA=RES12Q
              CLOSE(DT75)
        RETURN
        END 

         SUBROUTINE COEF12(AL12,AL12QSUBA)
*********************************************
         IMPLICIT COMPLEX*16 (X-Z)
         IMPLICIT REAL*8 (A-H,O-W)
         INTEGER DT75
         PARAMETER (DT75=75)
         DIMENSION FF12(3001)
         
         EXTERNAL FERM
         EXTERNAL ZGLOCAL
         EXTERNAL PIMGCOEF1
         EXTERNAL ZGCON
         EXTERNAL AELE12

         COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
         COMMON/OCC/TXFA,TDS,TXFC,TXC,TOT,TSOMA,FRIED,RF,RC
         COMMON/RF22/DELTA22
         COMMON/QD2/EF2,VQ2,EL2,TXF2,U2
         COMMON/QD1/EF1,VQ1,EL1,TXF1,U1

         OPEN (UNIT=DT75,FILE='AL12.dat',STATUS='UNKNOWN')
           DATA PI /3.14159265358979323D0/
               DELTA=0.01D0
               D=1.0D2*DELTA
               ETA=1.0D-4
               EPSIL=1.0D-10
c               ZN=ZGCON(AMU)
               GAM=4.0D0*CC*CC*D*D*VQ2*VQ2/(VQ1*VQ1)
c               GAM=DREAL(ZN)*DREAL(ZN)+DIMAG(ZN)*DIMAG(ZN)
c               GAM=1.0D0/GAM
c               GAM=1.0D0/DABS(DIMAG(ZN))
c               GAM=VQ2*VQ2/DELTA22
c               GAM=GAM*GAM
                       
               CC=1.0d0/PI
               TT=1.D0/BETA
                AI=-15.0D0*TT+AMU
               AF=15.0D0*TT+AMU
              NM=3001
C              BETA2=BETA*BETA
              AL120=0.D0
              RES12=0.0D0

              RES=0.0d0
              DO 70 KK=1,NM
              FF12(KK)=0.D0
70             CONTINUE
              CM1=(AF-AI)/(NM-1)
              DO 50 I=1,NM
              EW=AI+(I-1)*CM1
 
C              ZW=DCMPLX(EW,ETA) 
c              FER=1.0d0/(1.0d0+DEXP((EW-AMU)*BETA))

	      FER=FERM(EW)
              DFER=FER*(1.D0-FER)
              
C              CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
C     *ZGFXATSS,ZG00,ZM)
              ZAT=ZGLOCAL(EW)
*************************************************

            ATRANS=DIMAG(ZAT)*DIMAG(ZAT)+DREAL(ZAT)*DREAL(ZAT)
            ATRANS=ATRANS*GAM

c             IF (ATRANS.LE.(1.0D0)) THEN
c             ATRANS=ATRANS
c             ELSE
c             DAT=ATRANS-1.0D0
c             ATRANS=ATRANS-1.0D0*DAT
c             ENDIF
c             CONTINUE 

              FF12(I)=2.0d0*DFER*ATRANS*(EW-AMU)*BETA

              CALL SIMP(NM,CM1,FF12,RES12)
*****************************************
              AL120=AL120+RES12
              WRITE(DT75,18)EW,ATRANS,FF12(I)
50            CONTINUE
C MULTIPLIQUEI POR BETA O AL12
              AL12=AL120
              RES=QSUBA(AI,AF,EPSIL,NPTS,ICHECK,RELERR,AELE12)
              AL12QSUBA=RES
18            FORMAT(10F15.8)
              CLOSE(DT75)
        RETURN
        END

         SUBROUTINE COEF22(AL22,AL22QSUBA)
*********************************************
         IMPLICIT COMPLEX*16 (X-Z)
         IMPLICIT REAL*8 (A-H,O-W)
         INTEGER DT76
         PARAMETER (DT76=76)
         DIMENSION FF22(3001)
          
         EXTERNAL FERM
         EXTERNAL ZGLOCAL
         EXTERNAL PIMGCOEF1
         EXTERNAL ZGCON 
         EXTERNAL AELE22
        
         COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
         COMMON/OCC/TXFA,TDS,TXFC,TXC,TOT,TSOMA,FRIED,RF,RC
         COMMON/RF22/DELTA22
         COMMON/QD2/EF2,VQ2,EL2,TXF2,U2
         COMMON/QD1/EF1,VQ1,EL1,TXF1,U1
         
         OPEN (UNIT=DT76,FILE='AL22.dat',STATUS='UNKNOWN')
         DATA PI /3.14159265358979323D0/

               ETA=1.0D-4
               EPSIL=1.0D-10
               DELTA=0.01D0
               D=1.0D2*DELTA
               CC=1.0d0/PI
               TT=1.D0/BETA
c               ZN=ZGCON(AMU)
               GAM=4.0D0*CC*CC*D*D*VQ2*VQ2/(VQ1*VQ1)
c               GAM=DREAL(ZN)*DREAL(ZN)+DIMAG(ZN)*DIMAG(ZN)
c               GAM=1.0D0/GAM
c               GAM=1.0D0/DABS(DIMAG(ZN))
c               GAM=VQ2*VQ2/DELTA22
c               GAM=GAM*GAM

               AI=-15.0D0*TT+AMU
               AF=15.0D0*TT+AMU
              NM=3001
C              T=1.D0/BETA
C              BETA2=BETA*BETA
              AL220=0.D0
              RES22=0.0D0
              RES=0.0d0
              DO 70 KK=1,NM
              FF22(KK)=0.D0
70             CONTINUE
              CM1=(AF-AI)/(NM-1)
              DO 50 I=1,NM
              EW=AI+(I-1)*CM1
C	      ZW=DCMPLX(EW,ETA)
	      FER=FERM(EW)
c              FER=1.0d0/(1.0d0+DEXP((EW-AMU)*BETA))
	      DFER=FER*(1.D0-FER)
C              CALL GKONDO(ZW,ZGF,ZGC,ZGFC,ZG11,ZG33,ZG22,ZG44,ZG13,ZG31,
C     *ZGFXATSS,ZG00,ZM)              
****************************************************
              ZAT=ZGLOCAL(EW)

            ATRANS=DIMAG(ZAT)*DIMAG(ZAT)+DREAL(ZAT)*DREAL(ZAT)
            ATRANS=ATRANS*GAM

c             IF (ATRANS.LE.(1.0D0)) THEN
c             ATRANS=ATRANS
c             ELSE
c             DAT=ATRANS-1.0D0
c             ATRANS=ATRANS-1.0D0*DAT
c             ENDIF
c             CONTINUE 


              FF22(I)=2.0d0*DFER*(EW-AMU)*(EW-AMU)*ATRANS*BETA

              CALL SIMP(NM,CM1,FF22,RES22)
*****************************************
              AL220=AL220+RES22
              WRITE(DT76,18)EW,ATRANS,FF22(I)
50            CONTINUE
              AL22=AL220
              RES=QSUBA(AI,AF,EPSIL,NPTS,ICHECK,RELERR,AELE22)
              AL22QSUBA=RES
18            FORMAT(10F15.8)
              CLOSE(DT76)
        RETURN
        END
        

C*****************************************************************************************+
         FUNCTION XFER(ZW)
C********************************
          IMPLICIT REAL*8 (A-H,O-W)
          IMPLICIT COMPLEX*16 (X-Z)
          REAL*8 MU
          COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
          DATA BEXP/35.0d0/
          DATA ZUM/(1.D0,0.D0)/
          T=1.0d0/BETA
          Z=ZW-AMU
          ZEX=Z*BETA
          RX=DREAL(Z)
          RY=DIMAG(Z)
          IF(RY.NE.0.D0)GO TO 100
          XFER=DCMPLX(FERM(RX),0.D0)
          RETURN
100       ABSX=DABS(RX)
          IF(ABSX - T*BEXP) 250,200,200
200       IF(RX)210,220,230
210       XFER= (1.0D0,0.D0)
          RETURN
220       XFER=(0.5D0,0.D0)
          RETURN
230       XFER=(0.D0,0.D0)
          RETURN
250       CONTINUE
          XFER=ZUM/(ZUM+CDEXP(ZEX))
          RETURN
          END

C****************************************************************************************
          FUNCTION FERM(X)
C   *****************************
          IMPLICIT REAL*8(A-H,O-Z)
          COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
          T=1.D0/BETA
          BEXP=35.0d0
          ABSX=DABS(X)
          IF(ABSX - T*BEXP)250,200,200
200       IF(X-AMU)210,220,230
210       FERM=1.D0
          RETURN
220       FERM=0.5D0
          RETURN
230       FERM=0.D0
          RETURN
250       ARRG=(X)*BETA
          FERM=1.D0/(1.D0+DEXP(ARRG))
300       RETURN
          END
C********************************************************************************************

C******************************************************************************************
          FUNCTION FMU(DEQ)
C************************************
          IMPLICIT REAL*8 (A-H,O-W)
          IMPLICIT COMPLEX*16 (X-Z)
            
          COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
          COMMON/OCUP_EXAT/ANDD,ANSS,ANDDS,AVNSS,SOMA,ANDM13,ANDM24
          COMMON/OCC/TXFA,TDS,TXFC,TXC,TOT,TSOMA,FRIED,RF,RC
          COMMON/OCC2/TXFF,TXCC,TXUP,TXDOWN,TXDD,TXVAC,TSOMAF,TSOMA2
          COMMON/QD1/EF1,VQ1,EL1,TXF1,U1
          COMMON/QD2/EF2,VQ2,EL2,TXF2,U2

          V=VQ1
          V2=V*V
          EF=EF1
          COU=U1

          EQ=DEQ
          CALL EXATA_UFINITO
          CALL OCUP3
          
c          TXX=TXFF
c          CALL FASE(FASELO,FASECON,FASEQ1,FASEQ2)
c          FASEE=FASECON+5.0d-1
c          FASEE=DABS(FASEE)
c          FASEE=(DABS(FASEE-TXX)/TXX)

c          FRI2=DABS((FRIED/RF)-1.0D0)*1.0D5+FASEE*1.0D2
c          FMU=FRI2+DABS(1.0D0-DEQ*1.0D12)*1.0D2

C      AA=-1.0D0+DABS((FRIED/RF)-1.0D0)*1.0D12
C      AA=AA/(1.0D0+DABS(DEQ-AMU)*1.0d11)
C      BB=-1.0D0+DABS(AMU-DEQ)*1.0D11
C      BB=BB/(1.0D0+DABS((FRIED/RF)-1.0D0)*1.0d12)
      FMU=DABS(EL2-DEQ)*1.0D2+DABS((FRIED/RF)-1.0D0)*1.0d5
C      +DABS(AMU-DEQ)*1.0D8
C      FMU=0.0D0
C      FMU=AA+BB+FMU

          RETURN
          END

          FUNCTION FMU2(DEQ)
C************************************
          IMPLICIT REAL*8 (A-H,O-W)
          IMPLICIT COMPLEX*16 (X-Z)

c          EXTERNAL PIMGCOEF1
           
          COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU
          COMMON/OCUP_EXAT/ANDD,ANSS,ANDDS,AVNSS,SOMA,ANDM13,ANDM24
          COMMON/OCC/TXFA,TDS,TXFC,TXC,TOT,TSOMA,FRIED,RF,RC
          COMMON/OCC2/TXFF,TXCC,TXUP,TXDOWN,TXDD,TXVAC,TSOMAF,TSOMA2
          COMMON/RF22/DELTA22
          COMMON/QD2/EF2,VQ2,EL2,TXF2,U2
          COMMON/QD1/EF1,VQ1,EL1,TXF1,U1


          PI=1.0D0/CC

c          EF=EF1
c          V=VQ1
c          V2=VQ1*VQ1
c          COU=U1
c          EQ=EL1

c         FRIED2=FRIED*PI*1.0D-2
c         FRIED=FRIED2/(PI*DELTA22)

          EF=EF2
          V=VQ2
          V2=VQ2*VQ2
          COU=U2

          EQ=DEQ

          CALL EXATA_UFINITO
          CALL OCUP3
c          TXX=TXFF
c          CALL FASE(FASELO,FASECON,FASEQ1,FASEQ2)
c          FASEE=FASELO
c          FASEE=DABS(FASEE)
c          FASEE=(DABS(FASEE-TXX)/TXX)
C          FRI2=DABS((FRIED/RF)-1.0D0)*1.0D5
C          FRI2=FRI2-(TXF2/2.0D0)
C          FRI2=DABS(FRI2)
C      AA=-1.0D0+DABS((FRIED/RF)-1.0D0)*1.0D12
C      AA=AA/(1.0D0+DABS(DEQ-AMU)*1.0D11)
C      BB=-1.0D0+DABS(AMU-DEQ)*1.0D11
C      BB=BB/(1.0D0+DABS((FRIED/RF)-1.0D0)*1.0D12)
      FMU2=1.0d5*DABS((FRIED/RF)-1.0D0)+DABS(AMU-DEQ)*1.0d2
C      FMU2=AA+BB+FMU2
C      +BB+FMU2
C       FMU2=FMU2+FASEE*1.0D2
C     *    +DABS(AMU-DEQ)/DABS(FRIED-RF)
C     *     +1.0D2*(FRI2)
          RETURN
          END



C*********************************************************************************************
      Function Fmin(ax,bx,f,Tol)

c      Para funcion f(x) unimodular dentro del intervalo (ax,bx)
c      se encuentra el punto x=Fmin, con incertidumbre Tol,
c      donde la funcion acepta el valor minimo
c
c  Parametros de la entrada:
c        ax, bx -  los extremos del intervalo analizado;
c        f      -  subprograma-funcion, que calcula el valor
c                 de la funcion f(x) dentro (ax,bx)
c        Tol    - la precision requerida
c  Parametros de la salida:
c        Fmin   - la abscisa, que aproxima el punto, donde
c                 f(x) tiene el minimo
c  Metodo:
c     Se utiliza la combinacion del metodo de Fibonacci y
c     de interpolacion cuadratica sucesiva. El algoritmo
c     no puede ser menos rapido que el metodo de Fibonacci
c     y es mucho mas rapido si la funcion f(x) tiene la segunda
c     derivada continua.
c     Programa utiliza el minimo posible accesos al calculo
c     de la funcion, en cado paso se calcula el nuevo valor
c     no mas que una vez y nunca se calcula el nuevo valor,
c     si la distancia entre dos puntos sucesivos es menor que
c     macheps*abs(x)+tol/3, donde macheps es la precision
c     aritmetica del computador y x la abscisa del punto anterior.
c     Se utilizo el algoritmo, que esta descrito en el libro:
c     RICHARD BRENT. "ALGORITHMS FOR MINIMIZATION
c     WITHOUT DERIVATVES", Prentice -Hall,1973

       Implicit double precision (a-h,o-z)

       External f
c
c    El parametro C es el cuadrato del valor inversa
c            a seccion de oro
c
        C=0.5d0*(3.0d0-dsqrt(5.0d0))
c
c   El parametro Eps es igual a la raiz cuadrata de
c   la precision aritmetica del computador
       Eps=1.0d0
10     Eps=Eps/2.0d0
       Tol1=1.0d0+Eps
       If(tol1.gt.1.0d0) go to 10
       Eps=dsqrt(Eps)
c
c  Asignacion de los valores iniciales
c
       A=ax
       B=bx
       V=A+C*(B-A)
       W=V
       X=V
       E=0.0d0
       FX=f(X)
       FV=FX
       FW=FX
c
c  Aqui se inicia el ciclo principal
c
20     XM=0.5d0*(A+B)
       TOL1=Eps*Dabs(X)+tol/3.0d0
       TOL2=2.0d0*TOL1
c
c  Chequeo de terminacion del ciclo principal
c
       If(Dabs(X-XM).le.(TOL2-0.5d0*(B-A))) go to 90
c
c  � Es necesario la seccion de oro ?
c
       If(Abs(E).le.Tol1) go to 40
c
c   Construir la parabola
c
       R=(X-W)*(FX-FV)
       Q=(X-V)*(FX-FW)
       P=(X-V)*Q-(X-W)*R
       Q=2.0d0*(Q-R)
       If(Q.gt.0.0d0) P=-P
       Q=Dabs(Q)
       R=E
       E=D
c
c   � Parabola es aceptable ?
c
30     If(dabs(P).ge.dabs((0.5d0*Q*R))) go to 40
       If(P.ge.Q*(A-X)) go to 40
       If(P.ge.Q*(B-X)) go to 40
c
c   El paso de la interpolacion parabolica
c
       D=P/Q
       U=X+D
c
c   F no se permite calcular demasiado cerca al AX o BX
c
       If((U-A).lt.Tol2) D=Dsign(Tol1,XM-X)
       If((B-U).lt.Tol2) D=Dsign(Tol1,XM-X)
       Go to 50
c
c   El paso de seccion de oro
c
40     If(X.ge.XM) E=A-X
       If(X.lt.XM) E=B-X
       D=C*E
c
c   F no se permite calcular demasiado cerca al X
c
50     If(Dabs(D).ge.Tol1) U=X+D
       If(Dabs(D).lt.Tol1) U=X+Dsign(Tol1,D)
       FU=f(u)
c
c  Asignar valores nuevos a los parametros A,B,V,W y X
c
       If(FU.gt.FX) go to 60
       If(U.ge.X) A=X
       If(U.lt.X) B=X
       V=W
       FV=FW
       W=X
       FW=FX
       X=U
       FX=FU
       go to 20
60     if(U.lt.X) A=U
       If(U.ge.X) B=U
       If(FU.le.FW) go to 70
       If(W.eq.X) go to 70
       If(FU.le.FV) go to 80
       If(V.eq.X) go to 80
       If(V.eq.W) go to 80
       go to 20
70     V=W
       FV=FW
       W=U
       FW=FU
       go to 20
80     V=U
       FV=FU
       go to 20
c
c   Fin del ciclo principal
c
90     Fmin=X
       Return
       End
       
C ************************************************************************************************
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zig_mauro(cw,Dx,Dy,g00,g0n)
      implicit none
      integer l,j,i,npoles,teste,n
      parameter (n=3)
      double precision fj,fM,tM,wght,sq3,pi,Dx,Dy
      complex*16 cw,g00,g0n,argexp
      complex*16 k_pole,gp,zz,cacos,argpole,residue,aux_g00,cum


      sq3 = sqrt(3.d0)
      pi = 4.d0*datan(1.d0)
      cum = dcmplx(0.d0,1.d0)
      fM = float(n)
      g00 = dcmplx(0.d0,0.d0)
      g0n = dcmplx(0.d0,0.d0)
c
      do j=0,n
         wght = 1.d0
         if (j==0) wght = 0.5d0
         if (j==n) wght = 0.5d0
         fj = float(j)
         zz = argpole(cw,fM,fj)
         k_pole = 2.d0*cacos(zz)/sq3
         teste = 0
         if (dabs(dreal(k_pole))<=2.d0*pi/sq3) then
            teste = 1
            if (dimag(k_pole)>0.d0) then
               gp = k_pole
            else
               gp = -k_pole
            end if
         end if
         if (teste==1) then
            aux_g00 = residue(cw,gp,fM,fj)
            g00 = g00 + aux_g00*wght
            argexp = cum*(gp*Dx + 2.d0*fj*pi*Dy/fM)
            g0n = g0n + cdexp(argexp)*aux_g00*wght
         end if
      end do
c
      return
      end
c
c -----------------------------------------
c
      function cacos(z)
      implicit none
      complex*16 cacos,z
c
      cacos = -2.d0*dcmplx(0.d0,1.d0)*cdlog(cdsqrt(0.5d0*(1.d0+z)) +
     &        dcmplx(0.d0,1.d0)*cdsqrt(0.5d0*(1.d0-z)) )
      return
      end
c
c -----------------------------------------
c
      function residue(cw,k_pole,fm,fj)
      implicit none
      double precision fm,fj,t,sq3,pi
C      DOUBLE PRECISION  EF,EQ,AMU,V,V2,BETA,CC,D
      complex*16 residue,cw,k_pole,cum,cseno
C      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D,COU

      cum = dcmplx(0.d0,1.d0)
      sq3 = dsqrt(3.d0)
      t = 1.d0
      pi = 4.d0*datan(1.d0)
      residue = cum*cw/(4.d0*fM*t*t*cseno(sq3*k_pole*0.5d0)*
     &                 cos(fj*pi/fM))
c
      return
      end
c
c -----------------------------------------
c
      function argpole(cw,fM,fj)
      implicit none
      double precision fM,fj,t,pi
      complex*16 argpole,cw
c
      pi = 4.d0*datan(1.d0)
      t = 1.d0
c
      argpole = ((cw/t)*(cw/t) - 1.d0 -4.d0*(cos(fj*pi/fM))*
     &          (cos(fj*pi/fM)))/(4.d0*cos(fj*pi/fM))
      return
      end
c
c ------------------------------------------
c
      function cseno(zz)
      implicit none
      complex*16 cseno,zz,cum
      double precision zr,zi
c
      cum = dcmplx(0.d0,1.d0)
      zr = dreal(zz)
      zi = dimag(zz)
      cseno = dsin(zr)*dcosh(zi)+cum*dcos(zr)*dsinh(zi)
c
      return
      end
C*****************************************************************************

C*************************************************************************
C*************************************************************************

**********************************************************************
*********************************************************************
        SUBROUTINE QUAD(A,B,RESULT,K,EPSIL,NPTS,ICHECK,F)
        IMPLICIT REAL*8 (A-H,P-Z)
        DIMENSION FUNCT(127),P(381),RESULT(8)
CC      THIS SUBROUTINE ATTEMPTS TO CALCULATE THE INTEGRAL OF F(X)
CC      OVER THE INTERVAL  (A,B) WITH RELATIVE ERROR NOT
CC      EXCEEDING  EPSIL.
CC      THE RESULT IS OBTAINED USING A SEQUENCE OF 1, 3, 7, 15, 31, 63,
CC      127 AND 255 POINT INTERLACING FORMALAE (NO INTEGRAND
CC      EVALUATIONS ARE WASTED) OF RESPECTIVE DEGREE 1,5,11,23,
CC      47,95,191 AND 383. THE FORMULAE ARE BASED ON THE OPTIMAL
CC      EXTENSION OF THE 3-POINT GAUSS FORMULA. DETAILS  OF
CC      THE FORMULAE ARE GIVEN IN THE OPTIMUM ADDITION OF POINTS
CC      TO QUADRATURE FORMULAE  BY T. N. L. PATRERSON,MATHS.COMP.
C       VOL. 22,847-856,1968.
C                 ******INPUT *****
CC      A          LOWER LIMIT OF INTEGRATION.
C       B          UPPER LIMIT OF INTEGRATION.
C       EPSIL      RELATIVE ACCURANCY REQUIRED. WHEN THE RELATIVE
C                  DIFFERENCE OF TWO SUCCESSIVE FORMULAE DOES NOT
C                  EXCEED  *EPSL* THE LAST FORMULA COMPUTED IS TAKEN
C                  AS THE RESULT.
C       F          F(X) IS THE INTEGRAND.
C                 ******OUTPUT *****
C       RESULT     THIS ARRAY,WHICH SHOULD BE DECLARED TO HAVE AT
C                  LEAST 8 ELEMENTS,HOLDS THE RESULTS OBTAINED BY
C                  THE 1,3,7, ETC., POINT FORMULAE.THE NUMBER OF
C                  FORMULAE COMPUTED DEPENDS ON *EPSIL*
C       K           RESULT(K) HOLDS THE VALUE OF THE INTERVAL TO THE
C                  SPECIFIED RELATIVE ACCURANCY.
C       NPTS       NUMBER INTEGRAND EVALUATIONS.
C       ICHECK     ON  EXIT NORMALLY ICHECK =0. HOWEVER IF CONVERGENCE
C                  TO THE ACCURANCY REQUESTED IS NOT ACHIEVED ICHECK=1
C                  ON EXIT
CC      ABSCISSAE AND WEIGHTS OF QUADRATURE RULES ARE STACKED IN
C       ARRAY *P* IN THE ORDER IN WHICH THEY ARE NEEDED.
        DATA
     *P( 1),P( 2),P( 3),P( 4),P( 5),P( 6),P( 7),
     *P( 8),P( 9),P(10),P(11),P(12),P(13),P(14),
     *P(15),P(16),P(17),P(18),P(19),P(20),P(21),
     *P(22),P(23),P(24),P(25),P(26),P(27),P(28)/
     *0.77459666924148337704D+00,0.55555555555555555556D+00,
     *0.88888888888888888889D+00,0.26848808986833344073D+00,
     *0.96049126870802028342D+00,0.10465622602646726519D+00,
     *0.43424374934680255800D+00,0.40139741477596222291D+00,
     *0.45091653865847414235D+00,0.13441525524378422036D+00,
     *0.51603282997079739697D-01,0.20062852937698902103D+00,
     *0.99383196321275502221D+00,0.17001719629940260339D-01,
     *0.88845923287225699889D+00,0.92927195315124537686D-01,
     *0.62110294673722640294D+00,0.17151190913639138079D+00,
     *0.22338668642896688163D+00,0.21915685840158749640D+00,
     *0.22551049979820668739D+00,0.67207754295990703540D-01,
     *0.25807598096176653565D-01,0.10031427861179557877D+00,
     *0.84345657393211062463D-02,0.46462893261757986541D-01,
     *0.85755920049990351154D-01,0.10957842105592463824D+00/
        DATA
     *P(29),P(30),P(31),P(32),P(33),P(34),P(35),
     *P(36),P(37),P(38),P(39),P(40),P(41),P(42),
     *P(43),P(44),P(45),P(46),P(47),P(48),P(49),
     *P(50),P(51),P(52),P(53),P(54),P(55),P(56)/
     *0.99909812496766759766D+00,0.25447807915618744154D-02,
     *0.98153114955374010687D+00,0.16446049854387810934D-01,
     *0.92965485742974005667D+00,0.35957103307129322097D-01,
     *0.83672593816886873550D+00,0.56979509494123357412D-01,
     *0.70249620649152707861D+00,0.76879620499003531043D-01,
     *0.53131974364437562397D+00,0.93627109981264473617D-01,
     *0.33113539325797683309D+00,0.10566989358023480974D+00,

     *0.11248894313318662575D+00,0.11195687302095345688D+00,
     *0.11275525672076869161D+00,0.33603877148207730542D-01,
     *0.12903800100351265626D-01,0.50157139305899537414D-01,
     *0.42176304415588548391D-02,0.23231446639910269443D-01,
     *0.42877960025007734493D-01,0.54789210527962865032D-01,
     *0.12651565562300680114D-02,0.82230079572359296693D-02,
     *0.17978551568128270333D-01,0.28489754745833548613D-01/
        DATA
     *P(57),P(58),P(59),P(60),P(61),P(62),P(63),
     *P(64),P(65),P(66),P(67),P(68),P(69),P(70),
     *P(71),P(72),P(73),P(74),P(75),P(76),P(77),
     *P(78),P(79),P(80),P(81),P(82),P(83),P(84)/
     *0.38439810249455532039D-01,0.46813554990628012403D-01,
     *0.52834946790116519862D-01,0.55978436510476319408D-01,
     *0.99987288812035761194D+00,0.36322148184553065969D-03,
     *0.99720625937222195908D+00,0.25790497946856882724D-02,
     *0.98868475754742947994D+00,0.61155068221172463397D-02,
     *0.97218287474858179658D+00,0.10498246909621321898D-01,
     *0.94634285837340290515D+00,0.15406750466559497802D-01,
     *0.91037115695700429250D+00,0.20594233915912711149D-01,
     *0.86390793819369047715D+00,0.25869679327214746911D-01,
     *0.80694053195021761186D+00,0.31073551111687964880D-01,
     *0.73975604435269475868D+00,0.36064432780782572640D-01,
     *0.66290966002478059546D+00,0.40715510116944318934D-01,
     *0.57719571005204581484D+00,0.44914531653632197414D-01,
     *0.48361802694584102756D+00,0.48564330406673198716D-01/
        DATA
     *P(85),P(86),P(87),P(88),P(89),P(90),P(91),
     *P(92),P(93),P(94),P(95),P(96),P(97),P(98),
     *P(99),P(100),P(101),P(102),P(103),P(104),P(105),
     *P(106),P(107),P(108),P(109),P(110),P(111),P(112)/
     *0.38335932419873034692D 00,0.51583253952048458777D-01,
     *0.27774982202182431507D 00,0.53905499335266063927D-01,
     *0.16823525155220746498D 00,0.55481404356559363988D-01,
     *0.56344313046592789972D-01,0.56277699831254301273D-01,
     *0.56377628360384717388D-01,0.16801938574103865271D-01,
     *0.64519000501757369228D-02,0.25078569652949768707D-01,
     *0.21088152457266328793D-02,0.11615723319955134727D-01,
     *0.21438980012503867246D-01,0.27394605263981432516D-01,
     *0.63260731936263354422D-03,0.41115039786546930472D-02,
     *0.89892757840641357233D-02,0.14244877372916774306D-01,
     *0.19219905124727766019D-01,0.23406777495314006201D-01,
     *0.26417473395058259931D-01,0.27989218255238159704D-01,
     *0.18073956444538835782D-03,0.12895240826104173921D-02,
     *0.30577534101755311361D-02,0.52491234548088591251D-02/
        DATA
     *P(113),P(114),P(115),P(116),P(117),P(118),P(119),
     *P(120),P(121),P(122),P(123),P(124),P(125),P(126),
     *P(127),P(128),P(129),P(130),P(131),P(132),P(133),
     *P(134),P(135),P(136),P(137),P(138),P(139),P(140)/
     *0.77033752332797418482D-02,0.10297116957956355524D-01,
     *0.12934839663607373455D-01,0.15536775555843982440D-01,
     *0.18032216390391286320D-01,0.20357755058472159467D-01,
     *0.22457265826816098707D-01,0.24282165203336599358D-01,
     *0.25791626976024229388D-01,0.26952749667633031963D-01,
     *0.27740702178279681994D-01,0.28138849915627150636D-01,
     *0.99998243035489159858D 00,0.50536095207862517625D-04,
     *0.99959879967191068325D 00,0.37774664632698466027D-03,
     *0.99831663531840739253D 00,0.93836984854238150079D-03,
     *0.99572410469840718851D 00,0.16811428654214699063D-02,
     *0.99149572117810613240D 00,0.25687649437940203731D-02,
     *0.98537149959852037111D 00,0.35728927835172996494D-02,
     *0.97714151463970571416D 00,0.46710503721143217474D-02,
     *0.96663785155841656709D 00,0.58434498758356395076D-02/
        DATA
     *P(141),P(142),P(143),P(144),P(145),P(146),P(147),
     *P(148),P(149),P(150),P(151),P(152),P(153),P(154),
     *P(155),P(156),P(157),P(158),P(159),P(160),P(161),
     *P(162),P(163),P(164),P(165),P(166),P(167),P(168)/
     *0.95373000642576113641D 00,0.70724899954335554680D-02,
     *0.93832039777959288365D 00,0.83428387539681577056D-02,
     *0.92034002547001242073D 00,0.96411777297025366953D-02,
     *0.89974489977694003664D 00,0.10955733387837901648D-01,
     *0.87651341448470526974D 00,0.12275830560082770087D-01,
     *0.85064449476835027976D 00,0.13591571009765546790D-01,
     *0.82215625436498040737D 00,0.14893641664815182035D-01,
     *0.79108493379984836143D 00,0.16173218729577719942D-01,
     *0.75748396638051363793D 00,0.17421930159464173747D-01,
     *0.72142308537009891548D 00,0.18631848256138790186D-01,
     *0.68298743109107922809D 00,0.19795495048097499488D-01,
     *0.64227664250975951377D 00,0.20905851445812023852D-01,
     *0.59940393024224289297D 00,0.21956366305317824939D-01,
     *0.55449513263193254887D 00,0.22940964229387748761D-01/
        DATA
     *P(169),P(170),P(171),P(172),P(173),P(174),P(175),
     *P(176),P(177),P(178),P(179),P(180),P(181),P(182),
     *P(183),P(184),P(185),P(186),P(187),P(188),P(189),
     *P(190),P(191),P(192),P(193),P(194),P(195),P(196)/
     *0.50768775753371660215D 00,0.23854052106038540080D-01,
     *0.45913001198983233287D 00,0.24690524744487676909D-01,
     *0.40897982122988867241D 00,0.25445769965464765813D-01,
     *0.35740383783153215238D 00,0.26115673376706097680D-01,
     *0.30457644155671404334D 00,0.26696622927450359906D-01,
     *0.25067873030348317661D 00,0.27185513229624791819D-01,
     *0.19589750271110015392D 00,0.27579749566481873035D-01,
     *0.14042423315256017459D 00,0.27877251476613701609D-01,
     *0.84454040083710883710D-01,0.28076455793817246607D-01,
     *0.28184648949745694339D-01,0.28176319033016602131D-01,
     *0.28188814180192358694D-01,0.84009692870519326354D-02,
     *0.32259500250878684614D-02,0.12539284826474884353D-01,
     *0.10544076228633167722D-02,0.58078616599775673635D-02,
     *0.10719490006251933623D-01,0.13697302631990716258D-01/
        DATA
     *P(197),P(198),P(199),P(200),P(201),P(202),P(203),
     *P(204),P(205),P(206),P(207),P(208),P(209),P(210),
     *P(211),P(212),P(213),P(214),P(215),P(216),P(217),
     *P(218),P(219),P(220),P(221),P(222),P(223),P(224)/
     *0.31630366082226447689D-03,0.20557519893273465236D-02,
     *0.44946378920320678616D-02,0.71224386864583871532D-02,
     *0.96099525623638830097D-02,0.11703388747657003101D-01,
     *0.13208736697529129966D-01,0.13994609127619079852D-01,
     *0.90372734658751149261D-04,0.64476204130572477933D-03,
     *0.15288767050877655684D-02,0.26245617274044295626D-02,
     *0.38516876166398709241D-02,0.51485584789781777618D-02,
     *0.64674198318036867274D-02,0.77683877779219912200D-02,
     *0.90161081951956431600D-02,0.10178877529236079733D-01,
     *0.11228632913408049354D-01,0.12141082601668299679D-01,
     *0.12895813488012114694D-01,0.13476374833816515982D-01,
     *0.13870351089139840997D-01,0.14069424957813575318D-01,
     *0.25157870384280661489D-04,0.18887326450650491366D-03,
     *0.46918492424785040975D-03,0.84057143271072246365D-03/
        DATA
     *P(225),P(226),P(227),P(228),P(229),P(230),P(231),
     *P(232),P(233),P(234),P(235),P(236),P(237),P(238),
     *P(239),P(240),P(241),P(242),P(243),P(244),P(245),
     *P(246),P(247),P(248),P(249),P(250),P(251),P(252)/
     *0.12843824718970101768D-02,0.17864463917586498247D-02,
     *0.23355251860571608737D-02,0.29217249379178197538D-02,
     *0.35362449977167777340D-02,0.41714193769840788528D-02,
     *0.48205888648512683476D-02,0.54778666939189508240D-02,
     *0.61379152800413850435D-02,0.67957855048827733948D-02,
     *0.74468208324075910174D-02,0.80866093647888599710D-02,
     *0.87109650797320868736D-02,0.93159241280693950932D-02,
     *0.98977475240487497440D-02,0.10452925722906011926D-01,
     *0.10978183152658912470D-01,0.11470482114693874380D-01,
     *0.11927026053019270040D-01,0.12345262372243838455D-01,
     *0.12722884982732382906D-01,0.13057836688353048840D-01,
     *0.13348311463725179953D-01,0.13592756614812395910D-01,
     *0.13789874783240936517D-01,0.13938625738306850804D-01,
     *0.14038227896908623303D-01,0.14088159516508301065D-01/
        DATA
     *P(253),P(254),P(255),P(256),P(257),P(258),P(259),
     *P(260),P(261),P(262),P(263),P(264),P(265),P(266),
     *P(267),P(268),P(269),P(270),P(271),P(272),P(273),
     *P(274),P(275),P(276),P(277),P(278),P(279),P(280)/
     *0.99999759637974846462D 00,0.69379364324108267170D-05,
     *0.99994399620705437576D 00,0.53275293669780613125D-04,
     *0.99976049092443204733D 00,0.13575491094922871973D-03,
     *0.99938033802502358193D 00,0.24921240048299729402D-03,
     *0.99874561446809511470D 00,0.38974528447328229322D-03,
     *0.99780535449595727456D 00,0.55429531493037471492D-03,
     *0.99651414591489027385D 00,0.74028280424450333046D-03,
     *0.99483150280062100052D 00,0.94536151685852538246D-03,
     *0.99272134428278861533D 00,0.11674841174299594077D-02,
     *0.99015137040077015918D 00,0.14049079956551446427D-02,
     *0.98709252795403406719D 00,0.16561127281544526052D-02,
     *0.98351865757863272876D 00,0.19197129710138724125D-02,
     *0.97940628167086268381D 00,0.21944069253638388388D-02,
     *0.97473445975240266776D 00,0.24789582266575679307D-02/
        DATA
     *P(281),P(282),P(283),P(284),P(285),P(286),P(287),
     *P(288),P(289),P(290),P(291),P(292),P(293),P(294),
     *P(295),P(296),P(297),P(298),P(299),P(300),P(301),
     *P(302),P(303),P(304),P(305),P(306),P(307),P(308)/
     *0.96948465950245923177D 00,0.27721957645934509940D-02,
     *0.96364062156981213252D 00,0.30730184347025783234D-02,
     *0.95718821610986096274D 00,0.33803979910869203823D-02,
     *0.95011529752129487656D 00,0.36933779170256508183D-02,
     *0.94241156519108305981D 00,0.40110687240750233989D-02,
     *0.93406843615772578800D 00,0.43326409680929828545D-02,
     *0.92507893290707565236D 00,0.46573172997568547773D-02,
     *0.91543758715576504064D 00,0.49843645647655386012D-02,
     *0.90514035881326159519D 00,0.53130866051870565663D-02,
     *0.89418456833555902286D 00,0.56428181013844441585D-02,
     *0.88256884024734190684D 00,0.59729195655081658049D-02,
     *0.87029305554811390585D 00,0.63027734490857587172D-02,
     *0.85735831088623215653D 00,0.66317812429018878941D-02,
     *0.84376688267270860104D 00,0.69593614093904229394D-02/
        DATA
     *P(309),P(310),P(311),P(312),P(313),P(314),P(315),
     *P(316),P(317),P(318),P(319),P(320),P(321),P(322),
     *P(323),P(324),P(325),P(326),P(327),P(328),P(329),
     *P(330),P(331),P(332),P(333),P(334),P(335),P(336)/
     *0.82952219463740140018D 00,0.72849479805538070639D-02,
     *0.81462878765513741344D 00,0.76079896657190565832D-02,
     *0.79909229096084140180D 00,0.79279493342948491103D-02,
     *0.78291939411828301639D 00,0.82443037630328680306D-02,
     *0.76611781930376009072D 00,0.85565435613076896192D-02,
     *0.74869629361693660282D 00,0.88641732094824942641D-02,
     *0.73066452124218126133D 00,0.91667111635607884067D-02,
     *0.71203315536225203459D 00,0.94636899938300652943D-02,
     *0.69281376977911470289D 00,0.97546565363174114611D-02,
     *0.67301883023041847920D 00,0.10039172044056840798D-01,
     *0.65266166541001749610D 00,0.10316812330947621682D-01,
     *0.63175643771119423041D 00,0.10587167904885197931D-01,
     *0.61031811371518640016D 00,0.10849844089337314099D-01,
     *0.58836243444766254143D 00,0.11104461134006926537D-01/
        DATA
     *P(337),P(338),P(339),P(340),P(341),P(342),P(343),
     *P(344),P(345),P(346),P(347),P(348),P(349),P(350),
     *P(351),P(352),P(353),P(354),P(355),P(356),P(357),
     *P(358),P(359),P(360),P(361),P(362),P(363),P(364)/
     *0.56590588542365442262D 00,0.11350654315980596602D-01,
     *0.54296566649831149049D 00,0.11588074033043952568D-01,
     *0.51955966153745702199D 00,0.11816385890830235763D-01,
     *0.49570640791876146017D 00,0.12035270785279562630D-01,
     *0.47142506587165887693D 00,0.12244424981611985899D-01,
     *0.44673538766202847374D 00,0.12443560190714035263D-01,
     *0.42165768662616330006D 00,0.12632403643542078765D-01,
     *0.39621280605761593918D 00,0.12810698163877361967D-01,
     *0.37042208795007823014D 00,0.12978202239537399286D-01,
     *0.34430734159943802278D 00,0.13134690091960152836D-01,
     *0.31789081206847668318D 00,0.13279951743930530650D-01,
     *0.29119514851824668196D 00,0.13413793085110098513D-01,
     *0.26424337241092676194D 00,0.13536035934956213614D-01,
     *0.23705884558982972721D 00,0.13646518102571291428D-01/
        DATA
     *P(365),P(366),P(367),P(368),P(369),P(370),P(371),
     *P(372),P(373),P(374),P(375),P(376),P(377),P(378),
     *P(379),P(380),P(381)/
     *0.20966523824318119477D 00,0.13745093443001896632D-01,
     *0.18208649675925219825D 00,0.13831631909506428676D-01,
     *0.15434681148137810869D 00,0.13906019601325461264D-01,
     *0.12647058437230196685D 00,0.13968158806516938516D-01,
     *0.98482396598119202090D-01,0.14017968039456608810D-01,
     *0.70406976042855179063D-01,0.14055382072649964277D-01,
     *0.42269164765363603212D-01,0.14080351962553661325D-01,
     *0.14093886410782462614D-01,0.14092845069160408355D-01,
     *0.14094407090096179347D-01/
        ZERO=0.D0
        DOIS=2.D0
        ICHECK=0
C       CHECK FOR TRIVIAL CASE
        IF(A.EQ.B)GO TO 70
C       SCALE FACTORS
        SUM=(B+A)/DOIS
        DIFF=(B-A)/DOIS
C       1-POIN GAUSS
        FZERO=F(SUM)
        RESULT(1)=DOIS*FZERO*DIFF
        I=0
        IOLD=0
        INEW=1
        K=2
        ACUM=ZERO
        GO TO 30
10      IF(K.EQ.8)GO TO 50
        K=K+1
        ACUM=ZERO
C       CONTRIBUTION FROM FUNCTION VALUES ALREADY COMPUTED
        DO 20 J=1,IOLD
        I=I+1
        ACUM=ACUM+P(I)*FUNCT(J)
20      CONTINUE
C       CONTRIBUTION FROM NEW FUNCTION VALUES
30      IOLD=IOLD+INEW
        DO 40 J=INEW,IOLD
        I=I+1
        X=P(I)*DIFF
        FUNCT(J)=F(SUM+X)+F(SUM-X)
        I=I+1
        ACUM=ACUM+P(I)*FUNCT(J)
40      CONTINUE
        INEW=IOLD+1
        I=I+1
        RESULT(K)=(ACUM+P(I)*FZERO)*DIFF
C       CHECK FOR CONVERGENCE
        IF(DABS(RESULT(K)-RESULT(K-1))-EPSIL*DABS(RESULT(K)))60,
     *60,10
C       CONVERGENCE NOT ACHIEVED
50      ICHECK=1
C       NORMAL TERMINATION
60      NPTS=INEW+IOLD
        RETURN
C       TRIVIAL CASE
70      K=2
        RESULT(1)=ZERO
        RESULT(2)=ZERO
        NPTS=0
        RETURN
        END



        FUNCTION QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,F)
CC
CC
CC        THIS FUNCTION ROUTINE PERFORMS AUTOMATIC INTEGRATION
C       OVER A FINITE INTERVAL USING THE BASIC INTEGRATION
C       ALGORITHM QUAD TOGETHER WIHT, IF NECESSARY AN ADAPTIVE
C       SUBDIVISION PROCESS. IT IS GENERALLY MORE EFFICIENT THAN
C       THE NON-ADAPTIVE ALGORITHM QSUB BUT IS LIKILY TO BE LESS
C       RELIABLE(SEE COMP. J., 14,189,1971).
C         THE CALL TAKES THE FORM
C                QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,F)
C       AND CAUSES F(X) TO BE  INTEGRATED OER (A,B) WIHT RELATIVE
C       ERROR HOPEFULLY NOT  EXCCEEDING EPSIL. SHOULD QUAD CONVERGE
C       (ICHECK=0) THEN QSUBA WILL RETURN THE VALUE OBTAINED BY IT
C       OTHERWISE SUBDIVISION  WILL BE INVOKED AS A RESCUE
C       OPERATION IN AN ADAPTIVE MANNER.THE ARGUMENT RELERR GIVES
C       A CRUDE ESTIMATE OF THE ACTUAL  RELATIVE ERROR OBATAINED.
C       THE SUBDIVISION STRATEGY IS AS FOLLOWS
C       AT EACH STAGE OF THE PROCESS AN INTERVAL IS PRESENTED FOR
C       SUBDIVISION (INITIALLY THIS WILL BE THE WHOLE INTERVAL
C       (A,B)) THE INTERVAL IS HALVED AND QUAD APPLIED TO EACH
C       SUBINTERVAL. SHOULD QUAD FAIL ON THE FIRST SUBINTERVAL
C       THE SUBINTERVAL IS STACKED FOR FUTURE SUBDIVISION AND  THE
C       SECOND SUBINTERVAL IMMEDIATILY EXAMINED. SHOULD QUAD FAIL
C       ON THE SECOND SUBINTERVAL THE INTERVAL IS
C       IMMEDIATILY SUBDIVIDED AND THE WHOLE PROCESS REPEATED
C       EACH TIME A CONVERGED RESULT IS OBTAINED IT IS
C       ACCUMULATED AS THE PARTIAL VALUE OF THE INTEGRAL. WHEN
C       QUAD CONVERGES ON BOTH SUBINTERVALS THE INTERVAL LAST
C       STACKED IS CHOSEN NEXT FOR SUBDIVISION AND THE PROCESS
C       REPEATED. A SUBINTERVAL IS NOT EXAMINED AGAIN ONCE A
C       CONVERGED RESULT IS OBTAINED FOR IT SO THAT A  SPURIOUS
C       CONVERGENGE IS MORE LIKELY TO SLIP THROUGH THAN FOR THE
C       NON-ADAPTIVE ALGORITHM QSUB.
C       THE CONVERGE CRITERION OF QUAD IS SLIGHTLY RELAXED
C       IN THAT A PANEL IS DEEMED TO HAVE BEEN SUCCESSFULLY
C       INTEGRATED IF EITHER QUAD CONVERGES OR THE ESTIMATED
C       ABSOLUTE ERRROR  COMMITTED ON THIS PANEL DOES NOT EXCCEED
C       EPSIL TIMES THE ESTIMATED ABSOLUTE VALUE OF THE INTEGRAL
C       OVER (A,B). THIS RELAXATION IS TO TRY TO TAKES ACCOUNT OF
C       A COMMON SITUATION WHERE ONE PARTICULAR PANEL CAUSES
C       SPECIAL DIFFICULTY ,PERHAPS DUE TO A SINGULARITY OF
C       SOME TYPE. IN THIS CASE QUAD COULD OBTAIN NEARLY EXACT
C       ANSWERS ON ALL OTHER PANELS AND SO THE RELATIVE ERROR FOR
C       THE TOTAL INTEGRATION WOULD BE ALMOST ENTIRELY DUE TO THE
C       DELINQUENT PANEL. WITHOUT THIS CONDITION THE COMPUTATION
C       MIGHT CONTINUE DESPITE THE REQUESTED RELATIVE ERROR BEING
C       ACHIEVED. IF THIS RELAXED CONVERGENCE CRITERION IS APPLIED TOO
C       MANY TIMES,I.E. ICHECK = 2 IN MANY PANELS, THE RELATIVE
C       ERROR WOULD BE LARGER THAN EPSIL. IN THIS CASE CHECK
C        EPSIL  VS  RELERR .
C          THE OUTCOME OF INTEGRATION IS INDICATED BY ICHECK.
C          ICHECK=0  CONVERGENCE OBTAINED WITHOUT INVOKING SUB-
C                DIVISION. THIS WOULD CORRESPOND TO THE
C                DIRECT USE OF QUAD.
C          ICHECK=1   RESULT OBTAINED AFTER INVOKING SUBDIVISION.
C          ICHECK=2   AS FOR ICHECK=1,BUT AT SOME POINT THE
C                 RELAXED CONVERGENCE CRITERION WAS USED.
C                 THE RISK OF UNDERESTIMATING THE RELATIVE
C                 ERROR WILL BE INCREASED. IF NECESSARY,
C                 CONFIDENCE MAY BE RESTORED BY CHEKING
C                 EPSIL AND RELERR FOR A SERIOUS DISCREPANCY.
C          ICHECK NEGATIVE
C                 IF DURING THE SUBDIVISION PROCESS THE STACK
C                 OF DELINQUENT INTERVALS BECOMES FULL(IT IS
C                 PRESENTLY SET TO HOLD AT MOST 100 NUMBERS)
C                 A RESULT IS OBTAINED BY CONTINUING THE
C                 INTEGRATION IGNORING CONVERGENCE FAILURES
C                 WHICH CANNOT BE ACCOMODATED ON THE STACK.
C                 THIS OCCURRENCE IS FLAGGED BY RETURNING
C                 ICHECK WITH NEGATIVE SIGN.
C          THE RELIABILITY OF THE  ALGORITHM WILL DECREASE FOR LARGE
C          VALUES OF EPSIL.IT IS RECOMMENDED THAT EPSIL SHOULD
C          GENERALY BE LESS THAN ABOUT 0.001.
        IMPLICIT REAL*8 (A-H,P-Z)
        DIMENSION RESULT(8),STACK(100)
        EXTERNAL F
      ABCTR = 0.D0
      VALUE = QSUBC(A,B,EPSIL,NPTS,ICHECK,RELERR,F,ABCTR)
      QSUBA = VALUE
      RETURN
        END
        
      FUNCTION QSUBC(A,B,EPSIL,NPTS,ICHECK,RELERR,F,ABCTR)
C
C     ESTA SUBRUTINA TIENE UNA VARIABLE MAS:  ABCTR
C     SI  ABCTR = 0.D0  ES IGUAL A  QSUBA.
C     TIENE COMO OBJETO COMPARAR  EL RESULTADO DE
C     LA PRIMERA PASADA POR  QUAD  CON  EPSIL*ABCTR
C     CUANDO  DABS(QSUBC) .LT. EPSIL*DABS(ABCTR)  NO SUBDIVIDE.
C     DABS(QSUBC).LT.DABS(ABCTR)  ->  ESTIM = DABS(ABCTR)*EPSIL.
C     ESTA VARIANTE ES INTERESANTE QUANDO SE HA SUBDIVIDIDO LA
C     INTEGRAL EXTERNAMENTE A  QSUB2 , Y PUEDE PRETENDERSE CALCULAR
C     UNA DE LAS PARTES QUE ES RELATIVAMENTE MUY PEQUENA CON
C     DEMASIADA APROXIMACION.
C     HAY QUE TENER MUCHO CUIDADO CUANDO DIFERENTES PARTES TIENEN
C     SIGNO DIFERENTE, PUES EN ESE CASO SE COMPENSAN Y EL ERROR
C     RELATIVO PODRIA SER INADEQUADO. EN ESE CASO HAY QUE SEPARAR
C     LAS PARTES CON CUIDADO, PUES EL ORDEN DE SEPARACION AFECTARA
C     EL RESULTADO.
C
C
C
C
      IMPLICIT REAL*8 (A-H,P-Z)
        DIMENSION RESULT(8),STACK(100)
        EXTERNAL F
        DATA ISMAX/100/
      ABCTRQ = DABS(ABCTR)
      QSBMIN = EPSIL*ABCTRQ
        ZERO=0.D0
        CALL QUAD(A,B,RESULT,K,EPSIL,NPTS,ICHECK,F)
        QSUBC=RESULT(K)
      IF(ICHECK.NE.0)  GO TO 5
C     HAY CONVERGENCIA SIN SUBDIVISIONES ( ICHECK = 0 )  Y CALCULA  RELERR
C     PARA DEVOLVERLO A SUB. LLAMADA.
      RELERR = 0.D0
      IF(QSUBC.NE.ZERO)
     1 RELERR=DABS((RESULT(K)-RESULT(K-1))/QSUBC)
        RETURN
C
C*****CHECK IF SUBDIVISION IS NEEDED
C     SI ABS(QSUBC) < ABCTR*EPSIL  SUSPENDE SUBDIVISIONES, PUES
C     EL VALOR ESTIMADO DE LA INTEGRAL ES DEMASIADO PEQUENO.
    5 ABQSUBC = DABS(QSUBC)
      IF(ABQSUBC.LT.QSBMIN) RETURN
C       SUBDIVIDED
        ESTIM=DABS(QSUBC*EPSIL)
C     CUANDO  ABS(QSUBC) < ABCTR  USA  ESTIM = EPSIL*ABCTR EN LUGAR
C     DE  ESTIM = EPSIL*ABC(QSUBC) PARA APLICAR CRITERIO RELAJADO
C     DE CONVERGENCIA. ESTO PUEDE EVITAR CALCULOS INNECESARIOS.
      IF(ABQSUBC.LT.ABCTRQ) ESTIM = QSBMIN
C     USA  RELERR PARA ACUMULAR EL ERROR ABSOLUTO. SOLO AL FINAL
C     LO DIVIDE POR  QSUBC  PARA OBTENER EL VERDADERO  RELERR.
        RELERR=0.0D0
        QSUBC=0.0D0
        IS=1
        IC=1
        SUB1=A
        SUB3=B
10      SUB2=(SUB1+SUB3)*0.5D0
        CALL QUAD(SUB1,SUB2,RESULT,K,EPSIL,NF,ICHECK,F)
        NPTS=NPTS+NF
C     COMP ES EL VALOR ABSOLUTO DE LA DIFERENCIA ENTRE LAS DOS ULTIMAS
C     SUBDIVISIONES DE QUAD, QUE USA COMO VALOR ABSOLUTO Y ACUMULA
C     EN  RELERR
        COMP=DABS(RESULT(K)-RESULT(K-1))
        IF(ICHECK.EQ.0)GO TO 30
        IF(COMP.LE.ESTIM)GO TO 70
        IF(IS.GE.ISMAX)GO TO 20
C       STACK SUBINTERVAL (SUB1,SUB2) FOR FURTURE EXAMINATION
        STACK(IS)=SUB1
        IS=IS+1
        STACK(IS)=SUB2
        IS=IS+1
        GO TO 40
20      IC=-IABS(IC)
30      QSUBC=QSUBC+RESULT(K)
C     ACUMULA ERRORES ABSOLUTOS EN  RELERR
        RELERR=RELERR+COMP
40      CALL QUAD(SUB2,SUB3,RESULT,K,EPSIL,NF,ICHECK,F)
        NPTS=NPTS+NF
        COMP=DABS(RESULT(K)-RESULT(K-1))
        IF(ICHECK.EQ.0)GO TO 50
        IF(COMP.LE.ESTIM)GO TO 80
C       SUBDIVIDE INTEVAL (SUB2,SUB3)
        SUB1=SUB2
        GO TO 10
50      QSUBC=QSUBC+RESULT(K)
        RELERR=RELERR+COMP
        IF(IS.EQ.1) GO TO 60
C       SUBDIVIDE THE DELINQUENT INTERVAL LAST STACKED
        IS=IS-1
        SUB3=STACK(IS)
        IS=IS-1
        SUB1=STACK(IS)
        GO TO 10
C       SUBDIVISION RESULT
60      ICHECK=IC
      IF(QSUBC.NE.0.D0)  GO TO 65
C      SI EL ERROR ABSOLUTO ESTIMADO NO ES NULO, Y  QSUBC = 0.D0,
C     DEVUELVE ERROR ABSOLUTO EN  RELERR CON SIGNO CAMBIADO, Y
C     AVISA EN LA TELA.
      IF(RELERR.EQ.0.D0)  RETURN
      WRITE(*,*) ' ATENCION !  QSUBC = ',QSUBC,'  RELERR = ',RELERR
      RELERR = - RELERR
      WRITE(*,*) ' CAMBIA SIGNO DE  RELERR EN LA SALIDA! '
   65 RELERR=RELERR/DABS(QSUBC)
        RETURN
C       RELAXED CONVERGENCE
70      IC=ISIGN(2,IC)
        GO TO 30
80      IC=ISIGN(2,IC)
        GO TO 50
        END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C   IMSL ROUTINE NAME   - UERTST
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - VAX/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - PRINT A MESSAGE REFLECTING AN ERROR CONDITION
C
C   USAGE               - CALL UERTST (IER,NAME)
C
C   ARGUMENTS    IER    - ERROR PARAMETER. (INPUT)
C                           IER = I+J WHERE
C                             I = 128 IMPLIES TERMINAL ERROR MESSAGE,
C                             I =  64 IMPLIES WARNING WITH FIX MESSAGE,
C                             I =  32 IMPLIES WARNING MESSAGE.
C                             J = ERROR CODE RELEVANT TO CALLING
C                                 ROUTINE.
C                NAME   - A CHARACTER STRING OF LENGTH SIX PROVIDING
C                           THE NAME OF THE CALLING ROUTINE. (INPUT)
C
C   PRECISION/HARDWARE  - SINGLE/ALL
C
C   REQD. IMSL ROUTINES - UGETIO,USPKD
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE ERROR MESSAGE PRODUCED BY UERTST IS WRITTEN
C                TO THE STANDARD OUTPUT UNIT. THE OUTPUT UNIT
C                NUMBER CAN BE DETERMINED BY CALLING UGETIO AS
C                FOLLOWS..   CALL UGETIO(1,NIN,NOUT).
C                THE OUTPUT UNIT NUMBER CAN BE CHANGED BY CALLING
C                UGETIO AS FOLLOWS..
C                                NIN = 0
C                                NOUT = NEW OUTPUT UNIT NUMBER
C                                CALL UGETIO(3,NIN,NOUT)
C                SEE THE UGETIO DOCUMENT FOR MORE DETAILS.
C
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE UERTST (IER,NAME)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      INTEGER            NAME(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IEQ,IEQDF,IOUNIT,LEVEL,LEVOLD,NAMEQ(6),
     *                   NAMSET(6),NAMUPK(6),NIN,NMTB
      DATA               NAMSET/1HU,1HE,1HR,1HS,1HE,1HT/
      DATA               NAMEQ/6*1H /
      DATA               LEVEL/4/,IEQDF/0/,IEQ/1H=/
C                                  UNPACK NAME INTO NAMUPK
C                                  FIRST EXECUTABLE STATEMENT
      CALL USPKD (NAME,6,NAMUPK,NMTB)
C                                  GET OUTPUT UNIT NUMBER
      CALL UGETIO(1,NIN,IOUNIT)
C                                  CHECK IER
      IF (IER.GT.999) GO TO 25
      IF (IER.LT.-32) GO TO 55
      IF (IER.LE.128) GO TO 5
      IF (LEVEL.LT.1) GO TO 30
C                                  PRINT TERMINAL MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,35) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,35) IER,NAMUPK
      GO TO 30
    5 IF (IER.LE.64) GO TO 10
      IF (LEVEL.LT.2) GO TO 30
C                                  PRINT WARNING WITH FIX MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,40) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,40) IER,NAMUPK
      GO TO 30
   10 IF (IER.LE.32) GO TO 15
C                                  PRINT WARNING MESSAGE
      IF (LEVEL.LT.3) GO TO 30
      IF (IEQDF.EQ.1) WRITE(IOUNIT,45) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,45) IER,NAMUPK
      GO TO 30
   15 CONTINUE
C                                  CHECK FOR UERSET CALL
      DO 20 I=1,6
         IF (NAMUPK(I).NE.NAMSET(I)) GO TO 25
   20 CONTINUE
      LEVOLD = LEVEL
      LEVEL = IER
      IER = LEVOLD
      IF (LEVEL.LT.0) LEVEL = 4
      IF (LEVEL.GT.4) LEVEL = 4
      GO TO 30
   25 CONTINUE
      IF (LEVEL.LT.4) GO TO 30
C                                  PRINT NON-DEFINED MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,50) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,50) IER,NAMUPK
   30 IEQDF = 0
      RETURN
   35 FORMAT(19H *** TERMINAL ERROR,10X,7H(IER = ,I3,
     1       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
   40 FORMAT(27H *** WARNING WITH FIX ERROR,2X,7H(IER = ,I3,
     1       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
   45 FORMAT(18H *** WARNING ERROR,11X,7H(IER = ,I3,
     1       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
   50 FORMAT(20H *** UNDEFINED ERROR,9X,7H(IER = ,I5,
     1       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
C
C                                  SAVE P FOR P = R CASE
C                                    P IS THE PAGE NAMUPK
C                                    R IS THE ROUTINE NAMUPK
   55 IEQDF = 1
      DO 60 I=1,6
   60 NAMEQ(I) = NAMUPK(I)
   65 RETURN
      END
C   IMSL ROUTINE NAME   - UGETIO
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - VAX/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - TO RETRIEVE CURRENT VALUES AND TO SET NEW
C                           VALUES FOR INPUT AND OUTPUT UNIT
C                           IDENTIFIERS.
C
C   USAGE               - CALL UGETIO(IOPT,NIN,NOUT)
C
C   ARGUMENTS    IOPT   - OPTION PARAMETER. (INPUT)
C                           IF IOPT=1, THE CURRENT INPUT AND OUTPUT
C                           UNIT IDENTIFIER VALUES ARE RETURNED IN NIN
C                           AND NOUT, RESPECTIVELY.
C                           IF IOPT=2, THE INTERNAL VALUE OF NIN IS
C                           RESET FOR SUBSEQUENT USE.
C                           IF IOPT=3, THE INTERNAL VALUE OF NOUT IS
C                           RESET FOR SUBSEQUENT USE.
C                NIN    - INPUT UNIT IDENTIFIER.
C                           OUTPUT IF IOPT=1, INPUT IF IOPT=2.
C                NOUT   - OUTPUT UNIT IDENTIFIER.
C                           OUTPUT IF IOPT=1, INPUT IF IOPT=3.
C
C   PRECISION/HARDWARE  - SINGLE/ALL
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      EACH IMSL ROUTINE THAT PERFORMS INPUT AND/OR OUTPUT
C                OPERATIONS CALLS UGETIO TO OBTAIN THE CURRENT UNIT
C                IDENTIFIER VALUES. IF UGETIO IS CALLED WITH IOPT=2 OR
C                IOPT=3, NEW UNIT IDENTIFIER VALUES ARE ESTABLISHED.
C                SUBSEQUENT INPUT/OUTPUT IS PERFORMED ON THE NEW UNITS.
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE UGETIO(IOPT,NIN,NOUT)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,NIN,NOUT
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NIND,NOUTD
      DATA               NIND/1/,NOUTD/2/
C                                  FIRST EXECUTABLE STATEMENT
      IF (IOPT.EQ.3) GO TO 10
      IF (IOPT.EQ.2) GO TO 5
      IF (IOPT.NE.1) GO TO 9005
      NIN = NIND
      NOUT = NOUTD
      GO TO 9005
    5 NIND = NIN
      GO TO 9005
   10 NOUTD = NOUT
 9005 RETURN
      END
C   IMSL ROUTINE NAME   - USPKD
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - VAX/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES THAT HAVE
C                           CHARACTER STRING ARGUMENTS
C
C   USAGE               - CALL USPKD  (PACKED,NCHARS,UNPAKD,NCHMTB)
C
C   ARGUMENTS    PACKED - CHARACTER STRING TO BE UNPACKED.(INPUT)
C                NCHARS - LENGTH OF PACKED. (INPUT)  SEE REMARKS.
C                UNPAKD - INTEGER ARRAY TO RECEIVE THE UNPACKED
C                         REPRESENTATION OF THE STRING. (OUTPUT)
C                NCHMTB - NCHARS MINUS TRAILING BLANKS. (OUTPUT)
C
C   PRECISION/HARDWARE  - SINGLE/ALL
C
C   REQD. IMSL ROUTINES - NONE
C
C   REMARKS  1.  USPKD UNPACKS A CHARACTER STRING INTO AN INTEGER ARRAY
C                IN (A1) FORMAT.
C            2.  UP TO 129 CHARACTERS MAY BE USED.  ANY IN EXCESS OF
C                THAT ARE IGNORED.
C
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
      SUBROUTINE USPKD  (PACKED,NCHARS,UNPAKD,NCHMTB)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NC,NCHARS,NCHMTB
C
      INTEGER            UNPAKD(1),IBLANK
      INTEGER            PACKED(1)
      DATA               IBLANK /1H /
C                                  INITIALIZE NCHMTB
      NCHMTB = 0
C                                  RETURN IF NCHARS IS LE ZERO
      IF(NCHARS.LE.0) RETURN
C                                  SET NC=NUMBER OF CHARS TO BE DECODED
      NC = MIN0 (129,NCHARS)
C      DECODE (NC,150,PACKED) (UNPAKD(I),I=1,NC)
      WRITE(PACKED(1),150) (UNPAKD(I),I=1,NC)
  150 FORMAT (129A1)
C                                  CHECK UNPAKD ARRAY AND SET NCHMTB
C                                  BASED ON TRAILING BLANKS FOUND
      DO 200 N = 1,NC
         NN = NC - N + 1
         IF(UNPAKD(NN) .NE. IBLANK) GO TO 210
  200 CONTINUE
  210 NCHMTB = NN
      RETURN
      END
C       RUTINA MODIFICADA POR MARIO FOGLIO,
C
C   IMSL ROUTINE NAME   - ZRPOLY
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - VAX/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ZEROS OF A POLYNOMIAL WITH REAL
C                           COEFFICIENTS (JENKINS-TRAUB)
C
C   USAGE               - CALL ZRPOLY (A,NDEG,Z,IER)
C
C   ARGUMENTS    A      - INPUT REAL VECTOR OF LENGTH NDEG+1
C                           CONTAINING THE COEFFICIENTS IN ORDER OF
C                           DECREASING POWERS OF THE VARIABLE.
C                NDEG   - INPUT INTEGER DEGREE OF POLYNOMIAL.
C                           NDEG MUST BE GREATER THAN 0 AND LESS
C                           THAN 101.
C                Z      - OUTPUT COMPLEX VECTOR OF LENGTH NDEG
C                           CONTAINING THE COMPUTED ROOTS OF THE
C                           POLYNOMIAL.
C                         NOTE - THE ROUTINE TREATS Z AS A REAL VECTOR
C                           OF LENGTH 2*NDEG. AN APPROPRIATE
C                           EQUIVALENCE STATEMENT MAY BE REQUIRED.
C                           SEE DOCUMENT EXAMPLE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129, INDICATES THAT THE DEGREE OF THE
C                             POLYNOMIAL IS GREATER THAN 100 OR LESS
C                             THAN 1.
C                           IER=130, INDICATES THAT THE LEADING
C                             COEFFICIENT IS ZERO.
C                           IER=131, INDICATES THAT ZRPOLY FOUND FEWER
C                             THAN NDEG ZEROS. IF ONLY M ZEROS ARE
C                             FOUND, Z(J),J=M+1,...,NDEG ARE SET TO
C                             POSITIVE MACHINE INFINITY.
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO,ZRPQLB,ZRPQLC,ZRPQLD,ZRPQLE,
C                           ZRPQLF,ZRPQLG,ZRPQLH,ZRPQLI
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
        
        
        
        
***********************************************************************
***********************************************************************
C     IMSL ROUTINE NAME   - ZRPOLY
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - VAX/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ZEROS OF A POLYNOMIAL WITH REAL
C                           COEFFICIENTS (JENKINS-TRAUB)
C
C   USAGE               - CALL ZRPOLY (A,NDEG,Z,IER)
C
C   ARGUMENTS    A      - INPUT REAL VECTOR OF LENGTH NDEG+1
C                           CONTAINING THE COEFFICIENTS IN ORDER OF
C                           DECREASING POWERS OF THE VARIABLE.
C                NDEG   - INPUT INTEGER DEGREE OF POLYNOMIAL.
C                           NDEG MUST BE GREATER THAN 0 AND LESS
C                           THAN 101.
C                Z      - OUTPUT COMPLEX VECTOR OF LENGTH NDEG
C                           CONTAINING THE COMPUTED ROOTS OF THE
C                           POLYNOMIAL.
C                         NOTE - THE ROUTINE TREATS Z AS A REAL VECTOR
C                           OF LENGTH 2*NDEG. AN APPROPRIATE
C                           EQUIVALENCE STATEMENT MAY BE REQUIRED.
C                           SEE DOCUMENT EXAMPLE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129, INDICATES THAT THE DEGREE OF THE
C                             POLYNOMIAL IS GREATER THAN 100 OR LESS
C                             THAN 1.
C                           IER=130, INDICATES THAT THE LEADING
C                             COEFFICIENT IS ZERO.
C                           IER=131, INDICATES THAT ZRPOLY FOUND FEWER
C                             THAN NDEG ZEROS. IF ONLY M ZEROS ARE
C                             FOUND, Z(J),J=M+1,...,NDEG ARE SET TO
C                             POSITIVE MACHINE INFINITY.
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO,ZRPQLB,ZRPQLC,ZRPQLD,ZRPQLE,
C                           ZRPQLF,ZRPQLG,ZRPQLH,ZRPQLI
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPOLY (A,NDEG,Z,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NDEG,IER
      DOUBLE PRECISION   A(1),Z(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,NN,J,JJ,I,NM1,ICNT,N2,L,NZ,NPI
            REAL*8               ARE,RMRE,ETA
      REAL               RINFP,REPSP,RADIX,RLO,XX,YY,SINR,
     1                   COSR,RMAX,RMIN,X,SC,XM,FF,DX,DF,BND,XXX
CHANGE      REAL               ETA,RMRE,RINFP,REPSP,RADIX,RLO,XX,YY,SINR,
CHANGE     1                   COSR,RMAX,RMIN,X,SC,XM,FF,DX,DF,BND,XXX,ARE
      REAL               PT(101)
      DOUBLE PRECISION   TEMP(101),P(101),QP(101),RK(101),QK(101),
     1                   SVK(101)
      DOUBLE PRECISION   SR,SI,U,V,RA,RB,C,D,A1,A2,A3,
     1                   A6,A7,E,F,G,H,SZR,SZI,RLZR,RLZI,
     2                   T,AA,BB,CC,FACTOR,REPSR1,ZERO,ONE,FN
      LOGICAL            ZEROK
          INTEGER OUT,DAT,DT1,DAD
        PARAMETER(DAT=4,OUT=3,DT1=7,DAD=5)
      COMMON /ZRPQLJ/    P,QP,RK,QK,SVK,SR,SI,U,V,RA,RB,C,D,A1,A2,A3,A6,
     1                   A7,E,F,G,H,SZR,SZI,RLZR,RLZI,ETA,ARE,RMRE,N,NN
C                                  THE FOLLOWING STATEMENTS SET MACHINE
C                                    CONSTANTS USED IN VARIOUS PARTS OF
C                                    THE PROGRAM. THE MEANING OF THE
C                                    FOUR CONSTANTS ARE - REPSR1 THE
C                                    MAXIMUM RELATIVE REPRESENTATION
C                                    ERROR WHICH CAN BE DESCRIBED AS
C                                    THE SMALLEST POSITIVE FLOATING
C                                    POINT NUMBER SUCH THAT 1.+REPSR1 IS
C                                    GREATER THAN 1
C                                  RINFP THE LARGEST FLOATING-POINT
C                                    NUMBER
C                                  REPSP THE SMALLEST POSITIVE
C                                    FLOATING-POINT NUMBER IF THE
C                                    EXPONENT RANGE DIFFERS IN SINGLE
C                                    AND DOUBLE PRECISION THEN REPSP
C                                    AND RINFP SHOULD INDICATE THE
C                                    SMALLER RANGE
C                                  RADIX THE BASE OF THE FLOATING-POINT
C                                    NUMBER SYSTEM USED
C
C       Para evitar uso de data
C
        IF(INIC.EQ.1)  GO TO 101
        INIC = 1
C       Cambie valores por los indicados en el manual de FORTRAN 5
        RINFP = 3.402823E+38
CHANGEHOME        RINFP = HUGE(RINFP)
C        RINFP = HUGE(RINFP)
        REPSP = 1.175495E-38
CHANGEHOME        REPSP = TINY(REPSP)
C        REPSP = TINY(REPSP)
        RADIX = 2.0
CHANGEHOME        REPSR1 = 0.5960465d-7
C       En precision simple, para el SID 501 obtuve 0.596047E-7
C       Pero, hay que usar precision doble, y se tiene  0.111076513d-15
        REPSR1 = 0.11076513D-15
CHANGEHOME        REPSR1 = EPSILON(REPSR1)
C        REPSR1 = EPSILON(REPSR1)
        ZERO = 0.0D0
        ONE = 1.0D0
C        WRITE(OUT,1101) RINFP,REPSP,REPSR1
C 1101 FORMAT(' PUEDE ALTERAR RINFP = ',E9.2,/,'  REPSP = ',E12.5,/,
C     1  '  REPSR1 = ',D17.10)
C        READ(*,*) RINFP,REPSP,REPSR1
        WRITE(OUT,1102) RINFP,REPSP,REPSR1
 1102 FORMAT('  RINFP = ',E9.2,'  REPSP = ',E12.5,
     1  '  REPSR1 = ',D17.10)
  101   CONTINUE
C
C                                  ZRPOLY USES SINGLE PRECISION
C                                    CALCULATIONS FOR SCALING, BOUNDS
C                                    AND ERROR CALCULATIONS.
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (NDEG .GT. 100 .OR. NDEG .LT. 1) GO TO 165
      ETA = REPSR1
      ARE = ETA
      RMRE = ETA
      RLO = REPSP/ETA
C                                  INITIALIZATION OF CONSTANTS FOR
C                                    SHIFT ROTATION
      XX = .7071068
      YY = -XX
      SINR = .9975641
      COSR = -.06975647
      N = NDEG
      NN = N+1
C                                  ALGORITHM FAILS IF THE LEADING
C                                    COEFFICIENT IS ZERO.
      IF (A(1).NE.ZERO) GO TO 5
      IER = 130
      GO TO 9000
C                                  REMOVE THE ZEROS AT THE ORIGIN IF
C                                    ANY
    5 IF (A(NN).NE.ZERO) GO TO 10
      J = NDEG-N+1
      JJ = J+NDEG
      Z(J) = ZERO
      Z(JJ) = ZERO
      NN = NN-1
      N = N-1
      IF (NN.EQ.1) GO TO 9005
      GO TO 5
C                                  MAKE A COPY OF THE COEFFICIENTS
   10 DO 15 I=1,NN
         P(I) = A(I)
   15 CONTINUE
C                                  START THE ALGORITHM FOR ONE ZERO
   20 IF (N.GT.2) GO TO 30
      IF (N.LT.1) GO TO 9005
C                                  CALCULATE THE FINAL ZERO OR PAIR OF
C                                    ZEROS
      IF (N.EQ.2) GO TO 25
      Z(NDEG) = -P(2)/P(1)
      Z(NDEG+NDEG) = ZERO
      GO TO 145
   25 CALL ZRPQLI (P(1),P(2),P(3),Z(NDEG-1),Z(NDEG+NDEG-1),Z(NDEG),
     1   Z(NDEG+NDEG))
      GO TO 145
C                                  FIND LARGEST AND SMALLEST MODULI OF
C                                    COEFFICIENTS.
   30 RMAX = 0.
      RMIN = RINFP
      DO 35 I=1,NN
         X = ABS(SNGL(P(I)))
         IF (X.GT.RMAX) RMAX = X
         IF (X.NE.0..AND.X.LT.RMIN) RMIN = X
   35 CONTINUE
C                                  SCALE IF THERE ARE LARGE OR VERY
C                                    SMALL COEFFICIENTS COMPUTES A
C                                    SCALE FACTOR TO MULTIPLY THE
C                                    COEFFICIENTS OF THE POLYNOMIAL.
C                                    THE SCALING IS DONE TO AVOID
C                                    OVERFLOW AND TO AVOID UNDETECTED
C                                    UNDERFLOW INTERFERING WITH THE
C                                    CONVERGENCE CRITERION.
C                                  THE FACTOR IS A POWER OF THE BASE
      SC = RLO/RMIN
      IF (SC.GT.1.0) GO TO 40
      IF (RMAX.LT.10.) GO TO 55
      IF (SC.EQ.0.) SC = REPSP*RADIX*RADIX
      GO TO 45
   40 IF (RINFP/SC.LT.RMAX) GO TO 55
   45 L = ALOG(SC)/ALOG(RADIX)+.5
      IF (L .EQ. 0) GO TO 55
      FACTOR = DBLE(RADIX)**L
      DO 50 I=1,NN
   50 P(I) = FACTOR*P(I)
C                                  COMPUTE LOWER BOUND ON MODULI OF
C                                    ZEROS.
   55 DO 60 I=1,NN
   60 PT(I) = ABS(SNGL(P(I)))
      PT(NN) = -PT(NN)
C                                  COMPUTE UPPER ESTIMATE OF BOUND
      X = EXP((ALOG(-PT(NN))-ALOG(PT(1)))/N)
      IF (PT(N).EQ.0.) GO TO 65
C                                  IF NEWTON STEP AT THE ORIGIN IS
C                                    BETTER, USE IT.
      XM = -PT(NN)/PT(N)
      IF (XM.LT.X) X = XM
C                                  CHOP THE INTERVAL (0,X) UNTIL FF.LE.0
   65 XM = X*.1
      FF = PT(1)
      DO 70 I=2,NN
   70 FF = FF*XM+PT(I)
      IF (FF.LE.0.) GO TO 75
      X = XM
      GO TO 65
   75 DX = X
C                                  DO NEWTON ITERATION UNTIL X
C                                    CONVERGES TO TWO DECIMAL PLACES
   80 IF (ABS(DX/X).LE..005) GO TO 90
      FF = PT(1)
      DF = FF
      DO 85 I=2,N
         FF = FF*X+PT(I)
         DF = DF*X+FF
   85 CONTINUE
      FF = FF*X+PT(NN)
      DX = FF/DF
      X = X-DX
      GO TO 80
   90 BND = X
C                                  COMPUTE THE DERIVATIVE AS THE INTIAL
C                                    K POLYNOMIAL AND DO 5 STEPS WITH
C                                    NO SHIFT
      NM1 = N-1
      FN = ONE/N
      DO 95 I=2,N
   95 RK(I) = (NN-I)*P(I)*FN
      RK(1) = P(1)
      AA = P(NN)
      BB = P(N)
      ZEROK = RK(N).EQ.ZERO
      DO 115 JJ=1,5
         CC = RK(N)
         IF (ZEROK) GO TO 105
C                                  USE SCALED FORM OF RECURRENCE IF
C                                    VALUE OF K AT 0 IS NONZERO
         T = -AA/CC
         DO 100 I=1,NM1
            J = NN-I
            RK(J) = T*RK(J-1)+P(J)
  100    CONTINUE
         RK(1) = P(1)
         ZEROK = DABS(RK(N)).LE.DABS(BB)*ETA*10.
         GO TO 115
C                                  USE UNSCALED FORM OF RECURRENCE
  105    DO 110 I=1,NM1
            J = NN-I
            RK(J) = RK(J-1)
  110    CONTINUE
         RK(1) = ZERO
         ZEROK = RK(N).EQ.ZERO
  115 CONTINUE
C                                  SAVE K FOR RESTARTS WITH NEW SHIFTS
      DO 120 I=1,N
  120 TEMP(I) = RK(I)
C                                  LOOP TO SELECT THE QUADRATIC
C                                    CORRESPONDING TO EACH NEW SHIFT
      DO 140 ICNT=1,20
C                                  QUADRATIC CORRESPONDS TO A DOUBLE
C                                    SHIFT TO A NON-REAL POINT AND ITS
C                                    COMPLEX CONJUGATE. THE POINT HAS
C                                    MODULUS BND AND AMPLITUDE ROTATED
C                                    BY 94 DEGREES FROM THE PREVIOUS
C                                    SHIFT
         XXX = COSR*XX-SINR*YY
         YY = SINR*XX+COSR*YY
         XX = XXX
         SR = BND*XX
         SI = BND*YY
         U = -SR-SR
         V = BND*BND
C                                  SECOND STAGE CALCULATION, FIXED
C                                    QUADRATIC
         CALL ZRPQLB (20*ICNT,NZ)
         IF (NZ.EQ.0) GO TO 130
C                                  THE SECOND STAGE JUMPS DIRECTLY TO
C                                    ONE OF THE THIRD STAGE ITERATIONS
C                                    AND RETURNS HERE IF SUCCESSFUL.
C                                  DEFLATE THE POLYNOMIAL, STORE THE
C                                    ZERO OR ZEROS AND RETURN TO THE
C                                    MAIN ALGORITHM.
         J = NDEG-N+1
         JJ = J+NDEG
         Z(J) = SZR
         Z(JJ) = SZI
         NN = NN-NZ
         N = NN-1
         DO 125 I=1,NN
  125    P(I) = QP(I)
         IF (NZ.EQ.1) GO TO 20
         Z(J+1) = RLZR
         Z(JJ+1) = RLZI
         GO TO 20
C                                  IF THE ITERATION IS UNSUCCESSFUL
C                                    ANOTHER QUADRATIC IS CHOSEN AFTER
C                                    RESTORING K
  130    DO 135 I=1,N
  135    RK(I) = TEMP(I)
  140 CONTINUE
C                                  RETURN WITH FAILURE IF NO
C                                    CONVERGENCE WITH 20 SHIFTS
      IER = 131
C                                  CONVERT ZEROS (Z) IN COMPLEX FORM
  145 DO 150 I=1,NDEG
         NPI= NDEG+I
         P(I) = Z(NPI)
  150 CONTINUE
      N2 = NDEG+NDEG
      J = NDEG
      DO 155 I=1,NDEG
         Z(N2-1) = Z(J)
         Z(N2) = P(J)
         N2 = N2-2
         J = J-1
  155 CONTINUE
      IF (IER .EQ. 0) GO TO 9005
C                                  SET UNFOUND ROOTS TO MACHINE INFINITY
      N2 = 2*(NDEG-NN)+3
      DO 160 I=1,N
         Z(N2) = RINFP
         Z(N2+1) = RINFP
         N2 = N2+2
  160 CONTINUE
      GO TO 9000
  165 IER = 129
 9000 CONTINUE
      CALL UERTST (IER,6HZRPOLY)
 9005 RETURN
      END
C   IMSL ROUTINE NAME   - ZRPQLB
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - VAX/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZRPOLY
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - ZRPQLC,ZRPQLD,ZRPQLE,ZRPQLF,ZRPQLG,ZRPQLH,
C                           ZRPQLI
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPQLB (L2,NZ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            L2,NZ
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,NN,J,ITYPE,I,IFLAG
CHANGE      REAL               ARE,BETAS,BETAV,ETA,OSS,OTS,OTV,OVV,RMRE,SS,
CHANGE     1                   TS,TSS,TV,TVV,VV
      REAL*8               ARE,ETA,RMRE
      REAL               BETAS,BETAV,OSS,OTS,OTV,OVV,SS,
     1                   TS,TSS,TV,TVV,VV
      DOUBLE PRECISION   P(101),QP(101),RK(101),QK(101),SVK(101)
      DOUBLE PRECISION   SR,SI,U,V,RA,RB,C,D,A1,A2,A3,
     1                   A6,A7,E,F,G,H,SZR,SZI,RLZR,RLZI,
     2                   SVU,SVV,UI,VI,S,ZERO
      LOGICAL            VPASS,SPASS,VTRY,STRY
      COMMON /ZRPQLJ/    P,QP,RK,QK,SVK,SR,SI,U,V,RA,RB,C,D,A1,A2,A3,A6,
     1                   A7,E,F,G,H,SZR,SZI,RLZR,RLZI,ETA,ARE,RMRE,N,NN
C      DATA               ZERO/0.0D0/
        ZERO = 0.0D0
C                                  FIRST EXECUTABLE STATEMENT
      NZ = 0
C                                  COMPUTES UP TO L2 FIXED SHIFT
C                                    K-POLYNOMIALS, TESTING FOR
C                                    CONVERGENCE IN THE LINEAR OR
C                                    QUADRATIC CASE. INITIATES ONE OF
C                                    THE VARIABLE SHIFT ITERATIONS AND
C                                    RETURNS WITH THE NUMBER OF ZEROS
C                                    FOUND.
C                                  L2 - LIMIT OF FIXED SHIFT STEPS
C                                  NZ -NUMBER OF ZEROS FOUND
      BETAV = .25
      BETAS = .25
      OSS = SR
      OVV = V
C                                  EVALUATE POLYNOMIAL BY SYNTHETIC
C                                    DIVISION
      CALL ZRPQLH (NN,U,V,P,QP,RA,RB)
      CALL ZRPQLE (ITYPE)
      DO 40 J=1,L2
C                                  CALCULATE NEXT K POLYNOMIAL AND
C                                    ESTIMATE V
         CALL ZRPQLF (ITYPE)
         CALL ZRPQLE (ITYPE)
         CALL ZRPQLG (ITYPE,UI,VI)
         VV = VI
C                                  ESTIMATE S
         SS = 0.
         IF (RK(N).NE.ZERO) SS = -P(NN)/RK(N)
         TV = 1.
         TS = 1.
         IF (J.EQ.1.OR.ITYPE.EQ.3) GO TO 35
C                                  COMPUTE RELATIVE MEASURES OF
C                                    CONVERGENCE OF S AND V SEQUENCES
         IF (VV.NE.0.) TV = ABS((VV-OVV)/VV)
         IF (SS.NE.0.) TS = ABS((SS-OSS)/SS)
C                                  IF DECREASING, MULTIPLY TWO MOST
C                                    RECENT CONVERGENCE MEASURES
         TVV = 1.
         IF (TV.LT.OTV) TVV = TV*OTV
         TSS = 1.
         IF (TS.LT.OTS) TSS = TS*OTS
C                                  COMPARE WITH CONVERGENCE CRITERIA
         VPASS = TVV.LT.BETAV
         SPASS = TSS.LT.BETAS
         IF (.NOT.(SPASS.OR.VPASS)) GO TO 35
C                                  AT LEAST ONE SEQUENCE HAS PASSED THE
C                                    CONVERGENCE TEST. STORE VARIABLES
C                                    BEFORE ITERATING
         SVU = U
         SVV = V
         DO 5 I=1,N
    5    SVK(I) = RK(I)
         S = SS
C                                  CHOOSE ITERATION ACCORDING TO THE
C                                    FASTEST CONVERGING SEQUENCE
         VTRY = .FALSE.
         STRY = .FALSE.
         IF (SPASS.AND.((.NOT.VPASS).OR.TSS.LT.TVV)) GO TO 20
   10    CALL ZRPQLC (UI,VI,NZ)
         IF (NZ.GT.0) RETURN
C                                  QUADRATIC ITERATION HAS FAILED. FLAG
C                                    THAT IT HAS BEEN TRIED AND
C                                    DECREASE THE CONVERGENCE
C                                    CRITERION.
         VTRY = .TRUE.
         BETAV = BETAV*.25
C                                  TRY LINEAR ITERATION IF IT HAS NOT
C                                    BEEN TRIED AND THE S SEQUENCE IS
C                                    CONVERGING
         IF (STRY.OR.(.NOT.SPASS)) GO TO 25
         DO 15 I=1,N
   15    RK(I) = SVK(I)
   20    CALL ZRPQLD (S,NZ,IFLAG)
         IF (NZ.GT.0) RETURN
C                                  LINEAR ITERATION HAS FAILED. FLAG
C                                    THAT IT HAS BEEN TRIED AND
C                                    DECREASE THE CONVERGENCE CRITERION
         STRY = .TRUE.
         BETAS = BETAS*.25
         IF (IFLAG.EQ.0) GO TO 25
C                                  IF LINEAR ITERATION SIGNALS AN
C                                    ALMOST DOUBLE REAL ZERO ATTEMPT
C                                    QUADRATIC INTERATION
         UI = -(S+S)
         VI = S*S
         GO TO 10
C                                  RESTORE VARIABLES
   25    U = SVU
         V = SVV
         DO 30 I=1,N
   30    RK(I) = SVK(I)
C                                  TRY QUADRATIC ITERATION IF IT HAS
C                                    NOT BEEN TRIED AND THE V SEQUENCE
C                                    IS CONVERGING
         IF (VPASS.AND.(.NOT.VTRY)) GO TO 10
C                                  RECOMPUTE QP AND SCALAR VALUES TO
C                                    CONTINUE THE SECOND STAGE
         CALL ZRPQLH (NN,U,V,P,QP,RA,RB)
         CALL ZRPQLE (ITYPE)
   35    OVV = VV
         OSS = SS
         OTV = TV
         OTS = TS
   40 CONTINUE
      RETURN
      END
C   IMSL ROUTINE NAME   - ZRPQLC
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - VAX/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZRPOLY
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - ZRPQLE,ZRPQLF,ZRPQLG,ZRPQLH,ZRPQLI
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPQLC (UU,VV,NZ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NZ
      DOUBLE PRECISION   UU,VV
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,NN,J,I,ITYPE
CHANGE      REAL               ARE,EE,ETA,OMP,RELSTP,RMP,RMRE,T,ZM
      REAL*8               ARE,EE,ETA,OMP,RELSTP,RMP,RMRE,T,ZM
      DOUBLE PRECISION   P(101),QP(101),RK(101),QK(101),SVK(101)
      DOUBLE PRECISION   SR,SI,U,V,RA,RB,C,D,A1,A2,A3,
     1                   A6,A7,E,F,G,H,SZR,SZI,RLZR,RLZI,
     2                   UI,VI,ZERO,PT01,ONE
      LOGICAL            TRIED
      COMMON /ZRPQLJ/    P,QP,RK,QK,SVK,SR,SI,U,V,RA,RB,C,D,A1,A2,A3,A6,
     1                   A7,E,F,G,H,SZR,SZI,RLZR,RLZI,ETA,ARE,RMRE,N,NN
CARAJO
CHANGE        REAL*8 DEE,DRMP
C      DATA               ZERO,PT01,ONE/0.0D0,0.01D0,1.0D0/
        ZERO = 0.0D0
        PT01 = 0.01D0
        ONE = 1.0D0
C                                  FIRST EXECUTABLE STATEMENT
      NZ = 0
C                                  VARIABLE-SHIFT K-POLYNOMIAL
C                                    ITERATION FOR A QUADRATIC FACTOR
C                                    CONVERGES ONLY IF THE ZEROS ARE
C                                    EQUIMODULAR OR NEARLY SO
C                                  UU,VV - COEFFICIENTS OF STARTING
C                                    QUADRATIC
C                                  NZ - NUMBER OF ZERO FOUND
      TRIED = .FALSE.
      U = UU
      V = VV
      J = 0
C                                  MAIN LOOP
    5 CALL ZRPQLI (ONE,U,V,SZR,SZI,RLZR,RLZI)
C                                  RETURN IF ROOTS OF THE QUADRATIC ARE
C                                    REAL AND NOT CLOSE TO MULTIPLE OR
C                                    NEARLY EQUAL AND OF OPPOSITE SIGN
      IF ( DABS(DABS(SZR)-DABS(RLZR)).GT.PT01*DABS(RLZR)) RETURN
C                                  EVALUATE POLYNOMIAL BY QUADRATIC
C                                    SYNTHETIC DIVISION
      CALL ZRPQLH (NN,U,V,P,QP,RA,RB)
      RMP = DABS(RA-SZR*RB)+DABS(SZI*RB)
C                                  COMPUTE A RIGOROUS BOUND ON THE
C                                    ROUNDING ERROR IN EVALUTING P
CHANGE      ZM = SQRT(ABS(SNGL(V)))
CHANGE      EE = 2.*ABS(SNGL(QP(1)))
      ZM = DSQRT(DABS(V))
      EE = 2.D0*DABS(QP(1))
      T = -SZR*RB
      DO 10 I=2,N
CHANGE   10 EE = EE*ZM+ABS(SNGL(QP(I)))
CHANGE      EE = EE*ZM+ABS(SNGL(RA)+T)
CHANGE      EE = (5.*RMRE+4.*ARE)*EE-(5.*RMRE+2.*ARE)*(ABS(SNGL(RA)+T)+
CHANGE     1     ABS(SNGL(RB))*ZM)+2.*ARE*ABS(T)
   10 EE = EE*ZM+DABS(QP(I))
      EE = EE*ZM+DABS(RA+T)
      EE = (5.D0*RMRE+4.D0*ARE)*EE-(5.D0*RMRE+2.D0*ARE)*(DABS(RA+T)+
     1     DABS(RB)*ZM)+2.D0*ARE*DABS(T)
C                                  ITERATION HAS CONVERGED SUFFICIENTLY
C                                    IF THE POLYNOMIAL VALUE IS LESS
C                                    THAN 20 TIMES THIS BOUND
CARAJO
C        DRMP = DFLOAT(RMP)
C        DEE = DFLOAT(EE)
      IF (RMP.GT.20.D0*EE) GO TO 15
C
      NZ = 2
      RETURN
   15 J = J+1
C                                  STOP ITERATION AFTER 20 STEPS
      IF (J.GT.20) RETURN
      IF (J.LT.2) GO TO 25
      IF (RELSTP.GT.1.D-2.OR.RMP.LT.OMP.OR.TRIED) GO TO 25
C                                  A CLUSTER APPEARS TO BE STALLING THE
C                                    CONVERGENCE. FIVE FIXED SHIFT
C                                    STEPS ARE TAKEN WITH A U,V CLOSE
C                                    TO THE CLUSTER
      IF (RELSTP.LT.ETA) RELSTP = ETA
      RELSTP = SQRT(RELSTP)
      U = U-U*RELSTP
      V = V+V*RELSTP
      CALL ZRPQLH (NN,U,V,P,QP,RA,RB)
      DO 20 I=1,5
         CALL ZRPQLE (ITYPE)
         CALL ZRPQLF (ITYPE)
   20 CONTINUE
      TRIED = .TRUE.
      J = 0
   25 OMP = RMP
C                                  CALCULATE NEXT K POLYNOMIAL AND NEW
C                                    U AND V
      CALL ZRPQLE (ITYPE)
      CALL ZRPQLF (ITYPE)
      CALL ZRPQLE (ITYPE)
      CALL ZRPQLG (ITYPE,UI,VI)
C                                  IF VI IS ZERO THE ITERATION IS NOT
C                                    CONVERGING
      IF (VI.EQ.ZERO) RETURN
      RELSTP = DABS((VI-V)/VI)
      U = UI
      V = VI
      GO TO 5
      END
C   IMSL ROUTINE NAME   - ZRPQLD
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - VAX/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZRPOLY
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPQLD (SSS,NZ,IFLAG)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NZ,IFLAG
      DOUBLE PRECISION   SSS
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,NN,J,I
CHANGE      REAL               ARE,EE,ETA,OMP,RMP,RMS,RMRE
      REAL*8               ARE,EE,ETA,OMP,RMP,RMS,RMRE
      DOUBLE PRECISION   P(101),QP(101),RK(101),QK(101),SVK(101)
      DOUBLE PRECISION   SR,SI,U,V,RA,RB,C,D,A1,A2,A3,
     1                   A6,A7,E,F,G,H,SZR,SZI,RLZR,RLZI,
     2                   PV,RKV,T,S,ZERO,PT001
      COMMON /ZRPQLJ/    P,QP,RK,QK,SVK,SR,SI,U,V,RA,RB,C,D,A1,A2,A3,A6,
     1                   A7,E,F,G,H,SZR,SZI,RLZR,RLZI,ETA,ARE,RMRE,N,NN
CARAJO
CHANGE        REAL*8  DEE,DRMP,DRMRE
C      DATA               ZERO/0.0D0/,PT001/0.001D0/
        ZERO = 0.0D0
        PT001 = 0.001D0
C                                  VARIABLE-SHIFT H POLYNOMIAL
C                                    ITERATION FOR A REAL ZERO SSS -
C                                    STARTING ITERATE
C                                  NZ - NUMBER OF ZERO FOUND
C                                  IFLAG - FLAG TO INDICATE A PAIR OF
C                                    ZEROS NEAR REAL AXIS
C                                  FIRST EXECUTABLE STATEMENT
      NZ = 0
      S = SSS
      IFLAG = 0
      J = 0
C                                  MAIN LOOP
    5 PV = P(1)
C                                  EVALUATE P AT S
      QP(1) = PV
      DO 10 I=2,NN
         PV = PV*S+P(I)
         QP(I) = PV
   10 CONTINUE
      RMP = DABS(PV)
C                                  COMPUTE A RIGOROUS BOUND ON THE
C                                    ERROR IN EVALUATING P
      RMS = DABS(S)
CHANGE      EE = (RMRE/(ARE+RMRE))*ABS(SNGL(QP(1)))
      EE = (RMRE/(ARE+RMRE))*DABS(QP(1))
      DO 15 I=2,NN
CHANGE   15 EE = EE*RMS+ABS(SNGL(QP(I)))
   15 EE = EE*RMS+DABS(QP(I))
C                                  ITERATION HAS CONVERGED SUFFICIENTLY
C                                    IF THE POLYNOMIAL VALUE IS LESS
C                                    THAN 20 TIMES THIS BOUND
CARAJO
C        DRMP = DFLOAT(RMP)
C        DEE = DFLOAT(EE)
C        DRMRE = DFLOAT(RMRE)
CHANGE      IF (DRMP.GT.20.D0*((ARE+RMRE)*DEE-DRMRE*DRMP)) GO TO 20
      IF (RMP.GT.20.D0*((ARE+RMRE)*EE-RMRE*RMP)) GO TO 20
C
      NZ = 1
      SZR = S
      SZI = ZERO
      RETURN
   20 J = J+1
C                                  STOP ITERATION AFTER 10 STEPS
      IF (J.GT.10) RETURN
      IF (J.LT.2) GO TO 25
      IF (DABS(T).GT.PT001*DABS(S-T).OR.RMP.LE.OMP) GO TO 25
C                                  A CLUSTER OF ZEROS NEAR THE REAL
C                                    AXIS HAS BEEN ENCOUNTERED RETURN
C                                    WITH IFLAG SET TO INITIATE A
C                                    QUADRATIC ITERATION
      IFLAG = 1
      SSS = S
      RETURN
C                                  RETURN IF THE POLYNOMIAL VALUE HAS
C                                    INCREASED SIGNIFICANTLY
   25 OMP = RMP
C                                  COMPUTE T, THE NEXT POLYNOMIAL, AND
C                                    THE NEW ITERATE
      RKV = RK(1)
      QK(1) = RKV
      DO 30 I=2,N
         RKV = RKV*S+RK(I)
         QK(I) = RKV
   30 CONTINUE
      IF (DABS(RKV).LE.DABS(RK(N))*10.D0*ETA) GO TO 40
C                                  USE THE SCALED FORM OF THE
C                                    RECURRENCE IF THE VALUE OF K AT S
C                                    IS NONZERO
      T = -PV/RKV
      RK(1) = QP(1)
      DO 35 I=2,N
   35 RK(I) = T*QK(I-1)+QP(I)
      GO TO 50
C                                  USE UNSCALED FORM
   40 RK(1) = ZERO
      DO 45 I=2,N
   45 RK(I) = QK(I-1)
   50 RKV = RK(1)
      DO 55 I=2,N
   55 RKV = RKV*S+RK(I)
      T = ZERO
      IF (DABS(RKV).GT.DABS(RK(N))*10.D0*ETA) T = -PV/RKV
      S = S+T
      GO TO 5
      END
C   IMSL ROUTINE NAME   - ZRPQLE
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - VAX/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZRPOLY
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - ZRPQLH
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPQLE (ITYPE)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            ITYPE
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,NN
CARAJO
CHANGE      REAL               ARE,ETA,RMRE
      REAL*8               ARE,ETA,RMRE
      DOUBLE PRECISION   P(101),QP(101),RK(101),QK(101),SVK(101)
      DOUBLE PRECISION   SR,SI,U,V,RA,RB,C,D,A1,A2,A3,
     1                   A6,A7,E,F,G,H,SZR,SZI,RLZR,RLZI
      COMMON /ZRPQLJ/    P,QP,RK,QK,SVK,SR,SI,U,V,RA,RB,C,D,A1,A2,A3,A6,
     1                   A7,E,F,G,H,SZR,SZI,RLZR,RLZI,ETA,ARE,RMRE,N,NN
C                                  THIS ROUTINE CALCULATES SCALAR
C                                    QUANTITIES USED TO COMPUTE THE
C                                    NEXT K POLYNOMIAL AND NEW
C                                    ESTIMATES OF THE QUADRATIC
C                                    COEFFICIENTS
C                                  ITYPE - INTEGER VARIABLE SET HERE
C                                    INDICATING HOW THE CALCULATIONS
C                                    ARE NORMALIZED TO AVOID OVERFLOW
C                                  SYNTHETIC DIVISION OF K BY THE
C                                    QUADRATIC 1,U,V
C                                  FIRST EXECUTABLE STATEMENT
      CALL ZRPQLH (N,U,V,RK,QK,C,D)
      IF (DABS(C).GT.DABS(RK(N))*100.*ETA) GO TO 5
      IF (DABS(D).GT.DABS(RK(N-1))*100.*ETA) GO TO 5
      ITYPE = 3
C                                  TYPE=3 INDICATES THE QUADRATIC IS
C                                    ALMOST A FACTOR OF K
      RETURN
    5 IF (DABS(D).LT.DABS(C)) GO TO 10
      ITYPE = 2
C                                  TYPE=2 INDICATES THAT ALL FORMULAS
C                                    ARE DIVIDED BY D
      E = RA/D
      F = C/D
      G = U*RB
      H = V*RB
      A3 = (RA+G)*E+H*(RB/D)
      A1 = RB*F-RA
      A7 = (F+U)*RA+H
      RETURN
   10 ITYPE = 1
C                                  TYPE=1 INDICATES THAT ALL FORMULAS
C                                    ARE DIVIDED BY C
      E = RA/C
      F = D/C
      G = U*E
      H = V*RB
      A3 = RA*E+(H/C+G)*RB
      A1 = RB-RA*(D/C)
      A7 = RA+G*D+H*F
      RETURN
      END
C   IMSL ROUTINE NAME   - ZRPQLF
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - VAX/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZRPOLY
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPQLF (ITYPE)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            ITYPE
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,NN,I
CHANGE      REAL               ARE,ETA,RMRE
      REAL*8               ARE,ETA,RMRE
      DOUBLE PRECISION   P(101),QP(101),RK(101),QK(101),SVK(101)
      DOUBLE PRECISION   SR,SI,U,V,RA,RB,C,D,A1,A2,A3,
     1                   A6,A7,E,F,G,H,SZR,SZI,RLZR,RLZI,TEMP,ZERO
      COMMON /ZRPQLJ/    P,QP,RK,QK,SVK,SR,SI,U,V,RA,RB,C,D,A1,A2,A3,A6,
     1                   A7,E,F,G,H,SZR,SZI,RLZR,RLZI,ETA,ARE,RMRE,N,NN
C      DATA               ZERO/0.0D0/
        ZERO = 0.0D0
C                                  COMPUTES THE NEXT K POLYNOMIALS
C                                    USING SCALARS COMPUTED IN ZRPQLE
C                                  FIRST EXECUTABLE STATEMENT
      IF (ITYPE.EQ.3) GO TO 20
      TEMP = RA
      IF (ITYPE.EQ.1) TEMP = RB
      IF (DABS(A1).GT.DABS(TEMP)*ETA*10.) GO TO 10
C                                  IF A1 IS NEARLY ZERO THEN USE A
C                                    SPECIAL FORM OF THE RECURRENCE
      RK(1) = ZERO
      RK(2) = -A7*QP(1)
      DO 5 I=3,N
    5 RK(I) = A3*QK(I-2)-A7*QP(I-1)
      RETURN
C                                  USE SCALED FORM OF THE RECURRENCE
   10 A7 = A7/A1
      A3 = A3/A1
      RK(1) = QP(1)
      RK(2) = QP(2)-A7*QP(1)
      DO 15 I=3,N
   15 RK(I) = A3*QK(I-2)-A7*QP(I-1)+QP(I)
      RETURN
C                                  USE UNSCALED FORM OF THE RECURRENCE
C                                    IF TYPE IS 3
   20 RK(1) = ZERO
      RK(2) = ZERO
      DO 25 I=3,N
   25 RK(I) = QK(I-2)
      RETURN
      END
C   IMSL ROUTINE NAME   - ZRPQLG
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - VAX/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZRPOLY
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPQLG (ITYPE,UU,VV)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            ITYPE
      DOUBLE PRECISION   UU,VV
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,NN
CHANGE      REAL               ARE,ETA,RMRE
      REAL*8               ARE,ETA,RMRE
      DOUBLE PRECISION   P(101),QP(101),RK(101),QK(101),SVK(101)
      DOUBLE PRECISION   SR,SI,U,V,RA,RB,C,D,A1,A2,A3,
     1                   A6,A7,E,F,G,H,SZR,SZI,RLZR,RLZI,
     2                   A4,A5,B1,B2,C1,C2,C3,C4,TEMP,ZERO
      COMMON /ZRPQLJ/    P,QP,RK,QK,SVK,SR,SI,U,V,RA,RB,C,D,A1,A2,A3,A6,
     1                   A7,E,F,G,H,SZR,SZI,RLZR,RLZI,ETA,ARE,RMRE,N,NN
C      DATA               ZERO/0.0D0/
        ZERO = 0.0D0
C                                  COMPUTE NEW ESTIMATES OF THE
C                                    QUADRATIC COEFFICIENTS USING THE
C                                    SCALARS COMPUTED IN ZRPQLE
C                                  USE FORMULAS APPROPRIATE TO SETTING
C                                    OF TYPE.
C                                  FIRST EXECUTABLE STATEMENT
      IF (ITYPE.EQ.3) GO TO 15
      IF (ITYPE.EQ.2) GO TO 5
      A4 = RA+U*RB+H*F
      A5 = C+(U+V*F)*D
      GO TO 10
    5 A4 = (RA+G)*F+H
      A5 = (F+U)*C+V*D
C                                  EVALUATE NEW QUADRATIC COEFFICIENTS.
C
   10 B1 = -RK(N)/P(NN)
      B2 = -(RK(N-1)+B1*P(N))/P(NN)
      C1 = V*B2*A1
      C2 = B1*A7
      C3 = B1*B1*A3
      C4 = C1-C2-C3
      TEMP = A5+B1*A4-C4
      IF (TEMP.EQ.ZERO) GO TO 15
      UU = U-(U*(C3+C2)+V*(B1*A1+B2*A7))/TEMP
      VV = V*(1+C4/TEMP)
      RETURN
C                                  IF TYPE=3 THE QUADRATIC IS ZEROED
   15 UU = ZERO
      VV = ZERO
      RETURN
      END
C   IMSL ROUTINE NAME   - ZRPQLH
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - VAX/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZRPOLY
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPQLH (NN,U,V,P,Q,RA,RB)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NN
      DOUBLE PRECISION   P(NN),Q(NN),U,V,RA,RB
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      DOUBLE PRECISION   C
C                                  DIVIDES P BY THE QUADRATIC 1,U,V
C                                    PLACING THE QUOTIENT IN Q AND THE
C                                    REMAINDER IN A,B
C                                  FIRST EXECUTABLE STATEMENT
      RB = P(1)
      Q(1) = RB
      RA = P(2)-U*RB
      Q(2) = RA
      DO 5 I=3,NN
         C = P(I)-U*RA-V*RB
         Q(I) = C
         RB = RA
         RA = C
    5 CONTINUE
      RETURN
      END
C   IMSL ROUTINE NAME   - ZRPQLI
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - VAX/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZRPOLY
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPQLI (RA,B1,C,SR,SI,RLR,RLI)
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   RA,B1,C,SR,SI,RLR,RLI
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   RB,D,E,ZERO,ONE,TWO
C      DATA               ZERO,ONE,TWO/0.0D0,1.0D0,2.0D0/
        ZERO = 0.0D0
        ONE = 1.0D0
        TWO = 2.0D0
C                                  CALCULATE THE ZEROS OF THE QUADRATIC
C                                    A*Z**2 + B1*Z + C. THE QUADRATIC
C                                    FORMULA, MODIFIED TO AVOID
C                                    OVERFLOW, IS USED TO FIND THE
C                                    LARGER ZERO IF THE ZEROS ARE REAL
C                                    AND BOTH ZEROS ARE COMPLEX.
C                                  THE SMALLER REAL ZERO IS FOUND
C                                    DIRECTLY FROM THE PRODUCT OF THE
C                                    ZEROS C/A
C                                  FIRST EXECUTABLE STATEMENT
      IF (RA.NE.ZERO) GO TO 10
      SR = ZERO
      IF (B1.NE.ZERO) SR = -C/B1
      RLR = ZERO
    5 SI = ZERO
      RLI = ZERO
      RETURN
   10 IF (C.NE.ZERO) GO TO 15
      SR = ZERO
      RLR = -B1/RA
      GO TO 5
C                                  COMPUTE DISCRIMINANT AVOIDING
C                                    OVERFLOW
   15 RB = B1/TWO
      IF (DABS(RB).LT.DABS(C)) GO TO 20
      E = ONE-(RA/RB)*(C/RB)
      D = DSQRT(DABS(E))*DABS(RB)
      GO TO 25
   20 E = RA
      IF (C.LT.ZERO) E = -RA
      E = RB*(RB/DABS(C))-E
      D = DSQRT(DABS(E))*DSQRT(DABS(C))
   25 IF (E.LT.ZERO) GO TO 30
C                                  REAL ZEROS
      IF (RB.GE.ZERO) D = -D
      RLR = (-RB+D)/RA
      SR = ZERO
      IF (RLR.NE.ZERO) SR = (C/RLR)/RA
      GO TO 5
C                                  COMPLEX CONJUGATE ZEROS
   30 SR = -RB/RA
      RLR = SR
      SI = DABS(D/RA)
      RLI = -SI
      RETURN
      END






       
