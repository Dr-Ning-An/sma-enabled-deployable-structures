c ********************************************************************
c ********************************************************************
c
c     SMA_UM:  Shape Memory Alloy - Unified Model: 
c     User Material Subroutine for Thermomechanical 
c     Constitutive Model of Shape Memory Alloys
c
c     Version 4.1
c     February, 2003
c
c ********************************************************************
c
c     APPLICATIONS
c
c        3-D 
c        2-D plane strain and generalized plane strain
c        1-D problems
c
c ********************************************************************
c
c     AUTHORS
c
c        Dimitris C. Lagoudas
c        Zhonghe Bo
c        Muhammad A. Qidwai
c        Pavlin B. Entchev
c
c        Department of Aerospace Engineering
c        Texas A&M University
c        3141 TAMU
c        College Station, TX 77843-3141
c
c ********************************************************************
c ********************************************************************
c
c            
c     The following three sets of lines are written in the way required 
c     by ABAQUS 
c

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1RPL,DDSDDT,DRPLDE,DRPLDT,
     2STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KPST,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (MAXMTP=40000)

      CHARACTER*8 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      REAL*8 DSM,DSA,DT,ST,SD,TIME1,DTIME1,DA1,DA2,SDDT,RPLE,
     *       yumma,tumma,TEMPAL, props, props_um
      
      COMMON /DATA01/DSA(6,6),DSM(6,6),DA1(6,6),DA2(3),TEMPAL(MAXMTP)
      COMMON /DATA02/MAXMTP1,NELMTP,IMTP,ITER,KEYN,ISOLVE,KFLAG
      COMMON /DATA03/TIME1,DTIME1,TIME0,NUMMIN
      COMMON /PUNGA/ imli,MODEL
      
      DIMENSION DT(6,6),SD(6,2),ST(6),SDDT(6),RPLE(6),yumma(6)  
      Dimension PropsUM(40), StatevUM(90)
                                    

c  Here we initialize PropsUM using Props
c
      NPropsUM = 40
      Call Clear(PropsUM, 40)
      PropsUM(1) = Props(1)
      PropsUM(3) = Props(2)
      PropsUM(5) = Props(3)
      PropsUM(6) = Props(4)
      PropsUM(7) = Props(5)
      Do i=1, 12
         PropsUM(10+i) = Props(5+i)
      End Do
      Do i=1, 4
         PropsUM(15+i) = Props(10+i) - 273.0d0
      end do
      Do i=1, 6
         PropsUM(30+i) = Props(17+i)
      End Do
	PropsUM(38) = Props(24)

      model = Int(Props(2)+0.1d0)
      if(model.eq.1) then
         PropsUM(23) = dlog(0.01d0)/(Props(13)-Props(14))
         PropsUM(24) = dlog(0.01d0)/(Props(11)-Props(12))
      else if(model.eq.2) then
         PropsUM(23) = -Props(16)*(Props(14)-Props(13))
         PropsUM(24) = -Props(17)*(Props(11)-Props(12))
      else
         PropsUM(23) = 3.14159265358979/(Props(14)-Props(13))
         PropsUM(24) = 3.14159265358979/(Props(11)-Props(12))
      end if


c     initiating data only for the 1st material point of the 1st iter. 
c     of the first increment     

      IF(KFLAG.NE.920716)THEN
          KFLAG=920716          
          imli=0
          IMTP=0    ! initializing the counter for material points
          MAXMTP1=MAXMTP
          NELMTP=INT(PROPSUM(7)+0.1)  ! given # of material points
          
          IF(NELMTP.GT.MAXMTP)THEN  ! incase max. # of mat. pts. is low
              WRITE(6,4000)MAXMTP,NELMTP
              WRITE(*,4000)MAXMTP,NELMTP
              STOP
          ENDIF       

          MODEL=INT(PROPSUM(3)+0.1D0)

          IF (MODEL.NE.1.AND.MODEL.NE.2.AND.MODEL.NE.3) THEN
              WRITE(6,*) 'THE MODEL NUMBER IS WRONG IN THE INPUT'
              WRITE(*,*) 'THE MODEL NUMBER IS WRONG IN THE INPUT'
              stop 'Wrong model number'
          ENDIF    
          
      CALL CLEAR(TEMPAL, MAXMTP)
          
          IF (NDI.EQ.1) THEN    ! 1-D case
              CALL FORMDD_1(PROPSUM,NPROPSUM)
          ELSEIF (NDI.EQ.2.OR.NDI.EQ.3) THEN  ! 2-D/3-D case  
              CALL FORMDD(PROPSUM,DSM,DSA,DA1,DA2,NPROPSUM) ! calc. mat.ten.
          ENDIF


          TIME1=-999999999.D0 ! initializing time and inc. of time
          DTIME1=TIME1               
      ENDIF
      

      IOPT=INT(PROPSUM(2)+0.1) ! info. about plane stress/strain 
      
c     count the material point number

      IMTP=IMTP+1
      IF(IMTP.GT.NELMTP) IMTP=1 ! beginning of a new inc.
      REALTEMP=TEMP-273.0D0 
      TIME0=TIME(2)

      if(imtp.eq.1) imli=imli+1

      IF (TIME(2).EQ.0.0D0) THEN
         CALL CLEAR(STATEV,90+10*NUMMIN)

         STATEV(1)=PROPSUM(1)
         STATEV(2)=PROPSUM(6)
         STATEV(17)=0.0D0
         STATEV(26)=0.0D0

         TEMPAL(IMTP)=REALTEMP

         IF (NDI.EQ.2.OR.NDI.EQ.3) THEN
            CALL MTASSIGN(PROPSUM(31),STATEV(5),6,1)
            
            DO 5 I=1,6
 5             STATEV(19+I)=PROPSUM(10)*PROPSUM(20)

         ELSE IF (NDI.EQ.1) THEN
            STATEV(5)=PROPSUM(31)
            
            STATEV(20)=PROPSUM(10)*PROPSUM(20)   
            
         ENDIF
      ENDIF                        

c     check if this is the first iteration cycle of an increment

      ISOLVE=3 ! if the same increment
      KEYN=0
      ITER=2           
      
      IF(TIME1.NE.TIME(2))THEN  ! if new increment
          ITER=1
          KEYN=1
          TIME1=TIME(2)
          DTIME1=DTIME
      ELSEIF(DTIME.NE.DTIME1)THEN  ! if time increment changes
          ITER=1
          KEYN=2
          TIME1=TIME(2)
          DTIME1=DTIME
      ENDIF
                                             
      
      IF (NDI.EQ.2.OR.NDI.EQ.3) THEN  ! 2-d/3-d case     
                                             
c     converting the input quantities into second order tensorial form

         CALL STD3_2(ST,SD,STRESS,STRAN,DSTRAN,STATEV,NDI,NSHR,
     *        NSTATV,NTENS)



c     calculate stress increment for given strain and temperature
c     increments, and update the tangent stiffness matrix

         CALL SMA(ST,SD,STATEV,SDDT,RPLE,RPL,DRPLDT,PROPSUM,DT,REALTEMP,
     *           DTEMP,DSA,DSM,DA1,DA2,TEMPAL(IMTP),DTIME,NSTATV,
     *           NPROPSUM,yumma,tumma,NOEL)

c     convert 3-D quantities into proper dimensions

         CALL STF3_2(DT,DDSDDE,STRESS,ST,DDSDDT,DRPLDE,SDDT,RPLE,NDI,
     *        NSHR,NTENS)


         STATEV(48)=DSQRT(2.0D0)/2.0D0*DSQRT((STRESS(1)-STRESS(2))**2
     *              +(STRESS(1)-STRESS(3))**2+(STRESS(2)-STRESS(3))**2)
         STATEV(49)=DSQRT(2.0D0)/3.0D0*DSQRT((STRAN(1)+DSTRAN(1)
     *              -STRAN(2)-DSTRAN(2))**2+(STRAN(1)+DSTRAN(1)
     *              -STRAN(3)-DSTRAN(3))**2+(STRAN(2)+DSTRAN(2)
     *              -STRAN(3)-DSTRAN(3))**2)
         STATEV(50)=(STRESS(1)+STRESS(2)+STRESS(3))/3.0D0

      ELSEIF (NDI.EQ.1) THEN  ! 1-d case     
      
c     assigning the input quantities to local variable
   
          CALL STD_1(STRESS(1),STRAN(1),DSTRAN(1),DFGRD0,DROT(1,1),
     *               STATEV,PROPSUM,NSTATV,NPROPSUM,ST_1,SD_1,DSD_1)          
                                                                     
c     calculate stress increment for given strain and temperature
c     increments, and update the tangent stiffness matrix             

          CALL SMA_1(ST_1,SD_1,DSD_1,STATEV,STATEV(91),SDDT_1,RPLE_1,
     *           RPL,DRPLDT,PROPSUM,DT_1,REALTEMP,DTEMP,
     *           TEMPAL(IMTP),DTIME,NSTATV,NPROPSUM,NUMMIN)
                                                                 
c     assigning local output to global variables

          CALL STF_1(DT_1,DDSDDE(1,1),STRESS(1),ST_1,DDSDDT(1),DRPLDE(1)
     *           ,SDDT_1,RPLE_1,DFGRD1,RPL,DRPLDT,STATEV,PROPSUM,
     *               NSTATV,NPROPSUM)

      ENDIF                                                    
      
 4000 FORMAT(//2X
     *' THE REAL MATERIAL POINTS ARE MORE THAN THE MAXIUM MATERIAL ',
     *' POINTS ALLOWED!'//2X,
     *' PLEASE ENLARGE THE PARAMETER, MAXMTP, IN SUBROUTINE, UMAT TO ',
     *' ENSURE THAT IT IS LARGER THAN THE REAL MATERIAL POINTS.'//2X,
     *' THE MAXIUM MATERIAL POINTS ALLOWED NOW . . . . (MAXMTP)=',I8/2X,
     *' THE REAL MATERIAL POINTS INPUTTED. . . . . . . (NELMTP)=',I8)
  10  FORMAT(1X, I2, 7E12.5)
  11  FORMAT(1X, /)
      
      
      RETURN
      END 
      
c     ******************************************************************
c     *************************** END OF UMAT **************************
c     ******************************************************************

c     ****************************************************************** 
c     *********************** 3-D CASE STARTS HERE *********************
c     ******************************************************************

c     this subroutine is to perform an integration scheme to find the
c     real stresses under the given strain input, and then to update the
c     state variables, for a 2-d/3-d case.       
       
      SUBROUTINE SMA(ST0,SD,STATEV,SDDT,RPLE,RPL,DRPLDT,PROPS,DT,TEMP,
     *               DTEMP,DSA,DSM,DA1,DA2,TEMPINI,DTIME,NSTATV,
     *               NPROPS,yumma,tumma,NOEL)

      IMPLICIT REAL*8 (A-H,O-Z) 
      
      COMMON /DATA02/MAXMTP,NELMTP,IMTP,ITER,KEYN,ISOLVE,KFLAG 
      COMMON /PUNGA/ imli,MODEL
      
      DIMENSION ST0(6),SD(6,2),STATEV(NSTATV),PROPS(NPROPS),DT(6,6),
     * DSA(6,6),DSM(6,6),ST(6),SDT(6),RL(6),DA1(6,6),DA2(3),
     * SDDT(6),RPLE(6),yumma(6),FOLD(6),FNEW(6),DEV_ST(6),
     * STRN(6),RLNEW(6)


c ********************* A R R A Y S ********************

c     ST0(1) : stress vector from MAIN 
c     ST(6)  : local stress vector 
c     DD(6)  : local total strain tensor 
c     RL(6)  : lamda vector 
c     SP(6)  : deviatoric stress vector 
c     SDT(6) : transformation strain vector 
c     SDDT(6):thermal stifness vector 
c     RPLE(6) :

      
      STR_FLG=0.00001D0
      IOPT=INT(PROPS(2)+0.1) ! info about plane stress/strain
      PSI=STATEV(2) ! martensitic volume fraction
      IPHASE=INT(STATEV(1)+0.1) ! info. about no yield/for./rev. trans.
      SEF0 =STATEV(13) ! von-mises stress
      ALPH=PROPS(20)            
      BETA=PROPS(27)                
      GAMA=PROPS(28)

	IRULE = INT(PROPS(38)+0.1)
      
      CALL MTASSIGN(STATEV(5),SDT,1,6) ! assigning state variables to

      NRL=INT(STATEV(17)+0.1D0)

      CALL MTASSIGN(STATEV(20),RL,1,6) ! local variables

      TEMP0= TEMPINI ! initial temperature, To
      DTEMP0=TEMP+DTEMP-TEMP0 ! start temp.-To
      TEMP1=TEMP+DTEMP ! present temperature
      ALPHA=ALFA(PROPS(14),PROPS(15),PSI) ! coeff. of thermal
                                          ! expansion


c calculate the elastic predictor if no time inc. change

      CALL ELASTF(PSI,DSM,DSA,DT)

      DO 15 I=1,3
         ST(I)=ST0(I) ! use local stress vector
         I3=I+3
         ST(I3)=ST0(I3)+DT(I3,I3)*SD(I3,2) ! update stress
         DO 15 J=1,3
 15         ST(I)=ST(I)+DT(I,J)*(SD(J,2)-ALPHA*DTEMP)

      IF(STRS_DOT(ST,ST).LT.STR_FLG) THEN
         DO 17 I=1,6
 17         RLNEW(I)=PROPS(10)*PROPS(20)         
      ELSE
         CALL GETRL(RLNEW,ST,SDT,ALPH,BETA,GAMA,IPHASE, IRULE)
      ENDIF

c check yield of the material point

      CALL CLEAR(FOLD,6)
      CALL CLEAR(FNEW,6)

      DO 20 I=1,3
         I3=I+3
         FOLD(I3)=DA1(I3,I3)*ST0(I3)
         FNEW(I3)=DA1(I3,I3)*ST(I3)
         DO 20 J=1,3
            FOLD(I)=FOLD(I)+DA1(I,J)*ST0(J)
 20         FNEW(I)=FNEW(I)+DA1(I,J)*ST(J)
      
      DOLD=DOT(ST0,RL,6)+0.5D0*DOT(FOLD,ST0,6)+DOT(ST0,DA2,3)
     * *(TEMP-TEMP0)+PROPS(22)*(TEMP-PROPS(16))


      DNEW=DOT(ST,RLNEW,6)+0.5D0*DOT(FNEW,ST,6)+DOT(ST,DA2,3)
     * *(TEMP1-TEMP0)+PROPS(22)*(TEMP1-PROPS(16))

      DECIDE=DNEW-DOLD
      

      IF(DECIDE.GT.0.0D0)THEN
         IPHASE=1
         IFYD=1
      ELSEIF(DECIDE.LT.0.0D0)THEN
         IPHASE=2
         IFYD=1
      ELSE
         IFYD=0
         GOTO 70
      ENDIF

c     if yielding happened, return mapping integration scheme will 
c     be performed to find the real stress increment under the given 
c     strain and temperature increments and current state input.

      CALL BEINT1(ST0,SDT,PSI,SD,TEMP,DTEMP,TEMP0,PROPS,DSA,DSM,
     *            IPHASE,DA1,DA2,DT,SDDT,RPLE,RPL,DRPLDT,DTIME,NPROPS,
     *            yumma,tumma,RL,NRL,IFYD,NOEL)

c     update stresses, and all state vaiables

 70   STATEV(1)=DFLOAT(IPHASE)  ! whether forward/reverse trans.
      STATEV(2)=PSI
      STATEV(3)=DFLOAT(IFYD)  
      
      CALL MTASSIGN(SDT,STATEV(5),6,1) ! assign transform. strain 
      
      STATEV(11)=SD(2,1)+SD(2,2) ! 33-component of strain
      STATEV(12)=ST0(2)  ! 33-component of stress
      
      CALL DEV_STRS(ST0,DEV_ST)    
      
      SEF=STRS_DOT(DEV_ST,DEV_ST) ! von-mises stress
      STATEV(13)=DSQRT(1.5D0*SEF)
      STATEV(14)=DSQRT(2.D0/3.D0*STRN_DOT(SDT,SDT))! inelastic fiber strain
      STATEV(15)=TEMP1  ! assigning temp. 
      
c     get lamda(RL), effective stress(U)      
                                          
      STATEV(16)=DOT(RL,ST0,6)/PROPS(20) ! modified effective stress
      STATEV(17)=DFLOAT(NRL)
      STATEV(26)=DTEMP0

      CALL MTASSIGN(RL,STATEV(20),6,1)

      DO 80 I=1,6
 80      STRN(I)=SD(I,1)+SD(I,2)
     

      RETURN
      END   
      
c     ******************************************************************

c     this subroutine is to integrate the stress, inelastic strain, and
c     internal state variable, increment for gaven strain and temperature
c     increments by using Return Mapping Method (elastic predictor-plastic
c     corrector scheme), for a 2-d/3-d case.

      SUBROUTINE BEINT1(ST0,SDT0,PSI0,SD,TEMP,DTEMP,TEMP0,PROPS,DSA,DSM,
     *                  IPHASE,DA1,DA2,DT,SDDT,RPLE,RPL,DRPLDT,DTIME,
     *                  NPROPS,yumma,tumma,RL,NRL,IFYD,NOEL)

      IMPLICIT REAL*8 (A-H,O-Z)
      
      COMMON /DATA02/MAXMTP1,NELMTP,IMTP,ITER,KEYN,ISOLVE,KFLAG
      COMMON /PUNGA/ imli,MODEL
      
      DIMENSION SD(6,2),ST0(6),SDT0(6),PROPS(NPROPS),DSA(6,6),DSM(6,6),
     *          DT(6,6),ST(6),SDT(6),DD1(6),DD(6),Q(6),RL(6),R(6),
     *          DA1(6,6),DA2(3),SDDT(6),RPLE(6),yumma(6),F(6) 
     

c     ********************* A  R   R   A   Y   S *********************** 

c     SDT0(6)     : transformation strain, actually SDT(6) from sma
c     DD1(6)      : local total strain vector 
c     DD(6)       : local strain inc. vector
c     Q(6)        : see papers
c     R(6)        : see papers
c     D(6,6)      : elastic stiffness matrix
c     ST(6)       : local stress vector
c     SDT(6)      : local transformation strain vector
c     RL(6)       : lamda vector, see papers
c     ES1(6)      : elastic strain at previous increment
c     Q1(6)       : see papers
              
      TOL=PROPS(5)              ! criteria of convergence
      DTMP=DTEMP                ! temp. inc. for each inc. within iter.
      TMP=TEMP                  ! starting temp.
      PSI=PSI0                  ! martensitic volume frac.
      ALPH=PROPS(20)            
      BETA=PROPS(27)                
      GAMA=PROPS(28)
	
	IRULE = INT(PROPS(38)+0.1)

      DO 10 I=1,6           ! assigining stress, strain and strain inc. 
          ST(I)=ST0(I)      ! values from sma to local vectors
          SDT(I)=SDT0(I)
          DD1(I)=SD(I,1)
 10       DD(I)=SD(I,2)

      FLAG1F=1.0D0
      FLAG1R=0.0D0
      FLAG2R=0.001D0
      FLAG2F=0.99D0
      
      if(psi.le.0.0001D0) then
        iphase=1
      end if


c     to implement the return mapping integration scheme, calculate 
c     elastic predictor

      TMP=TMP+DTMP          ! present temp.
      DTMP0=TMP-TEMP0           ! present temp-To
      ALPHA=ALFA(PROPS(14),PROPS(15),PSI)
      
      CALL ELASTF(PSI,DSM,DSA,DT) !get elas. stiff. mat.     
      CALL ELATTF(DT,ALPHA,SDDT) 
      

      DO 20 I=1,6               ! updating strain and stress
         DD1(I)=DD1(I)+DD(I)
         DO 20 J=1,6
 20         ST(I)=ST(I)+DT(I,J)*DD(J)
 
      DO 30 I=1,3     ! adding the thermal part
 30      ST(I)=ST(I)+SDDT(I)*DTMP   
         

c     iteration for plastic corrector

      DO 50 LOOP=1,100

c     calculate rl and the value of the yield function
         
      CALL GETRL(RL,ST,SDT,ALPH,BETA,GAMA,IPHASE, IRULE)


         CALL YDFUN(FYD,ST,SDT,RL,TMP,DTMP0,PSI,PROPS,IPHASE,
     *              DA1,DA2,NPROPS,F) 

         IF(LOOP.EQ.1)THEN
            IF(FYD.LE.0.0D0)THEN
               IFYD=0
               GOTO 9040  
            ELSEIF(IPHASE.EQ.1.AND.PSI.GE.(FLAG1F-0.0001).AND.
     *              DABS(FYD).GT.0.0D0)THEN
               IFYD=0
               GOTO 9040
            ELSEIF(IPHASE.EQ.2.AND.PSI.LE.(FLAG1R+0.0001).AND.
     *              DABS(FYD).GT.0.0D0)THEN
               IFYD=0
               GOTO 9040
            ENDIF
         ENDIF

         
c     get rl, and b, also see papers/scheme              

         CALL GETQ(Q,DA1,DA2,DTMP0,RL,ST)

         CALL GETB(B,PSI,Q,PROPS,DT,NPROPS,IPHASE)

c     update PSI, SDT, and ST     
     
         DPSI=-FYD/B            ! inc. of mart. vol. frac.
         psii=psi
         PSI=PSI+DPSI           ! updating mart. vol. frac.
                                     
         if(psi.gt.FLAG1F)then
            if(iphase.eq.1)then
               psi=FLAG1F
               dpsi=psi-psii
            elseif(iphase.eq.2.and.psii.lt.FLAG1F)then
               psi=FLAG2R
               dpsi=psi-psii
            endif
         elseif(psi.lt.FLAG1R)then
            if(iphase.eq.1.and.psii.gt.FLAG1R)then
               psi=FLAG2F
               dpsi=psi-psii
            elseif(iphase.eq.2)then
               psi=FLAG1R
               dpsi=psi-psii
            endif
         endif


         DO 40 I=1,6            ! see scheme
 40         SDT(I)=SDT(I)+RL(I)*DPSI 

         CALL ELASTF(PSI,DSM,DSA,DT)
         ALPHA=ALFA(PROPS(14),PROPS(15),PSI)
         
         CALL CLEAR(ST,6)

c         pause"is the problem here?"

         DO 45 I=1,3
            I3=I+3
            ST(I3)=DT(I3,I3)*(DD1(I3)-(SDT(I3)-PROPS(30+I3)))
            DO 45 J=1,3
 45            ST(I)=ST(I)+DT(I,J)*(DD1(J)-ALPHA*DTMP0
     *                -(SDT(J)-PROPS(30+J)))
     

c     check convergence

         IF(DABS(DPSI).LT.TOL) GOTO 9020

 50   CONTINUE

      WRITE(*,*)'ITERATION FAILS TO CONVERGE! Del_xi=',DPSI,'  xi=', PSI

 9020 CONTINUE
  100 CONTINUE  

c     update tangent stifness matrix      

      if(iphase.eq.1.and.psi.lt.FLAG1F
     *     .or.iphase.eq.2.and.psi.gt.FLAG1R)then
      CALL TANSTF(DT,SDDT,PSI,ST,SDT,DTMP0,DSA,DSM,PROPS,IPHASE,DA1,
     *            DA2,RL,NPROPS)
      endif
      
      if(psi.le.0.0001D0) then
        iphase=1
      end if

    
      GOTO 9040   ! return

 9030 CONTINUE  
 

c     if PSI is out of bound, forward Euler integration scheme is applied. 

      IF(PSI.LT.1.0E-3.AND.IPHASE.EQ.2.OR.PSI.GT.99.99E-2.AND.
     *     IPHASE.EQ.1)THEN
         NSS=100
         DTMP=DTEMP/DFLOAT(NSS)
         DTMP0=TEMP-TEMP0
         PSI=PSI0
         DO 110 I=1,6           ! assigning local variables with global values
            ST(I)=ST0(I)
            SDT(I)=SDT0(I)
            DD1(I)=SD(I,1)
 110        DD(I)=SD(I,2)/DFLOAT(NSS) 
            
         CALL ELASTF(PSI,DSM,DSA,DT)    
 
         DO 200 N=1,NSS 
            DTMP1=DTMP0+DTMP 
            
            CALL GETRL(RL,ST,SDT,H,B2,IPHASE,IRULE)
            CALL GETQ(Q,DA1,DA2,DTMP1,RL,ST)
            CALL GETB(B,PSI,Q,PROPS,DT,NPROPS,IPHASE)
            CALL GETRS(R,S,Q,DT,DA2,ST,PROPS,NPROPS,IPHASE)
            
            DPSI=-(DOT(R,DD,6)+S*DTMP)/B
            PSI=PSI+DPSI      
            
            IF(PSI.GT.1.0D-5.AND.PSI.LT.0.999999D0)THEN ! if not near the 
               DTMP0=DTMP0+DTMP ! finish
               DO 120 I=1,6
                  DD1(I)=DD1(I)+DD(I) ! strain
 120              SDT(I)=SDT(I)+RL(I)*DPSI ! transformation strain
                  
               ALPHA=ALFA(PROPS(14),PROPS(15),PSI)
               
               CALL ELASTF(PSI,DSM,DSA,DT) ! elas. stiff. mat. 
               CALL ELATTF(DT,ALPHA,SDDT)
               CALL CLEAR(ST,6)
               
               DO 130 I=1,3     ! stress
                  I3=I+3
                  ST(I3)=DT(I3,I3)*(DD1(I3)-(SDT(I3)-PROPS(30+I3)))
                  DO 130 J=1,3
 130                 ST(I)=ST(I)+DT(I,J)*(DD1(J)-ALPHA*DTMP0
     *                     -(SDT(J)-PROPS(30+J)))
            ELSE       ! near the finish and is done elastic
               DO 140 I=1,6
 140              DD(I)=SD(I,1)+SD(I,2)-DD1(I) ! strain inc.
                     
               PSI=PSI-DPSI ! why
               
               ALPHA=ALFA(PROPS(14),PROPS(15),PSI)
               
               CALL ELASTF(PSI,DSM,DSA,DT) ! elastic stiffness
               CALL ELATTF(DT,ALPHA,SDDT)
               
               DO 150 I=1,3     ! stress
                  I3=I+3
                  ST(I3)=ST(I3)+DT(I3,I3)*DD(I3)
                  ST(I )=ST(I )-(DT(I,1)+DT(I,2)+DT(I,3))*ALPHA
     *                   *(TEMP+DTEMP-DTMP0-TEMP0)
                  DO 150 J=1,3
 150                 ST(I)=ST(I)+DT(I,J)*DD(J)
  
               GOTO 9040
            ENDIF      
 200     CONTINUE
      ELSE
         WRITE(*  ,*)'IPHASE-PSI WRONG!  IPHASE,PSI=',IPHASE,PSI
         WRITE(MSG,*)'IPHASE-PSI WRONG!  IPHASE,PSI=',IPHASE,PSI
      ENDIF

 9040 CONTINUE 

      IF(INT(PROPS(8)+0.1D0).EQ.1.AND.IFYD.EQ.1)THEN
         if(iphase.eq.1.and.psi.lt.(FLAG1F-0.0001d0)
     *        .or.iphase.eq.2.and.psi.gt.(FLAG1R+0.0001d0))then         
            T=TEMP+DTEMP        ! current temperature
            DTEMP0=TEMP+DTEMP-TEMP0 ! current temperature-reference temp.
            DPSI=PSI-PSI0
 
c     get lamda(RL), and then Q      
                                          
            CALL GETRL(RL,ST,SDT,ALPH,BETA,GAMA,IPHASE,IRULE) ! model 1,3
            CALL GETQ(Q,DA1,DA2,DTEMP0,RL,ST)
      
c     getting the variables necessary for coupled heat transfer analysis

            CALL CALCRPL(RPL,RPLE,DRPLDT,PROPS,IPHASE,T,SD,ST,st0,DT,
     *                   DA2,SDDT,Q,DTEMP,DTIME,DPSI,PSI,NPROPS,yumma,
     *                   tumma)
         else
            CALL CLEAR(RPLE,6)
            
            RPL=0.0D0
            DRPLDT=0.0D0
         endif
      ELSE
         CALL CLEAR(RPLE,6)
         
         RPL=0.0D0
         DRPLDT=0.0D0          
      ENDIF

c     assign local state variables to global ones
 
      PSI0=PSI
      DO 90 I=1,6
         SDT0(I)=SDT(I)         ! transformation strain
 90      ST0(I)=ST(I)           ! stress

      
      RETURN
      END   
      
c     ******************************************************************

c     this subroutine is to get the tangent STIFFNESS matrix DT, and
c     thermal STIFFNESS coefficients, ALS.

      SUBROUTINE TANSTF(DT,SDDT,PSI,ST,SDT,DTEMP0,DSA,DSM,PROPS,
     *                  IPHASE,DA1,DA2,RL,NPROPS)

      IMPLICIT REAL*8 (A-H,O-Z)  
      
      COMMON /PUNGA/ imli,MODEL

      DIMENSION DT(6,6),SDDT(6),ST(6),SDT(6),DSM(6,6),DSA(6,6),TEMP(6),
     *          PROPS(NPROPS),R(6),Q(6),RL(6),DA1(6,6),DA2(3)   
     
     
c     *********************  A    R   R   A   Y   S ********************

c     R(6)        : see scheme
c     Q(6)        : see scheme
c     RL(6)       : lamda vector          

      CALL CLEAR(SDDT,6)
      CALL CLEAR(TEMP,6)
      
      ALPH=PROPS(20)
      BETA=PROPS(27)    
      GAMA=PROPS(28)

      ALPHA=ALFA(PROPS(14),PROPS(15),PSI) ! coeff. of ther. exp.

c     get q, r, s, and b, also see scheme                        
                       
      CALL GETQ(Q,DA1,DA2,DTEMP0,RL,ST)
      CALL ELASTF(PSI,DSM,DSA,DT)
      CALL GETRS(R,DS,Q,DT,DA2,ST,PROPS,NPROPS,IPHASE)
      CALL GETB(B,PSI,Q,PROPS,DT,NPROPS,IPHASE)


      DO 10 I=1,3               ! thermal stiffness vector
         I3=I+3
         SDDT(I+3)=DT(I3,I3)*Q(I3)*DS/B
         DO 10 J=1,3
 10         SDDT(I)=SDDT(I)+DT(I,J)*(Q(J)*DS/B-ALPHA)

      DO 15 I=1,3
         I3=I+3
         TEMP(I3)=DT(I3,I3)*Q(I3)
         DO 15 J=1,3
 15         TEMP(I)=TEMP(I)+DT(I,J)*Q(J)                

      DO 20 I=1,6
         DO 20 J=1,6
 20         DT(I,J)=DT(I,J)+TEMP(I)*R(J)/B ! tangent stiffness matrix
 
      RETURN
      END   
      
c     ******************************************************************

c     this subroutine calculates q for a 2-d/3-d case
c     also see scheme.           

      SUBROUTINE GETQ(Q,DA1,DA2,DTEMP0,RL,ST) 
     
      IMPLICIT REAL*8 (A-H,O-Z)
      
      DIMENSION Q(6),ST(6),DA1(6,6),DA2(3),RL(6)
     
  
      CALL CLEAR(Q,6)   

      DO 10 I=1,3
      I3=I+3
      Q(I3)=DA1(I3,I3)*ST(I3)+RL(I3)
          Q(I)=DA2(I)*DTEMP0+RL(I)
      DO 10 J=1,3
 10       Q(I)=Q(I)+DA1(I,J)*ST(J)


         
      RETURN   
      END   

c     ************************************************************
      
c     this subroutine calculates b for a 2-d/3-d case
c     also see scheme.           

      SUBROUTINE GETB(B,PSI0,Q,PROPS,D,NPROPS,IPHASE)
    
      IMPLICIT REAL*8 (A-H,O-Z)
    
      COMMON /PUNGA/ imli,MODEL

      DIMENSION Q(6),PROPS(NPROPS),D(6,6),TEMP(6)
    
      PSI=PSI0  

      IF(PSI.GT.0.99993) PSI=0.99993
      IF(PSI.LT.0.00007) PSI=0.00007

      CALL CLEAR(TEMP,6)    

      DO 10 I=1,3
      I3=I+3
      TEMP(I3)=Q(I3)*D(I3,I3)
      DO 10 J=1,3   
 10       TEMP(I)=TEMP(I)+Q(J)*D(J,I)       
     
  
c     see scheme again                       
                       
      IF(IPHASE.EQ.1)THEN  ! forward transforamtion 
          IF (MODEL.EQ.1) THEN  ! model 1
              B=-DOT(Q,TEMP,6)-PROPS(22)/PROPS(24)/(1.0D0-PSI)
          ELSEIF (MODEL.EQ.2) THEN  ! model 2
              B=-DOT(Q,TEMP,6)-(2.0D0*PROPS(26)*PSI+PROPS(24))
          ELSEIF (MODEL.EQ.3) THEN  ! model 3
              B=-DOT(Q,TEMP,6)+2.0*PROPS(22)/PROPS(24)
     *          /DSQRT(1.0D0-(2.0*PSI-1.0D0)**2)               
          ENDIF        
      ELSEIF(IPHASE.EQ.2)THEN  ! reverse transformation 
          IF (MODEL.EQ.1) THEN  ! model 1
              B=DOT(Q,TEMP,6)-PROPS(21)/PROPS(23)/PSI
          ELSEIF (MODEL.EQ.2) THEN  ! model 2
              B=DOT(Q,TEMP,6)+(2.0D0*PROPS(25)*PSI+PROPS(23))    
          ELSEIF (MODEL.EQ.3) THEN  ! model 3
              B=DOT(Q,TEMP,6)-2.0*PROPS(21)/PROPS(23)
     *          /DSQRT(1.0D0-(2.0*PSI-1.0D0)**2)    
          ENDIF  
      ENDIF


      RETURN
      END


c***********************************************************************

c     this subroutine calculates r and s for a 2-d/3-d case
c     also see scheme. 


      SUBROUTINE GETRS(R,DS,Q,D,DA2,ST,PROPS,NPROPS,IPHASE)
    
      IMPLICIT REAL*8 (A-H,O-Z)
    
      DIMENSION R(6),Q(6),D(6,6),DA2(3),ST(6),PROPS(NPROPS)
    
     
      CALL CLEAR(R,6)   
 
      DO 20 I=1,3
      I3=I+3
      R(I3)=Q(I3)*D(I3,I3)
      IF(IPHASE.EQ.2) R(I3)=-R(I3)
      DO 10 J=1,3
 10       R(I)=R(I)+Q(J)*D(J,I)
 20   IF(IPHASE.EQ.2) R(I)=-R(I)        

        ALPHA=ALFA(PROPS(14),PROPS(15),PSI)
    
        
c     see scheme again                       
                       
      IF(IPHASE.EQ.1)THEN  ! forward transforamtion 
      DS=DOT(DA2,ST,3)+PROPS(22)-(R(1)+R(2)+R(3))*ALPHA
      ELSEIF(IPHASE.EQ.2)THEN  ! reverse transformation 
      DS=-DOT(DA2,ST,3)-PROPS(21)-(R(1)+R(2)+R(3))*ALPHA
      ENDIF

      RETURN
      END


c***********************************************************************

c     this subroutine calculates lamda(RL) and effevtice stress(U) for a 
c     2-d/3-d case
     
      SUBROUTINE GETRL(RL,ST,SDT,ALPH,BETA,GAMA,IPHASE, IRULE)
      
      IMPLICIT REAL*8 (A-H,O-Z)              
      
      COMMON /PUNGA/ imli,MODEL

      DIMENSION RL(6),ST(6),SDT(6),DST(6),DSTINV(3,3),DSTTEMP(6) 
      
      
c     *********************** A   R   R   A   Y   s ********************

c     SP(6)       : deviatoric stress vector

      IF((IPHASE.EQ.2).AND.(IRULE.NE.2))THEN
         DEN=DSQRT(2.0D0/3.0D0*STRN_DOT(SDT,SDT))

         IF(DEN.EQ.0.0D0)THEN
            DO 10 I=1,6
 10            RL(I)=ALPH
         ELSE
            DEN=ALPH/DEN
         
            DO 20 I=1,6
 20            RL(I)=DEN*SDT(I)
         ENDIF            
      ELSE
         CALL DEV_STRS(ST,DST)

         DJ2=0.5D0*STRS_DOT(DST,DST)

         IF(DJ2.EQ.0.D0)THEN    ! anomaly
            DO 30 I=1,6
 30            RL(I)=1.0D0
         ELSE
            DEN=ALPH*(3.0d0/2.0D0)/DSQRT(3.0D0*DJ2)
         
            DO 40 I=1,3
               RL(I)=DEN*DST(I)+GAMA
 40            RL(I+3)=2.0D0*DEN*DST(I+3)
         ENDIF
      ENDIF

      RETURN
      END

c     ******************************************************************  

c     calculates the value of the yield function for a material point of 
c     a 2-d/3-d case.

      SUBROUTINE YDFUN(FYD,ST,SDT,RL,TEMP1,DTEMP1,PSI0,PROPS,IPHASE,
     *                 DA1,DA2,NPROPS,F)   
      
      IMPLICIT REAL*8 (A-H,O-Z)          
      
      COMMON /PUNGA/ imli,MODEL

      DIMENSION ST(6),SDT(6),PROPS(NPROPS),RL(6),DA1(6,6),
     *      DA2(3),F(6) 
      
c     ********************* A   R   R   A   Y   S **********************

c     Rl(6)       : local lamda vector      

      FYD=0.D0  ! yield value of the function
      PSI=PSI0  ! martensitic volume fraction

      IF(PSI.GT.0.99993) PSI=0.99993
      IF(PSI.LT.0.00007) PSI=0.00007

      ALPH=PROPS(20)  ! max. transformation strain
      BETA=PROPS(27)
      GAMA=PROPS(28)
      
c     calculate lamda vector and effective stress      
                                                      
      CALL CLEAR(F,6)

      DO 10 I=1,3
      I3=I+3
          F(I3)=DA1(I3,I3)*ST(I3)
          DO 10 J=1,3
  10          F(I)=F(I)+DA1(I,J)*ST(J) 
  
      IF(IPHASE.EQ.1)THEN   ! forward transforamtion
          IF (MODEL.EQ.1) THEN  ! model 1
              FYD=DOT(ST,RL,6)+(1.0D0/2.0D0)*DOT(F,ST,6)+DOT(DA2,ST,3)
     *            *DTEMP1+PROPS(22)*TEMP1+PROPS(22)/PROPS(24)
     *        *DLOG(1.D0-PSI)-PROPS(22)*PROPS(16)
          ELSEIF (MODEL.EQ.2) THEN   ! model 2 
              FYD=DOT(ST,RL,6)+(1.0D0/2.0D0)*DOT(F,ST,6)+DOT(DA2,ST,3)
     *            *DTEMP1+PROPS(22)*(TEMP1-PROPS(16))-(PROPS(26)*PSI
     *            +PROPS(24))*PSI
          ELSEIF (MODEL.EQ.3) THEN   ! model 3 
              FYD=DOT(ST,RL,6)+(1.0D0/2.0D0)*DOT(F,ST,6)+DOT(DA2,ST,3)
     *            *DTEMP1+PROPS(22)*TEMP1-PROPS(22)/PROPS(24)
     *            *(DACOS(2.0D0*PSI-1.0D0)-22.0D0/7.0D0)-PROPS(22)
     *        *PROPS(16)     
          ENDIF               
      ELSEIF(IPHASE.EQ.2)THEN  ! reverse transformation
          IF (MODEL.EQ.1) THEN  ! model 1
              FYD=-DOT(ST,RL,6)-(1.0D0/2.0D0)*DOT(F,ST,6)-DOT(DA2,ST,3)
     *            *DTEMP1-PROPS(21)*TEMP1-PROPS(21)/PROPS(23)
     *        *DLOG(PSI)+PROPS(21)*PROPS(18)
        ELSEIF (MODEL.EQ.2) THEN   ! model 2 
              FYD=-DOT(ST,RL,6)-(1.0D0/2.0D0)*DOT(F,ST,6)-DOT(DA2,ST,3)
     *            *DTEMP1-PROPS(21)*(TEMP1-PROPS(19))+(PROPS(25)*PSI
     *            +PROPS(23))*PSI
          ELSEIF (MODEL.EQ.3) THEN   ! model 3
              FYD=-DOT(ST,RL,6)-(1.0D0/2.0D0)*DOT(F,ST,6)-DOT(DA2,ST,3)
     *            *DTEMP1-PROPS(21)*TEMP1+PROPS(21)/PROPS(23)
     *        *(DACOS(2.0D0*PSI-1.0D0)-22.0D0/7.0D0)+PROPS(21)
     *        *PROPS(19) 
          ENDIF
      ENDIF
                                                               
      RETURN
      END


c     ******************************************************************

c     this subroutine finds the elastic stiffness matrix

      SUBROUTINE ELASTF(PSI,DSM,DSA,D) 
      
      IMPLICIT REAL*8 (A-H,O-Z)  
      
      DIMENSION DSM(6,6),DSA(6,6),D(6,6),TEMP(3,3)


      DO 10 I=1,6
          DO 10 J=1,6
 10           D(I,J)=DSA(I,J)+PSI*(DSM(I,J)-DSA(I,J))


      DO 20 I=1,3
      DO 20 J=1,3
 20       TEMP(I,J)=D(I,J)
    
      CALL GET_INV_DET(TEMP,GARB,0)

      DO 30 I=1,3
      DO 30 J=1,3
 30       D(I,J)=TEMP(I,J)
    
      D(4,4)=1.0D0/D(4,4)           
      D(5,5)=1.0D0/D(5,5)           
      D(6,6)=1.0D0/D(6,6)           

            
      RETURN
      END

c     ******************************************************************

c     this subroutine finds the elastic thermal stiffness matrix

      SUBROUTINE ELATTF(D,ALPHA,TST)
    
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION D(6,6),TST(6)

      CALL CLEAR(TST,6)

      DO 10 I=1,3
 10   TST(I)=-ALPHA*(D(I,1)+D(I,2)+D(I,3))
 
      RETURN
      END
    
     
c***********************************************************************        

c     this subroutine is used to calculate the variables associated with
c     the coupled mechanical-thermal analysis

      SUBROUTINE CALCRPL(RPL,RPLE,DRPLDT,PROPS,IPHASE,T,SD,ST,st0,DT,
     *                   DA2,SDDT,Q,DTEMP,DTIME,DPSI,PSI,NPROPS,yumma,
     *                   tumma)

      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON /PUNGA/ imli,MODEL

      DIMENSION RPLE(6),PROPS(NPROPS),ST(6),SD(6,2),DT(6,6),SDDT(6),
     *          G(6),DA2(3),Q(6),DRPE(6),yumma(6),st0(6)


      CALL CLEAR(DRPE,6)    
      RPL=0.0D0
      TEMP=T

c     calling the subroutine to calculate the a and g vectors beside h

      CALL RPLEXT(G,HT,DA2,TEMP,ST,IPHASE,PROPS,Q,DPSI,RPLE,
     *            DRPLDT,DT,PSI,SDDT,NPROPS,st0,rpl)

      DO 10 I=1,6
         RPLE(I)=RPLE(I)/DTIME  
         DO 10 J=1,6
 10         DRPE(I)=DRPE(I)+G(J)*DT(J,I)

      DRPT=DOT(G,SDDT,6)+HT 


      DRPLDT=DRPLDT/DTIME   
      RPL=RPL/DTIME 

      do 40 i=1,6
 40      yumma(i)=drpe(i)

      tumma=drpt   


      RETURN
      END

      
c     ******************************************************************


      SUBROUTINE RPLEXT(G,HT,DA2,T,ST,IPHASE,PROPS,Q,DPSI,RPLE,
     *              DRPLDT,DT,PSI,SDDT,NPROPS,st0,rpl)

      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON /PUNGA/ imli,MODEL

      DIMENSION A(6),G(6),DA2(3),ST(6),PROPS(NPROPS),Q(6),RPLE(6),
     *        DT(6,6),SDDT(6),st0(6),dst(6)         

c   only 2nd model now  

      CALL CLEAR(RPLE,6)    
      CALL CLEAR(A,6)
      CALL CLEAR(G,6)

      T=T+273.0D0
      DALPH=DA2(1)
      ALPHA=ALFA(PROPS(14),PROPS(15),PSI)

      IF (IPHASE.EQ.1) THEN   ! forward transformation
         CALL MTASSIGN(Q,A,6,1) 
         B=DOT(DA2,ST,3)+PROPS(22)
         C=-(2.0D0*PROPS(26)*PSI+PROPS(24))
         DMU2=1.0D0/6.0D0*(PROPS(25)-PROPS(26))+1.0D0/4.0D0
     *        *(PROPS(23)-PROPS(24))
         FYD=-1.0D0/2.0D0*PROPS(22)*(PROPS(19)-PROPS(16))
     *        -DMU2
         HELLO=FYD-T*DALPH*(ST(1)+ST(2)+ST(3))-PROPS(22)*T
         YELLO=PROPS(22)*T-FYD
         DELLO=PROPS(22)
      ELSEIF(IPHASE.EQ.2) THEN  ! reverse transformation
         DO 15 I=1,6
 15         A(I)=-Q(I)      
         B=-DOT(DA2,ST,3)-PROPS(21)
         C=2.0D0*PROPS(25)*PSI+PROPS(23)
         DMU2=1.0D0/6.0D0*(PROPS(25)-PROPS(26))+1.0D0/4.0D0
     *        *(PROPS(23)-PROPS(24))
         FYD=1.0D0/2.0D0*PROPS(21)*(PROPS(19)-PROPS(16))
     *        +DMU2
         HELLO=FYD-T*DALPH*(ST(1)+ST(2)+ST(3))-PROPS(21)*T
         YELLO=PROPS(21)*T-FYD
         DELLO=PROPS(21)
      ENDIF     
         
      DO 20 I=1,6
         DO 20 J=1,6
 20         RPLE(I)=RPLE(I)+A(J)*DT(J,I)
               
      DRPLDT=YELLO/C*(DOT(A,SDDT,6)+B)-DELLO*DPSI               

      DO 25 I=1,3
         RPLE(I)=RPLE(I)*YELLO/C    
         G(I)=-(T*ALPHA+A(I)*HELLO/C)
         RPLE(I+3)=RPLE(I+3)*YELLO/C    
 25      G(I+3)=-A(I+3)*HELLO/C

      HT=-B*HELLO/C

      do 26 i=1,6
 26      dst(i)=st(i)-st0(i)

      rpl=hello*dpsi-t*alpha*(dst(1)+dst(2)+dst(3))

      RETURN    
      END 


c     ******************************************************************

c     this subroutine is to calculate the deviatoric stress tensor, SP
c     and the effective stress, SEF.

      SUBROUTINE DEV_STRS(ST,SP)

      IMPLICIT REAL*8 (A-H,O-Z)
      
      DIMENSION ST(6),SP(6)

      S1=(ST(1)+ST(2)+ST(3))/3.D0  ! calculating diagonal elements
      SP(1)=ST(1)-S1
      SP(2)=ST(2)-S1
      SP(3)=ST(3)-S1
      
      DO 10 I=1,3
 10       SP(I+3)=ST(I+3) ! calculating off-diagonal elements

      
      RETURN
      END

c     ******************************************************************

c     this function calculates the product of strain components

      FUNCTION STRN_DOT(Q,S) 
      
      IMPLICIT REAL*8 (A-H,O-Z)
      
      DIMENSION Q(6),S(6)

      D=0.D0
      
      DO 10 I=1,3
 10       D=D+Q(I)*S(I)+0.5D0*Q(I+3)*S(I+3) 
 
      STRN_DOT=D 
      
      
      RETURN
      END   
      
c     ******************************************************************

c     this function calculates the product of stress components      

      FUNCTION STRS_DOT(Q,S)    
      
      IMPLICIT REAL*8 (A-H,O-Z)
      
      DIMENSION Q(6),S(6)

      D=0.D0 
      
      DO 10 I=1,3
 10       D=D+Q(I)*S(I)+2.0D0*Q(I+3)*S(I+3)   
 
      STRS_DOT=D    
      
      
      RETURN
      END   
      

c     ******************************************************************      

c   this subroutine inverts a given matrix and return it  
        

      SUBROUTINE GET_INV_DET(XS,DET,NFLAG)
                                      
      IMPLICIT REAL*8 (A-H,O-Z)                                         
                 
      DIMENSION XS(3,3),A(3,3) 
    
      DO 10 I=1,3                                                       
      I1=I+1                                                            
      I2=I+2                                                            
      IF(I1.GT.3) I1=I1-3                                               
      IF(I2.GT.3) I2=I2-3                                               
      DO 10 J=1,3                                                       
          J1=J+1                                                            
          J2=J+2                                                            
          IF(J1.GT.3) J1=J1-3                                               
          IF(J2.GT.3) J2=J2-3                                               
   10         A(I,J)=XS(I1,J1)*XS(I2,J2)-XS(I1,J2)*XS(I2,J1)
                      
      DET=0.D0 
                                                             
      DO 20 I=1,3                                                       
   20     DET=DET+A(1,I)*XS(1,I) 
   
      NE=1  

      if (NFLAG.EQ.0.AND.det.le.0.0d0) then
      write(*,*) 'error in the jacobain, DETERMINANT', det
      pause'error in the jacobian' 
      endif         
                                           
      DO 30  I=1,3                                                      
      DO 30 J=1,3                                                       
   30         XS(I,J)=A(J,I)/DET 
                                                  
 
           
      RETURN                                                            
      END           
    

c***********************************************************************
     
c     this subroutine is to initiate the material data for a 2-3/3-d case.      

      SUBROUTINE FORMDD(PROPS,DSM,DSA,DA1,DA2,NPROPS)

      IMPLICIT REAL*8 (A-H,O-Z)  
      
      COMMON /PUNGA/ imli,MODEL

      DIMENSION PROPS(NPROPS),DSM(6,6),DSA(6,6),DA1(6,6),DA2(3)

      CALL CLEAR(DSA,36) ! initializing 
      CALL CLEAR(DSM,36)
      CALL CLEAR(DA1,36) ! initializing
      CALL CLEAR(DA2,3)  ! initializing
      
      V=PROPS(13)  ! poisson's ratio
      EA=PROPS(11) ! young's modulii
      EM=PROPS(12) 
      
c     working on martenstic material matrix, see any mechanics book for
c     isotropic material matrix      


      DO 30 I=1,3  ! diagonal terms
      I3=I+3
      DSA(I,I)=1.0D0/EA
      DSM(I,I)=1.0D0/EM
      DA1(I,I)=DSM(I,I)-DSA(I,I)

      DSA(I3,I3)=2.0D0*(1+V)/EA
      DSM(I3,I3)=2.0D0*(1+V)/EM
 30   DA1(I3,I3)=DSM(I3,I3)-DSA(I3,I3)
                             
C         off-diagonal terms                             
                             
      DSA(1,2)=-V/EA
      DSM(1,2)=-V/EM
      DA1(1,2)=DSM(1,2)-DSA(1,2)

      DSA(1,3)=DSA(1,2)
      DSA(2,1)=DSA(1,2)
      DSA(2,3)=DSA(1,2)
      DSA(3,1)=DSA(1,2)
      DSA(3,2)=DSA(1,2)

      DSM(1,3)=DSM(1,2)
      DSM(2,1)=DSM(1,2)
      DSM(2,3)=DSM(1,2)
      DSM(3,1)=DSM(1,2)
      DSM(3,2)=DSM(1,2)

      DA1(1,3)=DA1(1,2)
      DA1(2,1)=DA1(1,2)
      DA1(2,3)=DA1(1,2)
      DA1(3,1)=DA1(1,2)
      DA1(3,2)=DA1(1,2)
                              
      DO 40 I=1,3
 40   DA2(I)=PROPS(15)-PROPS(14)       
     
     
      RETURN
      END   
      
c     ******************************************************************

c     this subroutine is to transform 3/2-D stress and strain tensors 
c     into 1-D ones.      

      SUBROUTINE STD3_2(ST,SD,STRESS,STRAIN,DSTRAN,STATEV,NDI,NSHR,
     *                NSTATV,NTENS)

      IMPLICIT REAL*8 (A-H,O-Z) 
      
      DIMENSION ST(6),SD(6,2),STRESS(NTENS),STRAIN(NTENS),DSTRAN(NTENS),
     *      STATEV(NSTATV)
                                               
      CALL CLEAR(ST,6) ! initializing
      CALL CLEAR(SD,12)  
                                               
      IF(NDI.EQ.3)THEN  ! if 3 hydrostatic comp.
          IF (NSHR.EQ.3) THEN  
              DO 10 I=1,6
                  ST(I)=STRESS(I)     ! stress
                  SD(I,1)=STRAIN(I)   ! strain
 10               SD(I,2)=DSTRAN(I)   ! inc. of strain
          ELSEIF (NSHR.EQ.1) THEN     ! plane strain
                  ST(1)=STRESS(1)     ! stress
                  ST(2)=STRESS(2)
                  ST(3)=STRESS(3)
                  ST(4)=STRESS(4)
                  SD(1,1)=STRAIN(1)   ! strain
                  SD(2,1)=STRAIN(2)
                  SD(3,1)=STRAIN(3)
                  SD(4,1)=STRAIN(4)
                  SD(1,2)=DSTRAN(1)   ! inc. of strain
                  SD(2,2)=DSTRAN(2)
                  SD(3,2)=DSTRAN(3)  
                  SD(4,2)=DSTRAN(4)
          ENDIF        
      ELSEIF(NDI.EQ.2)THEN    ! plane stress
          ST(1)=STRESS(1)     ! stress
          ST(2)=STRESS(2)
          ST(4)=STRESS(3)
          SD(1,1)=STRAIN(1)   ! strain
          SD(2,1)=STRAIN(2)
          SD(4,1)=STRAIN(3)
          SD(1,2)=DSTRAN(1)   ! inc. of strain
          SD(2,2)=DSTRAN(2)
          SD(4,2)=DSTRAN(3) 
          
          IF(NDI.EQ.2)THEN   ! plane stress
              SD(3,1)=STATEV(11)      ! 33-comp. of strain
              STATEV(11)=STATEV(11)+SD(3,2)
          ENDIF   
          
      ELSE
          WRITE(6,*)'THIS CASE IS NOT 1-, 2-, OR 3-D CASE.'
      ENDIF
      
      
      RETURN
      END   
      
c     ******************************************************************

c     this subroutine is to convert 1-D stiffness matrix and stress
c     tensor to 3/2-D ones.

      SUBROUTINE STF3_2(D3,D2,STRESS,ST,DDSDDT,DRPLDE,SDDT,RPLE,NDI,
     *                  NSHR,NTENS)

      IMPLICIT REAL*8 (A-H,O-Z) 
      
      DIMENSION D3(6,6),D2(NTENS,NTENS),STRESS(NTENS),ST(6),
     *          DDSDDT(NTENS),DRPLDE(NTENS),SDDT(6),RPLE(6)

c     assigning local strain and strain inc. components to global 
c     variables      
      
      IF(NDI.EQ.3)THEN  ! if 3 comp. of hydrostatic strain
          IF (NSHR.EQ.3) THEN
              CALL MTASSIGN(D3,D2,6,6)  
          ELSEIF (NSHR.EQ.1) THEN  ! plane strain             
              D2(1,1)=D3(1,1)
              D2(1,2)=D3(1,2)
              D2(1,3)=D3(1,3)
              D2(1,4)=D3(1,4)
              D2(2,1)=D3(2,1)
              D2(2,2)=D3(2,2)
              D2(2,3)=D3(2,3)
              D2(2,4)=D3(2,4)
              D2(3,1)=D3(3,1)
              D2(3,2)=D3(3,2)
              D2(3,3)=D3(3,3)
              D2(3,4)=D3(3,4)
              D2(4,1)=D3(4,1)
              D2(4,2)=D3(4,2)
              D2(4,3)=D3(4,3)
              D2(4,4)=D3(4,4)
          ENDIF              
      ELSEIF(NDI.EQ.2)THEN        ! plane stress
              D2(1,1)=D3(1,1)-D3(1,3)*D3(3,1)/D3(3,3)
              D2(2,2)=D3(2,2)-D3(2,3)*D3(3,2)/D3(3,3)
              D2(1,2)=D3(1,2)-D3(1,3)*D3(3,2)/D3(3,3)
              D2(2,1)=D3(2,1)-D3(2,3)*D3(3,1)/D3(3,3)
              D2(1,3)=D3(1,6)-D3(1,3)*D3(3,6)/D3(3,3)
              D2(2,3)=D3(2,6)-D3(2,3)*D3(3,6)/D3(3,3)
              D2(3,1)=D3(6,1)-D3(6,3)*D3(3,1)/D3(3,3)
              D2(3,2)=D3(6,2)-D3(6,3)*D3(3,2)/D3(3,3)
              D2(3,3)=D3(6,6)-D3(6,3)*D3(3,6)/D3(3,3)
      ELSE
          WRITE(6,*)'THIS CASE IS NOT 1-, 2-, OR 3-D CASE.'
      ENDIF 
                   
c     doing the same thing for stress                   
                   
      IF(NDI.EQ.3)THEN 
          IF (NSHR.EQ.3) THEN
              DO 30 I=1,6
                  DDSDDT(I)=SDDT(I)
                  DRPLDE(I)=RPLE(I)
 30               STRESS(I)=ST(I) 
          ELSEIF (NSHR.EQ.1) THEN
              DO 40 I=1,4
                 STRESS(I)=ST(I)
                 DDSDDT(I)=SDDT(I)
 40              DRPLDE(I)=RPLE(I)
          ENDIF              
      ELSEIF(NDI.EQ.2)THEN
          STRESS(1)=ST(1)
          STRESS(2)=ST(2)
          STRESS(3)=ST(4)
          DDSDDT(1)=SDDT(1)
          DDSDDT(2)=SDDT(2)
          DDSDDT(3)=SDDT(4)
          DRPLDE(1)=RPLE(1)
          DRPLDE(2)=RPLE(2)
          DRPLDE(3)=RPLE(4)
      ENDIF 
            
            
      RETURN
      END  
      
c     ******************************************************************


c     this subroutine assigns matrix a's values to matrix b      

      SUBROUTINE MTASSIGN(A,B,M,N)
      
      IMPLICIT REAL*8 (A-H,O-Z)  
      
      DIMENSION A(M,N),B(M,N)

      DO 10 I=1,M
          DO 10 J=1,N
 10           B(I,J)=A(I,J) 
 
      
      RETURN
      END 
      
c     ******************************************************************      

c     this subroutine initializes a real matrix to zero

      SUBROUTINE CLEAR(A,N)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      
      DIMENSION A(N)
      
      DO 10 I=1,N
 10       A(I)=0.D0 
 
 
      RETURN
      END   
      
c     ******************************************************************

c     this subroutine initializes an integer matrix to zero      

      SUBROUTINE ICLEAR(IA,N)   
      
      DIMENSION IA(N)  
      
      DO 10 I=1,N
 10       IA(I)=0 
 
 
      RETURN
      END   
      
c     ******************************************************************

c     this function calculates the dot product of two vectors      

      FUNCTION DOT(A,B,N)

      IMPLICIT REAL*8 (A-H,O-Z) 
      
      DIMENSION A(N),B(N)      
      
      D=0.D0    
      
      DO 10 I=1,N
 10       D=D+A(I)*B(I)
 
      DOT=D 
      
      
      RETURN
      END   
      

c     ******************************************************************
c     *********************** 1-D CASE STARTS HERE *********************
c     ******************************************************************

c     this subroutine is to perform an integration scheme to find the
c     real stresses under the given strain input, and then to update the
c     state variables, for a 1-d case.

      SUBROUTINE SMA_1(ST0_1,SD_1,DSD_1,STATEV,STATEML,SDDT_1,RPLE_1,
     *               RPL,DRPLDT,PROPS,DT_1,TEMP,DTEMP,TEMPINI,
     *                 DTIME,NSTATV,NPROPS,NUMMIN) 
     
      IMPLICIT REAL*8 (A-H,O-Z)
      
      COMMON /DATA02/MAXMTP,NELMTP,IMTP,ITER,KEYN,ISOLVE,KFLAG
      COMMON /PUNGA/ imli,MODEL
      
      DIMENSION STATEV(NSTATV),PROPS(NPROPS), STATEML(10,NUMMIN)
      
      PSI=STATEV(2)             ! martensitic volume fraction
      IPHASE=INT(STATEV(1)+0.1) ! info. about no yield/for./rev. trans. 
      SDT_1=STATEV(5)           ! assigning state variabbles to 
      RL_1=STATEV(20)


      TEMP0=TEMPINI ! initial temperature, To
      DTEMP0=TEMP+DTEMP-TEMP0 ! start temp.-To
      TEMP1=TEMP+DTEMP  ! present temperature
      ALPHA=ALFA(PROPS(14),PROPS(15),PSI) ! coeff. of thermal
                                               ! expansion 

c     here a previous inc. has not converged and reduction to a smaller
c     time increment has been performed, only partial analysis will be 
c     done due to futility, and a tangent or a elastic stiffness  
c     is returned.

      IF(KEYN.EQ.2)THEN  ! time increment has changed
         IF(STATEV(3).GT.0.D0)THEN ! if yielded
            CALL TANSTF_1(DT_1,SDDT_1,PSI,ST0_1,SDT_1,DTEMP0,PROPS,
     *                    STATEML,IPHASE,NPROPS,SDT_R,PSI_R,
     *                    MINORLOOP,NUMMIN) 
            ST0_1=ST0_1+SDDT_1*DTEMP ! update stress, no 
         ELSE                   ! no yield
            CALL ELASTF_1(PSI,PROPS,DT_1,NPROPS) ! cal. elast. stiff.
            CALL ELATTF_1(DT_1,ALPHA,SDDT_1)
            ST0_1=ST0_1+SDDT_1*DTEMP ! update stress.
         ENDIF
         
         RETURN                 ! return for next increment
          
      ENDIF

c     calculate the elastic predictor if no time inc. change

      CALL ELASTF_1(PSI,PROPS,DT_1,NPROPS)

      ST_1=DT_1*((SD_1+DSD_1)-ALPHA*DTEMP0-(SDT_1-PROPS(31)))     

c     check yield of the material point

      DELTA_A1=1.0D0/PROPS(12)-1.0D0/PROPS(11)
      DELTA_ALFA=PROPS(15)-PROPS(14)
      
      DOLD=ST0_1*RL_1+0.5D0*DELTA_A1*ST0_1**2+ST0_1*DELTA_ALFA
     *      *(TEMP-TEMP0)+PROPS(22)*(TEMP-PROPS(16))
      DNEW=ST_1*RL_1+0.5D0*DELTA_A1*ST_1**2+ST_1*DELTA_ALFA
     *      *(TEMP1-TEMP0)+PROPS(22)*(TEMP1-PROPS(16))

      DECIDE=DNEW-DOLD

      IF(DECIDE.GT.0.0D0)THEN
         IPHASE=1
         IFYD=1
      ELSEIF(DECIDE.LT.0.0D0)THEN
         IPHASE=2
         IFYD=1
      ELSE
         IFYD=0 
         GOTO 70
      ENDIF

         
c     if yielding happened, return mapping integration scheme will
c     be performed to find the real stress increment under the given
c     strain and temperature increments and current state input.

      CALL BEINT1_1(ST0_1,SDT_1,PSI,SD_1,DSD_1,TEMP,DTEMP,TEMP0,
     *              PROPS,STATEML,IPHASE,DT_1,SDDT_1,RPLE_1,RPL,
     *          DRPLDT,DTIME,NPROPS,SDT_R,PSI_R,HCUR,BACK,
     *              DRAG,YN1,RL_1,EP_1,APSI,MINORLOOP,NUMMIN,IFYD)

c     update stresses, and all state vaiables

 70   STATEV(1)=DFLOAT(IPHASE)  ! whether forward/reverse trans.
      STATEV(2)=PSI
      STATEV(3)=DFLOAT(IFYD)  

      STATEV(5)=SDT_1 ! assign transform. strain
      STATEV(11)=SD_1+DSD_1 ! 11-component of strain
      STATEV(12)=ST0_1  ! 11-component of stress               
      STATEV(14)=SDT_1
      STATEV(15)=TEMP1  ! assigning temp. 

      STATEV(16)=RL_1*ST0_1/PROPS(20) ! modified effective stress


      STATEV(26)=DTEMP0      

         
      STATEV(20)=RL_1
     
      RETURN
      END 
      
c     ******************************************************************   
      
c     this subroutine is to integrate the stress, inelastic strain, and
c     internal state variable, increment for gaven strain and temperature
c     increments by using Return Mapping Method (elastic predictor-plastic
c     corrector scheme), for a 1-d case.

      SUBROUTINE BEINT1_1(ST0_1,SDT0_1,PSI0,SD_1,DSD_1,TEMP,DTEMP,TEMP0,
     *                    PROPS,STATE,IPHASE,DT_1,SDDT_1,RPLE_1,RPL,
     *                  DRPLDT,DTIME,NPROPS,SDT_R,PSI_R,HCUR0,BACK0,
     *                    DRAG0,YN10,RL_10,EP0_1,APSI,MLOOP,NMLOOP,IFYD)

      IMPLICIT REAL*8 (A-H,O-Z)
      
      COMMON /DATA02/MAXMTP,NELMTP,IMTP,ITER,KEYN,ISOLVE,KFLAG 
      COMMON /PUNGA/ imli,MODEL
      
      DIMENSION PROPS(NPROPS),STATE(10,NMLOOP)
              
      TOL=PROPS(5)          ! criteria of convergence
      DTMP=DTEMP                ! temp. inc. for each inc. within iter.
      TMP=TEMP              ! starting temp.
      PSI=PSI0              ! martensitic volume frac.
      H=PROPS(20)
      
      ST_1=ST0_1      ! assigining stress, strain and strain inc.
      SDT_1=SDT0_1    ! values from sma to local vectors
      DD1_1=SD_1
      DD_1=DSD_1

      FLAG1F=1.0D0
      FLAG1R=0.0D0
 
      FLAG2R=0.001D0
      FLAG2F=0.99D0

c     to implement the return mapping integration scheme, calculate 
c     elastic predictor


      TMP=TMP+DTMP              ! present temp.
      DTMP0=TMP-TEMP0           ! present temp-To
      
      ALPHA=ALFA(PROPS(14),PROPS(15),PSI)
      
      CALL ELASTF_1(PSI,PROPS,DT_1,NPROPS) ! elas. stiff. mat.     
      CALL ELATTF_1(DT_1,ALPHA,SDDT_1)
      
      DD1_1=DD1_1+DD_1          ! updating strain and stress
      
      ST_1=DT_1*(DD1_1-ALPHA*DTMP0-(SDT_1-PROPS(31)))


c     iteration for plastic corrector

      DO 10 LOOP=1,35
         
         
         CALL GETRL_1(RL_1,ST_1,SDT_1,H,B2,IPHASE,ES_1,SDT_R,PSI_R,
     *                BACK)

c     calculate the value of the yield function

         CALL YDFUN_1(FYD,ST_1,ES_1,SDT_1,RL_1,TMP,DTMP0,PSI,PROPS,
     *                STATE,IPHASE,NPROPS,HCUR,SDT_R,PSI_R,BACK,DRAG,
     *                MLOOP,NMLOOP) 

         IF(LOOP.EQ.1)THEN
            IF(FYD.LE.0.0D0)THEN
               IFYD=0
               GOTO 50  
            ELSEIF(IPHASE.EQ.1.AND.PSI.GE.(FLAG1F-0.0001).AND.
     *              DABS(FYD).GT.0.0D0)THEN
               IFYD=0
               GOTO50
            ELSEIF(IPHASE.EQ.2.AND.PSI.LE.(FLAG1R+0.0001).AND.
     *              DABS(FYD).GT.0.0D0)THEN
               IFYD=0
               GOTO50
            ENDIF
         ENDIF

c     get rl, and b, also see papers/scheme              
            
         CALL GETQ_1(Q_1,DTMP0,RL_1,ST_1,PROPS,NPROPS)
         CALL GETB_1(B_1,PSI,Q_1,PROPS,STATE,DT_1,NPROPS,IPHASE,HCUR,
     *               RL_1,SDT_1,MLOOP,NMLOOP)
    
c     update mart. vol. frac. (PSI), trans. strain (SDT), and stress (ST)     
     
         DPSI=-FYD/B_1          ! inc. of mart. vol. frac.
         psii=psi
         PSI=PSI+DPSI           ! updating mart. vol. frac.
                                               
         if(psi.gt.FLAG1F)then
            if(iphase.eq.1)then
               psi=FLAG1F
               dpsi=psi-psii
            elseif(iphase.eq.2.and.psii.lt.FLAG1F)then
               psi=FLAG2R
               dpsi=psi-psii
            endif
         elseif(psi.lt.FLAG1R)then
            if(iphase.eq.1.and.psii.gt.FLAG1R)then
               psi=FLAG2F
               dpsi=psi-psii
            elseif(iphase.eq.2)then
               psi=FLAG1R
               dpsi=psi-psii
            endif
         endif
         
         SDT_1=SDT_1+RL_1*DPSI  ! updating trans. strain

         CALL ELASTF_1(PSI,PROPS,DT_1,NPROPS) ! elas. stiff. mat.     
         
         ALPHA=ALFA(PROPS(14),PROPS(15),PSI)
         ST_1=DT_1*(DD1_1-ALPHA*DTMP0-(SDT_1-PROPS(31)))
              

c     check convergence

         IF(DABS(DPSI).LT.TOL) GOTO 20

 10   CONTINUE

      WRITE(*,*)'ITERATION FAILS TO CONVERGE!',DPSI

 20   CONTINUE

c     update tangent stifness       

      if(iphase.eq.1.and.psi.lt.FLAG1F
     *     .or.iphase.eq.2.and.psi.gt.FLAG1R)then
         CALL TANSTF_1(DT_1,SDDT_1,PSI,ST_1,SDT_1,DTMP0,PROPS,STATE,
     *                 IPHASE,NPROPS,SDT_R,PSI_R,MLOOP,NMLOOP)
      endif
 
c     calling the subroutine to calculate the heat related terms

 50   CONTINUE 

      DPSI=PSI-PSI0

      IF(MODEL.EQ.2)THEN
         T=TEMP+DTEMP
         DTEMP0=TEMP+DTEMP-TEMP0
   
         CALL GETRL_1(RL_1,ST_1,SDT_1,H,B2,IPHASE,ES_1,SDT_R,PSI_R,BACK)
         CALL GETQ_1(Q_1,DTMP0,RL_1,ST_1,PROPS,NPROPS)

         CALL CALCRPL_1(RPL,RPLE_1,DRPLDT,PROPS,IPHASE,T,DSD_1,ST_1,
     *                   DT_1,SDDT_1,Q_1,DTEMP,DTIME,DPSI,NPROPS)
      ENDIF   
         
c     assign local state variables to global ones

      ST0_1=ST_1                ! stress 

      IF(IFYD.EQ.1)THEN
         PSI0=PSI
         SDT0_1=SDT_1           ! transformation strain

      ENDIF  
      
      RETURN
      END   
      
c     ******************************************************************
           
c     this subroutine is to get the tangent STIFFNESS  DT_1, and
c     thermal STIFFNESS coefficient, ALS_1.

      SUBROUTINE TANSTF_1(DT_1,ALS_1,PSI,ST_1,SDT_1,DTEMP0,PROPS,STATE,
     *                    IPHASE,NPROPS,SDT_R,PSI_R,MLOOP,NMLOOP)

      IMPLICIT REAL*8 (A-H,O-Z)  

      COMMON /PUNGA/ imli,MODEL
      
      DIMENSION PROPS(NPROPS),STATE(10,NMLOOP)

      H=PROPS(20)
     

      B2=PROPS(25)

      DALPH=PROPS(15)-PROPS(14)

c     get q, r, rl, s, and b, also see scheme                        
                       
      CALL GETRL_1(RL_1,ST_1,SDT_1,H,B2,IPHASE,ES_1,SDT_R,PSI_R,BACK)
      CALL GETQ_1(Q_1,DTEMP0,RL_1,ST_1,PROPS,NPROPS)
      CALL ELASTF_1(PSI,PROPS,DT_1,NPROPS) ! elastic stiff.
      CALL GETB_1(B_1,PSI,Q_1,PROPS,STATE,DT_1,NPROPS,IPHASE,H,RL_1,
     *            SDT_1,MLOOP,NMLOOP)
      CALL GETRS_1(R_1,S_1,Q_1,DT_1,DALPH,ST_1,PROPS,NPROPS,IPHASE)

      ALPHA=ALFA(PROPS(14),PROPS(15),PSI) ! coeff. of ther. exp.

      ALS_1=DT_1*(Q_1*S_1/B_1-ALPHA)  ! thermal stiffness vector
      DT_1=DT_1+DT_1*Q_1*R_1/B_1      ! tangent stifness vector

 
      RETURN
      END   
      
c     ******************************************************************


c     this subroutine calculates q for a 1-d case
c     also see scheme.           

      SUBROUTINE GETQ_1(Q_1,DTEMP0,RL_1,ST_1,PROPS,NPROPS) 
     
      IMPLICIT REAL*8 (A-H,O-Z)
      
      DIMENSION PROPS(NPROPS)
     
      Q_1=(1.0D0/PROPS(12)-1.0D0/PROPS(11))*ST_1+(PROPS(15)-PROPS(14))
     *     *DTEMP0+RL_1
        
         
      RETURN   
      END   


c*****************************************************************
      
c     this subroutine calculates b for a 1-d case
c     also see scheme.           

      SUBROUTINE GETB_1(B_1,PSI0,Q_1,PROPS,STATE,D_1,NPROPS,IPHASE,
     *                  HCUR,RL_1,SDT_1,MLOOP,NMLOOP)
    
      IMPLICIT REAL*8 (A-H,O-Z)
    
      COMMON /PUNGA/ imli,MODEL

      DIMENSION PROPS(NPROPS),STATE(10,NMLOOP) 
    
      PSI=PSI0  

      IF(PSI.GT.0.99993) PSI=0.99993
      IF(PSI.LT.0.00007) PSI=0.00007

      H=PROPS(20)
      
c     see scheme again                       
                       
      IF(IPHASE.EQ.1)THEN  ! forward transforamtion 
          IF (MODEL.EQ.1) THEN  ! model 1
              B_1=-Q_1*D_1*Q_1-PROPS(22)/PROPS(24)/(1.0D0-PSI)
          ELSEIF (MODEL.EQ.2) THEN  ! model 2
              B_1=-Q_1*D_1*Q_1-PROPS(24)
          ELSEIF (MODEL.EQ.3) THEN  ! model 3
              B_1=-Q_1*D_1*Q_1+2.0*PROPS(22)/PROPS(24)
     *          /DSQRT(1.0D0-(2.0*PSI-1.0D0)**2)     
          ENDIF        
      ELSEIF(IPHASE.EQ.2)THEN  ! reverse transformation 
          IF (MODEL.EQ.1) THEN  ! model 1
              B_1=Q_1*D_1*Q_1-PROPS(21)/PROPS(23)/PSI
          ELSEIF (MODEL.EQ.2) THEN  ! model 2
              B_1=Q_1*D_1*Q_1+PROPS(23)    
          ELSEIF (MODEL.EQ.3) THEN  ! model 3
              B_1=Q_1*D_1*Q_1-2.0*PROPS(21)/PROPS(23)
     *          /DSQRT(1.0D0-(2.0*PSI-1.0D0)**2)    
          ENDIF  
      ENDIF


      RETURN
      END


c***********************************************************************

c     this subroutine calculates r and s for a 2-d/3-d case
c     also see scheme. 


      SUBROUTINE GETRS_1(R_1,S_1,Q_1,D_1,DALPH,ST_1,PROPS,NPROPS,IPHASE)
    
      IMPLICIT REAL*8 (A-H,O-Z)
    
      COMMON /PUNGA/ imli,MODEL

      DIMENSION PROPS(NPROPS)
    
      H=PROPS(20)
 
      R_1=Q_1*D_1 
      IF(IPHASE.EQ.2) R_1=-R_1      

      ALPHA=ALFA(PROPS(14),PROPS(15),PSI)
        
c     see scheme again                       
                       
      IF(IPHASE.EQ.1)THEN  ! forward transforamtion 
      S_1=DALPH*ST_1+PROPS(22)-R_1*ALPHA
      ELSEIF(IPHASE.EQ.2)THEN  ! reverse transformation 
      S_1=-DALPH*ST_1-PROPS(21)-R_1*ALPHA
      ENDIF


      RETURN
      END


c     ******************************************************************

c     this subroutine calculates  effevtice stress (es_1)
                 
      SUBROUTINE GETRL_1(RL_1,ST_1,SDT_1,H,B2,IPHASE,ES_1,SDT_R,PSI_R,
     *                   BACK)
      
      IMPLICIT REAL*8 (A-H,O-Z)              
      
      COMMON /PUNGA/ imli,MODEL

      IF(INT(MODEL+0.1D0).EQ.4)THEN
         ES_1=ST_1+BACK
      ELSE   
         ES_1=ST_1-B2*SDT_1     ! B2 is the hardening factor
      ENDIF   

      IF(IPHASE.EQ.1)THEN
         IF(ES_1.EQ.0.0D0)THEN
            RL_1=1.0D0
         ELSEIF(ES_1.NE.0.0D0)THEN
            RL_1=H*ES_1/DABS(ES_1)
         ENDIF 
      ELSEIF(IPHASE.EQ.2)THEN
         IF(SDT_1.EQ.0.0D0)THEN
            RL_1=1.0D0
         ELSEIF(SDT_1.NE.0.0D0)THEN
            RL_1=H*SDT_1/DABS(SDT_1)
         ENDIF   
      ENDIF
        

      RETURN
      END

c     ******************************************************************

c     calculates the value of the yield function for a material point of 
c     a 1-d case.      
                  
      SUBROUTINE YDFUN_1(FYD,ST_1,ES_1,SDT_1,RL_1,TEMP1,DTEMP1,PSI0,
     *               PROPS,STATE,IPHASE,NPROPS,HCUR,SDT_R,PSI_R,BACK,
     *                 DRAG,MLOOP,NMLOOP)   
      
      IMPLICIT REAL*8 (A-H,O-Z)          

      COMMON /PUNGA/ imli,MODEL
      
      DIMENSION PROPS(NPROPS),STATE(10,NMLOOP) 

      FYD=0.D0  ! yield value of the function
      PSI=PSI0  ! martensitic volume fraction

      IF(PSI.GT.0.99993) PSI=0.99993
      IF(PSI.LT.0.00007) PSI=0.00007

      H=PROPS(20)               ! max. transformation strain
        
     
c     calculate  effective stress      

      DELTA_A1=1.0D0/PROPS(12)-1.0D0/PROPS(11)
      DELTA_ALFA=PROPS(15)-PROPS(14)

      IF(IPHASE.EQ.1)THEN   ! forward transforamtion
          IF (MODEL.EQ.1) THEN  ! model 1
              FYD=ES_1*RL_1+(1.0D0/2.0D0)*DELTA_A1*ST_1**2+DELTA_ALFA
     *            *ST_1*DTEMP1+PROPS(22)*TEMP1+PROPS(22)/PROPS(24)
     *          *DLOG(1.D0-PSI)-PROPS(22)*PROPS(16)
          ELSEIF (MODEL.EQ.2) THEN   ! model 2     
              FYD=ES_1*RL_1+(1.0D0/2.0D0)*DELTA_A1*ST_1**2
     *            +DELTA_ALFA*ST_1*DTEMP1+PROPS(22)*TEMP1-PROPS(24)
     *            *PSI-PROPS(22)*PROPS(16)
          ELSEIF (MODEL.EQ.3) THEN   ! model 3 
              FYD=ES_1*RL_1+(1.0D0/2.0D0)*DELTA_A1*ST_1**2
     *            +DELTA_ALFA*ST_1*DTEMP1+PROPS(22)*TEMP1-PROPS(22)
     *            /PROPS(24)*(DACOS(2.0D0*PSI-1.0D0)-22.0D0/7.0D0)
     *          -PROPS(22)*PROPS(16)         
          ENDIF               
      ELSEIF(IPHASE.EQ.2)THEN  ! reverse transformation
          IF (MODEL.EQ.1) THEN  ! model 1
              FYD=-ES_1*RL_1-(1.0D0/2.0D0)*DELTA_A1*ST_1**2
     *            -DELTA_ALFA*ST_1*DTEMP1-PROPS(21)*TEMP1-PROPS(21)
     *          /PROPS(23)*DLOG(PSI)+PROPS(21)*PROPS(18)
          ELSEIF (MODEL.EQ.2) THEN   ! model 2  
              FYD=-ES_1*RL_1-(1.0D0/2.0D0)*DELTA_A1*ST_1**2
     *            -DELTA_ALFA*ST_1*DTEMP1-PROPS(21)*TEMP1+PROPS(23)
     *            *PSI+PROPS(21)*PROPS(19)
          ELSEIF (MODEL.EQ.3) THEN   ! model 3
              FYD=-ES_1*RL_1-(1.0D0/2.0D0)*DELTA_A1*ST_1**2
     *            -DELTA_ALFA*ST_1*DTEMP1-PROPS(21)*TEMP1+PROPS(21)
     *          /PROPS(23)*(DACOS(2.0D0*PSI-1.0D0)-22.0D0/7.0D0)
     *          +PROPS(21)*PROPS(19)             
           ENDIF
      ENDIF

      RETURN
      END

c     ******************************************************************

c     this subroutine finds the elastic stiffness
 
      SUBROUTINE ELASTF_1(PSI,PROPS,D_1,NPROPS) 
      
      IMPLICIT REAL*8 (A-H,O-Z)  
      
      DIMENSION PROPS(NPROPS)

      D_1=1.0D0/PROPS(11)+PSI*(1.0D0/PROPS(12)-1.0D0/PROPS(11)) 
      D_1=1.0D0/D_1
            
          
      RETURN
      END

c     ******************************************************************

c     this subroutine finds the elastic thermal stiffness matrix

      SUBROUTINE ELATTF_1(D_1,ALPHA,TST_1)
    
      IMPLICIT REAL*8 (A-H,O-Z)

      TST_1=-D_1*ALPHA  
 
      RETURN
      END
    
     
c     ******************************************************************

c     this subroutine is used to calculate the variables associated with
c     the coupled mechanical-thermal analysis

      SUBROUTINE CALCRPL_1(RPL,RPLE_1,DRPLDT,PROPS,IPHASE,T,DSD_1,
     *                 ST_1,DT_1,SDDT_1,Q_1,DTEMP,DTIME,DPSI,NPROPS)

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION PROPS(NPROPS)

      TEMP=T

c     calling the subroutine to calculate the a and g vectors beside h

      CALL RPLEXT_1(G_1,HT_1,TEMP,ST_1,IPHASE,PROPS,Q_1,DPSI,RPLE_1,
     *          DRPLDT,DT_1,SDDT_1,NPROPS)

      RPLE_1=RPLE_1/DTIME   
      DRPLDT=DRPLDT/DTIME   

      DRPE_1=G_1*DT_1
      DRPT=G_1*SDDT_1+HT_1  
      RPL=DRPE_1*DSD_1+DRPT*DTEMP
      RPL=RPL/DTIME 


      RETURN
      END

      
c     ******************************************************************


      SUBROUTINE RPLEXT_1(G_1,HT_1,T,ST_1,IPHASE,PROPS,Q_1,DPSI,RPLE_1,
     *            DRPLDT,DT_1,SDDT_1,NPROPS)

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION PROPS(NPROPS)

c   only 2nd model now  

      H=PROPS(20)
      T=T+273.0D0
      DALPH=PROPS(15)-PROPS(14)
      ALPHA=ALFA(PROPS(14),PROPS(15),PSI)

      IF (IPHASE.EQ.1) THEN   ! forward transformation  
        A=Q_1
        B=DALPH*ST_1+PROPS(22)
        C=-PROPS(24)
        FYD=-1.0D0/2.0D0*PROPS(22)*(PROPS(19)-PROPS(16))
     *           -1.0D0/4.0D0*PROPS(22)*(PROPS(16)-PROPS(17)-PROPS(19)
     *            +PROPS(18))
        HELLO=FYD-T*DALPH*ST_1-PROPS(22)*T
        YELLO=PROPS(22)*T-FYD
        DELLO=PROPS(22)
      ELSEIF(IPHASE.EQ.2) THEN  ! reverse transformation
        A=-Q_1      
        B=-DALPH*ST_1-PROPS(21)
        C=PROPS(23)
        FYD=1.0D0/2.0D0*PROPS(21)*(PROPS(19)-PROPS(16))
     *          +1.0D0/4.0D0*PROPS(21)*(PROPS(16)-PROPS(17)-PROPS(19)
     *           +PROPS(18))
        HELLO=FYD-T*DALPH*ST_1-PROPS(21)*T
        YELLO=PROPS(21)*T-FYD
        DELLO=PROPS(21)
      ENDIF     

      RPLE_1=A*DT_1*YELLO/C
      DRPLDT=YELLO/C*(A*SDDT_1+B)-DELLO*DPSI                
    
      G_1=-(T*ALPHA+A*HELLO/C)
      HT_1=-B*HELLO/C

  
      RETURN    
      END 


c     ******************************************************************
                                                            
c     this subroutine is to ouput the material data for a 1-d case.
      
      SUBROUTINE FORMDD_1(PROPS,NPROPS)

      IMPLICIT REAL*8 (A-H,O-Z)  
      
      COMMON /PUNGA/ imli,MODEL

      DIMENSION PROPS(NPROPS)
    
     
      RETURN
      END   
      
c     ******************************************************************
                                                        
c     this subroutine assigns global values to local values
                                                        
      SUBROUTINE STD_1(STRESS,STRAN,DSTRAN,F0,DROT,STATEV,PROPS,
     *                 NSTATV,NPROPS,ST_1,SD_1,DSD_1)
      
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION F0(3,3),STATEV(NSTATV),PROPS(NPROPS)
      
      IF(INT(PROPS(37)+0.1D0).EQ.0)THEN
         ST_1=STRESS
         SD_1=STRAN
         DSD_1=DSTRAN
      ELSE IF(INT(PROPS(37)+0.1D0).EQ.1)THEN
         PSI=STATEV(2)
         V=PROPS(13)
         ET=STATEV(5)
         DTEMP0=STATEV(26)
         ALPHA=ALFA(PROPS(14),PROPS(15),PSI)
         F0(2,2)=-V*(F0(1,1)**2-1.0D0)+ET*(2.0D0*V-1.0D0)+2.0D0
     *            *ALPHA*DTEMP0*(V+1.0D0)+1.0D0
         F0(2,2)=DSQRT(F0(2,2))
         F0(3,3)=F0(2,2)
         DET=F0(1,1)*F0(2,2)*F0(3,3)

         IF(DROT.EQ.0.0d0) DROT=1.0D0 ! 1-d patch

         ST_1=DET*(STRESS/F0(1,1)**2)/DROT**2
         SD_1=0.5D0*(F0(1,1)**2-1.0D0)
         R=DROT*STATEV(61)
         DSD_1=R*(DSTRAN+0.5D0*DSTRAN**2+(DSTRAN**3/6.0D0))*R
         STATEV(61)=R
      ENDIF   
      
      RETURN
      END
      
c     ******************************************************************      

c     this subroutine assigns local values to global values
 
      SUBROUTINE STF_1(DT_1,DDSDDE,STRESS,ST_1,DDSDDT,DRPLDE
     *           ,SDDT_1,RPLE_1,F1,RPL,DRPLDT,STATEV,PROPS,
     *               NSTATV,NPROPS)
      
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION PROPS(NPROPS),STATEV(NSTATV),F1(3,3)
      
      IF(INT(PROPS(37)+0.1D0).EQ.0)THEN
         STRESS=ST_1
         DDSDDE=DT_1
         DDSDDT=SDDT_1
         DRPLDE=RPLE_1
      ELSE IF(INT(PROPS(37)+0.1D0).EQ.1)THEN
         PSI=STATEV(2)
         V=PROPS(13)
         ET=STATEV(5)
         DTEMP0=STATEV(26)
         ALPHA=ALFA(PROPS(14),PROPS(15),PSI)
         F1(2,2)=-V*(F1(1,1)**2-1.0D0)+ET*(2.0D0*V-1.0D0)+2.0D0
     *            *ALPHA*DTEMP0*(V+1.0D0)+1.0D0
         F1(2,2)=DSQRT(F1(2,2))
         F1(3,3)=F1(2,2)
         DET=F1(1,1)*F1(2,2)*F1(3,3)   
         STRESS=ST_1*F1(1,1)**2/DET
         DDSDDE=F1(1,1)**2*DT_1*F1(1,1)**2/DET
         DDSDDT=F1(1,1)**2*SDDT_1/DET
         DRPLDE=RPLE_1*F1(1,1)**2
         DRPLDT=DRPLDT/DET
         RPL=RPL/DET
      ENDIF   
         
      RETURN
      END 

c     ******************************************************************      

c     this function calculates the present thermal expansion coefficient

      FUNCTION ALFA(ALPHAA,ALPHAM,PSI)
    
      IMPLICIT REAL*8 (A-H,O-Z)
    
      ALFA=ALPHAA+PSI*(ALPHAM-ALPHAA)
    
      RETURN
      END 

c     ******************************************************************
c     ******************************************************************
