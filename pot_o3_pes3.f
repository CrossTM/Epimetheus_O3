
      SUBROUTINE potv(V,R1,R2,xcos)

      IMPLICIT DOUBLE PRECISION (A-H,O-Y),logical(z)
      COMMON /MASS/ XMASS(3),G1,G2,xmassr(3)
      DATA X1/1.0D0/,X0/0.0D0/,TINY/9.0D-15/,X2/2.0D0/

      pi=dacos(-1.d0)

      IF (G1 .EQ. X0) THEN
         Q1 = R1
         Q2 = R2
         THETA = ACOS(XCOS)
      ELSE IF (G2 .EQ. X0) THEN
         XX = R1 * G1
         YY = R1 * (X1 - G1)
         IF (R2 .EQ. X0 .OR. XCOS .GE. (X1 - TINY)) THEN
            Q1 = ABS(XX - R2)
            Q2 = (YY + R2)
            COST = -X1
         ELSE IF (XCOS .LE. (TINY - X1)) THEN
            Q1 = (XX + R2)
            Q2 = ABS(YY + R2)
            COST = X1
         ELSE
            Q1 = SQRT(XX*XX + R2*R2 - X2*XX*R2*XCOS)
            Q2 = SQRT(YY*YY + R2*R2 + X2*YY*R2*XCOS)
            COST = (Q1**2 + Q2**2 - R1**2) / (X2 * Q1 * Q2)
         ENDIF
         THETA = ACOS(COST)
      ELSE
         F1= X1/G1
         F2= X1/G2
         F12= X1 - F1*F2
         P1= R1*(X1-F1)/(G2*F12)
         P2= R2*(X1-F2)/(G1*F12)
         S1= R1-P1
         S2= R2-P2
         Q1= SQRT(P1*P1 + S2*S2 + X2*P1*S2*XCOS)/(X1-G1)
         Q2= SQRT(P2*P2 + S1*S1 + X2*P2*S1*XCOS)/(X1-G2)
         Q3= SQRT(P1*P1 + P2*P2 - X2*P1*P2*XCOS)
         COST = (Q1*Q1 + Q2*Q2 - Q3*Q3)/(X2*Q1*Q2)
         THETA = ACOS(COST)
      ENDIF

      att1=Q1*0.5291772d0
      att2=Q2*0.5291772d0
        
      call poten(vp,att1,att2,THETA)
      v=vp/219474.624d0

      end

!  wifin from ai 10.2Dr
      subroutine poten(v,r1,r2,th)

      implicit none
      double precision  r1,r2,th,v
      double precision  reoo,thetae,b1,alpha
      double precision  phh2,xs1,xs2,xst,rr1,rr2,xep1,xep2,xep3,xep4
      double precision  v0,vp,voo,roo
      
           reoo=1.28200d0
           thetae=116.8800d0
           b1=2.15d0
           phh2=3.20164303995d0

      thetae=thetae*.314159265358979312d01*.00555555555555555555d0

      xs1=(r1+r2)*0.5d0-reoo
      xs2=(r1-r2)*0.5d0
      xst=dcos(th)-dcos(thetae)

      rr1=r1-reoo
      rr2=r2-reoo

      alpha=3.31d0
      xep1=(dexp(-2.d0*alpha*rr1)-2.d0*dexp(-alpha*rr1)+1.d0)*11500.d0

      xep2=(dexp(-2.d0*alpha*rr2)-2.d0*dexp(-alpha*rr2)+1.d0)*11500.d0

      xep4=dexp(-10.d0*(th-thetae))
      xep3=dexp(-b1*((r1-reoo)**2+(r2-reoo)**2))
      roo=dsqrt(r1**2+r2**2-2.d0*r1*r2*dcos(th))
      voo=0.800D+06*dexp(-phh2*roo)

      v0 = -0.66435358784700D+03*xs1**0*xs2**0*xst**0
      vp = -0.25325058268700D+04*xs1**0*xs2**0*xst**1
     *+0.11120062755000D+05*xs1**1*xs2**0*xst**0
     *+0.61765305206700D+05*xs1**0*xs2**0*xst**2
     *-0.45780288856700D+05*xs1**0*xs2**2*xst**0
     *-0.46596921143500D+05*xs1**1*xs2**0*xst**1
     *+0.99721915480700D+05*xs1**2*xs2**0*xst**0
     *+0.46015252747700D+04*xs1**0*xs2**0*xst**3
     *+0.77881395192400D+05*xs1**0*xs2**2*xst**1
     *-0.18982404064800D+06*xs1**1*xs2**0*xst**2
     *+0.16660920043300D+06*xs1**1*xs2**2*xst**0
     *+0.15734360487300D+06*xs1**2*xs2**0*xst**1
     *-0.79763939269300D+05*xs1**3*xs2**0*xst**0
     *+0.43555606631600D+05*xs1**0*xs2**0*xst**4
     *+0.15974290450200D+06*xs1**0*xs2**2*xst**2
     *-0.76800022571700D+06*xs1**0*xs2**4*xst**0
     *-0.39029977392200D+05*xs1**1*xs2**0*xst**3
     *+0.21320955750600D+06*xs1**1*xs2**2*xst**1
     *+0.61354058659800D+06*xs1**2*xs2**0*xst**2
     *-0.13416067343700D+07*xs1**2*xs2**2*xst**0
     *-0.34758398470700D+06*xs1**3*xs2**0*xst**1
     *+0.39314461542800D+06*xs1**4*xs2**0*xst**0
     *-0.63063261766500D+05*xs1**0*xs2**0*xst**5
     *+0.28778210638700D+06*xs1**0*xs2**2*xst**3
     *-0.70629161067600D+06*xs1**0*xs2**4*xst**1
     *-0.15657092877900D+05*xs1**1*xs2**0*xst**4
     *-0.98320680679800D+06*xs1**1*xs2**2*xst**2
     *+0.19905188575100D+07*xs1**1*xs2**4*xst**0
     *-0.29833663403000D+06*xs1**2*xs2**0*xst**3
     *-0.35683419758700D+06*xs1**2*xs2**2*xst**1
     *+0.31204774303900D+06*xs1**3*xs2**0*xst**2
     *+0.92516264398000D+06*xs1**3*xs2**2*xst**0
     *-0.65465044716300D+06*xs1**4*xs2**0*xst**1
     *+0.14801095319700D+06*xs1**5*xs2**0*xst**0
     *+0.10989083107400D+06*xs1**0*xs2**0*xst**6
     *+0.18813319298700D+07*xs1**0*xs2**6*xst**0

       v=v0+vp*xep3+voo+xep1+xep2+xep4
       if(th*180/3.14<50.d0) v=30000.d0

      end subroutine poten
