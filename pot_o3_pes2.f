
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

      v0 = -0.66671958940200D+03*xs1**0*xs2**0*xst**0
      vp = -0.24199738812500D+04*xs1**0*xs2**0*xst**1
     *+0.10722773838100D+05*xs1**1*xs2**0*xst**0
     *+0.62096281361000D+05*xs1**0*xs2**0*xst**2
     *-0.44608770711800D+05*xs1**0*xs2**2*xst**0
     *-0.46601019614100D+05*xs1**1*xs2**0*xst**1
     *+0.10069243428400D+06*xs1**2*xs2**0*xst**0
     *+0.44931212309300D+04*xs1**0*xs2**0*xst**3
     *+0.72068647625700D+05*xs1**0*xs2**2*xst**1
     *-0.18358442305500D+06*xs1**1*xs2**0*xst**2
     *+0.15542891601900D+06*xs1**1*xs2**2*xst**0
     *+0.13697760022600D+06*xs1**2*xs2**0*xst**1
     *-0.77951001090000D+05*xs1**3*xs2**0*xst**0
     *+0.38539454203600D+05*xs1**0*xs2**0*xst**4
     *+0.12736400889900D+06*xs1**0*xs2**2*xst**2
     *-0.74830417713800D+06*xs1**0*xs2**4*xst**0
     *-0.54952308462100D+05*xs1**1*xs2**0*xst**3
     *+0.28771443254800D+06*xs1**1*xs2**2*xst**1
     *+0.69191197175200D+06*xs1**2*xs2**0*xst**2
     *-0.12140625865000D+07*xs1**2*xs2**2*xst**0
     *-0.36968132337800D+06*xs1**3*xs2**0*xst**1
     *+0.39540262607300D+06*xs1**4*xs2**0*xst**0
     *-0.73147161076900D+05*xs1**0*xs2**0*xst**5
     *+0.41837210967700D+06*xs1**0*xs2**2*xst**3
     *-0.60979169135100D+06*xs1**0*xs2**4*xst**1
     *-0.15657092877900D+05*xs1**1*xs2**0*xst**4
     *-0.13151291140200D+07*xs1**1*xs2**2*xst**2
     *+0.17207206957400D+07*xs1**1*xs2**4*xst**0
     *-0.29671711706300D+06*xs1**2*xs2**0*xst**3
     *-0.62129931201900D+06*xs1**2*xs2**2*xst**1
     *+0.56321588223200D+06*xs1**3*xs2**0*xst**2
     *+0.82389893488000D+06*xs1**3*xs2**2*xst**0
     *-0.39802228339200D+06*xs1**4*xs2**0*xst**1
     *-0.12649250581000D+06*xs1**5*xs2**0*xst**0
     *+0.77019440555400D+05*xs1**0*xs2**0*xst**6
     *+0.18874620717300D+07*xs1**0*xs2**6*xst**0

       v=v0+vp*xep3+voo+xep1+xep2+xep4
       if(th*180/3.14<50.d0) v=30000.d0

      end subroutine poten
