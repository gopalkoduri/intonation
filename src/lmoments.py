"""
##LMOMENT PYTHON LIBRARY:

This file contains the lmoments.f library created by:
      J. R. M. HOSKING                                                   
      IBM RESEARCH DIVISION                                              
      T. J. WATSON RESEARCH CENTER                                       
      YORKTOWN HEIGHTS                                                   
      NEW YORK 10598, U.S.A.      
      AUGUST    1996

The base Fortran code is copyright of the IBM Corperation, and the licensing
information is shown below:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IBM software disclaimer

LMOMENTS: Fortran routines for use with the method of L-moments
Permission to use, copy, modify and distribute this software for any purpose
and without fee is hereby granted, provided that this copyright and permission
notice appear on all copies of the software. The name of the IBM Corporation
may not be used in any advertising or publicity pertaining to the use of the
software. IBM makes no warranty or representations about the suitability of the
software for any purpose. It is provided "AS IS" without any express or implied
warranty, including the implied warranties of merchantability, fitness for a
particular purpose and non-infringement. IBM shall not be liable for any direct,
indirect, special or consequential damages resulting from the loss of use,
data or projects, whether in an action of contract or tort, arising out of or
in connection with the use or performance of this software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

All methodologies in this library are repeated from this source file verbatim,
with the exception of the samlmu() function.  This was redesigned to take
advantage of pythons language, as direct FORTRAN conversion of the method
did not translate speed.

This library was designed to use L-moments to predict optimal parameters
for a number of distributions.  Distributions supported in this file are
listed below, with their distribution suffix:
    *Exponential (EXP)
    *Gamma (GAM)
    *Generalised Extreme Value (GEV)
    *Generalised Logistic (GLO)
    *Generalised Normal (GNO)
    *Generalised Pareto (GPA)
    *Gumbel (GUM)
    *Kappa (KAP)
    *Normal (NOR)
    *Pearson III (PE3)
    *Wakeby (WAK)
    *Weibull (WEI)

The primary function in this file is the samlmu(x,nmom) function, which takes
an input dataset x and  input of the number of moments to produce the log
moments of that dataset.

For Instance, given a list "Data", if 5 l-moments are needed, the function
would be called by lmom.samlmu(Data,5)

In this file contains four different functions for using each distribution.
Each function can be called by the prefix FUN with the suffix DIS.

*PEL: (x,nmom):
      Parameter Estimates.  This takes the L-Moments calculated by samlmu()
      and predicts the parameter estimates for that function.
      
      EXAMPLE: Find Wakeby distribution that best fits dataset DATA:

      import lmom
      para = lmom.pelwak(lmom.samlmu(DATA,5))

*QUA: (f,para)
      Quartile Estimates.  This takes the parameter estimates for a
      distribution, and a given quartile value to calculate the quartile for the
      given function.

      EXAMPLE: Find the Upper Quartile (75%) of the Kappa distribution that
      best fits dataset DATA:

      import lmom
      para = lmom.pelkap(lmom.samlmu(DATA,5))
      UQ = lmom.quakap(0.75,para)

*LMR: (para,nmom):
      L-Moment Ratios.  This takes the parameter estimates for a distribution
      and calculates nmom L-Moment ratios.

      EXAMPLE: Find 4 lmoment ratios for the Gumbel distribution that
      best fits dataset DATA:

      import lmom
      para = lmom.pelgum(lmom.samlmu(DATA,5))
      LMR = lmom.lmrgum(para,4)

*CDF: (x,para):
      Cumulative Distribution Function.  This takes the parameter estimates
      for a distribution and calculates the quartile for a given value x.

      EXAMPLE: Find the quartile of the datapoint 6.4 for the Weibull
      Distribution that best fits the dataset DATA:

      import lmom
      para = lmom.pelwei(lmom.samlmu(DATA,5))
      quartile = lmom.quawei(6.4,para)

Python Translation conducted by:
    Sam Gillespie
    Numerical Analyst
    C&R Consulting
    Townsville Australia
    June 2013

For more information, or to report bugs, contact:
    sam@candrconsulting.com.au

Licensing for Python Translation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Copyright (c) 2013 C&R Consulting.
All rights reserved.

Redistribution and use in source and binary forms are permitted
provided that the above copyright notice and this paragraph are
duplicated in all such forms and that any documentation,
advertising materials, and other materials related to such
distribution and use acknowledge that the software was developed
by the C&R Consulting.  The name of the
C&R Consulting may not be used to endorse or promote products derived
from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""

import scipy as sp
import scipy.special as spsp

################################################################
##L-MOMENT CALCULATION FUNCTION samlmu
################################################################
def comb(N,k,exact=1):
    if exact:
        if (k > N) or (N < 0) or (k < 0):
            return 0
        val = 1
        for j in xrange(min(k, N-k)):
            val = (val*(N-j))//(j+1)
        return val
    else:
        from scipy import special
        k,N = sp.asarray(k), sp.asarray(N)
        lgam = special.gammaln
        cond = (k <= N) & (N >= 0) & (k >= 0)
        sv = special.errprint(0)
        vals = sp.exp(lgam(N+1) - lgam(N-k+1) - lgam(k+1))
        sv = special.errprint(sv)
        return sp.where(cond, vals, 0.0)


def samlmu(x,nmom=5):
    x = sorted(x)
    n = len(x)   
    ##Calculate first order
    ##Pretty efficient, no loops
    coefl1 = 1.0/comb(n,1)
    suml1 = sum(x)
    l1 = coefl1*suml1

    if nmom == 1:
        ret = l1
        return(ret)

    ##Calculate Second order

    #comb terms appear elsewhere, this will decrease calc time
    #for nmom > 2, and shouldn't decrease time for nmom == 2
    #comb1 = comb(i-1,1)
    #comb2 = comb(n-i,1)
    comb1 = []
    comb2 = []
    for i in range(1,n+1):
        comb1.append(comb(i-1,1))
        comb2.append(comb(n-i,1))
    
    coefl2 = 0.5 * 1.0/comb(n,2)
    xtrans = []
    for i in range(1,n+1):
        coeftemp = comb1[i-1]-comb2[i-1]
        xtrans.append(coeftemp*x[i-1])
    
    l2 = coefl2 * sum(xtrans)

    if nmom  ==2:
        ret = [l1,l2]
        return(ret)

    ##Calculate Third order
    #comb terms appear elsewhere, this will decrease calc time
    #for nmom > 2, and shouldn't decrease time for nmom == 2
    #comb3 = comb(i-1,2)
    #comb4 = comb(n-i,2)
    comb3 = []
    comb4 = []
    for i in range(1,n+1):
        comb3.append(comb(i-1,2))
        comb4.append(comb(n-i,2))
    
    coefl3 = 1.0/3 * 1.0/comb(n,3)
    xtrans = []
    for i in range(1,n+1):
        coeftemp = (comb3[i-1]-
                    2*comb1[i-1]*comb2[i-1] +
                    comb4[i-1])
        xtrans.append(coeftemp*x[i-1])

    l3 = coefl3 *sum(xtrans) /l2

    if nmom  ==3:
        ret = [l1,l2,l3]
        return(ret)

    ##Calculate Fourth order
    #comb5 = comb(i-1,3)
    #comb6 = comb(n-i,3)
    comb5 = []
    comb6 = []
    for i in range(1,n+1):
        comb5.append(comb(i-1,3))
        comb6.append(comb(n-i,3))
    
    coefl4 = 1.0/4 * 1.0/comb(n,4)
    xtrans = []
    for i in range(1,n+1):
        coeftemp = (comb5[i-1]-
                    3*comb3[i-1]*comb2[i-1] +
                    3*comb1[i-1]*comb4[i-1] -
                    comb6[i-1])
        xtrans.append(coeftemp*x[i-1])

    l4 = coefl4 *sum(xtrans)/l2

    if nmom  ==4:
        ret = [l1,l2,l3,l4]
        return(ret)

    ##Calculate Fifth order
    coefl5 = 1.0/5 * 1.0/comb(n,5)
    xtrans = []
    for i in range(1,n+1):
        coeftemp = (comb(i-1,4)-
                    4*comb5[i-1]*comb2[i-1] +
                    6*comb3[i-1]*comb4[i-1] -
                    4*comb1[i-1]*comb6[i-1] +
                    comb(n-i,4))
        xtrans.append(coeftemp*x[i-1])

    l5 = coefl5 *sum(xtrans)/l2

    if nmom ==5:
        ret = [l1,l2,l3,l4,l5]
        return(ret)

#######################################################
#CDF FUNCTIONS
#######################################################

def cdfexp(x,para):
    U = para[0]
    A = para[1]
    if A <= 0:
        cdfexp = 0
        print("Parameters Invalid")
        return(cdfexp)
    else:
        Y = (x-U)/A
        if U <= 0:
            cdfexp = 0
            print("Parameters Invalid")
            return(cdfexp)
        else:
            cdfexp = 1-sp.exp(-Y)
            if cdfexp >= 0:
                return(cdfexp)
            else:
                return(0)
            
#############################################################
            
def cdfgam(x,para):
    CDFGAM=0
    Alpha=para[0]
    Beta=para[1]
    if Alpha <= 0 or Beta <= 0:
        print("Parameters Invalid")
        return
    if x <= 0:
        print("x Parameter Invalid")
        return
    CDFGAM = spsp.gammainc(Alpha,x/Beta)
    return(CDFGAM)

#############################################################

def cdfgev(x,para):
    SMALL = 1e-15
    U=para[0]
    A=para[1]
    G=para[2]
    if A <= 0:
        print("Parameters Invalid")
        return
    Y = (x-U)/A
    if G==0:
        CDFGEV = sp.exp(-sp.exp(-Y))
    else:
        Arg = 1-G*Y
        if Arg > SMALL:
            Y = -sp.log(Arg)/G
            CDFGEV = sp.exp(-sp.exp(-Y))
        elif G<0:
            CDFGEV = 0
        else:
            CDFGEV = 1

    return(CDFGEV)

#############################################################

def cdfglo(x,para):
    SMALL = 1e-15
    U=para[0]
    A=para[1]
    G=para[2]
    if A <= 0:
        print("Parameters Invalid")
        return
    Y = (x-U)/A
    if G==0:
        CDFGLO=1/(1+sp.exp(-Y))
    else:
        Arg = 1-G*Y
        if Arg > SMALL:
            Y = -sp.log(Arg)/G
            CDFGLO=1/(1+sp.exp(-Y))
        elif G<0:
            CDFGLO = 0
        else: 
            CDFGLO = 1
            
    return(CDFGLO)

#############################################################

def cdfgno(x,para):
    SMALL = 1e-15
    U=para[0]
    A=para[1]
    G=para[2]
    if A <= 0:
        print("Parameters Invalid")
        return
    Y = (x-U)/A
    if G==0:
        CDFGNO = 0.5+0.5*sp.erg(Y*1/sp.sqrt(2))
    else:
        Arg = 1-G*Y

        if Arg > SMALL:
            Y = -sp.log(Arg)/G
            CDFGNO = 0.5+0.5*spsp.erf(Y*1/sp.sqrt(2))
        elif G<0:
            CDFGNO = 0
        else:
            CDFGNO = 1

    return(CDFGNO)

#############################################################

def cdfgpa(x,para):
    SMALL = 1e-15
    U=para[0]
    A=para[1]
    G=para[2]
    CDFGPA = 0
    if A <= 0:
        print("Parameters Invalid")
        return
    Y = (x-U)/A
    if Y <= 0:
        print("Parameters Invalid")
        return
    if G==0:
        CDFGPA=1-sp.exp(-Y)
    else:
        Arg = 1-G*Y
        if Arg > SMALL:
            Y = -sp.log(Arg)/G
            CDFGPA=1-sp.exp(-Y)
        else:
            CDFGPA = 1

    return(CDFGPA)

#############################################################

def cdfgum(x,para):
    U = para[0]
    A = para[1]
    if A <= 0:
        print("Parameters Invalid")
        return
    else:
        Y = (x-U)/A
        CDFGUM = sp.exp(-sp.exp(-Y))
        return(CDFGUM)

#############################################################
    
def cdfkap(x,para):
    SMALL = 1e-15
    U = para[0]
    A = para[1]
    G = para[2]
    H = para[3]

    if A <= 0:
        print("Invalid Parameters")
        return
    Y = (x-U)/A
    if G == 0:
        pass
    else:
        ARG = 1-G*Y
        if ARG > SMALL:
            pass
        else:
            if G < 0:
                CDFKAP = 0
            if G > 0:
                CDFKAP = 1
            return(CDFKAP)

        Y = -sp.log(ARG)/G
    Y = sp.exp(-Y)
    if H == 0:
        CDFKAP = sp.exp(-Y)
    else:
        ARG = 1-H*Y
        if ARG > SMALL:
            Y = -sp.log(ARG)/H
            CDFKAP = sp.exp(-Y)
            return(CDFKAP)
        else:
            CDFKAP = 0
            return(CDFKAP)

#############################################################
        
def cdfnor(x,para):
    if para[1] < 0:
        print("Invalid Parameters")
    cdfnor = 0.5+0.5*spsp.erf((x-para[0])/para[1]*1.0/sp.sqrt(2))
    return(cdfnor)

#############################################################

def cdfpe3(x,para):
    SMALL = 1e-6
    CDFPE3 = 0
    if para[1]<= 0:
        print("Parameters Invalid")
        return
    else:
        Gamma = para[2]
        if abs(Gamma) <= SMALL:
            Z = (x-para[0])/para[1]
            CDFPE3 = 0.5+0.5*spsp.erf(Z*1/sp.sqrt(2))
            return(CDFPE3)
        else:
            Alpha = 4/(Gamma**2)
            Z = 2*(x-para[0])/(para[1]*Gamma)+Alpha
            if Z > 0:
                CDFPE3 = spsp.gammainc(Alpha,Z)
            if Gamma < 0:
                CDFPE3 = 1-CDFPE3
            return(CDFPE3)

        
#############################################################

        
def cdfwak(x,para):
    
    EPS = 1e-8
    MAXIT = 20
    ZINCMX =3
    ZMULT = 0.2
    UFL = -170
    XI = para[0]
    A = para[1]
    B = para[2]
    C = para[3]
    D = para[4]

    if B+D <= 0 and (B!=0 or C!=0 or D!= 0):
        print("Invalid Parameters")
        return
    if A == 0 and B!= 0:
        print("Invalid Parameters")
        return
    if C == 0 and D != 0:
        print("Invalid Parameters")
        return
    if C < 0 or A+C < 0:
        print("Invalid Parameters")
        return
    if A == 0 and C == 0:
        print("Invalid Parameters")
        return

    CDFWAK = 0
    if x <= XI:
        return(CDFWAK)

    #Test for special cases
    if B == 0 and C == 0 and D == 0:
        
        Z = (x-XI)/A
        CDFWAK = 1
        if -Z >= UFL:
            CDFWAK = 1-sp.exp(-Z)
        return(CDFWAK)
    
    if C == 0:
        CDFWAK = 1
        if x >= (XI+A/B):
            return(CDFWAK)
        Z = -sp.log(1-(x-XI)*B/A)/B
        if -Z >= UFL:
            CDFWAK = 1-sp.exp(-Z)
        return(CDFWAK)

    
    if A == 0:
        Z = sp.log(1+(x-XI)*D/C)/D
        if -Z >= UFL:
            CDFWAK = 1-sp.exp(-Z)
        return(CDFWAK)


    CDFWAK=1
    if D <0 and x >= (XI+A/B-C/D):
        return(CDFWAK)

    Z=0.7
    if x < quawak.quawak(0.1,para):
        Z = 0
    if x < quawak.quawak(0.99,para):
        pass
    else:
        if D < 0:
            Z = sp.log((x-XI-A/B)*D/C+1)/D
        if D == 0:
            Z = (x-XI-A/B)/C
        if D > 0:
            Z = sp.log((x-XI)*D/C+1)/D
            
    for IT in range(1,MAXIT+1):
        EB = 0
        BZ = -B*Z
        if BZ >= UFL:
            EB = sp.exp(BZ)
        GB = Z
        
        if abs(B)>EPS:
            GB = (1-EB)/B
        ED = sp.exp(D*Z)
        GD = -Z

        if abs(D)>EPS:
            GD = (1-ED)/D

        XEST =XI +A*GB-C*GD
        FUNC = x-XEST
        DERIV1 = A*EB+C*ED
        DERIV2 = -A*B*EB+C*D*ED
        TEMP = DERIV1+0.5*FUNC*DERIV2/DERIV1

        if TEMP <= 0:
            TEMP = DERIV1
        ZINC = FUNC/TEMP

        if ZINC > ZINCMX:
            ZINC = ZINCMX

        ZNEW = Z+ZINC

        if ZNEW <= 0:
            Z = Z*ZMULT
        else:
            Z = ZNEW
            if abs(ZINC) <= EPS:
                CDFWAK = 1
                if -Z >= UFL:
                    CDFWAK = 1-sp.exp(-Z)
                return(CDFWAK)

#############################################################
            
def cdfwei(x,para):
    U = para[0]
    A = para[1]
    G = para[2]

    if len(para) != 3:
        print("Invalid number of parameters")
        return
    elif para[1] <= 0 or para[2] <= 0:
        print("Invalid Parameters")
        return
    else:
        cdfwei = 1-sp.exp(-((x-para[0])/para[1])**para[2])
        return(cdfwei)
    
#############################################################
#LMR FUNCTIONS
#############################################################

def lmrexp(para,nmom):
    A=para[1]
    if A <= 0:
        print("Invalid Parameters")
        return
    if nmom > 20:
        print("Parameter nmom too large")
    xmom = []
    xmom.append(para[0]+A)
    if nmom == 1:
        return(xmom)

    xmom.append(0.5*A)
    if nmom ==2:
        return(xmom)

    for i in range(3,nmom+1):
        xmom.append(2/float(i*(i-1)))

    return(xmom)
    
#############################################################

def lmrgam(para,nmom):
    A0 = 0.32573501
    [A1,A2,A3] = [0.16869150, 0.078327243,-0.0029120539]
    [B1,B2] = [0.46697102, 0.24255406]
    C0 = 0.12260172
    [C1,C2,C3] = [0.053730130, 0.043384378, 0.011101277]
    [D1,D2]    = [0.18324466, 0.20166036]
    [E1,E2,E3] = [2.3807576, 1.5931792, 0.11618371]
    [F1,F2,F3] = [5.1533299, 7.1425260, 1.9745056]
    [G1,G2,G3] = [2.1235833, 4.1670213, 3.1925299]
    [H1,H2,H3] = [9.0551443, 26.649995, 26.193668]

    Alpha = para[0]
    Beta = para[1]
    if Alpha <= 0 or Beta <= 0:
        print("Invalid Parameters")
        return
    if nmom > 4:
        print("Parameter nmom too large")
        return
    
    xmom = []
    xmom.append(Alpha*Beta)
    if nmom == 1:
        return(xmom)

    xmom.append(Beta*1/sp.sqrt(sp.pi)*sp.exp(spsp.gammaln(Alpha+0.5)-spsp.gammaln(Alpha)))
    if nmom == 2:
        return(xmom)

    if Alpha < 1:
        Z= Alpha
        xmom.append((((E3*Z+E2)*Z+E1)*Z+1)/(((F3*Z+F2)*Z+F1)*Z+1))
        if nmom == 3:
            return(xmom)
        xmom.append((((C3*Z+C2)*Z+C1)*Z+C0)/((D2*Z+D1)*Z+1))
        if nmom == 4:
            return(xmom)
    else:
        Z=1/Alpha
        xmom.append(sp.sqrt(Z)*(((A3*Z+A2)*Z+A1)*Z+A0)/((B2*Z+B1)*Z+1))
        if nmom == 3:
            return(xmom)
        
        xmom.append((((C3*Z+C2)*Z+C1)*Z+C0)/((D2*Z+D1)*Z+1))
        if nmom == 4:
            return(xmom)

#############################################################

def lmrgev(para,nmom):

    ZMOM=[0.577215664901532861, 0.693147180559945309,
        0.169925001442312363,0.150374992788438185,
        0.558683500577583138e-1,0.581100239999710876e-1,
        0.276242584297309125e-1,0.305563766579053126e-1,
        0.164650282258328802e-1,0.187846624298170912e-1,
        0.109328215063027148e-1,0.126973126676329530e-1,
        0.778982818057231804e-2,0.914836179621999726e-2,
        0.583332389328363588e-2,0.690104287590348154e-2,
        0.453267970180679549e-2,0.538916811326595459e-2,
        0.362407767772368790e-2,0.432387608605538096e-2]
    SMALL = 1e-6
    U = para[0]
    A = para[1]
    G = para[2]
    if A<= 0 or G <= -1:
        print("Invalid Parameters")
        return
    if nmom > 20:
        print("Parameter nmom too large")
        return

    if abs(G)>SMALL:
        GAM = sp.exp(spsp.gammaln(1+G))
        xmom = [U+A*(1-GAM)/G]
        if nmom == 1:
            return(xmom)

        XX2 = 1-2**(-G)
        xmom.append(A*XX2*GAM/G)
        if nmom == 2:
            return(xmom)
 
        Z0=1
        for j in range(2,nmom):
            DJ=j+1
            BETA = (1-DJ**(-G))/XX2
            Z0 = Z0*(4*DJ-6)/DJ
            Z = Z0*3*(DJ-1)/(DJ+1)
            SUM = Z0*BETA-Z
            if j == 2:
                xmom.append(SUM)
            else:
                for i in range(1,j-1):
                    DI = i+1
                    Z = Z*(DI+DI+1)*(DJ-DI)/((DI+DI-1)*(DJ+DI))
                    SUM = SUM-Z*xmom[i+1]
                xmom.append(SUM)
        return(xmom)
    
    else:
        xmom = [U]
        if nmom == 1:
            return(xmom)

        xmom.append(A*ZMOM[1])
        if nmom == 2:
            return(xmom)

        for i in range(2,nmom):
            xmom.append(zmom[i-1])

        return(xmom)
  
#############################################################
    
def lmrglo(para,nmom):
    SMALL = 1e-4
    C1 = sp.pi**2/6
    C2 = 7*sp.pi**4/360


    Z = [[0],[0]]
    Z.append([1])
    Z.append([0.166666666666666667,  0.833333333333333333])
    Z.append([0.416666666666666667,  0.583333333333333333])
    Z.append([0.666666666666666667e-1,  0.583333333333333333,
              0.350000000000000000])
    Z.append([0.233333333333333333,  0.583333333333333333,
              0.183333333333333333])

    Z.append([0.357142857142857143e-1,  0.420833333333333333,
              0.458333333333333333,  0.851190476190476190e-1])

    Z.append([0.150992063492063492,  0.515625000000000000,
              0.297916666666666667,  0.354662698412698413e-1])

    Z.append([0.222222222222222222e-1,  0.318893298059964727,
              0.479976851851851852,  0.165509259259259259,
              0.133983686067019400e-1])

    Z.append([0.106507936507936508,  0.447663139329805996,
              0.360810185185185185,  0.803902116402116402e-1,
              0.462852733686067019e-2])

    Z.append([0.151515151515151515e-1,  0.251316137566137566,
              0.469695216049382716,  0.227650462962962963,
              0.347139550264550265e-1,  0.147271324354657688e-2])

    Z.append([0.795695045695045695e-1,  0.389765946502057613,
              0.392917309670781893,  0.123813106261022928,
              0.134998713991769547e-1,  0.434261597456041900e-3])

    Z.append([0.109890109890109890e-1,  0.204132996632996633,
              0.447736625514403292,  0.273053442827748383,
              0.591917438271604938e-1,  0.477687757201646091e-2,
              0.119302636663747775e-3])

    Z.append([0.619345205059490774e-1,  0.342031759392870504,
              0.407013705173427396,  0.162189192806752331,
              0.252492100235155791e-1,  0.155093427662872107e-2,
              0.306778208563922850e-4])

    Z.append([0.833333333333333333e-2,  0.169768364902293474,
              0.422191282868366202,  0.305427172894620811,
              0.840827939972285210e-1,  0.972435791446208113e-2,
              0.465280282988616322e-3,  0.741380670696146887e-5])

    Z.append([0.497166028416028416e-1,  0.302765838589871328,
              0.410473300089185506,  0.194839026503251764,
              0.386598063704648526e-1,  0.341399407642897226e-2,
              0.129741617371825705e-3,  0.168991182291033482e-5])
             
    Z.append([0.653594771241830065e-2,  0.143874847595085690,
              0.396432853710259464,  0.328084180720899471,
              0.107971393165194318,  0.159653369932077769e-1,
              0.110127737569143819e-2,  0.337982364582066963e-4,
              0.364490785333601627e-6])

    Z.append([0.408784570549276431e-1,  0.270244290725441519,
              0.407599524514551521,  0.222111426489320008,
              0.528463884629533398e-1,  0.598298239272872761e-2,
              0.328593965565898436e-3,  0.826179113422830354e-5,
              0.746033771150646605e-7])

    Z.append([0.526315789473684211e-2,  0.123817655753054913,
              0.371859291444794917,  0.343568747670189607,
              0.130198662812524058,  0.231474364899477023e-1,
              0.205192519479869981e-2,  0.912058258107571930e-4,
              0.190238611643414884e-5,  0.145280260697757497e-7])

    U = para[0]
    A = para[1]
    G = para[2]

    if A <= 0 or abs(G) >= 1:
        print("Invalid Parameters")
        return

    if nmom > 20:
        print("Parameter nmom too large")
        return
    GG = G*G
    ALAM1 = -G*(C1+GG*C2)
    ALAM2 = 1+GG*(C1+GG*C2)
    if abs(G) > SMALL:
        ALAM2=G*sp.pi/sp.sin(G*sp.pi)
        ALAM1=(1-ALAM2)/G

    xmom = [U+A*ALAM1]
    if nmom == 1:
        return(xmom)
             
    xmom.append(A*ALAM2)
    if nmom == 2:
        return(xmom)

    for M in range(3,nmom+1):
        kmax = M/2
        SUMM=Z[M-1][kmax-1]
        for K in range(kmax-1,0,-1):
            SUMM = SUMM*GG+Z[M-1][K-1]
        if M != M/2*2:
            SUMM = -G*SUMM
        xmom.append(SUMM)

    return(xmom)

#############################################################

def lmrgno(para,nmom):

    ZMOM = [0,   0.564189583547756287, 0,   0.122601719540890947,
            0,   0.436611538950024944e-1,0, 0.218431360332508776e-1,
            0,   0.129635015801507746e-1,0, 0.852962124191705402e-2,
            0,   0.601389015179323333e-2,0, 0.445558258647650150e-2,
            0,   0.342643243578076985e-2,0, 0.271267963048139365e-2]


    RRT2 = 1/sp.sqrt(2)
    RRTPI = 1/sp.sqrt(sp.pi)
    
    RANGE = 5
    EPS = 1e-8
    MAXIT = 10

    U = para[0]
    A = para[1]
    G = para[2]
    if A <= 0:
        print("Invalid Parameters")
        return
    if nmom > 20:
        print("Parameter nmom too large")
        return

    

    if abs(G)<=EPS:
        xmom = [U]
        if nmom == 1:
            return(xmom)

        xmom.append(A*ZMOM[1])
        if nmom == 2:
            return(xmom)

        for i in range(3,nmom+1):
            xmom.append(zmom[i-1])

        return(xmom)


    EGG = sp.exp(0.5*G**2)
    ALAM1 = (1-EGG)/G
    xmom = [U+A*ALAM1]
    if nmom == 1:
        return(xmom)
    
    ALAM2=EGG*spsp.erf(0.5*G)/G
    xmom.append(A*ALAM2)
    if nmom == 2:
        return(xmom)
  
    CC=-G*RRT2
    XMIN=CC-RANGE
    XMAX=CC+RANGE
    SUMM = [0]*nmom
    
    N=16
    XINC=(XMAX-XMIN)/N

    for i in range(1,N):
        X = XMIN+i*XINC
        E = sp.exp(-((X-CC)**2))
        D = spsp.erf(X)
        P1 = 1
        P = D
        for M in range(3,nmom+1):
            C1=M+M-3
            C2=M-2
            C3=M-1
            P2=P1
            P1=P
            P=(C1*D*P1-C2*P2)/C3
            SUMM[M-1] = SUMM[M-1]+E*P

    EST = []
    for i in SUMM:
        EST.append(i*XINC)


    for IT in range(1,MAXIT+1):
        ESTX = EST
        N=N*2
        XINC=(XMAX-XMIN)/N
        for i in range(1,N-1,2):
            X = XMIN+i*XINC
            E = sp.exp(-((X-CC)**2))
            D = spsp.erf(X)
            P1 = 1
            P = D
            for M in range(3,nmom+1):
                C1=M+M-3
                C2=M-2
                C3=M-1
                P2=P1
                P1=P
                P=(C1*D*P1-C2*P2)/C3
                SUMM[M-1] = SUMM[M-1]+E*P

        NOTCGD = 0
        for M in range(nmom,2,-1):
            EST[M-1] = SUMM[M-1]*XINC
            if abs(EST[M-1]-ESTX[M-1]) > EPS*abs(EST[M-1]):
                NOTCGD = M

        if NOTCGD == 0:
            CONST = -sp.exp(CC**2)*RRTPI/(ALAM2*G)
            
            for M in range(3,nmom+1):
                xmom.append(CONST*EST[M-1])
            return(xmom)
        else:
            print("Did Not Converge")
            return
        
#############################################################
                
def lmrgpa(para,nmom):
    U = para[0]
    A = para[1]
    G = para[2]
    if A <=0 or G < -1:
        print("Invalid Parameters")
        return
    if nmom > 20:
        print("Parameter nmom too large")
        return

    Y = 1/(1+G)
    xmom = [U+A*Y]
    if nmom == 1:
        return(xmom)
    
    Y = Y/(2+G)
    xmom.append(A*Y)
    if nmom == 2:
        return(xmom)
    
    Y = 1
    for i in range(3,nmom+1):
        AM = i-2
        Y = Y*(AM-G)/(i+G)
        xmom.append(Y)
    return(xmom)

#############################################################

def lmrgum(para,nmom):
    ZMOM = [0.577215664901532861,  0.693147180559945309,
     0.169925001442312363,  0.150374992788438185,
     0.0558683500577583138,  0.0581100239999710876,
     0.0276242584297309125,  0.0305563766579053126,
     0.0164650282258328802,  0.0187846624298170912,
     0.0109328215063027148,  0.0126973126676329530,
     0.00778982818057231804,  0.00914836179621999726,
     0.00583332389328363588,  0.00690104287590348154,
     0.00453267970180679549,  0.00538916811326595459,
     0.00362407767772368790,  0.00432387608605538096]

    A = para[1]
    if A <=0:
        print("Invalid Parameters")
        return
    if nmom >20:
        print("Parameter nmom too large")
        return
    xmom = [para[0]+A*ZMOM[0]]
    if nmom == 1:
        return(xmom)
    xmom.append(A*ZMOM[1])
    if nmom == 2:
        return(xmom)

    for i in range(2,nmom):
        xmom.append(ZMOM[i])
    return(xmom)

#############################################################

def lmrkap(para,nmom):
    EU = 0.577215664901532861
    SMALL = 1e-8
    OFL = 170
    U = para[0]
    A = para[1]
    G = para[2]
    H = para[3]

    if A <= 0 or G <= -1: 
        print("Invalid Parameters")
        return
    if H < 0 and (G*H)<= -1:
        print("Invalid Parameters")
        return
    if nmom > 20:
        print("Parameter nmom too large")
        return

    DLGAM = spsp.gammaln(1+G)
    ICASE = 1
    if H > 0:
        ICASE = 3
    elif abs(H) < SMALL:
        ICASE = 2
    elif G == 0:
        ICASE = ICASE+3

    if ICASE == 1:
        Beta = []
        for IR in range(1,nmom+1):
            ARG = DLGAM + spsp.gammaln(-IR/H-G) - spsp.gammaln(-IR/H)-G*sp.log(-H)
            if abs(ARG) > OFL:
                print("Calculation of L-Moments Failed")
                return
            Beta.append(sp.exp(ARG))


    elif ICASE == 2:
        Beta = []
        for IR in range(1,nmom+1):
            Beta.append(sp.exp(DLGAM-G*sp.log(IR))*(1-0.5*H*G*(1+G)/IR))

    elif ICASE == 3:
        Beta = []
        for IR in range(1,nmom+1):
            ARG = DLGAM+ spsp.gammaln(1+IR/H)-spsp.gammaln(1+G+IR/H)-G*sp.log(H)
            if abs(ARG) > OFL:
                print("Calculation of L-Moments Failed")
                return
            Beta.append(sp.exp(ARG))
            
    elif ICASE == 4:
        Beta = []
        for IR in range(1,nmom+1):
            Beta.append(EU+sp.log(-H)+spsp.psi(-IR/H))
            
    elif ICASE == 5:
        Beta = []
        for IR in range(1,nmom+1):
            Beta.append(EU+sp.log(IR))

    elif ICASE == 6:
        Beta = []
        for IR in range(1,nmom+1):
            Beta.append(EU+sp.log(H)+spsp.psi(1+IR/H))

    if G == 0:
        xmom = [U+A*Beta[0]]
    else:
        xmom = [U+A*(1-Beta[0])/G]

    if nmom == 1:
        return(xmom)

    ALAM2 = Beta[1]-Beta[0]
    if G == 0:
        xmom.append(A*ALAM2)
    else:
        xmom.append(A*ALAM2/(-G))

    if nmom == 2:
        return(xmom)

    Z0 = 1
    for j in range(3,nmom+1):
        Z0 = Z0*(4.0*j-6)/j
        Z = 3*Z0*(j-1)/(j+1)
        SUMM = Z0*(Beta[j-1]-Beta[0])/ALAM2 - Z
        if j == 3:
            xmom.append(SUMM)
        else:
            for i in range(2,j-1):
                Z = Z*(i+i+1)*(j-i)/((i+i-1)*(j+i))
                SUMM = SUMM - Z*xmom[i]
            xmom.append(SUMM)
    return(xmom)

#############################################################

def lmrnor(para,nmom):

    ZMOM =[0, 0.564189583547756287, 0,   0.122601719540890947,
        0,  0.0436611538950024944, 0,   0.0218431360332508776,
        0,  0.0129635015801507746, 0,   0.00852962124191705402,
        0,  0.00601389015179323333, 0,   0.00445558258647650150,
        0,  0.00342643243578076985, 0,   0.00271267963048139365]

    if para[1] <= 0:
        print("Invalid Parameters")
        return
    if nmom > 20:
        print("Parameter nmom too large")
        return
    xmom = [para[0]]
    if nmom == 1:
        return(xmom)

    xmom.append(para[1]*ZMOM[1])
    if nmom == 2:
        return(xmom)

    for M in range(2,nmom):
        xmom.append(ZMOM[M])

    return(xmom)

#############################################################

def lmrpe3(para,nmom):
    SMALL = 1e-6
    CONST = 1/sp.sqrt(sp.pi)
    A0 = 0.32573501
    [A1,A2,A3] = [0.16869150, 0.078327243,-0.0029120539]
    [B1,B2] = [0.46697102,0.24255406]
    C0 = 0.12260172
    [C1,C2,C3] = 0.053730130, 0.043384378, 0.011101277
    [D1,D2] = [0.18324466, 0.20166036]
    [E1,E2,E3] = [2.3807576, 1.5931792, 0.11618371]
    [F1,F2,F3] = [5.1533299, 7.1425260, 1.9745056]
    [G1,G2,G3] = [2.1235833, 4.1670213, 3.1925299]
    [H1,H2,H3] = [9.0551443, 26.649995, 26.193668]

    SD = para[1]
    if SD <= 0:
        print("Invalid Parameters")
        return
    if nmom > 4:
        print("Parameter nmom too large")
        return

    xmom = [para[0]]
    if nmom == 1:
        return(xmom)

    Gamma = para[2]
    if abs(Gamma) < SMALL:
        xmom = [para[0]]
        if nmom == 1:
            return(xmom)

        xmom.append(CONST*Para[1])
        if nmom == 2:
            return(xmom)

        xmom.append(0)
        if nmom == 3:
            return(xmom)

        xmom.append(C0)
        return(xmom)
    else:
        Alpha = 4/(Gamma*Gamma)
        Beta = abs(0.5*SD*Gamma)
        ALAM2 = CONST*sp.exp(spsp.gammaln(Alpha+0.5)-spsp.gammaln(Alpha))
        xmom.append(ALAM2*Beta)
        if nmom == 2:
            return(xmom)

        if Alpha < 1:
            Z = Alpha
            xmom.append((((E3*Z+E2)*Z+E1)*Z+1)/(((F3*Z+F2)*Z+F1)*Z+1))
            if Gamma<0:
                xmom[2] = -xmom[2]
            if nmom == 3:
                return(xmom)

            xmom.append((((G3*Z+G2)*Z+G1)*Z+1)/(((H3*Z+H2)*Z+H1)*Z+1))
            return(xmom)

        else:
            Z = 1.0/Alpha
            xmom.append(sp.sqrt(Z)*(((A3*Z+A2)*Z+A1)*Z+A0)/((B2*Z+B1)*Z+1))
            if Gamma < 0:
                xmom[2] = -xmom[2]

            if nmom == 3:
                return(xmom)

            xmom.append((((C3*Z+C2)*Z+C1)*Z+C0)/((D2*Z+D1)*Z+1))
            return(xmom)


#############################################################

def lmrwak(para,nmom):
    [XI,A,B,C,D]=para
    fail = 0
    if D >= 1:
        fail = 1
    if (B+D)<= 0 and (B!= 0 or C != 0 or D!=0):
        fail = 1
    if A == 0 and B != 0:
        fail = 1
    if C == 0 and D != 0:
        fail = 1
    if C < 0:
        fail = 1
    if (A+C) < 0:
        fail = 1
    if A == 0 and C == 0:
        fail = 1
    if nmom >= 20:
        fail = 2

    if fail == 1:
        print("Invalid Parameters")
        return
    if fail == 2:
        print("Parameter nmom too large")
        return

    Y=A/(1+B)
    Z=C/(1-D)
    xmom = []
    xmom.append(XI+Y+Z)
    if nmom == 1:
        return
    
    Y=Y/(2+B)
    Z=Z/(2-D)
    ALAM2=Y+Z
    xmom.append(ALAM2)
    if nmom == 2:
        return

    for i in range(2,nmom):
        AM=i+1
        Y=Y*(AM-2-B)/(AM+B)
        Z=Z*(AM-2+D)/(AM-D)
        xmom.append((Y+Z)/ALAM2)

    return(xmom)

#############################################################
def lmrwei(para,nmom):
    if len(para) != 3:
        print("Invalid number of parameters")
        return
    if para[1] <= 0 or para[2] <= 0:
        print("Invalid Parameters")
        return
    
    xmom = lmrgev.lmrgev([0,para[1]/para[2],1/para[2]],nmom)
    xmom[0] = para[0]+para[1] - xmom[0]
    xmom[2] = -xmom[2]
    return(xmom)


#############################################################
###PEL FUNCTIONS
#############################################################


def pelexp(xmom):
    if xmom[1] <= 0:
        print("L-Moments Invalid")
        return
    else:
        para = [xmom[0]-2*xmom[1],2*xmom[1]]
        return(para)

#############################################################

def pelgam(xmom):
    A1 = -0.3080
    A2 = -0.05812
    A3 = 0.01765
    B1 = 0.7213
    B2 = -0.5947
    B3 = -2.1817
    B4 = 1.2113
    
    if xmom[0] <= xmom[1] or xmom[1]<= 0:
        print("L-Moments Invalid")
        return
    CV = xmom[1]/xmom[0]
    if CV >= 0.5:
        T = 1-CV
        ALPHA =T*(B1+T*B2)/(1+T*(B3+T*B4))
    else:
        T=sp.pi*CV**2
        ALPHA=(1+A1*T)/(T*(1+T*(A2+T*A3)))
        
    para = [ALPHA,xmom[0]/ALPHA]
    return(para)

#############################################################

def pelgev(xmom):
    SMALL = 1e-5
    eps = 1e-6
    maxit = 20
    EU =0.57721566
    DL2 = sp.log(2)
    DL3 = sp.log(3)
    A0 =  0.28377530
    A1 = -1.21096399
    A2 = -2.50728214
    A3 = -1.13455566
    A4 = -0.07138022
    B1 =  2.06189696 
    B2 =  1.31912239 
    B3 =  0.25077104
    C1 =  1.59921491
    C2 = -0.48832213
    C3 =  0.01573152
    D1 = -0.64363929
    D2 =  0.08985247

    T3 = xmom[2]
    if xmom[1]<= 0 or abs(T3)>= 1:
        print("L-Moments Invalid")
        return
    if T3<= 0:
        G=(A0+T3*(A1+T3*(A2+T3*(A3+T3*A4))))/(1+T3*(B1+T3*(B2+T3*B3)))
        if T3>= -0.8:
            para3 = G
            GAM = sp.exp(sp.special.gammaln(1+G))
            para2=xmom[1]*G/(GAM*(1-2**(-G)))
            para1=xmom[0]-para2*(1-GAM)/G
            para = [para1,para2,para3]
            return(para)

        if T3 <= -0.97:
            G = 1-sp.log(1+T3)/DL2
            
        T0=(T3+3)*0.5
        for IT in range(1,maxit):
            X2=2**(-G)
            X3=3**(-G)
            XX2=1-X2
            XX3=1-X3
            T=XX3/XX2
            DERIV=(XX2*X3*DL3-XX3*X2*DL2)/(XX2**2)
            GOLD=G
            G=G-(T-T0)/DERIV
            if abs(G-GOLD) <= eps*G:
                para3 = G
                GAM = sp.exp(sp.special.gammaln(1+G))
                para2=xmom[1]*G/(GAM*(1-2**(-G)))
                para1=xmom[0]-para2*(1-GAM)/G
                para = [para1,para2,para3]
                return(para)
            
        print("Iteration has not converged")

    Z=1-T3
    G=(-1+Z*(C1+Z*(C2+Z*C3)))/(1+Z*(D1+Z*D2))
    if abs(G)<SMALL:
        para2 = xmom[1]/DL2
        para1 = xmom[0]-EU*para2
        para = [para1,para2,0]
        return(para)
    else:
        para3 = G
        GAM = sp.exp(sp.special.gammaln(1+G))
        para2=xmom[1]*G/(GAM*(1-2**(-G)))
        para1=xmom[0]-para2*(1-GAM)/G
        para = [para1,para2,para3]
        return(para)
 
#############################################################

def pelglo(xmom):
    SMALL = 1e-6
    
    G=-xmom[2]
    if xmom[1]<= 0 or abs(G)>= 1:
        print("L-Moments Invalid")
        return

    if abs(G)<= SMALL:
        para = [xmom[0],xmom[1],0]
        return(para)

    GG = G*sp.pi/sp.sin(G*sp.pi)
    A = xmom[1]/GG
    para1 = xmom[0]-A*(1-GG)/G
    para = [para1,A,G]
    return(para)

#############################################################

def pelgno(xmom):
    A0 =  0.20466534e+01
    A1 = -0.36544371e+01
    A2 =  0.18396733e+01
    A3 = -0.20360244e+00
    B1 = -0.20182173e+01
    B2 =  0.12420401e+01
    B3 = -0.21741801e+00
    SMALL = 1e-8

    T3=xmom[2]
    if xmom[1] <= 0 or abs(T3) >= 1:
        print("L-Moments Invalid")
        return
    if abs(T3)>= 0.95:
        para = [0,-1,0]
        return(para)

    if abs(T3)<= SMALL:
        para =[xmom[0],xmom[1]*sp.sqrt(sp.pi),0] 

    TT=T3**2
    G=-T3*(A0+TT*(A1+TT*(A2+TT*A3)))/(1+TT*(B1+TT*(B2+TT*B3)))
    E=sp.exp(0.5*G**2)
    A=xmom[1]*G/(E*sp.special.erf(0.5*G))
    U=xmom[0]+A*(E-1)/G
    para = [U,A,G]
    return(para)

#############################################################

def pelgpa(xmom):
    T3=xmom[2]
    if xmom[1]<= 0:
        print("L-Moments Invalid")
        return
    if abs(T3)>= 1:
        print("L-Moments Invalid")
        return

    G=(1-3*T3)/(1+T3)
    
    PARA3=G
    PARA2=(1+G)*(2+G)*xmom[1]
    PARA1=xmom[0]-PARA2/(1+G)
    para = [PARA1,PARA2,PARA3]
    return(para)

#############################################################

def pelgum(xmom):
    EU = 0.577215664901532861
    if xmom[1] <= 0:
        print("L-Moments Invalid")
        return
    else:
        para2 = xmom[1]/sp.log(2)
        para1 = xmom[0]-EU*para2
        para = [para1, para2]
        return(para)

#############################################################

def pelkap(xmom):
    EPS = 1e-6
    MAXIT = 20
    MAXSR = 10
    HSTART = 1.001
    BIG = 10
    OFLEXP = 170
    OFLGAM = 53

    T3 = xmom[2]
    T4 = xmom[3]
    para = [0]*4
    if xmom[1] <= 0:
        print("L-Moments Invalid")
        return
    if abs(T3) >= 1 or abs(T4) >=  1:
        print("L-Moments Invalid")
        return

    if T4 <= (5*T3*T3-1)/4:
        print("L-Moments Invalid")
        return

    if T4 >= (5*T3*T3+1)/6:
        print("L-Moments Invalid")
        return

    G = (1-3*T3)/(1+T3)
    H = HSTART
    Z = G+H*0.725
    Xdist = BIG

    #Newton-Raphson Iteration
    for it in range(1,MAXIT+1):
        for i in range(1,MAXSR+1):
            if G > OFLGAM:
                print("Failed to converge")
                return
            if H > 0:
                U1 = sp.exp(spsp.gammaln(1/H)-spsp.gammaln(1/H+1+G))
                U2 = sp.exp(spsp.gammaln(2/H)-spsp.gammaln(2/H+1+G))
                U3 = sp.exp(spsp.gammaln(3/H)-spsp.gammaln(3/H+1+G))
                U4 = sp.exp(spsp.gammaln(4/H)-spsp.gammaln(4/H+1+G))
            else:
                U1 = sp.exp(spsp.gammaln(-1/H-G)-spsp.gammaln(-1/H+1))
                U2 = sp.exp(spsp.gammaln(-2/H-G)-spsp.gammaln(-2/H+1))
                U3 = sp.exp(spsp.gammaln(-3/H-G)-spsp.gammaln(-3/H+1))
                U4 = sp.exp(spsp.gammaln(-4/H-G)-spsp.gammaln(-4/H+1))

            ALAM2 =  U1-2*U2
            ALAM3 = -U1+6*U2-6*U3
            ALAM4 =  U1-12*U2+30*U3-20*U4
            if ALAM2 == 0:
                print("Failed to Converge")
                return
            TAU3 = ALAM3/ALAM2
            TAU4 = ALAM4/ALAM2
            E1 = TAU3-T3
            E2 = TAU4-T4

            DIST = max(abs(E1),abs(E2))
            if DIST < Xdist:
                Success = 1
                break
            else:
                DEL1 = 0.5*DEL1
                DEL2 = 0.5*DEL2
                G = XG-DEL1
                H = XH-DEL2
                
        if Success == 0:
            print("Failed to converge")
            return

        #Test for convergence
        if DIST < EPS:
            para[3]=H
            para[2]=G
            TEMP = spsp.gammaln(1+G)
            if TEMP > OFLEXP:
                print("Failed to converge")
                return
            GAM = sp.exp(TEMP)
            TEMP = (1+G)*sp.log(abs(H))
            if  TEMP > OFLEXP:
                print("Failed to converge")
                return

            HH = sp.exp(TEMP)
            para[1] = xmom[1]*G*HH/(ALAM2*GAM)
            para[0] = xmom[0]-para[1]/G*(1-GAM*U1/HH)
            return(para)
        else:
            XG=G
            XH=H
            XZ=Z
            Xdist=DIST
            RHH=1/(H**2)
            if H > 0:
                U1G=-U1*spsp.psi(1/H+1+G)
                U2G=-U2*spsp.psi(2/H+1+G)
                U3G=-U3*spsp.psi(3/H+1+G)
                U4G=-U4*spsp.psi(4/H+1+G)
                U1H=  RHH*(-U1G-U1*spsp.psi(1/H))
                U2H=2*RHH*(-U2G-U2*spsp.psi(2/H))
                U3H=3*RHH*(-U3G-U3*spsp.psi(3/H))
                U4H=4*RHH*(-U4G-U4*spsp.psi(4/H))
            else:
                U1G=-U1*spsp.psi(-1/H-G)
                U2G=-U2*spsp.psi(-2/H-G)
                U3G=-U3*spsp.psi(-3/H-G)
                U4G=-U4*spsp.psi(-4/H-G)
                U1H=  RHH*(-U1G-U1*spsp.psi(-1/H+1))
                U2H=2*RHH*(-U2G-U2*spsp.psi(-2/H+1))
                U3H=3*RHH*(-U3G-U3*spsp.psi(-3/H+1))
                U4H=4*RHH*(-U4G-U4*spsp.psi(-4/H+1))

            DL2G=U1G-2*U2G
            DL2H=U1H-2*U2H
            DL3G=-U1G+6*U2G-6*U3G
            DL3H=-U1H+6*U2H-6*U3H
            DL4G=U1G-12*U2G+30*U3G-20*U4G
            DL4H=U1H-12*U2H+30*U3H-20*U4H
            D11=(DL3G-TAU3*DL2G)/ALAM2
            D12=(DL3H-TAU3*DL2H)/ALAM2
            D21=(DL4G-TAU4*DL2G)/ALAM2
            D22=(DL4H-TAU4*DL2H)/ALAM2
            DET=D11*D22-D12*D21
            H11= D22/DET
            H12=-D12/DET
            H21=-D21/DET
            H22= D11/DET
            DEL1=E1*H11+E2*H12
            DEL2=E1*H21+E2*H22
            
##          TAKE NEXT N-R STEP
            G=XG-DEL1
            H=XH-DEL2
            Z=G+H*0.725

##          REDUCE STEP IF G AND H ARE OUTSIDE THE PARAMETER SPACE

            FACTOR=1
            if G <= -1:
                FACTOR = 0.8*(XG+1)/DEL1
            if H <= -1:
                FACTOR = min(FACTOR,0.8*(XH+1)/DEL2)
            if Z <= -1:
                FACTOR = min(FACTOR,0.8*(XZ+1)/(XZ-Z))
            if H <= 0 and G*H<= -1:
                FACTOR = min(FACTOR,0.8*(XG*XH+1)/(XG*XH-G*H))

            if FACTOR == 1:
                pass
            else:
                DEL1 = DEL1*FACTOR
                DEL2 = DEL2*FACTOR
                G = XG-DEL1
                H = XH-DEL2
                Z = G+H*0.725
    
#############################################################

def pelnor(xmom):
    if xmom[1] <= 0:
        print("L-Moments Invalid")
        return
    else:
        para = [xmom[0],xmom[1]*sp.sqrt(sp.pi)]
        return(para)

#############################################################
    
def pelpe3(xmom):
    Small = 1e-6
    #Constants used in Minimax Approx:

    C1 = 0.2906
    C2 = 0.1882
    C3 = 0.0442
    D1 = 0.36067
    D2 = -0.59567
    D3 = 0.25361
    D4 = -2.78861
    D5 = 2.56096
    D6 = -0.77045

    T3=abs(xmom[2])
    if xmom[1] <= 0 or T3 >= 1:
        para = [0]*3
        print("L-Moments Invalid")
        return(para)

    if T3<= Small:
        para = []
        para.append(xmom[0])
        para.append(xmom[1]*sp.sqrt(sp.pi))
        para.append(0)
        return(para)

    if T3 >= (1.0/3):
        T = 1-T3
        Alpha = T*(D1+T*(D2+T*D3))/(1+T*(D4+T*(D5+TD6)))
    else:
        T=3*sp.pi*T3*T3
        Alpha=(1+C1*T)/(T*(1+T*(C2+T*C3)))
                    
    RTALPH=sp.sqrt(Alpha)
    BETA=sp.sqrt(sp.pi)*xmom[1]*sp.exp(sp.special.gammaln(Alpha)-sp.special.gammaln(Alpha+0.5))
    para = []
    para.append(xmom[0])
    para.append(BETA*RTALPH)
    para.append(2/RTALPH)
    if xmom[2] < 0:
        para[2]=-para[2]

    return(para)

#############################################################

def pelwak(xmom):
    iFail = 0
    FitPareto = 0
    tryxiiszero = 0
    if abs(xmom[1]) <= 0:
        iFail=3
    if abs(xmom[2]) > 1:        
        iFail=3
    if abs(xmom[3]) > 1:
        iFail=3
    if abs(xmom[4]) > 1:
        iFail=3

    if iFail ==3:
        print("L-Moments Invalid")
        para = [0]*5
        return para
         
    iFail = 0

#CALCULATE THE L-MOMENTS (LAMBDA'S)
    alam1 = xmom[0]
    alam2 = xmom[1]
    alam3 = xmom[2]*alam2
    alam4 = xmom[3]*alam2
    alam5 = xmom[4]*alam2
    
#ESTIMATE N1,N2,N3,C1,C2,C3 WHEN XI.NE.0

    XN1= 3*alam2-25*alam3 +32*alam4
    XN2=-3*alam2 +5*alam3  +8*alam4
    XN3= 3*alam2 +5*alam3  +2*alam4
    XC1= 7*alam2-85*alam3+203*alam4-125*alam5
    XC2=-7*alam2+25*alam3  +7*alam4 -25*alam5
    XC3= 7*alam2 +5*alam3  -7*alam4  -5*alam5

#Estimate B and D
    
    XA=XN2*XC3-XC2*XN3
    XB=XN1*XC3-XC1*XN3
    XC=XN1*XC2-XC1*XN2
    Disc=XB*XB-4*XA*XC
    
    tryxiiszero = 0
    if Disc < 0:
        tryxiiszero = 1
    else:
        Disc=sp.sqrt(Disc)
        ROOT1=0.5*(-XB+Disc)/XA
        ROOT2=0.5*(-XB-Disc)/XA
        B= max(ROOT1,ROOT2)
        D=-min(ROOT1,ROOT2)
        if D >= 1:
            tryxiiszero = 1
        else:
            A=(1+B)*(2+B)*(3+B)/(4*(B+D))*((1+D)*alam2-(3-D)*alam3)
            C=-(1-D)*(2-D)*(3-D)/(4*(B+D))*((1-B)*alam2-(3+B)*alam3)
            XI=alam1-A/(1+B)-C/(1-D)
            success = 0
            if C >= 0 and (A+C)>= 0:
               success = 1

##         CAN'T FIND VALID ESTIMATES FOR XI UNRESTRICTED, SO TRY XI=0
##         ESTIMATE B AND D FOR XI=0
    if tryxiiszero == 1:
        iFail=1
        XI=0
        ZN1=4*alam1-11*alam2+9*alam3
        ZN2=-alam2+3*alam3
        ZN3=alam2+alam3
        ZC1=10*alam1-29*alam2+35*alam3-16*alam4
        ZC2=-alam2+5*alam3-4*alam4
        ZC3=alam2-alam4
        ZA=ZN2*ZC3-ZC2*ZN3
        ZB=ZN1*ZC3-ZC1*ZN3
        ZC=ZN1*ZC2-ZC1*ZN2
        Disc=ZB*ZB-4*ZA*ZC
        FitPareto = 0
        if Disc < 0:
            FitPareto = 1
        else:
            Disc=sp.sqrt(Disc)
            ROOT1=0.5*(-ZB+Disc)/ZA
            ROOT2=0.5*(-ZB-Disc)/ZA
            B= max(ROOT1,ROOT2)
            D=-min(ROOT1,ROOT2)
            if D >= 1:
                FitPareto = 1
            else:
                ##         ESTIMATE A AND C
                A= (1+B)*(2+B)/(B+D)*(alam1-(2-D)*alam2)
                C=-(1-D)*(2-D)/(B+D)*(alam1-(2+B)*alam2)
                if C >= 0 and (A+C) >= 0:
                    success = 1
               

    if FitPareto == 1:
        iFail=2
        D=-(1-3*xmom[2])/(1+xmom[2])
        C=(1-D)*(2-D)*xmom[1]
        B=0
        A=0
        XI=xmom[0]-C/(1-D)
        if D > 0:
            success = 1
        else:
            A=C
            B=-D
            C=0
            D=0


    if success == 1:
        para = []
        para.append(XI)
        para.append(A)
        para.append(B)
        para.append(C)
        para.append(D)
        return(para)

#############################################################
def pelwei(lmom):
    if len(lmom) < 3:
        print("Insufficient L-Moments: Need 3")
        return
    if lmom[1] <= 0 or lmom[2] >= 1:
        print("L-Moments Invalid")
        return
    pg = pelgev.pelgev([-lmom[0],lmom[1],-lmom[2]])
    delta = 1/pg[2]
    beta = pg[1]/pg[2]
    out = [-pg[0]-beta,beta,delta]
    return(out)

#############################################################
##QUARTILE FUNCTIONS
#############################################################

def quaexp(F,para):
    U = para[0]
    A = para[1]
    if A <= 0:
        print("Parameters Invalid")
        return
    if F <= 0 or F >= 1:
        print("F Value Invalid")
        return

    QUAEXP = U-A*sp.log(1-F)
    return(QUAEXP)

#############################################################

def quagam(F,para):
    EPS = 1e-10
    maxit = 30
    QUAGAM = 0
    Alpha = para[0]
    Beta = para[1]
    if Alpha <= 0 or Beta <= 0:
        print("Parameters Invalid")
        return
    if F<=0 or F>= 1:
        print("F Value Invalid")
        return
    
    AM1 = Alpha - 1
    if AM1 != 0:
        DLOGG = spsp.gammaln(Alpha)
        if AM1 <= 0:
            Root = sp.exp((sp.log(Alpha*F)+DLOGG)/ALPHA)
        else:
            Root = Alpha*(1-1/(9*Alpha) + quastn.quastn(F)/sp.sqrt(9*Alpha))**3

        if Root <= 0.01*Alpha:
            Root = sp.exp((sp.log(Alpha*F)+DLOGG)/ALPHA)

        for it in range(1,maxit+1):
            FUNC = spsp.gammainc(Alpha,Root)-F
            RINC = FUNC*sp.exp(DLOGG+Root-AM1*sp.log(Root))
            Root = Root-RINC
            if abs(FUNC) <= EPS:
                QUAGAM = Root*Beta
                return(QUAGAM)
    else:
        QUAGAM = -sp.log(1-F)*Beta
        return(QUAGAM)

    print("Result failed to converge")
    return

#############################################################

def quagev(F,para):
    U = para[0]
    A = para[1]
    G = para[2]
    if A <= 0:
        print("Parameters Invalid")
        return
    if F <= 0 or F >= 1:
        if F == 0 and G < 0:
            QUAGEV = U+A/G
        elif F == 1 and G > 0:
            QUAGEV = U+A/G
        else:
            print("F Value Invalid")
            return

        
        print("F Value Invalid")
        return
    else:
        Y = -sp.log(-sp.log(F))
        if G != 0:
            Y = (1-sp.exp(-G*Y))/G
        QUAGEV = U+A*Y
        return(QUAGEV)

#############################################################

def quaglo(F,para):
    U = para[0]
    A = para[1]
    G = para[2]
    if A <= 0:
        print("Invalid Parameters")
        return
    if F <= 0 or F >= 1:
        if F == 0 and G < 0:
            QUAGLO = U+A/G
            return(QUAGLO)
        elif F == 1 and G > 0:
            QUAGLO = U+A/G
            return(QUAGLO)
        else:
            print("F Value Invalid")
            return

    Y = sp.log(F/(1-F))
    if G != 0:
        Y = (1-sp.exp(-G*Y))/G
    QUAGLO = U+A*Y
    return(QUAGLO)

#############################################################

def quagno(F,para):
    U = para[0]
    A = para[1]
    G = para[2]
    if A <= 0:
        print("Invalid Parameters")
        return
    if F <= 0 or F >= 1:
        if F == 0 and G < 0:
            QUAGNO = U+A/G
            return(QUAGNO)
        elif F == 1 and G > 0:
            QUAGNO = U+A/G
            return(QUAGNO)
        else:
            print("F Value Invalid")
            return

    Y = quastn.quastn(F)
    if G != 0:
        Y = (1-sp.exp(-G*Y))/G
    QUAGNO = U+A*Y
    return(QUAGNO)    
 
#############################################################

def quagpa(F,para):
    U = para[0]
    A = para[1]
    G = para[2]
    if A <= 0:
        print("Invalid parameters")
        return
    if F <= 0 or F >= 1:
        if F == 0:
            QUAGPA = U
            return(QUAGPA)
        elif F == 1 and G > 0:
            QUAGPA = U + A/G
            return(QUAGPA)
        else:
            print("F Value Invalid")
            return

    Y = -sp.log(1-F)
    if G !=0:
        Y = (1-sp.exp(-G*Y))/G
    QUAGPA = U+A*Y
    return(QUAGPA)

#############################################################

def quagum(F,para):

    U = para[0]
    A = para[1]
    
    if A <= 0:
        print("Parameters Invalid")
        return
    if F <= 0 or F >= 1:
        print("F Value Invalid")
        return
    QUAGUM = U-A*sp.log(-sp.log(F))
    return(QUAGUM)

#############################################################

def quakap(F,para):
    U = para[0]
    A = para[1]
    G = para[2]
    H = para[3]
    if A <= 0:
        print("Invalid Parameters")
        return
    if F <= 0 or F>= 1:
        if F==0:
            if H<=0 and G < 0:
                QUAKAP = U+A/G
            if H<= 0 and G>= 0:
                print("F Value Invalid")
                return
            if H > 0 and G!= 0:
                QUAKAP = U+A/G*(1-H**(-G))
            if H > 0 and G == 0:
                QUAKAP = U+A*sp.log(H)

            return(QUAKAP)
        
        if F == 1:
            if G <= 0:
                print("F Value Invalid")
                return
            else:
                QUAKAP = U+A/G
                return(QUAKAP)

    else:
        Y = -sp.log(F)
        if H!=0:
            Y = (1-sp.exp(-H*Y))/H
            
        Y = -sp.log(Y)
        if G!= 0:
            Y = (1-sp.exp(-G*Y))/G
        QUAKAP = U+A*Y
        return(QUAKAP)

#############################################################

def quanor(F,para):

    if para[1] <= 0:
        print("Parameters Invalid")
        return
    if F <= 0 or F >= 1:
        print("F Value Invalid")
        return
    QUANOR = para[0]+para[1]*quastn.quastn(F)
    return(QUANOR)

#############################################################

def quape3(F,para):
    SMALL = 1e-6

    if para[1]<= 0:
        print("Paremters Invalid")
        return
    Gamma = para[2]
    if F <= 0 or F >= 1:
        if F == 0 and Gamma >0:
            QUAPE3 = para[0]-2*para[1]/Gamma
            return(QUAPE3)
        elif F == 1 and Gamma < 0:
            QUAPE3 = para[0]-2*para[1]/Gamma
            return(QUAPE3)
        else:
            print("F Value Invalid")
            return


    if abs(Gamma) < SMALL:
        QUAPE3 = para[0] + para[1]*quastn(F)
        return(QUAPE3)

    Alpha = 4/(Gamma*Gamma)
    Beta = abs(0.5*para[1]*Gamma)
    par = [Alpha,Beta]
    if Gamma > 0:
        QUAPE3 = para[0]-Alpha*Beta+quagam.quagam(F,par)
    if Gamma < 0:
        QUAPE3 = para[0]+Alpha*Beta-quagam.quagam(1-F,par)
    return(QUAPE3)

#############################################################

def quastn(F):
    split1 = 0.425
    split2 = 5
    const1 = 0.180625
    const2 = 1.6
    [A0,A1,A2,A3,A4,A5,A6,A7,B1,B2,B3,B4,B5,B6,B7] = [0.338713287279636661e1,
     0.133141667891784377e3,  0.197159095030655144e4,
     0.137316937655094611e5,  0.459219539315498715e5,
     0.672657709270087009e5,  0.334305755835881281e5,
     0.250908092873012267e4,  0.423133307016009113e2,
     0.687187007492057908e3,  0.539419602142475111e4,
     0.212137943015865959e5,  0.393078958000927106e5,
     0.287290857357219427e5,  0.522649527885285456e4]

    [C0,C1,C2,C3,C4,C5,C6,C7,D1,D2,D3,D4,D5,D6,D7] = [0.142343711074968358e1,
     0.463033784615654530e1,  0.576949722146069141e1,
     0.364784832476320461e1,  0.127045825245236838e1,
     0.241780725177450612e0,  0.227238449892691846e-1,
     0.774545014278341408e-3,  0.205319162663775882e1,
     0.167638483018380385e1,  0.689767334985100005e0,
     0.148103976427480075e0,  0.151986665636164572e-1,
     0.547593808499534495e-3,  0.105075007164441684e-8]
                                                      
    [E0,E1,E2,E3,E4,E5,E6,E7,F1,F2,F3,F4,F5,F6,F7] = [0.665790464350110378e1,
     0.546378491116411437e1,  0.178482653991729133e1,
     0.296560571828504891e0,  0.265321895265761230e-1,
     0.124266094738807844e-2,  0.271155556874348758e-4,
     0.201033439929228813e-6,  0.599832206555887938e0,
     0.136929880922735805e0,  0.148753612908506149e-1,
     0.786869131145613259e-3,  0.184631831751005468e-4,
     0.142151175831644589e-6,  0.204426310338993979e-14]

    Q = F-0.5
    if abs(Q) > split1:
        R=F
        if Q >= 0:
            R = 1-F
        if R <= 0:
            print("F Value Invalid")
        R = sp.sqrt(-sp.log(R))
        if R > split2:
            R = R - split2
            QUASTN=((((((((E7*R+E6)*R+E5)*R+E4)*R+E3)*R+E2)*R+E1)*R+E0)/
            (((((((F7*R+F6)*R+F5)*R+F4)*R+F3)*R+F2)*R+F1)*R+1))
            if Q < 0:
                QUASTN = -QUASTN
            return(QUASTN)
        else:
            R=R-const2
            QUASTN=((((((((C7*R+C6)*R+C5)*R+C4)*R+C3)*R+C2)*R+C1)*R+C0)/
            (((((((D7*R+D6)*R+D5)*R+D4)*R+D3)*R+D2)*R+D1)*R+1))
            if Q < 0:
                QUASTN = -QUASTN
            return(QUASTN)

    else:
        R = const1-Q*Q
        QUASTN = Q*((((((((A7*R+A6)*R+A5)*R+A4)*R+A3)*R+A2)*R+A1)*R+A0)/
        (((((((B7*R+B6)*R+B5)*R+B4)*R+B3)*R+B2)*R+B1)*R+1))
        return(QUASTN)
    
#############################################################

def quawak(F,para):
    ufl = -170
  
    XI = para[0]
    A = para[1]
    B = para[2]
    C = para[3]
    D = para[4]
    
    fail = 0
    if (B+D) <= 0 and (B != 0 or C != 0 or D!= 0):
        fail = 1
    if A == 0 and B != 0:
        fail = 1
    if C == 0 and D != 0:
        fail = 1
    if C <0 or (A+C)< 0:
        fail = 1
    if A == 0 and C == 0:
        fail = 1

    if fail == 1:
        print("Parameters Invalid")
        return
    
    if F<= 0 or F>= 1:
        if F == 0:
            QUAWAK = XI
        elif F == 1:
            if D > 0:
                fail = 1
            if D < 0:
                QUAWAK = XI + A/B - C/D
            if D == 0 and C > 0:
                fail = 1
            if D == 0 and C == 0 and B == 0:
                fail = 1
            if D == 0 and C == 0 and B >0:
                QUAWAK = XI+A/B

            if fail == 1:
                print("Function Failed")
            else:
                return(QUAWAK)


    Z=-sp.log(1-F)
    Y1 = Z
    if B == 0:
        Y2 = Z
        if D !=0:
            Y2 = (1-sp.exp(D*Y2))/(-D)
        QUAWAK = XI+A*Y1+C*Y2
        return(QUAWAK)
    else:
        TEMP = -B*Z
        if TEMP < ufl:
            Y1 = 1/B
        if TEMP >= ufl:
            Y1 = (1-sp.exp(TEMP))/B

        Y2 = Z
        if D !=0:
            Y2 = (1-sp.exp(D*Y2))/(-D)
        QUAWAK = XI+A*Y1+C*Y2
        return(QUAWAK)

#############################################################

def quawei(f,para):
    if len(para) != 3:
        print("Invalid number of parameters")
    if para[1] <= 0 or para[2] <= 0:
        print("Invalid Parameters")
    if f <0 or f > 1:
        print("F Value Invalid")

    return(para[0]+para[1]*((-sp.log(1-f))**(1/para[2])))


        
      
     
