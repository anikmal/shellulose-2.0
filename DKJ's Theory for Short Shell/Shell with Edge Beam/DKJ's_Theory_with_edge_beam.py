import sympy as sp
import scipy as syp
import numpy as np
from math import * 
import matplotlib.pyplot as plt
l=83.25#float(input("enter span in meter"))
r=25#float(input("enter radius in meter"))
d=0.25#float(input("enter thickness in meter"))
o=35#float(input("semi circular angle in degree"))
dl=36#float(input("enter dead load intensity in kN/m^2"))
ll=14#float(input("enter live load intensity kN/m^2"))

#Dimension of edge beam
a=5#float(input("enter depth of edge beam in meter"))     
b=0.75#float(input("enter width of edge beam in meter"))
Wb=a*b*150##############################
I=b*a**3/12
W1=Wb*4/pi
a1=a/2
F=E*I
Ar=a*b
g=dl+ll
E=360*10**6#float(input("Modulus of Elasticity kN/m^2"))
Dn=E*d**3/12 #Dn=Flexural Rigidity of shell
print("Flexural Rigidity of shell" ,Dn)
o1=0
ln=pi*r/l
p=(12*(pi**4)*(r**6)/(l**4)/(d**2))**(0.125)
k=((pi**2)*(r**2)/(l**2)/(p**2))
print("Aas-Jakobsen Parameters" ,p  ,k) 
#roots calculation
alpha1=p*(((1+k*2**0.5)**2+1)**0.5+1+k*2**0.5)**0.5/8**0.25
beta1=p*(((1+k*2**0.5)**2+1)**0.5-(1+k*2**0.5))**0.5/8**0.25
alpha2=p*(((1-k*2**0.5)**2+1)**0.5-(1-k*2**0.5))**0.5/8**0.25
beta2=p*(((1-k*2**0.5)**2+1)**0.5+1-k*2**0.5)**0.5/8**0.25
print("Roots of the characterestic equation" , alpha1,beta1,alpha2,beta2)

def coeff( t, n ):       #'coeff' function used make successive differentiation easier of unknown co-efficients of solution
    A=sp.Symbol('A')
    B=sp.Symbol('B')
    C=sp.Symbol('C')
    D=sp.Symbol('D')
    X=0
    Y=0
    X1=0
    Y1=0
    for i in range(n):
        X=sp.simplify(-alpha1*A+beta1*B)
        Y=sp.simplify(-beta1*A-alpha1*B)
        X1=sp.simplify(-alpha2*C+beta2*D)
        Y1=sp.simplify(-beta2*C-alpha2*D)
        A=X
        B=Y
        C=X1
        D=Y1
    if t=='A':
        return(A)
    elif t=='B':
        return(B)
    elif t=='C':
        return(C)
    else:
        return(D)
    
#Equation for stress-resultants based on charecterestic roots at boundary condition

o1=0
x=sp.Symbol('x')
x=0
o21=2*radians(o)
#Equation of Mθ
Mql=sp.simplify(-Dn*((exp(-alpha1*o1)*(coeff('A',2)*cos(beta1*o1)+coeff('B',2)*sin(beta1*o1)))+exp(-alpha2*o1)*(coeff('C',2)*cos(beta2*o1)+coeff('D',2)*sin(beta2*o1)))*sp.cos(ln*x/r)/r**2)
Mqr=sp.simplify(-Dn*((exp(-alpha1*o21)*(coeff('A',2)*cos(beta1*o21)+coeff('B',2)*sin(beta1*o21)))+exp(-alpha2*o21)*(coeff('C',2)*cos(beta2*o21)+coeff('D',2)*sin(beta2*o21)))*sp.cos(ln*x/r)/r**2)
Mq=Mql+Mqr #Considering Edge disturbance
print("Mq=" ,Mq) 

#Equation of Mxθ
dMxql=sp.simplify(Dn*ln**2*(exp(-alpha1*o1)*(coeff('A',1)*cos(beta1*o1)+coeff('B',1)*sin(beta1*o1))+exp(-alpha2*o1)*(coeff('C',1)*cos(beta2*o1)+coeff('D',1)*sin(beta2*o1)))*sp.cos(ln*x/r)/r**3)
dMxqr=sp.simplify(Dn*ln**2*(exp(-alpha1*o21)*(coeff('A',1)*cos(beta1*o21)+coeff('B',1)*sin(beta1*o21))+exp(-alpha2*o21)*(coeff('C',1)*cos(beta2*o21)+coeff('D',1)*sin(beta2*o21)))*sp.cos(ln*x/r)/r**3)
dMxq=dMxql-dMxqr #Considering Edge disturbance

#Equation of Qθ
Qql=sp.simplify(-Dn*p**2*(exp(-alpha1*o1)*((coeff('A',1)-coeff('B',1))*cos(beta1*o1)+(coeff('A',1)+coeff('B',1))*sin(beta1*o1))-exp(-alpha2*o1)*((coeff('C',1)+coeff('D',1))*cos(beta2*o1)-(coeff('C',1)-coeff('D',1))*sin(beta2*o1)))*sp.cos(ln*x/r)/r**3/2**0.5)
Qqr=sp.simplify(-Dn*p**2*(exp(-alpha1*o21)*((coeff('A',1)-coeff('B',1))*cos(beta1*o21)+(coeff('A',1)+coeff('B',1))*sin(beta1*o21))-exp(-alpha2*o21)*((coeff('C',1)+coeff('D',1))*cos(beta2*o21)-(coeff('C',1)-coeff('D',1))*sin(beta2*o21)))*sp.cos(ln*x/r)/r**3/2**0.5)
Qq=Qql-Qqr+dMxq 
print('Qq=' ,Qq) #Considering Edge disturbance

#Equation of Nθ
Nql=sp.simplify(Dn*p**4*(exp(-alpha1*o1)*(coeff('A',0)*sin(beta1*o1)-coeff('B',0)*cos(beta1*o1))-exp(-alpha2*o1)*(coeff('C',0)*sin(beta2*o1)-coeff('D',0)*cos(beta2*o1)))*sp.cos(ln*x/r)/r**3)
Nqr=sp.simplify(Dn*p**4*(exp(-alpha1*o21)*(coeff('A',0)*sin(beta1*o21)-coeff('B',0)*cos(beta1*o21))-exp(-alpha2*o21)*(coeff('C',0)*sin(beta2*o21)-coeff('D',0)*cos(beta2*o21)))*sp.cos(ln*x/r)/r**3)
Nq=Nql+Nqr-4*g*r*cos(radians(o)-o1)*sp.cos(ln*x/r)/pi
print('Nq=',Nq)#Considering Edge disturbance

#Equation of Nxθ
Nxql=sp.simplify(-Dn*p**4*(exp(-alpha1*o1)*(coeff('A',1)*sin(beta1*o1)-coeff('B',1)*cos(beta1*o1))-exp(-alpha2*o1)*(coeff('C',1)*sin(beta2*o1)-coeff('D',1)*cos(beta2*o1)))*sp.sin(ln*l/2/r)/r**3/ln)
Nxqr=sp.simplify(-Dn*p**4*(exp(-alpha1*o21)*(coeff('A',1)*sin(beta1*o21)-coeff('B',1)*cos(beta1*o21))-exp(-alpha2*o21)*(coeff('C',1)*sin(beta2*o21)-coeff('D',1)*cos(beta2*o21)))*sp.sin(ln*l/2/r)/r**3/ln)
Nxq=Nxql-Nxqr+8*g*l*sin(radians(o)-o1)*sp.sin(ln*l/2/r)/pi**2
print('Nxq=',Nxq) #Considering Edge disturbance

#Equation of rotation 
rotationl=sp.simplify((exp(-alpha1*o1)*((coeff('A',1)*(1+p**2*2**0.5)+coeff('B',1)*(1-k*2**0.5))*cos(beta1*o1)+(coeff('A',1)*(1-k*2**0.5)-coeff('B',1)*(1+p**2*2**0.5))*sin(beta1*o1))+exp(alpha2*o1)*((coeff('C',1)*(1-p**2*2**0.5)-coeff('D',1)*(1+k*2**0.5))*cos(beta2*o1)+(coeff('C',1)*(1+k*2**0.5)+coeff('B',1)*(1-p**2*2**0.5))*sin(beta2*o1)))/(r*p**2*2**0.5)*sp.cos(ln*x/r))
rotationr=sp.simplify((exp(-alpha1*o21)*((coeff('A',1)*(1+p**2*2**0.5)+coeff('B',1)*(1-k*2**0.5))*cos(beta1*o21)+(coeff('A',1)*(1-k*2**0.5)-coeff('B',1)*(1+p**2*2**0.5))*sin(beta1*o21))+exp(alpha2*o21)*((coeff('C',1)*(1-p**2*2**0.5)-coeff('D',1)*(1+k*2**0.5))*cos(beta2*o21)+(coeff('C',1)*(1+k*2**0.5)+coeff('B',1)*(1-p**2*2**0.5))*sin(beta2*o21)))/(r*p**2*2**0.5)*sp.cos(ln*x/r))
rotation=rotationl-rotationr #Considering Edge disturbance

#Equation of horizontal displacement of edge of shell
Ul=sp.simplify(-Dn*p**6*(exp(-alpha1*o1)*((coeff('A',0)-coeff('B',0))*sin(beta1*o1)-(coeff('A',0)+coeff('B',0))*cos(beta1*o1))+exp(-alpha2*o1)*((coeff('C',0)+coeff('D',0))*sin(beta2*o1)-(coeff('D',0)-coeff('C',0))*cos(beta2*o1)))*sp.sin(ln*l/2/r)/(2**0.5*E*d*r**2*ln**3))
Ur=sp.simplify(-Dn*p**6*(exp(-alpha1*o21)*((coeff('A',0)-coeff('B',0))*sin(beta1*o21)-(coeff('A',0)+coeff('B',0))*cos(beta1*o21))+exp(-alpha2*o21)*((coeff('C',0)+coeff('D',0))*sin(beta2*o21)-(coeff('D',0)-coeff('C',0))*cos(beta2*o21)))*sp.sin(ln*l/2/r)/(2**0.5*E*d*r**2*ln**3))
U=Ul+Ur-8*g*l**3*cos(radians(o)-o1)*sp.sin(ln*l/2/r)/(pi**4*r*E*d)
print('Ushell=',U)
jd=pi*r*(E*d*r*(pi/l)**4+E*d**3/12/r**5) #a factor below membrane formula

#equation of vertical Displacement of edge of shell
Vl=sp.simplify(Dn*p**6*(exp(-alpha1*o1)*((coeff('A',1)+coeff('B',1))*cos(beta1*o1)-(coeff('A',1)-coeff('B',1))*sin(beta1*o1))-exp(-alpha2*o1)*((coeff('C',1)-coeff('D',1))*cos(beta2*o1)+(coeff('C',1)+coeff('D',1))*sin(beta2*o1)))*sp.cos(ln*x/r)/(2**0.5*E*d*r**2*ln**4))
Vr=sp.simplify(Dn*p**6*(exp(-alpha1*o21)*((coeff('A',1)+coeff('B',1))*cos(beta1*o21)-(coeff('A',1)-coeff('B',1))*sin(beta1*o21))-exp(-alpha2*o21)*((coeff('C',1)-coeff('D',1))*cos(beta2*o21)+(coeff('C',1)+coeff('D',1))*sin(beta2*o21)))*sp.cos(ln*x/r)/(2**0.5*E*d*r**2*ln**4))
V=Vl-Vr-8*g*sin(radians(o)-o1)*sp.cos(ln*x/r)/jd

#equation of deflection of edge beam due to dead weight 
Wl=sp.simplify((exp(-alpha1*o1)*(coeff('A',0)*cos(beta1*o1)+coeff('B',0)*sin(beta1*o1))+exp(-alpha2*o1)*(coeff('C',0)*cos(beta2*o1)+coeff('D',0)*sin(beta2*o1)))*sp.cos(ln*x/r))
Wr=sp.simplify((exp(-alpha1*o21)*(coeff('A',0)*cos(beta1*o21)+coeff('B',0)*sin(beta1*o21))+exp(-alpha2*o21)*(coeff('C',0)*cos(beta2*o21)+coeff('D',0)*sin(beta2*o21)))*sp.cos(ln*x/r))
W=Wl+Wr+8*g*cos(radians(o)-o1)*sp.cos(ln*x/r)/jd


ubeam=sp.simplify((l/pi)**2*Nxq/Ar/E+(l/pi)**2*a1**2*Nxq/F+(Nq*sin(radians(o))-Qq*cos(radians(o)))*(l/pi)**3*a1/F-W1*(l/pi)**3*a1/F)       
#print("Ubeam=",ubeam)

Verbeam=sp.simplify((Nq*sin(radians(o))-Qq*cos(radians(o)))*(l/pi)**4/F+Nxq*(l/pi)**3*a1/F-W1*(l/pi)**4/F)

H=sp.simplify(Nq*cos(radians(o))+Qq*sin(radians(o)))
#print('H=',H)
Ldef=sp.simplify(U-ubeam)
#print('Ldef=',Ldef)
Vdef=sp.simplify(V*sin(radians(o))-W*cos(radians(o))-Verbeam)
#print('Vdef=',Vdef)

#solution of the co-efficients based on the boundary condition
A=sp.Symbol('A')
B=sp.Symbol('B')
C=sp.Symbol('C')
D=sp.Symbol('D')
coefficients_solution=sp.solve([Mq,H,Ldef,Vdef],[A,B,C,D])
print(coefficients_solution)
o12=0
x=0
A1=coefficients_solution[A]
B1=coefficients_solution[B]
C1=coefficients_solution[C]
D1=coefficients_solution[D]
print(A1,B1,C1,D1)
def coeffn( t, n ):
    A=A1
    B=B1
    C=C1
    D=D1
    X=0
    Y=0
    X1=0
    Y1=0
    for i in range(n):
        X=(-alpha1*A+beta1*B)
        Y=(-beta1*A-alpha1*B)
        X1=(-alpha2*C+beta2*D)
        Y1=(-beta2*C-alpha2*D)
        A=X
        B=Y
        C=X1
        D=Y1
    if t=='A':
        return(A)
    elif t=='B':
        return(B)
    elif t=='C':
        return(C)
    else:
        return(D)
#print(-Dn*((exp(-alpha1*o12)*(coeffn('A',2)*cos(beta1*o12)+coeffn('B',2)*sin(beta1*o12)))+exp(-alpha2*o12)*(coeffn('C',2)*cos(beta2*o12)+coeffn('D',2)*sin(beta2*o12)))*cos(ln*x/r)/r**2)

#Calculation of the stress resultants' value at different location
Mq_a=np.zeros(31)                  #a stands for array
Mq_a1=np.zeros(31)                 #array to hold the value or to save the value 
Mq_a2=np.zeros(31)
Nq_a=np.zeros(31)
Nq_a1=np.zeros(31)
Nq_a2=np.zeros(31)
Nxq_a=np.zeros(31)
Nxq_a1=np.zeros(31)
Nxq_a2=np.zeros(31)
Mx_a=np.zeros(31)
Nx_a=np.zeros(31)
Nx_a1=np.zeros(31)
Nx_a2=np.zeros(31)
Mxq_a=np.zeros(31)
Qx_a=np.zeros(31)
Qq_a=np.zeros(31)
angle=np.zeros(31)   
#print('Mq=', Mq)

for i in range(0,31):
    o12=radians(o)*i/30
    Mq_a1[i]=-Dn*((exp(-alpha1*o12)*(coeffn('A',2)*cos(beta1*o12)+coeffn('B',2)*sin(beta1*o12)))+exp(-alpha2*o12)*(coeffn('C',2)*cos(beta2*o12)+coeffn('D',2)*sin(beta2*o12)))*cos(ln*x/r)/r**2
    
    Nq_a1[i]=(Dn*p**4*(exp(-alpha1*o12)*(coeffn('A',0)*sin(beta1*o12)-coeffn('B',0)*cos(beta1*o12))-exp(-alpha2*o12)*(coeffn('C',0)*sin(beta2*o12)-coeffn('D',0)*cos(beta2*o12)))*sp.cos(ln*x/r)/r**3)
    
    Nx_a1[i]=(-Dn*p**4*(exp(-alpha1*o12)*(coeffn('A',2)*sin(beta1*o12)-coeffn('B',2)*cos(beta1*o12))-exp(-alpha2*o12)*(coeffn('C',2)*sin(beta2*o12)-coeffn('D',2)*cos(beta2*o12)))*cos(ln*x/r)/r**3/ln**2)
    
    Nxq_a1[i]=(-Dn*p**4*(exp(-alpha1*o12)*(coeffn('A',1)*sin(beta1*o12)-coeffn('B',1)*cos(beta1*o12))-exp(-alpha2*o12)*(coeffn('C',1)*sin(beta2*o12)-coeffn('D',1)*cos(beta2*o12)))*sin(ln*l/2/r)/r**3/ln)
    angle[i]=o*i/30
for i in range(0,31):
    o12=radians(2*o-o*i/30)
    o123=radians(o)*i/30
    Mq_a2[i]=-Dn*((exp(-alpha1*o12)*(coeffn('A',2)*cos(beta1*o12)+coeffn('B',2)*sin(beta1*o12)))+exp(-alpha2*o12)*(coeffn('C',2)*cos(beta2*o12)+coeffn('D',2)*sin(beta2*o12)))*cos(ln*x/r)/r**2
    Mq_a[i]=Mq_a1[i]+Mq_a2[i]
    
    Nq_a2[i]=(Dn*p**4*(exp(-alpha1*o12)*(coeffn('A',0)*sin(beta1*o12)-coeffn('B',0)*cos(beta1*o12))-exp(-alpha2*o12)*(coeffn('C',0)*sin(beta2*o12)-coeffn('D',0)*cos(beta2*o12)))*sp.cos(ln*x/r)/r**3)
    Nq_a[i]=Nq_a1[i]+Nq_a2[i]-4*g*r*cos(radians(o)-o123)*cos(ln*x/r)/pi
    
    Nx_a2[i]=(-Dn*p**4*(exp(-alpha1*o12)*(coeffn('A',2)*sin(beta1*o12)-coeffn('B',2)*cos(beta1*o12))-exp(-alpha2*o12)*(coeffn('C',2)*sin(beta2*o12)-coeffn('D',2)*cos(beta2*o12)))*cos(ln*x/r)/r**3/ln**2)
    Nx_a[i]=Nx_a1[i]+Nx_a2[i]-8*g*l**2*cos(radians(o)-o123)*cos(ln*x/r)/pi**3/r
    
    Nxq_a2[i]=(-Dn*p**4*(exp(-alpha1*o12)*(coeffn('A',1)*sin(beta1*o12)-coeffn('B',1)*cos(beta1*o12))-exp(-alpha2*o12)*(coeffn('C',1)*sin(beta2*o12)-coeffn('D',1)*cos(beta2*o12)))*sin(ln*l/2/r)/r**3/ln)
    Nxq_a[i]=Nxq_a1[i]-Nxq_a2[i]+8*g*l*sin(radians(o)-o123)*sin(ln*l/2/r)/pi**2
    print(Nxq_a[i],angle[i])

#finals graphs
plt.figure(1)
plt.plot(Mq_a,angle)
plt.ylabel('Mq')
plt.figure(2)
plt.plot(Nq_a,angle,'k')
plt.ylabel('Nq')
plt.figure(3)
plt.plot(Nx_a,angle,'k')
plt.ylabel('Nx')
plt.figure(4)
plt.plot(Nxq_a,angle,'k')
plt.ylabel('Nxq')
plt.show()