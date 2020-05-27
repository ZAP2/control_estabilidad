import math
import numpy as np
from decimal import Decimal
from matplotlib import pyplot as plt 

#DATOS DEL LEARJET
m=float(6124)                                      			#masa
T=101000                                      			#Fuerza de empuje
#Dimensiones
d_MotorCDGx=3										#Distancia entre un motor y el centro de gravedad del avion respecto del eje X
d_MotorCDGy=10										#Distancia entre un motor y el centro de gravedad del avion respecto del eje Y
d_MotorCDGz=1										#Distancia entre un motor y el centro de gravedad del avion respecto del eje Z
S=22													#Superficie alar
b=10.4													#Envergadura del ala
#Coeficientes aerodinamicos:
Cl0=1												#Coeficientes de balanceo
Clbeta=-0.5
Cldelta_a=1
Cldelta_r=0.007

Cn0=0.06												#Coeficientes de cabeceo
Cnalfa=-0.7
Cndelta_e=-1.6

Cnbeta=0.12												#Coeficientes de guiñada
Cndelta_a=1
Cndelta_r=-0.074

Cd = 0.15
#Componentes del Tensor de Inercia (Jxy=Jyz=0):
Ix=0
Iy=0
Iz=0
Jxz=0

g=9.81
rho=1.225
#Aceleraciones:
uu=0													#Aceleración lineal en eje X
vv=0												#Aceleración lineal en eje Y
ww=0												#Aceleración lineal en eje Z
pp=0													#Aceleración angular en eje X
qq=0												#Aceleración angular en eje Y
rr=0

#Angulos de resbalamiento, balanceo y guiñada
beta=0
    #betaR = math.radians(beta)
fiM=0
    #fiM = math.radians(fiM)
fi=0
    #fiR = math.radians(fi)
#Angulos de alerones y timon
delta_a=0
    #delta_aR = math.radians(delta_a)
delta_r=0
    #delta_rR = math.radians(delta_r)
#Velocidad lineal (componente Y)
v=0
u=23
w=0
#Velocidad angular (componentes X y Z)
p=0
r=0
q=0

print( 0.5*rho*(u**2)*S*b)

def FProgLong(delta_eR):
    #Hallar valor de alfa:
    alfa_0 = 10*(math.pi/180)
    
    #condiciones iniciales
    iteracion = 1
    incremento = 1
    rel = 0.9

    #condiciones limitantes
    iteracionMax = 9000
    tolerancia = 1e-9
    
    while iteracion<iteracionMax and incremento>tolerancia:
        Mc = Iy*qq - (Iz-Iz)*p*r + Jxz*(p**2 - r**2)
        Mt = T*math.cos(alfa_0)*math.cos(beta)*d_MotorCDGy - T*math.sin(alfa_0)*d_MotorCDGy
        PresionDin= 0.5*rho*(u**2)*S*b
        
        alfa_i = (((Mc - Mt) / PresionDin) - Cn0 - Cndelta_e*delta_e) / Cnalfa

        incremento = abs((alfa_i - alfa_0)/alfa_0)
        iteracion = iteracion + 1
        alfa_iG = alfa_i*(180/math.pi)
        
        print(iteracion, "Control incremento: ", incremento, "\nControl alfa (en grados): ", alfa_iG)

        alfa_0 = rel*alfa_i + (1-rel)*alfa_0

    #Con este alfa_i, hallar valor de theta:
    Fres_x = float((m*(uu - r*v + q*w)))
    Fax = float(PresionDin*(Cd))
    Ftx = float(T*math.cos(alfa_i)*math.cos(beta))
    senTheta = (Fres_x - Fax - Ftx)/(-m*g)
    
    theta = math.asin(senTheta)

    thetaG = theta*(180/math.pi)
    print()
    print("Fuerza resultante: ", int(Fres_x)*1e-3, "[kN]", "\nFuerza Aerodinamica: ",
            int(Fax)*1e-3, "[kN]", "\nFuerza Thrust: ", int(Ftx)*1e-3, "[kN]", "\nPeso: ",
            int(m*g)*1e-3, "[kN]", "\nSeno Theta: ", senTheta)

    print()
    print("El valor final de alfa en grados es: ", alfa_iG)
    print("El valor de theta en grados es: ", thetaG)
    print()
    print("FIN DEL PROGRAMA") 
    print()

    
delta_e = float(input("Introduce el valor de delta_e en radianes: "))
FProgLong(delta_e)