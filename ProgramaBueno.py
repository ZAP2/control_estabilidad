repetir="si"

def SeleccionConfig(configuracion):
    if configuracion==1:
        m=28000                                      			#masa
        T=101000                                      			#Fuerza de empuje
        #Dimensiones
        d_MotorCDGx=10										#Distancia entre un motor y el centro de gravedad del avion respecto del eje X
        d_MotorCDGy=3										#Distancia entre un motor y el centro de gravedad del avion respecto del eje Y
        d_MotorCDGz=1										#Distancia entre un motor y el centro de gravedad del avion respecto del eje Z
        S=125													#Superficie alar
        b=20													#Envergadura del ala
        #Coeficientes aerodinamicos:
        Cl0=1												#Coeficientes de balanceo
        Clbeta=1
        Cldelta_a=1
        Cldelta_r=1
        Cn0=1												#Coeficientes de cabeceo
        Cnalfa=1
        Cndelta_e=1
        Cnbeta=1												#Coeficientes de guiñada
        Cndelta_a=1
        Cndelta_r=1
        #Componentes del Tensor de Inercia (Jxy=Jyz=0):
        Ix=0
        Iy=0
        Iz=0
        Jxz=0

        print("Masa: ", m)
        print("Dimensiones de ala: ", S, b)

    elif configuracion==2:
        m=45000                                      			#masa
        T=202000                                      			#Fuerza de empuje
        #Dimensiones
        d_MotorCDGx=10										#Distancia entre un motor y el centro de gravedad del avion respecto del eje X
        d_MotorCDGy=3										#Distancia entre un motor y el centro de gravedad del avion respecto del eje Y
        d_MotorCDGz=1										#Distancia entre un motor y el centro de gravedad del avion respecto del eje Z
        S=125													#Superficie alar
        b=20													#Envergadura del ala
        #Coeficientes aerodinamicos:
        Cl0=1												#Coeficientes de balanceo
        Clbeta=1
        Cldelta_a=1
        Cldelta_r=1
        Cn0=1												#Coeficientes de cabeceo
        Cnalfa=1
        Cndelta_e=1
        Cnbeta=1												#Coeficientes de guiñada
        Cndelta_a=1
        Cndelta_r=1
        #Componentes del Tensor de Inercia (Jxy=Jyz=0):
        Ix=0
        Iy=0
        Iz=0
        Jxz=0
        
        print("Masa: ", m)
        print("Dimensiones de ala: ", S, b)

    elif configuracion==3:
        m=32000                                      			#masa
        T=151000                                      			#Fuerza de empuje
        #Dimensiones
        d_MotorCDGx=10										#Distancia entre un motor y el centro de gravedad del avion respecto del eje X
        d_MotorCDGy=3										#Distancia entre un motor y el centro de gravedad del avion respecto del eje Y
        d_MotorCDGz=1										#Distancia entre un motor y el centro de gravedad del avion respecto del eje Z
        S=125													#Superficie alar
        b=20												#Envergadura del ala
        #Coeficientes aerodinamicos:
        Cl0=1												#Coeficientes de balanceo
        Clbeta=1
        Cldelta_a=1
        Cldelta_r=1
        Cn0=1												#Coeficientes de cabeceo
        Cnalfa=1
        Cndelta_e=1
        Cnbeta=1												#Coeficientes de guiñada
        Cndelta_a=1
        Cndelta_r=1
        #Componentes del Tensor de Inercia (Jxy=Jyz=0):
        Ix=0
        Iy=0
        Iz=0
        Jxz=0

        print("Masa: ", m)
        print("Dimensiones de ala: ", S, b)

    elif configuracion==4:
        m=29000                                      			#masa
        T=211000                                      			#Fuerza de empuje
        #Dimensiones
        d_MotorCDGx=10										#Distancia entre un motor y el centro de gravedad del avion respecto del eje X
        d_MotorCDGy=3										#Distancia entre un motor y el centro de gravedad del avion respecto del eje Y
        d_MotorCDGz=1										#Distancia entre un motor y el centro de gravedad del avion respecto del eje Z
        S=125													#Superficie alar
        b=20													#Envergadura del ala
        #Coeficientes aerodinamicos:
        Cl0=1												#Coeficientes de balanceo
        Clbeta=1
        Cldelta_a=1
        Cldelta_r=1
        Cn0=1												#Coeficientes de cabeceo
        Cnalfa=1
        Cndelta_e=1
        Cnbeta=1												#Coeficientes de guiñada
        Cndelta_a=1
        Cndelta_r=1
        #Componentes del Tensor de Inercia (Jxy=Jyz=0):
        Ix=0
        Iy=0
        Iz=0
        Jxz=0

        print("Masa: ", m)
        print("Dimensiones de ala: ", S, b)

def DensidadAire(altura):
    #Funcion densidad del aire
    rho_0=1.225											#densidad inicial [kg/m3]
    alfa_T=6.5e-3										#Ratio descenso temperatura [K/m]
    theta_0=288.15										#temperatura inicial [K]
    g=9.81												#gravedad [m/s2]
    Ra=287.05											#constante del aire [J/kg*K]
    rho=rho_0*((1-(alfa_T*altura)/(theta_0))**(((g)/(Ra*alfa_T))-1))
    print("La densidad es: ", rho, "kg/m3")
    
def AlturaCte(VueloConstante):
    if VueloConstante=="si":
        print("Rho es cte")
    elif VueloConstante=="no":
        print("Todavia no hemos programado esta parte")

def FuncionEstDin(EstDin):
    if EstDin=="e":
        #Aceleraciones:
        uu=0													#Aceleración lineal en eje X
        vv=0												#Aceleración lineal en eje Y
        ww=0												#Aceleración lineal en eje Z
        pp=0													#Aceleración angular en eje X
        qq=0												#Aceleración angular en eje Y
        rr=0
        
        print("Las aceleraciones lineales y angulares son nulas")


    elif EstDin=="d":
        print("Todavia no hemos programado esta parte")		

def FuncionLongLat(LongLat):
    if LongLat=="long":
        #Angulos de resbalamiento, balanceo y guiñada
        beta=0
        fiM=0
        fi=0
        #Angulos de alerones y timon
        delta_a=0
        delta_r=0
        #Velocidad lineal (componente Y)
        v=0
        #Velocidad angular (componentes X y Z)
        p=0
        r=0

        print("Los angulos de resbalamiento, balanceo, guiñada, alerones y timon, la velocidad lineal en Y y las velocidades angulares en X y Z son nulas")

    
    elif LongLat=="lat":
        print("Todavia no hemos programado esta parte")

while repetir=="si":
    seleccion=int(input("Selecciona una configuracion (1/2/3/4): ")) 
    while 1>seleccion>4:
        print("Esa configuracion no existe, debe estar entre 1 y 4")
        seleccion=int(input("Vuelve a seleccionar una configuracion: "))
    SeleccionConfig(seleccion)

    altura=float(input("Determina la altura de vuelo [m]: "))
    while altura<0 or altura>1000000:
        print("El avion no va a volar a esa altura.")
        altura=float(input("Vuelve a determinar una altura: "))
    DensidadAire(altura)
    VueloConstante=(input("El vuelo es de altura constante? (si/no): "))
    AlturaCte(VueloConstante)

    EstDin=(input("El vuelo es estatico o dinamico? (e/d): "))
    FuncionEstDin(EstDin)
    if EstDin=="e" :
        LongLat=(input("El vuelo es longitudinal o lateral-direccional? (long/lat): "))
        FuncionLongLat(LongLat)

    repetir=input("Quieres volver a ejecutar el programa? (si/no): ")

print("El programa ha finalizado")


