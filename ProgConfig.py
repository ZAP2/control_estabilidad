seleccion=input("Selecciona una configuracion: ")

def SeleccionConfig(configuracion):
    if configuracion=="configuracion1":
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

    elif configuracion=="configuracion2":
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

    elif configuracion=="configuracion3":
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

    elif configuracion=="configuracion4":
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


SeleccionConfig(seleccion)