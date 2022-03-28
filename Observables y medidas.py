import numpy as np
import math

def sumacplx(a, b): #Funcion para sumar numeros complejos
    real = a[0] + b[0]
    img = a[1] + b[1]
    return (real, img)

def multcplx(a, b): #Funcion para multiplicar numeros complejos
    real = (a[0] * b[0]) - (a[1] * b[1])
    img = (a[0] * b[1]) + (a[1] * b[0])
    return (real, img)

def conjugadocplx(a): #Funcion para realizar el conjugado de un numero complejo
    return (a[0], -1*a[1])

def modulocplx(a): #Funcion para calcular el modulo de un numero complejo
    return math.sqrt((a[0])**2 + (a[1])**2)

def productointerno(A, B): #Funcion para calcular el producto interno de vectores de numeros complejos
    for i in range(len(A)):
        A[i][0] = conjugadocplx(A[i][0])
    filas = len(A)
    matriz = [0 for i in range(filas)]
    for i in range(filas):
        matriz[i] = A[i][0]
    suma = (0,0)
    for i in range(len(matriz)):
        suma = sumacplx(suma, multcplx(matriz[i], B[i][0]))
    return suma

def normavector(v): #Funcion para hallar la norma de un vector
    w = [a[:] for a in v]
    resp = productointerno(v,w)
    return math.sqrt(resp[0])

def divComplex(c1, c2): #Funcion para dividir complejos
    return c1/c2

def probabilidadenpunto(posicion, vector): #Funcion para conocer la probabilidad en un punto
    return ((modulocplx(vector[posicion][0]))**2)/((normavector(vector))**2)

def probabilidadVectorAOtro(vector1, vector2): #Funcion para conocer la probabilidad de ir de un estado a otro
    multnorma = normavector(vector1) * normavector(vector2)
    prod = productointerno(vector2, vector1)
    return divComplex(complex(prod[0], prod[1]), multnorma)

def matrizhermitiana(A): #Funcion para conocer si una matriz es hermitiana
    filas = len(A)
    columnas = len(A[0])
    if filas != columnas:
        return False
    B = [[0 for i in range(filas)]for i in range(columnas)]
    for i in range(filas):
        for j in range(columnas):
            B[j][i] = A[i][j]
    for i in range(filas):
        for j in range(columnas):
            B[i][j] = conjugadocplx(B[i][j])
    for i in range(filas):
        for j in range(columnas):
            if B[i][j] != A[i][j]:
                return False
    return True

def productomatrices(A, B): #Funcion para realizar el producto entre matrices
    filas = len(A)
    columnas = len(A[0])
    matriz = []
    for i in range(filas):
        fila = []
        for j in range(filas):
            suma = (0,0)
            for k in range(columnas):
                suma = sumacplx(suma, multcplx(A[i][k], B[k][j]))
            fila += [suma]
        matriz += [fila]
    return matriz

def accionmatrizvectorComplex(A,v): #Funcion para realizar la accion de una matriz de numeros complejos sobre un vector de numeros complejos 
    filas = len(A)
    columnas = len(A[0])
    matriz = []
    for i in range(filas):
        suma = (0,0)
        for j in range(columnas):
            suma = sumacplx(suma, multcplx(A[i][j], v[j][0]))
        matriz.append([suma])
    return matriz

def valoresperado(matriz, vector): #Funcion para hallar el valor esperado de un observable y un vector ket
    return productointerno(accionmatrizvectorComplex(matriz, vector), vector)

def crearunitaria(n):#Funcion para crear una matriz unitaria con tamano nxn
    unitaria = []
    for i in range(n):
        fila = []
        for j in range(n):
            if j == i:
                fila += [(1,1)]
            else:
                fila += [(0, 0)]
        unitaria.append(fila)
    return unitaria

def inversamatriz(A): #Funcion par calcular la inversa de una matriz
    matriz = []
    for i in range(len(A)):
        fila = []
        for j in range(len(A[0])):
            fila = fila + [(-1*A[i][j][0],-1*A[i][j][1])]
        matriz = matriz + [fila]
    return matriz

def adicionmatrices(A, B): #Funcion para realizar la adicion de matrices
    filas = len(A)
    columnas = len(A[0])
    matriz = []
    for i in range(filas):
        fila = []
        for j in range(columnas):
            fila = fila + [sumacplx(A[i][j], B[i][j])]
        matriz = matriz + [fila]
    return matriz

def multescalarmatriz(c, A): #Funcion para multiplicar un escalar por una matriz
    matriz = []
    for i in range(len(A)):
        fila = []
        for j in range(len(A[0])):
            fila = fila + [(c*A[i][j][0],c*A[i][j][1])]
        matriz = matriz + [fila]
    return matriz

def mediaobservable(observable, vectorket): #Funcion para calcular la media de un observable
    valore = valoresperado(observable, vectorket)
    unitaria = crearunitaria(len(observable[0]))
    return adicionmatrices(observable, inversamatriz(multescalarmatriz(valore[0], unitaria)))

def traspuestacomplex(A): #Funcion para calcular la traspuesta de una matriz
    filas = len(A)
    columnas = len(A[0])
    matriz = [[0 for i in range(filas)]for i in range(columnas)]
    for i in range(filas):
        for j in range(columnas):
            matriz[j][i] = A[i][j]
    return matriz

def conjugadamatriz(A): #Funcion para calcular la conjugada de una matriz
    matriz = []
    for i in range(len(A)):
        fila = []
        for j in range(len(A[0])):
            fila = fila + [conjugadocplx(A[i][j])]
        matriz = matriz + [fila]
    return matriz

def adjuntamatriz(A): #Funcion para calcular la adjunta de una matriz
    return traspuestacomplex(conjugadamatriz(A))

def traceComplex(matriz): #Funcion para realizar la suma de la diagonal de una matriz
    suma = (0,0)
    for i in range(len(matriz)):
        suma += sumacplx(suma, matriz[i][i])
    return suma

def productointernomatriz(A, B): #Funcion para calcular el producto interno entre dos matricez
    tamano = len(A)
    suma = (0,0)
    A = adjuntamatriz(A)
    for i in range(tamano):
        for j in range(tamano):
            suma += sumacplx(suma, multcplx(A[i][j], B[i][j]))
    return suma

def varianza(observable, vectorket): #Funcion para calcular la varianza de un observable
    media = mediaobservable(observable, vectorket)
    return valoresperado(productomatrices(media, media),vectorket)

def main():
    '''
    posiciones = int(input("Numero de posiciones: "))
    vectorinicial = []
    for i in range(posiciones):
        valores = input("Escriba el valor de la amplitud, parte real e imaginaria separada por espacios: ")
        amplitud = tuple(int(x) for x in valores.split(" "))
        vectorinicial += [[amplitud]]
    posicion = int(input("Digite la posicion de la que desea conocer la probabilidad: "))
    vector1 = [x[:] for x in vectorinicial]
    print(probabilidadenpunto(posicion, vector1))
    otro = input("¿Desea calcular la probabilidad de transitar a otro vector? (Si/No): ") 
    if otro == "Si":
        vector2 = []
        for i in range(posiciones):
            valores = input("Escriba el valor de la amplitud, parte real e imaginaria separada por espacios: ")
            amplitud = tuple(int(x) for x in valores.split(" "))
            vector2 += [[amplitud]]
        print(probabilidadVectorAOtro(vector2, vectorinicial))
    '''
    tamano = int(input("tamaño de la matriz: "))
    observable = []
    for i in range(tamano):
        fila = []
        for j in range(tamano):
            valores = input("Escriba el valor de la amplitud en la posicion " + str(i) +","+ str(j) + " parte real e imaginaria separada por espacios: ")
            amplitud = tuple(int(x) for x in valores.split(" "))
            fila += [amplitud]
        observable += [fila]
    vectorket = []
    
    for i in range(tamano):
        valores = input("Escriba el valor de la amplitud del vector ket, parte real e imaginaria separada por espacios: ")
        amplitud = tuple(int(x) for x in valores.split(" "))
        vectorket += [[amplitud]]
    if matrizhermitiana(observable):
        print(mediaobservable(observable, vectorket))
        print(varianza(observable, vectorket))
    
main()