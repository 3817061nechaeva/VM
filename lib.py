import numpy
import random
import math

epsilon = 0.0001

def inputList(number, message):
    print(message)
    res = []
    while (len(res) < number):
        tmp = input().split()
        for elem in tmp:
            try:
                val = float(elem)
            except ValueError:
                continue
            else:
                if (len(res) < number):
                    res.append(val)
                else:
                    break
    return res

def KramerMethod(a, y):
    error1 = "Система линейных уравнений имеет множество бесконечных решений"
    error2 = "Система линейных уравнений не имеет решений"
    detC = []
    correct = 0
    n = len(a)
    x = numpy.array([0.0 for i in range(n)])
    for i in range(n):
        detC.append(GetDet(a, y, i))
    det = numpy.linalg.det(a)
    if det == 0:
        for i in range(n):
            correct += detC[i]
        if correct == 0:
            correct = 0
            for i in range(n):
                correct += y[i]
            if correct == 0:
                return error1
            else:
                return error2
        else:
            return error2
    else:
        for i in range(n):
            x[i] = round(detC[i] / det,3)
    return x

def GetDet(a, y, j):
    n = len(a)
    b = a.copy()
    for i in range(n):
        b.A[i][j] = y[i]
    det = numpy.linalg.det(b)
    return det


def Gaus(b,f):
    error1 = "Система линейных уравнений не имеет решения"
    error2 = "Система линейных уравнений имеет бесконечное множество решений"
    a = b.copy()
    y = f.copy()
    k, n = 0, len(a)
    t = 0
    x = numpy.array([0.0 for i in range(n)])
    while k < n:
        max = abs(a.A[k][k])
        index = k
        for i in range(k+1, n):
            if abs(a.A[i][k]) > max:
                max = abs(a.A[i][k])
                index = i
        if max == 0:
            if k == n-1:
                if y[k] !=0:
                    return error1
                else:
                    return error2
            else:
                for i in range(n):
                    t += a.A[k][i]
                if t == 0:
                    if y[k] != 0:
                        return error1
            k += 1
        else :
            if index != k:
                a.A[index], a.A[k] = a.A[k].copy(), a.A[index].copy()
                y[index], y[k] = y[k], y[index]
            for i in range(k,n):
                temp = a.A[i][k]
                if temp == 0:   continue
                for j in range(n):
                    a.A[i][j] = a.A[i][j] / temp
                y[i] = y[i] / temp
                if i == k:  continue
                for j in range(n):
                    a.A[i][j] = a.A[i][j] - a.A[k][j]
                y[i] = y[i] - y[k]
            k +=1
    for k in range(n-1, -1, -1):
        x[k] = round(y[k],3)
        for i in range(n):
            y[i] = round(y[i] - a.A[i][k] * x[k],3)
    return x

def Print(a,y):
    n = len(a)
    for i in range(n):
        for j in range(n):
            print(str(a.A[i][j])+" * x"+str(j), end="")
            if j < n-1:
                print(" +", end=" ")
        print(" = ", y[i], end="\n\n")

def simpleIterationMethod(mat, y):
    n = len(mat)
    precision = 1e-5
    if 0 in mat.diagonal().A[0]:
        return "Главная диагональ состоит из нулевых элементов"
    diag = mat.diagonal().A.copy()[0]
    A = mat.A.copy()
    f = y.copy()
    for i in range(n):
        temp = diag[i]
        f[i] /= temp
        for j in range(n):
            if i == j:
                A[i][j] = 0
            else:
                A[i][j] /= temp
    norm = numpy.linalg.norm(A, numpy.inf)
    if norm >= 1:
        return "Нет диагонального преобладания"
    xnew = f.copy()
    xold = f.copy()
    while True:
        xold = xnew
        xnew = f.copy()
        norm = 0
        for i in range(n):
            for j in range(n):
                xnew[i] -= A[i][j] * xold[j]
            if abs(xnew[i]-xold[i]) > norm:
                norm = abs(xnew[i]-xold[i])
        if norm <= precision:
            break
    for i in range(n):
        xnew[i] = round(xnew[i],3)
    return xnew

def LU(mat, f):
    if numpy.linalg.det(mat) == 0:
        return "Определитель данной матрицы равен 0, невозможно решить данным методом"
    A = numpy.matrix(mat.A.copy())
    n = len(A.A)
    L = numpy.matrix([[0.0]*n for i in range(n)])
    U = numpy.matrix([[0.0]*n for i in range(n)])
    x = numpy.array([0.0 for i in range(n)])
    y = numpy.array([0.0 for i in range(n)])
    for i in range(n):
        for j in range(i, n):
            L.A[j][i] = A.A[j][i] - sum(L.A[j][k]*U.A[k][i] for k in range(j))
            U.A[i][j] = (A.A[i][j] - sum(L.A[i][k]*U.A[k][j] for k in range(i)))/L.A[i][i]
    for i in range(n):
        y[i] = (f[i]-sum(L.A[i][k]*y[k] for k in range(i)))/L.A[i][i]
    for i in range(n-1, -1, -1):
        s = sum(U.A[i][k]*x[k] for k in range(n-1, i-1, -1))
        x[i] = round(y[i] - s,3)
    return x

def seidel(A, b, eps=0.001):
    n = len(A)
    x = numpy.array([.0 for i in range(n)])
    converge = False
    for i in range(n):
        if A.A[i][i] == 0:
            converge = True
            x = "Главная диагональ состоит из нулевых элементов"
            break
    while not converge:
        x_new = numpy.copy(x)
        for i in range(n):
            s1 = sum(A.A[i][j] * x_new[j] for j in range(i))
            s2 = sum(A.A[i][j] * x[j] for j in range(i + 1, n))
            x_new[i] = round((b[i] - s1 - s2) / A.A[i][i],3)
        if abs(x_new[0] - x[0]) > 100000:
            x = "Нет диагонального преобладания"
            break
        converge = math.sqrt(
            sum((x_new[i] - x[i]) ** 2 for i in range(n))) <= eps
        x = x_new
    return x


def check_symmetric(a, tol=1e-8):
    return numpy.all(numpy.abs(a-a.T) < tol)


def relax(A, b, eps=0.001, omega=1.5):
    n = len(A)
    x = numpy.array([.0 for i in range(n)])
    converge = False
    # if (not check_symmetric(A)):
    #     x = "The matrix is asymmetric"
    #     converge = True
    if not converge:
        for i in range(n):
            if A.A[i][i] == 0:
                converge = True
                x = "Главная диагональ состоит из нулевых элементов"
                break
    while not converge:
        x_new = numpy.copy(x)
        for i in range(n):
            s1 = sum(A.A[i][j] * x_new[j] * omega for j in range(i))
            s2 = sum(A.A[i][j] * x[j] * omega for j in range(i + 1, n))
            x_new[i] = (b[i] * omega - s1 - s2) / A.A[i][i] - x[i]*(omega - 1)
        if abs(x_new[0] - x[0]) > 100000:
            x = "Нет диагонального преобладания"
            break
        converge = math.sqrt(
            sum((x_new[i] - x[i]) ** 2 for i in range(n))) <= eps
        x = x_new
    for i in range(n):
        x[i] = round(x[i],3)   
    return x

def GausJordan(b,f):
    error1 = "Система линейных уравнений не имеет решения"
    error2 = "Система линейных уравнений имеет бесконечное множество решений"
    a = b.copy()
    y = f.copy()
    k, n = 0, len(a)
    t = 0
    x = numpy.array([0.0 for i in range(n)])
    while k < n:
        max = abs(a.A[k][k])
        index = k
        for i in range(k+1, n):
            if max == 0:
                max = abs(a.A[i][k])
                index = i
            else:
                break
        if max == 0:
            if k == n-1:
                if y[k] !=0:
                    return error1
                else:
                    return error2
            else:
                for i in range(n):
                    t += a.A[k][i]
                if t == 0:
                    if y[k] != 0:
                        return error1
            k += 1
        else :
            if index != k:
                a.A[index], a.A[k] = a.A[k].copy(), a.A[index].copy()
                y[index], y[k] = y[k], y[index]
            temp = a.A[k][k]
            for i in range(k,n):
                a.A[k][i] /= temp
            y[k] /= temp
            for i in range(n):
                temp = a.A[i][k]
                if i == k: 
                    continue
                for j in range(n):
                    a.A[i][j] -= temp * a.A[k][j]
                y[i] -= temp * y[k]
            k +=1
    for i in range(n):
        x[i] = round(y[i],3)
    return x