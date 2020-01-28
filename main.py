from lib import *
import timeit
import msvcrt

k = 0
print("1) Ввести матрицу\n2) Рандомная матрица\n3) Метод Гаусса\n4) Метод Kрамера\n5) Метод верхней релаксации\n6) Метод Зейделя\n7) Метод простых итераций\n8) Метод Гаусса-Жордана\n9) Метод LU-разложения\n10) Все методы решения\n11) Выход")
while k != 11:
    k = int(input("Введите действие "))
    if k == 1:
        n = int(input("Введите размер матрицы "))
        tmp = inputList(n * n, "Введите матрицу ")
        a = numpy.matrix([[tmp[i*n + j] for j in range(n)] for i in range(n)])
        y = numpy.array(inputList(n, "Введите свободный член "))
        Print(a,y)
    elif k == 2:
        n = int(input("Введите размер матрицы "))
        temp = [[0.0 for i in range(n)] for j in range(n)]
        for i in range(n):
            for j in range(n):
                temp[i][j] = round(random.uniform(0, 20),1)
            temp[i][i] = round(sum(temp[i]),1)
        a = numpy.matrix(temp)
        y = numpy.array([round(random.uniform(0, 30),1) for i in range(n)])
        Print(a,y)
    elif k == 3:
        x = Gaus(a,y)
        print("x = ", x)
    elif k == 4:
        x = KramerMethod(a,y)
        print("x = ", x)
    elif k == 5:
        x = relax(a,y)
        print("x = ", x)
    elif k == 6:
        x = seidel(a,y)
        print("x = ", x)
    elif k == 7:
        x = simpleIterationMethod(a,y)
        print("x = ", x)
    elif k == 8:
        x = GausJordan(a,y)
        print("x = ", x)
    elif k == 9:
        x = LU(a,y)
        print("x = ", x)
    elif k == 10:
        t = timeit.default_timer()
        x = Gaus(a,y)
        t = round(timeit.default_timer() - t,7)
        print("\nМетод Гаусса: ", x, "\nВремя работы метода = ", t, "\n")
        t = timeit.default_timer()
        x = KramerMethod(a,y)
        t = round(timeit.default_timer() - t,7)
        print("Метод Крамера: ", x, "\nВремя работы метода = ", t, "\n")
        t = timeit.default_timer()
        x = relax(a,y)
        t = round(timeit.default_timer() - t,7)
        print("Метод верхней релаксации: ", x, "\nВремя работы метода = ", t, "\n")
        t = timeit.default_timer()
        x = seidel(a,y)
        t = round(timeit.default_timer() - t,7)
        print("Метод Зейделя: ", x, "\nВремя работы метода = ", t, "\n")
        t = timeit.default_timer()
        x = simpleIterationMethod(a,y)
        t = round(timeit.default_timer() - t,7)
        print("Метод простых итераций: ", x, "\nВремя работы метода = ", t, "\n")
        t = timeit.default_timer()
        x = GausJordan(a,y)
        t = round(timeit.default_timer() - t,7)
        print("Метод Гаусса-Жордана ", x, "\nВремя работы метода = ", t, "\n")
        t = timeit.default_timer()
        x = LU(a,y)
        t = round(timeit.default_timer() - t,7)
        print("Метод LU-разложения ", x, "\nВремя работы метода = ", t, "\n")
msvcrt.getch()








