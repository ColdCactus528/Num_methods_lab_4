import matplotlib.pyplot as plt
import numpy as np
import math

def FunctionLab4(x):
    return x * math.atan(x) + x/2 * math.cos(x) - 3

def FunctionDerivative(x):
    return x/(x*x + 1) - 0.5*x*math.sin(x) + math.cos(x)*0.5 + math.atan(x)

def DichotomyMethod(a, b, epsilon):
    count = 0
    print("Итерация, вычисленная по формуле - " + str(math.floor(math.log((b-a) / epsilon, 2)) + 1))
    if (FunctionLab4(a) == 0):
        print("Итерация, на которой был получен корень - 1")
        return a

    elif (FunctionLab4(b) == 0):
        print("Итерация, на которой был получен корень - 1")
        return b

    while (b-a > epsilon):
        c = (b+a)/2
        count += 1

        if (FunctionLab4(a) * FunctionLab4(c) < 0):
            b = c

        if (FunctionLab4(b) * FunctionLab4(c) < 0):
            a = c

    print("Итерация, на которой был получен корень - " + str(count))
    return c

def NewtonMethod(a, b, epsilon):
    x = a
    x1 = b
    mass = []
    mass.append(x)
    while (abs(x-x1) > epsilon):
        x1 = x
        x = x - FunctionLab4(x) / FunctionDerivative(x)
        mass.append(x)

    print("Итерация, на которой был получен корень - " + str(len(mass)))
    return mass

def ChordMethod(a, b, epsilon):
    x0 = b
    fx0 = FunctionLab4(x0)
    x = a
    x1 = b
    mass = []
    mass.append(x)
    while (abs(x1 - x) > epsilon):
        x1 = x
        x = x - FunctionLab4(x) * (x - x0) / (FunctionLab4(x) - fx0)
        mass.append(x)

    print("Итерация, на которой был получен корень - " + str(len(mass)))
    return mass

def RateOfConvergenceMethod(n, a, b, Method):
    mass = Method(a, b, 1e-12)
    x1 = mass[n-2]
    x2 = mass[n-1]
    x3 = mass[n]
    x4 = mass[n+1]

    return math.log(abs((x4 - x3) / (x3-x2))) / math.log(abs((x3 - x2) / (x2-x1)))

epsilon1 = 1e-3
epsilon2 = 1e-6
epsilon3 = 1e-9
a = 0
b = 10
print("Для эпсилон = " + str(epsilon1))
print(str(DichotomyMethod(a, b, epsilon1)) + " - метод дихотомии\n")
print( str(NewtonMethod(a, b, epsilon1)[-1]) + " - метод Ньютона\n")
print(str(ChordMethod(a, b, epsilon1)[-1]) + " - метод хорд")
print("\nДля эпсилон = " + str(epsilon2))
print(str(DichotomyMethod(a, b, epsilon2)) + " - метод дихотомии\n")
print(str(NewtonMethod(a, b, epsilon2)[-1]) + " - метод Ньютона\n")
print(str(ChordMethod(a, b, epsilon2)[-1]) + " - метод хорд")
print("\nДля эпсилон = " + str(epsilon3))
print(str(DichotomyMethod(a, b, epsilon3)) + " - метод дихотомии\n")
print(str(NewtonMethod(a, b, epsilon3)[-1]) + " - метод Ньютона\n  ")
print(str(ChordMethod(a, b, epsilon3)[-1]) + " - метод хорд")
print("\n" + str(RateOfConvergenceMethod(10, a, b, ChordMethod)) + " - скорость сходимости метода хорд")
print(str(RateOfConvergenceMethod(6, a, b, NewtonMethod)) + "- скорость сходимости метода Ньютона")
