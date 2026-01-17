import numpy as np
import matplotlib.pyplot as plt
import math
import random

def f(x):
    return x**3*math.sin(x)

def E_build(x, y, n, m):
    E = np.empty((m, n+1))
    for i in range (m):
        for j in range(n+1):
            E[i][j]=x[i]**j
    #print(E)
    return E


def a_values (a, _x):
    res = 0
    for i in range (len(a)):
        res += a[i]*_x**i
    return res

def summa_kv_o_1(x, y, n):
    E = E_build(x, y, n, len(x))
    F = np.dot(np.transpose(E), y)
    E = np.matmul(np.transpose(E), E)
    a = np.linalg.solve(E, F)
    approx_values = []
    res = 0
    for i in range(len(x)):
        res += (y[i]-a_values(a, x[i]))**2
        approx_values.append(a_values(a, x[i]))

    plt.figure()
    plt.plot(x, approx_values, label = "Аппроксимация функции")
    plt.scatter(x, y, s=0.5, label = "Экспериментальные точки")
    plt.legend()
    plt.title('Нормальные полиномы')

    return res

def q(q, _x):
    res = 0
    for i in range (len(q)):
        res += q[i]*_x**i
    return res

def apha_build(x, q_1):
    s1 = 0
    s2 = 0
    for i in range (len(x)):
        s1 = s1+ x[i]*q(q_1, x[i])**2
        s2 = s2+ q(q_1, x[i])**2
    return -s1/s2

def beta_build(x, q_1, q_0):
    s1 = 0
    s2 = 0
    for i in range (len(x)):
        s1 += x[i]*q(q_1, x[i])*q(q_0, x[i])
        s2 += q(q_0, x[i])**2
    return -s1/s2

def Q_build(x, n):
    Q = np.zeros((n+1, n+1))
    Q[0][0] = 1
    Q[1][1] = 1
    Q[1][0] = -sum(x)/len(x)
    for i in range (2, n+1):
        alpha = apha_build(x, Q[i-1])
        beta = beta_build(x, Q[i-1], Q[i-2])
        for j in range (n+1):
            Q[i][j] += Q[i-1][j]*alpha
            Q[i][j] += Q[i-2][j]*beta
        for j in range (1, n+1):
            Q[i][j] += Q[i-1][j-1]
    return Q

def coeff_build(x, y, q_1):
    s1 = 0
    s2 = 0
    for i in range (len(x)):
        s1 += q(q_1, x[i])*y[i]
        s2 += q(q_1, x[i])**2
    return s1/s2

def ort_value(x, y, n, _x, Q, a):

    res = 0
    for i in range (n+1):
        res += q(Q[i], _x)*a[i]
    return res


def summa_kv_o_2(x, y, n):
    Q = Q_build(x, n)
    a = np.array([coeff_build(x, y, Q[i]) for i in range(n + 1)])
    approx_values = []
    res = 0
    for i in range(len(x)):
        res += (y[i] - ort_value(x, y, n, x[i], Q, a))**2
        approx_values.append(ort_value(x, y, n, x[i], Q, a))

    plt.figure()
    plt.plot(x, approx_values, label="Аппроксимация функции")
    plt.scatter(x, y, s=0.5, label="Экспериментальные точки")
    plt.legend()
    plt.title('Ортогональные полиномы')

    return res

x = [random.uniform(-1+(i-1)*2/55, -1+i*2/55) for i in range(1, 56)]
for j in range (2):
    for i in range(55):
        x.append(x[i])
x=np.array(x)
x.sort()
y = np.array([f(x[i])+random.uniform(-0.03, 0.03) for i in range(len(x))])
for i in range (1, 6):
    print(i, " ", summa_kv_o_1(x, y, i), " ", summa_kv_o_2(x, y, i))

plt.show()








