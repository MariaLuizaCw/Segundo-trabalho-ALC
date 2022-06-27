import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
class algcom:
    def __init__(self,passos ,t, m, c,k,a1,a2,a3,w1,w2,w3):
        self.passos = passos #delta t
        self.t = t #tempo total de integração 
        #valores dos parâmetros 
        self.m = m
        self.c = c
        self.k = k
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.w1 = w1
        self.w2 = w2
        self.w3 = w3
    
    def F(self,t):
        
        return self.a1*np.sin((self.w1*t)) + self.a2*np.sin((self.w2*t)) +self.a3*np.cos((self.w3*t))
        
    def f(self, t, x, x_linha):
        return (self.F(t) - self.c*x_linha - self.k*x)/self.m
    
    def Runge_Kutta(self):
        h =  self.passos
        qnt = int(self.t/self.passos) + 1
        x = [0]*qnt
        x_linha = [0]*qnt
        x_duas_linhas = [0]*qnt
        tempoAtual = [0]*qnt
        
        for k in range(qnt):
            if k > 0:
                tempoAtual[k] = k*h

                k1 = 0.5*h*self.f(tempoAtual[k-1], x[k-1],x_linha[k-1])
                Q = 0.5*h*(x_linha[k-1] + 0.5*k1)
                
                k2 = 0.5*h*self.f((tempoAtual[k-1] + 0.5*h ), (x[k-1]+ Q),(x_linha[k-1]+k1))
                
                k3 = 0.5*h*self.f((tempoAtual[k-1] + 0.5*h ), (x[k-1]+ Q),(x_linha[k-1]+k2))
                L = h*(x_linha[k-1] + k3)
                
                k4 = 0.5*h*self.f((tempoAtual[k-1] + h ), (x[k-1]+ L),(x_linha[k-1]+2*k3))
                
                x[k] = x[k-1] + h*(x_linha[k-1] + (k1+k2+k3)/3)
                x_linha[k] = x_linha[k-1] + (k1+2*k2+2*k3+k4)/3
                x_duas_linhas[k] = self.f(tempoAtual[k], x[k],x_linha[k])
            else:
                x[k]=0
                x_linha[k]=0 
                x_duas_linhas[k] = self.f(tempoAtual[k], x[k],x_linha[k])
        return  {"x":x, "x'": x_linha, "x''": x_duas_linhas, "tempo":tempoAtual}
    
    def output(self):
        ys = self.Runge_Kutta()
        y= ys['x']
        y_linha = ys["x'"]
        y_duas_linha = ys["x''"]
        tempo = ys["tempo"]
        plt.plot(tempo,y)
        plt.show()
        df = pd.DataFrame(list(zip(tempo,y,y_linha,y_duas_linha)), columns = ["tempo","deslocamento",'velocidade','aceleracao'])
        df.to_csv('output.csv', index=None)


passo = float(input('Qual é o passo de integração?'))
t = float(input('Qual o tempo total?'))
m = float(input('m?'))
c = float(input('c?'))
k = float(input('k?'))
a1 = float(input('a1?'))
a2= float(input('a2?'))
a3 = float(input('a3?'))
w1 = float(input('w1?'))
w2 = float(input('w2?'))
w3 = float(input('w3?'))
teste = algcom(passo,t,m,c,k,a1,a2,a3,w1,w2,w3) 

teste.output()
