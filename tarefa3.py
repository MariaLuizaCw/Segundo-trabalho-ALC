import numpy as np
import pandas as pd
import matplotlib as plt
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
        
        return self.a1*np.sin((np.pi*self.w1*t)/180) + self.a2*np.sin((np.pi*self.w2*t)/180) +self.a3*np.cos((np.pi*self.w3*t)/180)
        
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
                print("tempo Atual:")
                print(tempoAtual[k])
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
        plt.plot(tempo,y)
        y_linha = ys["x'"]
        y_duas_linha = ys["x''"]
        tempo = ys["tempo"]
        df = pd.DataFrame(list(zip(tempo,y,y_linha,y_duas_linha)), columns = ["tempo","deslocamento",'velocidade','aceleracao'])
        with open("output.txt", "w") as arquivo:
            arquivo.write(f"Dados lidos:\nPasso de integracao: {self.passos}\nTempo total de integracao: {self.t} segundos\nm:{self.m}\nc: {self.c}\nk: {self.k}\na1: {self.a1}\na2: {self.a2}\na3: {self.a3}\nw1: {self.w1}\nw2: {self.w2}\nw3: {self.w3}\n\nTabela:\n {df}\n")


#icod = int(input('Qual o código da operação?'))
#theta1 = float(input('Qual theta 1?'))
#theta2 = float(input('Qual theta 2?'))
#tolm = float(input('Qual tolerancia?'))
#self,passos ,t, m, c,k,a1,a2,a3,w1,w2,w3
teste = algcom(10,60,1,0.1,2,1,2,1.5,0.05,1,2) 
teste.output()
#print(teste.Runge_Kutta())