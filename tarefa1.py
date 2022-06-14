import numpy as np

class algcom:
    def __init__(self, icod, theta1, theta2, tolm):
        self.icod = icod
        self.theta1 = theta1
        self.theta2 = theta2
        self.tolm = tolm
        self.erro = 0

    def Jacobiana(self,x):
        # c2 = x[0], c3 = c3, c4 = x[2]
        c2 = x[0]
        c3 = x[1]
        c4 = x[2]
        deltaf1 = [2*c2,               4*c3,        12*c4]
        deltaf2 = [12*c2*c3 + 36*c3*c4, 24*(c3**2) + 6*c2**2 + 36*c2*c4 + 108*c4**2, 36*c2*c3 + 216*c3*c4]
        deltaf3 = [120*c2*c3**2 + 576*c4*c3**2 + 504*c2*c4**2 + 1296*c4**3 + 72*c4*c2**2 + 3, 240*c3**3 + 120*c3*c2**2 + 1152*c3*c4*c2 + 4464*c3*c4**2, 576*c2*c3**2 + 4464*c4*c3**2 + 504*c4*c2**2 + 3888*c2*c4**2 + 13392*c4**3 + 24*c2**3]
        return np.array([deltaf1,deltaf2,deltaf3])
    
    def Funcao(self, x):
        c2 = x[0]
        c3 = x[1]
        c4 = x[2] 
        f1 = 2*c3**2 + c2**2 + 6*c4**2 - 1
        f2 = 8*c3**3 + 6*c3*c2**2 + 36*c3*c2*c4 + 108*c3*c4**2 - self.theta1
        f3 = 60*c3**4 + 60*(c3**2)*c2**2+ 576*c2*c4*(c3**2) + 2232*(c3**2)*(c4**2)+ 252*(c4**2)*(c2**2) + 1296*c2*c4**3 + 3348*c4**4 + 24*c4*c2**3 +3*c2 - self.theta2
        return np.array([f1,f2,f3])
    
    def Newton(self):
        x = [1,0,0]
        x = np.array(x)
        maxiter = 1000
        while maxiter:
            J = self.Jacobiana(x)
            F = self.Funcao(x)
            try:
                deltax = - np.linalg.inv(J)@F
            except:
                self.erro = "The matrix J generated was a Sigunlar matrix"
                return 
            x = x + deltax
            tolk = np.linalg.norm(deltax)/ np.linalg.norm(x)
            maxiter -= 1
            if tolk < self.tolm:
                return x
        self.erro = "convergence not reached"
        return  
    def Broyden(self):
        x = [1,0,0]
        x = np.array(x)
        B = self.Jacobiana(x)
        maxiter = 1000
        while maxiter:
            J = B.copy()
            F = self.Funcao(x)
            try:
                deltax = - np.linalg.inv(J)@F
            except:
                self.erro = "The matrix B generated was a Sigunlar matrix"
                return
            x = x + deltax
            Y = self.Funcao(x) - F
            tolk = np.linalg.norm(deltax)/ np.linalg.norm(x)
            maxiter -= 1
            if tolk < self.tolm:
                return x
            else:
                B = B + ((Y-B@deltax)@(deltax.T))/((deltax.T)@deltax)
        self.erro = "convergence not reached"
        return  

    def output(self):
        if self.icod == 1:
            x = self.Newton()
        elif self.icod == 2:
            x = self.Broyden()
        with open("output.txt", "w") as arquivo:
            if self.erro == 0:
                arquivo.write(f"X: {x}\n")
            else:             
                arquivo.write(f"Erro: {self.erro}\n")

          


icod = int(input('Qual o código da operação?'))
theta1 = float(input('Qual theta 1?'))
theta2 = float(input('Qual theta 2?'))
tolm = float(input('Qual tolerancia?'))


    
teste = algcom(icod, theta1, theta2, tolm) 
teste.output()