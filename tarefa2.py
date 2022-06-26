import numpy as np
from gauss import tabela_gauss
from polinomial import polynomial_table


class algcom:
    def __init__(self, icod, cs, TOLm):
        self.icod = icod
        self.c1, self.c2, self.c3, self.c4 = cs
        self.result = ""
        self.tol = TOLm

    def f(self, x):
        return self.c1*np.e**(self.c2*x) + self.c3*x**(self.c4)
    def df(self, x):
        return self.c2*self.c1*np.e**(self.c2*x) + self.c4*self.c3*x**(self.c4-1)
    def bissecao(self, a, b, maxiter=1000):
        itercount = 0
        if self.f(a)*self.f(b) > 0 :
            self.output("Choose a different interval")
            return 
        while True:
            media  = (a + b)/2
            fb, fmedia = self.f(b), self.f(media)
            prod = fb * fmedia
            itercount += 1
            if abs(b-a) <= self.tol:
                self.output(f"Raiz: {media}")
                return
            elif maxiter == itercount:
                self.output("Method did not converge")
                return
            elif prod > 0:
                b = media
            elif prod < 0:
                a = media
        
        
    def Newton(self, x0, maxiter=1000): 
        while maxiter:
            x1 = x0 - self.f(x0)/self.df(x0)
            if abs(x1 - x0) <= self.tol:
                self.output(f"Raiz: {x1}")
                return 
            x0 = x1
            maxiter -= 1
        self.output("Method did not converge")
        return 

    def quadratura_gauss(self, n, a, b):
        L = b - a
        ws = tabela_gauss[n-2]["ws"]
        zs = tabela_gauss[n-2]["zs"]
        xs = [(a + b + L*zs[i])/2 for i in range(n)]
        res = L*sum([self.f(xs[i])*ws[i] for i in range(n)])/2
        self.output(f"A integral vale: {res}")
    def quadratura_polinomial(self, n, a, b):
        pol_table = polynomial_table(a, b, n)
        ws, zs = pol_table[0][n], pol_table[1][n]
        res = sum([self.f(zs[i])*ws[i] for i in range(n)])
        self.output(f"A integral vale: {res}")
            
            
    def derivada_frente(self, x, deltax):
        return (self.f(x + deltax) - self.f(x))/deltax
    def derivada_tras(self, x, deltax):
        return (self.f(x) - self.f(x - deltax))/deltax
    def derivada_central(self, x, deltax):
        return (self.f(x + deltax) - self.f(x - deltax))/(2*deltax)
    def derivar(self, qs, x, deltax):
        res = "O resultado da(s) derivada(s) foi: "
        if qs[0] == 's':
            res += f"{self.derivada_frente(x, deltax)}, "
        
        if qs[1] == 's':
            res += f"{self.derivada_tras(x, deltax)}, "
        
        if qs[2] == 's':
            res += f"{self.derivada_central(x, deltax)}, "
        self.output(res)
    def richard(self, x, deltax1, deltax2):
        p = 1
        q = deltax1/deltax2
        d1 = self.derivada_frente(x, deltax1)
        d2 = self.derivada_frente(x, deltax2)
        res = d1 + (d1-d2)/(q**(-p) - 1)
        teste.output(f"Derivada: {res}")


    def output(self, result):
        with open("output.txt", "w") as arquivo:
            arquivo.write(result)
                


icod = int(input('Qual o código da operação?'))
c1 = float(input('Qual o valor de c1?'))
c2 = float(input('Qual o  valor de c2?'))
c3 = float(input('Qual o valor de c3?'))
c4 = float(input('Qual o valor de c4?'))
teste = algcom(icod, [c1, c2, c3, c4], 0)
if icod == 1:
    tol = float(input('Qual tolerancia?'))
    teste.tol = tol
    met = input('Qual o método (n = newton, b = bissecao)?')
    if met == "b":
        a = float(input('Qual o valor de a?'))
        b = float(input('Qual o valor de b?'))
        teste.bissecao(a, b)
    elif met == "n":
        x0 = float(input('Qual o valor de inicial de x?'))
        teste.Newton(x0)
elif icod == 2:
    met = input('Qual o método (g = gauss, p = polinomial)? ')
    a = float(input('Qual o valor de a? '))
    b = float(input('Qual o valor de b? '))
    n = int(input('Quantos pontos de integração? '))
    if met == "g":
        teste.quadratura_gauss(n, a, b)
    if met == "p":
        teste.quadratura_polinomial(n, a, b)
elif icod == 3:
     q1 = input('Deseja derivar pra frente? (s/n) ')
     q2 = input('Deseja derivar pra trás? (s/n) ')
     q3 = input('Deseja derivada central? (s/n) ')
     x = float(input('Qual o valor de x? '))
     deltax = float(input('Qual o valor de deltax? '))

     teste.derivar(q1+q2+q3, x, deltax)
elif icod == 4:
    x = float(input('Qual o valor de x? '))
    deltax1 = float(input('Qual o primeiro deltax? '))
    deltax2 = float(input('Qual o segundo deltax? '))
    teste.richard(x, deltax1, deltax2)

   