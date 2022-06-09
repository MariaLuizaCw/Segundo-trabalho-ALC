from tkinter import N
import numpy as np

class algcom:
    def __init__(self, icod, theta1, theta2, tolm):
        self.icod = icod
        self.theta1 = theta1
        self.theta2 = theta2
        self.tolm = tolm

    def Jacobiana(x):
        # c2 = x[0], c3 = c3, c4 = x[2]
        c2 = x[0]
        c3 = x[1]
        c4 = x[2]
        deltaf1 = [2*c2,               4*c3,        12*c4]
        deltaf2 = [12*c2*c3 + 36*c3*c4, 24*(c3**2) + 6*c2**2 + 36*c2*c4 + 108*c4**2, 36*c2*c3 + 216*c3*c4]
    def Newton():

