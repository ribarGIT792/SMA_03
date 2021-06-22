

# Choleckio skaidos algoritmo taikymas lygciu sistemos sprendimui

import numpy as np
from PyFunkcijos import *

# textbox isvedimui
T=ScrollTextBox(100,20) # sukurti teksto isvedimo langa
T.insert(END,"# ********************* Choleckio algoritmas  ****************************")

A=np.matrix([[ 4,  3, -1,  1],
             [ 3,  9, -2, -2],
             [-1, -2, 11, -1],
             [ 1, -2, -1,  5]]).astype(np.float)
b =np.matrix([[12],[10],[-28],[16]]).astype(np.float)
n=(np.shape(A))[0]
SpausdintiMatrica(A,T,'A');SpausdintiMatrica(b,T,'b');SpausdintiMatrica(n,T,'n');
Aprad=np.matrix(A);bprad=np.matrix(b);

# Choleckio L'*L skaida
for i in range (0,n) :
    A[i,i]=sqrt(A[i,i]-np.sum(np.square(A[0:i,i])))
    for j in range (i+1,n) :
        A[i,j]=(A[i,j]-A[0:i,i].transpose()*A[0:i,j])/A[i,i];
    
SpausdintiMatrica(A,T,'A');

# 1-as atgalinis zingsnis, sprendziama  L'y=b, y->b
for i in range (0,n) :
    b[i,:]=(b[i,:]-A[0:i,i].transpose()*b[0:i])/A[i,i]

SpausdintiMatrica(b,T,'b');

#% 2-as atgalinis zingsnis   Ux=b,  x->b
for i in range (n-1,-1,-1) :
    b[i]=(b[i]-A[i,i+1:n]*b[i+1:n])/A[i,i]

SpausdintiMatrica(b,T,'sprendinys');
SpausdintiMatrica(Aprad*b-bprad,T,'liekana');

str1=input("Baigti darba? Press Enter \n")