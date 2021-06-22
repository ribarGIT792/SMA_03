# Pvz_3_06_Lygciu sistemos sprendimas taikant LU skaida

import numpy as np
from PyFunkcijos import *

# textbox isvedimui
T=ScrollTextBox(100,20) # sukurti teksto isvedimo langa
T.insert(END,"# ********************* LU algoritmas  ****************************")

# -------------------iseities duomnenys:
A=np.matrix([[1 , 1,  1,  1],
             [1, -1, -1,  1],
             [2,  1, -1,  2],
             [3,  1,  2, -1]]).astype(np.float)        # koeficientu matrica

#b=np.matrix([[2],[0],[9],[7]]).astype(np.float)
b=np.matrix([[2,1],[0,1],[9,3],[7,1]]).astype(np.float)
As=np.matrix(A);bs=np.matrix(b);  # koeficientru matrica ir laisvuju nariu vektorius issaugomi sprendinio patikrinimui
n=(np.shape(A))[0]   # lygciu skaicius nustatomas pagal ivesta matrica A
nb=(np.shape(b))[1]  # laisvuju nariu vektoriu skaicius nustatomas pagal ivesta matrica b
P=np.arange(0,n)
SpausdintiMatrica(A,T,"A");SpausdintiMatrica(b,T,"b");SpausdintiMatrica(n,T,"n");SpausdintiMatrica(nb,T,"nb")
SpausdintiMatrica(P,T,"P");

# tiesioginis etapas: 

for i in range (0,n-1):   # range pradeda 0 ir baigia n-2 (!)
    a=max(abs(A[i:n,i]));  iii=abs(A[i:n,i]).argmax()
    A[[i,i+iii],:]=A[[i+iii,i],:]  # sukeiciamos eilutes 
    P[[i,i+iii]]=P[[i+iii,i]]       # sukeiciami eiluciu numeriai 
    for j in range (i+1,n):  # range pradeda i+1 ir baigia n-1
        r=A[j,i]/A[i,i]
        A[j,i:n+nb]=A[j,i:n+nb]-A[i,i:n+nb]*r;
        A[j,i]=r;   

b=b[P,:]
SpausdintiMatrica(b,T,"b_pertvarkytas")

# 1-as atvirkstinis etapas, sprendziama  Ly=b, y->b
for i in range(1,n) :  
    b[i,:]=b[i,:]-A[i,0:i]*b[0:i,:]
# 2-as atvirkstinis etapas , sprendziama Ux=b, x->b
for i in range (n-1,-1,-1) : 
    b[i,:]=(b[i,:]-A[i,i+1:n]*b[i+1:n,:])/A[i,i]
    

T.insert(END,"sprendinys:");SpausdintiMatrica(b,T,'b');
T.insert(END,"------------ sprendinio patikrinimas ----------------");
liekana=As.dot(b)-bs;SpausdintiMatrica(liekana,T,'liekana');
SpausdintiMatrica(np.linalg.norm(liekana)/ np.linalg.norm(bs),T,"bendra santykine paklaida:")

str1=input("Baigti darba? Press Enter \n")