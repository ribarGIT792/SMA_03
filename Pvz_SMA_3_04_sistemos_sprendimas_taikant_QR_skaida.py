
import numpy as np
from PyFunkcijos import *
print("*******prdzia*********")
# textbox isvedimui
T=ScrollTextBox(100,20) # sukurti teksto isvedimo langa
T.insert(END,"# ********************* Gauso algoritmas  ****************************")

# -------------------iseities duomnenys:
A=np.matrix([[1 , 1,  1,  1],
             [1, -1, -1,  1],
             [2,  1, -1,  2],
             [3,  1,  2, -1]]).astype(np.float)        # koeficientu matrica
Ap=A   # bus naudojama patikrinimui
b=(np.matrix([2,0,9,7])).transpose().astype(np.float)   #laisvuju nariu vektorius-stulpelis
n=(np.shape(A))[0]   # lygciu skaicius nustatomas pagal ivesta matrica A
nb=(np.shape(b))[1]  # laisvuju nariu vektoriu skaicius nustatomas pagal ivesta matrica b

SpausdintiMatrica(A,T,'A');SpausdintiMatrica(b,T,'b');SpausdintiMatrica(n,T,'n');SpausdintiMatrica(nb,T,'nb');

# tiesioginis etapas(QR skaida):
Q=np.identity(n)
for i in range (0,n-1):
    z=A[i:n,i]
    zp=np.zeros(np.shape(z)); zp[0]=np.linalg.norm(z)
    SpausdintiMatrica(zp,T,'zp')
    omega=z-zp; omega=omega/np.linalg.norm(omega)
    Qi=np.identity(n-i)-2*omega*omega.transpose()
    A[i:n,:]=Qi.dot(A[i:n,:])
    SpausdintiMatrica(A,T,'A')
    Q=Q.dot(
            np.vstack(
             (
                np.hstack((np.identity(i),np.zeros(shape=(i,n-i)))),
                np.hstack((np.zeros(shape=(n-i,i)),Qi))
             )
                 )
           )
    SpausdintiMatrica(Q,T,'Q')
    SpausdintiMatrica(A,T,'A')

    # atgalinis etapas:
b1=Q.transpose().dot(b);
x=zeros(shape=(n,nb));
for i in range (n-1,-1,-1):    # range pradeda n-1 ir baigia 0 (trecias parametras yra zingsnis)
    x[i,:]=(b1[i,:]-A[i,i+1:n]*x[i+1:n,:])/A[i,i];

    SpausdintiMatrica(x,T,"x")

T.insert(END,"sprendinys:");SpausdintiMatrica(x,T,'x');
T.insert(END,"------------ sprendinio patikrinimas ----------------");
liekana=Ap.dot(x)-b1;SpausdintiMatrica(liekana,T,'liekana');
SpausdintiMatrica(np.linalg.norm(liekana)/ np.linalg.norm(x),T,"bendra santykine paklaida:")

str1=input("Baigti darba? Press Enter \n")


