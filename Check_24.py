import Fun as f
import numpy as np
from termcolor import colored
import itertools

n=24
D=2

A=[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18216.0, 0.0, 156492.0, 0.0, 1147608.0, 0.0, 3736557.0, 0.0, 6248088.0, 0.0, 4399164.0, 0.0, 1038312.0, 0.0, 32778.0]
for m in range(n+1):
	Ap.append(sum(D**(-m)*f.comb(n-j,n-m)*A[j] for j in range(m+1)))


#A=[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18216.0, -0.0, 156491.98, 0.11, 1147607.75, 0.25, 3736557.01, -0.28, 6248088.22, 0.1, 4399163.69, 0.26, 1038311.89, 0.03, 32778.0]
#A=[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18216.07, -0.18, 156490.54, 8.69, 1147589.06, 15.62, 3736568.76, -40.5, 6248125.67, -7.92, 4399148.96, 16.14, 1038304.47, 1.79, 32777.82]

Ap=[]

for m in range(n+1):
	Ap.append(sum(D**(-m)*f.comb(n-j,n-m)*A[j] for j in range(m+1)))

M=f.BetaMatrices(Ap,n)

MinEig_Beta=np.inf

for i in range(0,int(n/2)+1):
	
	Diag=np.linalg.eigvals(M[i])
	Diag.sort()
	
	MinEig=Diag[0]
	
	if MinEig_Beta >= MinEig:
		
		MinEig_Beta=MinEig 
		

print(MinEig_Beta)


for i in range(0,int(n/2)+1):
	
	Diag=np.linalg.eigvals(M[i])
	print(Diag)
