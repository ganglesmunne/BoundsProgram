import Fun as f
import numpy as np

#### Check Identity maps to Identity ####

n=5

Ap=[0]*(n+1)

Ap[0]=1


M=f.GammaTilde(Ap,n)
Id=np.eye(2**n)

print(np.allclose(Id,M))


M=f.BetaMatrices(Ap,n)

for k in range(0,int(n/2)+1):
	
	Id=np.eye(n-2*k+1)
	
	print(np.allclose(Id,M[k]))


#### Check minimum eigenvalue ####


Ap=np.random.rand(n+1)


M=f.GammaTilde(Ap,n)

Diag=np.linalg.eigvals(M)

Diag.sort
print(Diag)
