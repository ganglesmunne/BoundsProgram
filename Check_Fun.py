import Fun as f
import numpy as np
from termcolor import colored
import itertools


#############################################
#                                           #
####   Check Identity maps to Identity 	 ####
#                                           #
#############################################


#### Check Gamma tilde is the idenity ####

n=9

Ap=[0]*(n+1)

Ap[0]=1  # All the Ap are 0 except Ap[0] to have the idenity in the Gamma tilde matrix 

M=f.GammaTilde(Ap,n)

Id=np.eye(2**n)

if np.allclose(Id,M)==True:
	
	Id_gamma=1 # If the gamma tilde is the identity

else:
	
	Id_gamma=0 # If the gamma tilde is not the identity

#### Check Beta matrices are the idenity ####

M=f.BetaMatrices(Ap,n)

for k in range(0,int(n/2)+1):
	
	Id=np.eye(n-2*k+1)
	 
	if np.allclose(Id,M[k])==False:
		
		Id_beta=0 # Break if one of the beta matrices is not the identity
		break
		
	else:
		
		Id_beta=1

if Id_beta==1 and Id_gamma==1:
	
	print("Id_gamma <-> Id_beta:", colored('True', 'green', attrs=['bold'])) 
	
else:
	
	print("Id_gamma <-> Id_beta: ", colored('False', 'red', attrs=['bold']))



#############################################
#                                           #
### Minum eigenvalue of the Beta matrices ###
#                                           #
#############################################


n=9

#### Check minimum eigenvalue ####

Ap=np.random.rand(n+1) # Random selection of Ap

### Minimum eigenvalue of Gamma tilde ###

M=f.GammaTilde(Ap,n)

Diag=np.linalg.eigvals(M)

Diag.sort()

MinEig_Gamma=Diag[0] 

M=f.BetaMatrices(Ap,n)

MinEig_Beta=np.inf

for i in range(0,int(n/2)+1):
	
	Diag=np.linalg.eigvals(M[i])
	Diag.sort()
	
	MinEig=Diag[0]
	
	if MinEig_Beta >= MinEig:
		
		MinEig_Beta=MinEig 

if np.round(MinEig_Beta,9) == np.round(MinEig_Gamma,9):
	
	print("Min eig gamma = Min eig beta: ", colored('True', 'green', attrs=['bold']))
	
else:
	
	print("Min eig gamma = Min eig beta: ", colored('False', 'red', attrs=['bold']))



#############################################
#                                           #
###      Check if Betas are symetric      ###
#                                           #
#############################################


n=31
k=0

s=0

for i,j,t in itertools.product(range(n+1),range(n+1),range(n+1)):
	if f.Beta(i,j,k,t,n) !=  f.Beta(j,i,k,t,n):
		s=1
		break

if s==0:

	print("Symetric Beta: ", colored('True', 'green', attrs=['bold']))

else:

	print("Symetric Beta: ", colored('False', 'red', attrs=['bold']))


#############################################
#                                           #
###   Check mapping Gamma -> Gamma_tilde  ###
#                                           #
#############################################

n=6

AS=np.random.rand(2**n)

Gamma,Gamma_tilde,Aj=f.Gamma(AS,n,True,True)

Gamma_tilde2=np.array(f.GammaTilde(Aj,n))

if np.allclose(Gamma_tilde,Gamma_tilde2)==True:

	print("Same Gamma_tilde: ", colored('True', 'green', attrs=['bold']))

else:

	print("Same Gamma_tilde: ", colored('False', 'red', attrs=['bold']))


