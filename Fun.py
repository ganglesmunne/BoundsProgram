import numpy as np
import itertools
import picos as pic
import cvxopt as cvx
from scipy.special import comb as comb_s


def Swap_multiplication(A,B):

	C = tuple(set(A)^set(B))

	return C 
        
def K_fun(m,k,n):

	Fun=0

	for q in range(0,m+1):
		Fun+= comb(n-k,m-q)*comb(k,q)*(-1)**q
	return Fun 

def comb(n,k):

	value=comb_s(n,k,exact=True)

	return value


def Beta(i,j,k,t,n):

	beta=0

	for u in range(n+1):

		beta+= comb(u,t)*comb(n-2*k,u-k)*comb(n-k-u,i-u)*comb(n-k-u,j-u)*(-1)**(u-t)

	return beta

def BetaMatrices(Ap,n):

	BMatrix=[]

	for k in range(0,int(n/2)+1):

		Matrix=[[0]*(n-2*k+1) for _ in range(n-2*k+1)]
		
		for i, j, t in itertools.product(range(k,n-k+1), range(k,n-k+1), range(n+1)):

			l=i+j-2*t

			if l>n or t>i or t>j:

				Matrix[i-k][j-k]+=0

			else:
				const=(1/comb(n,l))*Beta(i,j,k,t,n)*(comb(n-2*k,i-k)*comb(n-2*k,j-k))**(-1/2)
				Matrix[i-k][j-k]+=const*Ap[l]

		BMatrix.append(Matrix)

	return BMatrix


def GammaTilde(Ap,n):
	
	Generators=list(range(0,n))
	
	Elements=[]
	
	for i in range(n+1):
		Elements += (list(itertools.combinations(Generators, i)))
	
	Elements_R={new_label:old_label for old_label, new_label in enumerate(Elements)}

	GTilde=[[0]*(2**n) for _ in range(2**n)]
	
	for S,T in itertools.product(Elements,Elements):
		
		i=len(S)
		j=len(T)
		t=len(tuple(set(S) & set(T)))
		
		l=i+j-2*t
		
		GTilde[Elements_R[S]][Elements_R[T]]=(1/comb(n,l))*Ap[l]

	return GTilde


def Gamma(AS,n,aver,weigths):
	
	T=2**n
		
	F=np.zeros((T, T))

	F[0]=AS

	Generators=list(range(0, n))

	Elements=[]

	for i in range(n+1):
		Elements += (list(itertools.combinations(Generators, i)))
	   
	Elements_R={new_label:old_label for old_label, new_label in enumerate(Elements)}

	rho=[]

	length=0
	
	Rep=[]

	for i in range(T):
		for j in range(T):
			C=Swap_multiplication(Elements[i],Elements[j])
			Rep.append([C,i,j])  

	Rep=sorted(Rep)

	for rep in Rep:
		if rep==Rep[0]:  

			rep_old=rep.copy()

		else:

			if set(rep[0])==set(rep_old[0]):
				
				F[rep[1],rep[2]]=F[rep_old[1],rep_old[2]]
			else:
				rep_old=rep.copy()
		
	if aver==True:
		
		F_Tilde=np.zeros([2**n,2**n])
		
		for X in Elements:
			for Y in Elements:

				i=len(X)
				j=len(Y)
				t=len(tuple(set(X) & set(Y)))
				l=i+j-2*t
				
				suma=0

				for S in Elements:
					for T in Elements:

						Inter=tuple(set(S) & set(T))

						if len(S)==i and len(T)==j and len(Inter)==t:

							suma+=F[Elements_R[S],Elements_R[T]]

				F_Tilde[Elements_R[X]][Elements_R[Y]]=1/(comb(n,l)*comb(n-l,t)*comb(l,j-t))*suma	
	else:
		F_Tilde=None

	if weigths==True:
		
		Aj=np.zeros(n+1)
		
		
		Aj[0]=AS[0]
		
		for X in Elements:
			i=len(X)
			
			if i==0: 
				Aj[0]=AS[Elements_R[X]]
			elif i == ip:
				Aj[i]+=AS[Elements_R[X]]
			else:
				Aj[i]=AS[Elements_R[X]]
			ip=i
	else:
		Aj=None
		
	return F,F_Tilde,Aj

	

