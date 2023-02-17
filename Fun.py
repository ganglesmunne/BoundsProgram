import numpy as np
import math
import itertools
import picos as pic
import cvxopt as cvx

def K_fun(m,k,n):

	Fun=0

	for q in range(0,m+1):
		Fun+= math.comb(n-k,m-q)*math.comb(k,q)*(-1)**q
	return Fun 

def comb_g(n,k):

	if n>=0 and k>=0:
		value=math.comb(n,k)
	else:
		value=0

	return value

def Beta(i,j,k,t,n):

	beta=0

	for u in range(n+1):

		beta+= comb_g(u,t)*comb_g(n-2*k,u-k)*comb_g(n-k-u,i-u)*comb_g(n-k-u,j-u)*(-1)**(u-t)

	return beta

def BetaMatrices(Ap,n):

	BMatrix=[]

	for k in range(0,int(n/2)+1):

		Matrix=[[0]*(n-2*k+1) for _ in range(n-2*k+1)]

		for i in range(k,n-k+1):
			for j in range(k,n-k+1):
				for t in range(n+1):

					l=i+j-2*t

					if l>n or t>i or t>j:

						Matrix[i-k][j-k]+=0

					else:

						Matrix[i-k][j-k]+=(1/math.comb(n,l))*Ap[l]*Beta(i,j,k,t,n)*(math.comb(n-2*k,i-k)*math.comb(n-2*k,j-k))**(-0.5)

		BMatrix.append(Matrix)

	return BMatrix


def GammaTilde(Ap,n):
	
	Generators=list(range(0,n))
	
	Elements=[]
	
	for i in range(n+1):
		Elements += (list(itertools.combinations(Generators, i)))
	
	Elements_R={new_label:old_label for old_label, new_label in enumerate(Elements)}

	GTilde=[[0]*(2**n) for _ in range(2**n)]
	
	for S in Elements:
		for T in Elements:

			i=len(S)
			j=len(T)
			t=len(tuple(set(S) & set(T)))
			
			l=i+j-2*t
			
			GTilde[Elements_R[S]][Elements_R[T]]=(1/math.comb(n,l))*Ap[l]

	return GTilde


