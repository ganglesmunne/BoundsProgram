import Fun as f
import numpy as np
from termcolor import colored
import itertools
from sigfig import round

#### Check Beta coeficient ####
n=31
k=0

"""
for i,j,t in itertools.product(range(n+1),range(n+1),range(n+1)):
    if round(f.Beta(i,j,k,t,n),sigfigs=8) !=  round(f.Beta(j,i,k,t,n),sigfigs=8):
        print(i,j)
        print(f.Beta(i,j,k,t,n))
        print(f.Beta(j,i,k,t,n))
        print("\n")
"""

i=15 
j=22

for t in range(n+1):
    if f.Beta(i,j,k,t,n) !=  f.Beta(j,i,k,t,n):
        print(f.Beta(i,j,k,t,n))
        print(f.Beta(j,i,k,t,n))
        print("\n")