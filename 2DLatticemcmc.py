# Code by Maria Davila with help from Dr. Rogers

#Prints out density profiles for lattices:
    # Constant number of particles
    # Increasing in y while x remains constant
    # Two solid wall

import numpy as np
rand = np.random # use numpy's random number module

import matplotlib.pyplot as plt

# Note: modifies ind in-place!
def lattice(N,r,c,ind):
    lat = np.zeros((r,c), np.int) # make these integers
    for i in range(c):
        lat[0,i] = 1
        lat[(r-1),i] = 1
    
    while np.sum(lat) < N:
        i = rand.randint(0, r)
        j = rand.randint(0, c)
        if lat[i,j] == 1:
            continue
        else:
            lat[i,j] = 1
            ind.append((i,j))
    ind.sort()
    return lat

def count_nbrs(r,c,lat):
    x,y = lat.shape
    nbrs = lat[r-1,c]*(r!=0) + lat[(r+1)%x,c]*(r!=x-1) + lat[r,c-1] + lat[r,(c+1)%y]
    return nbrs

def energycalc(ind, lat):
    nbrs = 0
    for i,j in ind:
        nbrs += count_nbrs(i,j,lat)
    return nbrs * 0.5 * beta

# Note: modifies ind, lat, and E in-place.
# returns 1 if successful and 0 otherwise (for accounting purposes)
def swap(ind,lat,E):
    n    = rand.randint(0, len(ind)) # choose existing particle
    i, j = ind[n] # position of existing particle

    while True: # find a clear location
        r = rand.randint(0, lat.shape[0]) # choose new location
        c = rand.randint(0, lat.shape[1])
        if lat[r,c] == 0: # clear.
            break
        
    old_nbrs = count_nbrs(i,j,lat)
    lat[i,j] = 0 # temporarily remove particle
    new_nbrs = count_nbrs(r,c,lat)

    #print (i,j,old_nbrs), "->", (r,c,new_nbrs)
    dE = beta * (new_nbrs - old_nbrs)
    if np.exp(-dE) > rand.random(): # accept with min(1, exp(-beta dE))
        ind[n] = (r,c)
        lat[r,c] = 1
        E[0] += dE
        return 1

    lat[i,j] = 1 # put particle back
    return 0

def mcmh(ind, lat, steps):
    E = [ energycalc(ind, lat) ]
    stat = np.zeros(lat.shape) # make these floats
    acc = 0
    for i in range(steps):
        for j in range(5): # skip
            acc += swap(ind, lat, E)
        stat += lat
    print ("%d / %d swaps accepted."%(acc, steps*5))
    return stat / float(steps)

p_row =[]
length = []

for h in range(8):
    r,c= 5+h,5
    particles = 10
    immovable = 2*c
    N = particles + immovable
    beta = -0.5
    ind = []
    lat = lattice(N,r,c,ind)
    stat = mcmh(ind, lat, 1000)
    print(stat)
    p_row.append(sum(stat[-2]/5))
    length.append(r-2)