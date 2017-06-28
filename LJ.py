#!/usr/bin/env python

import numpy as np
import sys
from math import pow
import matplotlib.pyplot as plt
import matplotlib.colors as color
import time
import cPickle as pickle

# usage: LJ.py <pdb file> <mol2 file> <pqr file>
# output: data.pickle of pairwise energy profile matrix and plots of heat map and energy values vs. atom number

atoms= ['C','H','N','O','P','S']
k_standard = 8.9875517873681764e9
angstrom_m = 10e10
electron = 1.60217646e-19
avogadro = 6.022e23
start_time = time.time()

def parsePDB():
	atomlist = list()
	coordlist = np.array([0,0,0])
	pdb = open(sys.argv[1],'r')
	for line in pdb:
		string = line.split()
		if string[0] == 'ATOM':
			atomlist.append(string[-1])
			coordlist = np.vstack((coordlist,np.array([float(string[6]),float(string[7]),float(string[8])])))
	coordlist = np.delete(coordlist,0,0)
	mol2 = open(sys.argv[2],'r')
	bonds = [[] for i in range(len(atomlist))]
	nex = False
	for line in mol2:
		if line[:13] == '@<TRIPOS>BOND':
			nex = True
		elif nex:
			string = line.split()
			bonds[int(string[1])-1].append(int(string[2])-1)
			bonds[int(string[2])-1].append(int(string[1])-1)
	pqr = open(sys.argv[3],'r')
	charges=list()
	atomlist2=list()
	for line in pqr:
		string = line.split()
		if string[0] == 'ATOM':
			charges.append(float(string[-2]))
			atomlist2.append(string[
	return coordlist, atomlist, bonds, charges

def iterMol():
	energyVals = list()
	coordlist, atomlist, bonds, charges = parsePDB()
	A, B = computeConstants()
	pairwise_mtx = np.zeros((len(atomlist),len(atomlist)))
	mu = float()
	m = int()
	totalLJ = float()
	totalCoul = float()
	for i in range(len(atomlist)):
		counter = float()	
		for j in [x for x in range(len(atomlist)) if x != i]:
			if (j in bonds[i]) or (i in bonds[j]):
				pass
			else:
				temp,Coul = atom_atomPotential(atomlist[i],atomlist[j],coordlist[i,:],coordlist[j,:],A,B,charges[i],charges[j])
				totalLJ += temp
				totalCoul += Coul
				mu += temp
				m += 1
				pairwise_mtx[i,j] = temp
				counter += temp
		energyVals.append(counter)
	mu = mu/m
	print('LJ=='+str(totalLJ)+' Coul=='+str(totalCoul))
	var = float()
	for i in range(len(atomlist)):
		for j in range(len(atomlist)):
			var += pow((pairwise_mtx[i,j]-mu),2)
	var = var/m
	std = pow(var,0.5)
	return pairwise_mtx,energyVals,mu,std

def computeConstants():
	length = len(atoms)
	r_eqm_X = np.array([2.00, 1.00, 1.75, 1.60, 2.10, 2.00])
	eps_X = np.array([0.15, 0.02,0.16,0.20,0.20,0.20])
	r_eqm_X = np.repeat(r_eqm_X,length).reshape(length,length).swapaxes(0,1)
	eps_X = np.repeat(eps_X,length).reshape(length,length).swapaxes(0,1)
	eps_XY = np.sqrt(eps_X*eps_X.T)
	r_eqm_XY = (r_eqm_X + r_eqm_X.T)/2.0
	A = eps_XY*(r_eqm_XY**12)
	B = 2 * eps_XY * (r_eqm_XY**6)
	return A, B

def atom_atomPotential(atom1,atom2,coord1,coord2,A,B,charge1,charge2):
	for num, val in enumerate(atoms):
		if val == atom1:
			atom1 = num
		if val == atom2:
			atom2 = num
	dist = coord1-coord2
	dist = np.dot(dist,dist)  #distance squared
	rootdist = pow(dist,0.5)  #distance
	k = k_standard*angstrom_m*electron**2 #adjust Coulomb constant for units
	LJ = A[atom1,atom2]/pow(dist,6) - B[atom1,atom2]/pow(dist,3)
	Coulomb = k*(charge1*charge2)/rootdist
	Coulomb = (Coulomb/1000)*avogadro
	if Coulomb > 350:
		print(Coulomb,rootdist,charge1,charge2)
	#total = LJ + Coulomb
	return LJ,Coulomb

def plot():
	pairwise_mtx, energyVals, mu, std = iterMol()
	pickle.dump((pairwise_mtx,mu,std),open('data.pickle','wb'))
	pairwise_mtx,mu,std=pickle.load(open('data.pickle','rb'))
	print("Program ran for %s seconds" %(time.time() - start_time))
	f = plt.figure(1)
	plt.plot(list(range(1,201)),energyVals[:200])
	f.show()
	g = plt.figure(2)
	plt.imshow(np.log(-pairwise_mtx+1e-16),cmap=plt.get_cmap('hot'))
	plt.colorbar()
	plt.xlabel('Atom number')
	plt.ylabel('Atom number')
	plt.title('Lennard Jones Pairwise Potential for PDB:1AB4')
	g.show()
	plt.show()
	h = plt.figure(3)
	plt.plot(list(range(len(energyVals))),energyVals)
	h.show()
	raw_input()
	
if __name__ == '__main__':
	plot()
