#%% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataScience.changeDirOnImportExport setting
# ms-python.python added
import os
try:
	os.chdir(os.path.join(os.getcwd(), '..\\..'))
	print(os.getcwd())
except:
	pass

#%%
import matplotlib
import numpy as np
get_ipython().run_line_magic('matplotlib', 'notebook')
import matplotlib.pyplot as plt
from scipy.spatial import distance
import random


#%%
A=0.000001
B=0.00001
d_min=0.01
d_max=1
imax=20
shift_d = 0.00000001
#initialize the parameters


#%%
def potential(q,p):
    #define the potential between atom q and p
   return A/(distance.euclidean(q,p)**12) - B/(distance.euclidean(p,q)**6)

class Atom:
    #class atom contain a atom's coordinate and force
    def __init__(self,coords):
        self.coords = coords
        self.force = [0, 0]


#%%
class seed():
    
    def __init__(self, size):
        self.size = size
        self.atoms = []
        return
    
    def generate_atoms(self):
        #generate specific atoms randomly with rantional distance between each other
        
        self.atoms = []
        site0 = Atom([0,0])
        self.atoms.append(site0)
        
        while True:
            site = [random.uniform(-5,5),random.uniform(-5,5)] 
            for atom in self.atoms:
                if (distance.euclidean(site,atom.coords) > d_min) and (distance.euclidean(site, atom.coords) < d_max):
                    atom = Atom(site)
                    self.atoms.append(atom)
                    break
            if len(self.atoms) > self.size :
                break
                
                    
    def print(self):
        for atom in self.atoms:
            print(atom.coords)
            
    
    def cal_V(self):
        #calculate the sum of potential of those atoms
        
        V_sum = 0
        for atom in self.atoms:
            for atom1 in self.atoms:
                if atom1.coords != atom.coords:
                    V_sum = potential(atom.coords,atom1.coords) + V_sum
        V_sum = V_sum/2
        return V_sum
    
    def cal_force(self):
        #calculate each atom's force and store it in class Atom
        
        for atom0 in self.atoms:
            fx, fy = 0, 0
            for atom1 in self.atoms:
                if atom1.coords != atom0.coords:
                    d = distance.euclidean(atom0.coords, atom1.coords)
                    f_sum = -12*A/(d**13) + 6*B/(d**7)
                    x = atom1.coords[0] - atom0.coords[0]
                    y = atom1.coords[1] - atom0.coords[1]
                    fx = f_sum*np.cos(x/d) + fx
                    fy = f_sum*np.sin(y/d) + fy
            atom0.force[0] = fx
            atom0.force[1] = fy
            
            
    def rebuild(self):
        
        #shift those atoms to find lower potential , the distance that each atom moves depends on its force plus a tiny coefficient 
        
        for atom in self.atoms:
            atom.coords[0] = atom.force[0]*shift_d + atom.coords[0]
            atom.coords[1] = atom.force[1]*shift_d + atom.coords[1]
            
        return


#%%

test = []
while len(test) < 10:
    a = seed(20)
    a.generate_atoms()
    while True :
        V_sum = a.cal_V()
        a.cal_force()
        a.rebuild()
        V_previous = V_sum
        V_sum = a.cal_V()
        if abs(V_sum - V_previous) < 0.1 :
            for atom in a.atoms:
                if distance.euclidean(atom.coords, [0,0]) > 6:
                    break
            test.append(V_sum)
            break

print(test)


#%%
b = seed(20)
b.generate_atoms()
b_V = b.cal_V()
x1,y1 = [], []
x2,y2 = [], []

for atom in b.atoms:
        atom.x = atom.coords[0]
        x1.append(atom.x)
        atom.y = atom.coords[1]
        y1.append(atom.y)
        
while True :
        V_sum = b.cal_V()
        b.cal_force()
        b.rebuild()
        V_previous = V_sum
        V_sum = b.cal_V()
        if (V_previous - V_sum) < 0.1 :
            break
            
for atom in b.atoms:
        atom.x = atom.coords[0]
        x2.append(atom.x)
        atom.y = atom.coords[1]
        y2.append(atom.y)
            
fig = plt.figure()
plt.scatter(x1,y1,color='black')
plt.scatter(x2,y2,color='red')
plt.show()


#%%
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.set(title="before")
ax2.set(title="after")
ax1.scatter(x1,y1,color='black')
ax2.scatter(x2,y2,color='red')
a_V = b.cal_V()
print(b_V)
print(a_V)


#%%
a_V = b.cal_V()
print(b_V)
print(a_V)


#%%



