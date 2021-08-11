# Day 2 Homework

Put your homework for day 2 here!

import matplotlib.pyplot as plt

def calculate_LJ(r_ij):
    r6_term = math.pow(1/r_ij, 6)
    r12_term = math.pow(r6_term, 2)
    
    pairwise_energy = 4 * (r12_term - r6_term)
    return pairwise_energy

for i in range(1, 51):
    r = i * 0.1
    E = calculate_LJ(r)
    #print(E)
    if E> 0:
        plt.scatter(r, math.log(E), c='red')
    else:
        plt.scatter(r, E, c='blue')
plt.xlabel('r')
plt.ylabel('E_lj/InE')
plt.show()

E_list = []
for i in range(1, 51):
    r = i * 0.1
    E = calculate_LJ(r)
    print(E)
    E_list.append(E)
