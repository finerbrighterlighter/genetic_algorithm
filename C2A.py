#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 18:05:53 2020

@author: hteza
"""
import numpy as np
import random as rnd
import bisect
import matplotlib.pyplot as plt
np.set_printoptions(suppress=True)

# Two seperate functions are used
# for population generation
# for passing through every stages of the gA process

def population (input,number_of_chrom,min_bound,max_bound):
    # creating population
    gene = len(input)
    no_chrom = number_of_chrom
    pop_size=(no_chrom,gene)
    population=np.random.randint(low=min_bound,high=max_bound,size=pop_size)
    return population

def GA(input, population, min_bound, max_bound, prob_of_co, prob_of_mu):
    # fitness funtion
    fitness = np.sum(population*input, axis=1)
    best_outputs.append(np.max(fitness))
    # ranking by fitness ( here maximum is optimum )
    rank=np.argsort(fitness)[::-1]
    pop_sort=population[rank]
    # building roulette wheel
    # full wheel
    full_range=sum(abs(fitness))
    # each section
    chance=abs(fitness)/full_range
    chance_sort=chance[rank]
    section=np.cumsum(chance_sort)
    # picking partners for every chromosome
    # turning the wheel
    partner=[]
    for i in chance_sort:
        pick=rnd.uniform(0,1)
        pair=bisect.bisect_right(section,pick)
        partner.append(pair)
    partner=np.ravel(partner)
    partner=partner.reshape(-1,1)
    pop_sort_pair=np.append(pop_sort,partner,axis=1)
    # chance of crossing over
    # assigning a random value to every pair
    # if it is less than prob of crossover
    # it is randomized between 0 to 1
    pc = prob_of_co
    parents=np.empty((1,7))
    for i in range(0,pop_sort_pair.shape[0],1):
        chance_for_co = rnd.uniform(0,1)
        if chance_for_co>=pc:
            p=pop_sort_pair[i] 
            p=p.reshape(-1,1)
            p=np.transpose(p)
            parents=np.concatenate((parents,p),axis=0)
    parents=np.delete(parents,(0),axis=0)
    # crossing over
    # point of crossover is randomly decided from 1 to 5
    # two chromsomes are crossed over at that point
    children=np.empty((1,6))
    for i in range(0,parents.shape[0],1):
        p1= parents[i, :-1]
        p2= pop_sort_pair[int(parents[i, -1]), :-1]
        k = rnd.randint(0,5)
        c1=np.append(p1[0:k],p2[k:6])
        c2=np.append(p2[0:k],p1[k:6])
        c1=c1.reshape(-1,1)
        c1=np.transpose(c1)
        c2=c2.reshape(-1,1)
        c2=np.transpose(c2)
        c3=np.concatenate((c1, c2), axis=0)
        children=np.concatenate((children,c3),axis=0)
        i=i+1
    children=np.delete(children,(0),axis=0)
    # old parents and new children
    family = np.concatenate((pop_sort,children),axis=0)
    # mutation
    # assigning a random value to every gene
    # if it is less than prob of mutation
    # it is randomized between 0 to 1
    pm = prob_of_mu
    final_family=family
    for i in range(0,final_family.shape[0],1):
        for j in range(0,final_family.shape[1],1):
            chance_for_mut = rnd.uniform(0,1)
            if chance_for_mut<pm :
                final_family[i,j]=np.random.randint(min_bound,max_bound)
            else:
                final_family[i,j]=final_family[i,j] 
    return final_family

##-----------------------------------------------------------##
    
# C2
 
input=[4,-2,7,-5,11,1] 
input_01  = pow(input[0], 2)
input_02  = pow(input[1], 3)

input[0] = input_01
input[1] = input_02   

best_outputs=[]
num_generation = 5

population = population(input,50,-100,100)       
    
for generation in range(num_generation):
    print("Generation: ",generation)
    population = GA(input,population,-100, 100, 0.8, 0.1)

fitness = np.sum(population*input, axis=1)
best_match_idx=np.where(fitness==np.max(fitness))

print("Optimum Y score : ", population[best_match_idx,:])
print("Best solution fitness : ",fitness[best_match_idx])

plt.plot(best_outputs)
plt.xlabel("iteration")
plt.ylabel("fitness")
plt.show()

##-----------------------------------------------------------##