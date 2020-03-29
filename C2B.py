#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 19:25:29 2020

@author: hteza
"""
import numpy as np
from matplotlib import pyplot as plt


def population (input,number_of_chrom,min_bound,max_bound):
    # creating population
    gene = len(input)
    no_chrom = number_of_chrom
    pop_size=(no_chrom,gene)
    population=np.random.randint(low=min_bound,high=max_bound,size=pop_size)
    return population

def GA(input, population, number_of_chrom, number_of_mu, min_bound, max_bound, prob_of_co, prob_of_mu):
    gene = len(input)
    population_size=(number_of_chrom,gene)
    # fitness funtion
    fitness = np.sum(population*input, axis=1)
    best_outputs.append(np.max(fitness))
    # parents
    # hlaf the population becomes parents
    num_parents = int(number_of_chrom / 2)
    parents = np.empty((num_parents, population.shape[1]))
    for parent_num in range(num_parents):
        max_fitness_idx = np.where(fitness == np.max(fitness))
        max_fitness_idx = max_fitness_idx[0][0]
        parents[parent_num, :] = population[max_fitness_idx, :]
        fitness[max_fitness_idx] = -999999999
    # offspring
    offspring_size=(population_size[0]-parents.shape[0],gene)
    offspring = np.empty(offspring_size)
    # always crossover at the middle
    crossover_point = np.uint8(offspring_size[1]/2)
    for k in range(offspring_size[0]): 
        # Random probability of crossover for each offspring
        prop_ran = np.random.uniform(low=0, high=1)
        # first half of first parent
        parent1_idx = k%parents.shape[0]
        # second half of second parent
        parent2_idx = (k+1)%parents.shape[0]
    #crossing over
        if prop_ran < prob_of_co: 
            offspring[k, 0:crossover_point] = parents[parent1_idx, 0:crossover_point]
            offspring[k, crossover_point:] = parents[parent2_idx, crossover_point:]
    # mutation
    # number of genes mutating = number of new children
    mutations_counter = np.uint8(offspring.shape[1] / number_of_mu)
    for idx in range(offspring.shape[0]):
        gene_idx = mutations_counter - 1
        for mutation_num in range(number_of_mu):
            #Random probability of mutation
            prop_ran = np.random.uniform(low=0, high=1)     
            if prop_ran < prob_of_mu: 
                random_value = np.random.uniform(-1.0, 1.0, 1)                
                offspring[idx, gene_idx] = offspring[idx, gene_idx] + random_value
            gene_idx = gene_idx + mutations_counter
            population[0:parents.shape[0], :] = parents
            population[parents.shape[0]:, :] = offspring
            return population
                
        
##-----------------------------------------------------------##
    
# C2
 
input=[4,-2,7,-5,11,1] 
input_01  = pow(input[0], 2)
input_02  = pow(input[1], 3)

input[0] = input_01
input[1] = input_02
   

best_outputs=[]
num_generation = 1000

population = population(input,50,-150,150)       
    
for generation in range(num_generation):
    print("Generation: ",generation)
    population =  GA(input, population, 50, 1, -150, 150, 0.8, 0.1)

fitness = np.sum(population*input, axis=1)

print("Optimum Y score : ", np.max(fitness)) 
print("Optimum Weight : ",population[np.argmax(fitness),:])


plt.plot(best_outputs)
plt.xlabel("iteration")
plt.ylabel("fitness")
plt.show()

##-----------------------------------------------------------##