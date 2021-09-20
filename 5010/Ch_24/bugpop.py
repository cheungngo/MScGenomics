"""Personal Module for the generation of bug populations"""

import random
import statistics

class bug:
    """class to generate individual bugs"""
    
    def __init__(self):
        """initializing genome for a bug"""
        self.bases = ["A", "T", "C", "G"]
        self.genome = [random.choice(self.bases) for i in range(100)]
        
    def mutate_base(self):
        """mutate a genome randomly"""
        base_no = random.randint(0,99)
        self.genome[base_no] = random.choice(self.bases)
        
    def set_base(self, index, base):
        """artificially mutate a genome"""
        self.genome[index] = base
        
    def get_fitness(self):
        """get a arbituary fitness score"""
        self.fitness = 0
        fitness = 0
        for i in range(100):
            if self.genome[i] == "G" or self.genome[i] == "C":
                fitness += 1
        genome_str = "".join(self.genome)
        if "AAA" in genome_str:
            fitness += 5
        self.fitness = fitness
        
class population:
    """Class to manage a population of bugs"""
    
    def __init__(self):
        """Generating 50 bugs by default"""
        self.bug_list = [bug() for i in range(50)]
    
    def create_offspring(self):
        """Create offspring based on existing genome"""
        new_pop = list()
        oldbug = [i for i in self.bug_list]
        for i in self.bug_list:
            i.mutate_base()
        newbug = [i for i in self.bug_list]
        new_pop = oldbug + newbug
        self.bug_list = new_pop
    
    def cull(self):
        """Keeping the fittest bugs"""
        new_pop = list()
        fitness = list()
        for i in self.bug_list:
            i.get_fitness()
            fitness.append(i.fitness)
        halfnumber = sorted(fitness, reverse = True)[len(fitness)//2]
        for i in range(len(fitness)):
            if fitness[i] > halfnumber:
                new_pop.append(self.bug_list[i])
        for i in range(len(fitness)):
            if fitness[i] == halfnumber:
                if len(new_pop) < len(fitness)//2:
                    new_pop.append(self.bug_list[i])
                else:
                    break
        self.bug_list = new_pop
        
    def mean_fitness(self):
        """Calculating the mean fitness"""
        fitness = list()
        for i in self.bug_list:
            i.get_fitness()
            fitness.append(i.fitness)
        return(statistics.mean(fitness))