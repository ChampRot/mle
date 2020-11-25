import random
from math import exp

import matplotlib.pyplot as pyplot

def random_bits(len):
    """Returns a list containing a random combination
    of 0s and 1s"""

    return [random.randint(0, 1) for x in range(len)]


def init_population(p, len):
    """Returns a list of random bit-lists"""

    return [random_bits(len) for x in range(p)]


def gen_volumes(len, max_v):
    """Return a list of positive random floats smaller max_v"""

    return [random.random() * max_v for x in range(len)]


def calculate_volume(hypo, volumes):
    """What more can I say?"""

    return sum(map(lambda b, v: b * v, hypo, volumes))


def calculate_fitness(hypo, volumes, c):
    """What more can I say?"""

    return exp(-c * ((100 - calculate_volume(hypo, volumes)) ** 2))


def create_probability_function(population, short_fittness_function):
    """Calculates the sum of the fittness of all the elements in population
    and returns a function using it.
    Expects a fittness functin which takes only one hypothesis as argument."""

    sum_fitness = sum(map(short_fittness_function, population))
    return lambda hypo: short_fittness_function(hypo) / sum_fitness


def select_index(population, probability_function):
    """Selects the index of an element in population
    by using probability_function"""

    index = random.randrange(0, len(population))
    rand_num = random.random()
    sum = 0
    while sum < rand_num:
        index = (index + 1) % len(population)
        sum += probability_function(population[index])

    return index


def selection(population, amount, probability_function):
    """Select amount elements from population with probability
    given by probability_function"""

    return [population[select_index(population, probability_function)] for x in range(amount)]


def select_fittest(population, amount, fittness_function):
    tmp = sorted(zip(population, map(fittness_function, population)), key=lambda x : x[1])
    
    return [tmp.pop()[0] for x in range(amount)]   


def crossover(bit_mother, bit_father):
    """Create one sibling with a single random crossover-point 
    for the given bit-lists"""

    cross_index = random.randrange(0, len(bit_mother))
    return bit_mother[:cross_index] + bit_father[cross_index:]


def mutate(gene):
    """Swaps a single Bit in the given list of bits"""

    i = random.randrange(0, len(gene))
    # => True is 1 and False is 0 for some reason
    gene[i] = not gene[i]


def evolution(p, r, m, run_ct, gene_length, short_fittness_function):
    """Try and find some good genes"""
    
    population = init_population(p, gene_length)
    stats = []
    for i in range(run_ct):
        p_function = create_probability_function(population, short_fittness_function)
        
        next_gen = selection(population, round((1 - r) * p), p_function)
        # next_gen = select_fittest(population, round((1 - r) * p), short_fittness_function)

        parents = selection(population, round(r * p), p_function)
        for j in range(1, len(parents), 2):
            next_gen += [crossover(parents[j], parents[j - 1])]
            next_gen += [next_gen[len(next_gen) - 1].copy()]
        # wenn es eine ungerade anzahl an eltern gibt glaube ich
        if len(next_gen) < len(population):
            next_gen += [crossover(parents[random.randrange(0, len(parents))], 
                                  parents[random.randrange(0, len(parents))])]

        for j in range(round(m * p)):
            mutate(next_gen[random.randrange(0, len(next_gen))])

        population = next_gen
        fittest = 0
        best_gene = None
        for gene in population:
            if short_fittness_function(gene) > fittest:
              fittest = short_fittness_function(gene) 
              best_gene = gene
        stats.append(fittest)
        print(fittest)
    return (fittest, best_gene, stats)


c = 0.01

volumes = gen_volumes(100, 10)
r = evolution(10, 0.4, 0.3, 50, 100, lambda gene : calculate_fitness(gene, volumes, c))
print(calculate_volume(r[1], volumes))

pyplot.plot(list(range(50)), r[2], )
pyplot.show()
