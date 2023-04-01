import pandas as pd

def epistatic_genome_frequency(length, number_of_mutations,neutral_prob,epistasis_rate):
    if(number_of_mutations < 2):
        raise Exception('number_of_mutations must be at least 2')
    elif(number_of_mutations%2!=0):
        raise Exception('number_of_mutations must be even')

    epistatic_prob = neutral_prob * epistasis_rate
    product = 1
    for k in range(0,number_of_mutations//2):
        N_k = (length-2*k) * (epistatic_prob) * (length-2*k+1) * (1-neutral_prob)
        product *= N_k
    
    return product

def aa_fitness():
    path_to_data = "../data/aa_fitness.csv"

    df = pd.read_csv(path_to_data)
    s_rows = df[df['gene'] == 'S']

    positive_fitness = s_rows[s_rows['fitness'] >= 0]

    positive_fitness_frequency = len(positive_fitness) / len(s_rows)

    print(f"S gene aa length: {len(s_rows)}")
    print(f"Rows with positive fitness: {len(positive_fitness)}")
    print(f"Positive fitness frequency: {positive_fitness_frequency}")

def qc():
    length = 200
    number_of_mutations = [2,4,6,8,10,12,14,16,18,20]
    neutral_prob = .8523
    epistasis_rate = .00001

    sum = 0

    for n in number_of_mutations:
        freq = epistatic_genome_frequency(length, n,neutral_prob,epistasis_rate)
        sum += freq
        path_to_data = "../data/epistatic_genome_frequency.csv"
        with open(path_to_data, 'a') as f:
            print(str(n) + ',' + str(freq))
            f.write(str(n) + ',' + str(freq))
            f.write("\n")

    print(f"Sum of frequencies: {sum}")


if __name__ == "__main__":
    qc()
    #aa_fitness()