# stage-1
Translate DNA to protein
team = [ 
{ 
"name": "Alabi Saheed", 
"slack_username": "@Alabi", 
"email": "alabi@example.com", "hobby": "Reading", 
"country": "Nigeria", 
"discipline": "Software Engineering", "preferred_language": "Python" }, 
{ 
"name": "John Doe", 
"slack_username": "@JohnD", 
"email": "johnd@example.com", "hobby": "Cycling", 
"country": "USA", 
"discipline": "Data Science", 
"preferred_language": "R" 
}, 
{ 
"name": "Jane Smith", 
"slack_username": "@JaneS", 
"email": "janes@example.com", "hobby": "Painting", 
"country": "Canada", 
"discipline": "Cybersecurity", 
"preferred_language": "C++" 
}, 
{ 
"name": "Carlos Mendoza", 
"slack_username": "@CarlosM", "email": "carlosm@example.com", "hobby": "Hiking", 
"country": "Mexico", 
"discipline": "Embedded Systems", "preferred_language": "C" 
}, 
{ 
"name": "Aisha Khan", 
"slack_username": "@AishaK", "email": "aishak@example.com", "hobby": "Writing", 
"country": "India", 
"discipline": "Web Development",
"preferred_language": "JavaScript" 
} 
] 
print("Team Members:\n") 
for member in team: 
print(f"Name: {member['name']}") 
print(f"Slack Username: {member['slack_username']}") print(f"Email: {member['email']}") 
print(f"Hobby: {member['hobby']}") 
print(f"Country: {member['country']}") 
print(f"Discipline: {member['discipline']}") 
print(f"Preferred Language: {member['preferred_language']}") print("-" * 40) 
Result: 
Team Members: 
Name: Alabi Saheed 
Slack Username: @Alabi 
Email: alabi@example.com 
Hobby: Reading 
Country: Nigeria 
Discipline: Software Engineering 
Preferred Language: Python 
---------------------------------------- 
Name: John Doe 
Slack Username: @JohnD 
Email: johnd@example.com 
Hobby: Cycling 
Country: USA 
Discipline: Data Science 
Preferred Language: R 
---------------------------------------- 
Name: Jane Smith 
Slack Username: @JaneS 
Email: janes@example.com 
Hobby: Painting 
Country: Canada 
Discipline: Cybersecurity 
Preferred Language: C++ 
---------------------------------------- 
Name: Carlos Mendoza
Slack Username: @CarlosM Email: carlosm@example.com Hobby: Hiking 
Country: Mexico 
Discipline: Embedded Systems Preferred Language: C 
---------------------------------------- Name: Aisha Khan 
Slack Username: @AishaK Email: aishak@example.com Hobby: Writing 
Country: India 
Discipline: Web Development Preferred Language: JavaScript ----------------------------------------






Stage zero;

team = [
{
"name": "Alabi Saheed",
"slack_username": "@Alabi",
"email": "alabi@example.com", "hobby": "Reading",
"country": "Nigeria",
"discipline": "Software Engineering", "preferred_language": "Python" },
{
"name": "John Doe",
"slack_username": "@JohnD",
"email": "johnd@example.com", "hobby": "Cycling",
"country": "USA",
"discipline": "Data Science",
"preferred_language": "R"
},
{
"name": "Jane Smith",
"slack_username": "@JaneS",
"email": "janes@example.com", "hobby": "Painting",
"country": "Canada",
"discipline": "Cybersecurity",
"preferred_language": "C++"
},
{
"name": "Carlos Mendoza",
"slack_username": "@CarlosM", "email": "carlosm@example.com", "hobby": "Hiking",
"country": "Mexico",
"discipline": "Embedded Systems", "preferred_language": "C"
},
{
"name": "Aisha Khan",
"slack_username": "@AishaK", "email": "aishak@example.com", "hobby": "Writing",
"country": "India",
"discipline": "Web Development",
"preferred_language": "JavaScript"
}
]
print(team)
15/02/25  
print("Hello, World!")
comment is basically any line of code or text you do not want to run.
Comments can be used to explain Python code.
Comments can be used to make the code more readable.
Comments can be used to prevent execution when testing code.
#the next line prints Hello, HackBio!
print("Hello, HackBio!") #this is the code line itself


A comment does not have to be text that explains the code, it can also be used to prevent Python from executing code:
#print("Hello, World!")
print("Cheers to learning bioinformatics!")

Python does not really have a syntax for multiline comments.
However, you can use the multi-line string literal in your code to place your comment. Further uses of this will be discussed in future sections
"""
This is a comment
written in
more than just one line
"""
16/02/2025
Write a function for translating DNA to protein
Write a function that simulates and generates a logistic population growth curve. Your function should include 2 extra parameters that randomize the length of the lag phase and the exponential phase [See population curve here] . Most living populations follow a logistic population growth. Therefore, your growth curve can be: Population Size vs Time, Cell density vs Time, OD vs Time, CFU vs Time, etc
Using your function, generate a dataframe with 100 different growth curves
Write a function for determining the time to reach 80% of the maximum growth; usually the carrying capacity
Finally, write a function for calculating the hamming distance between your Slack username and twitter/X handle (synthesize if you don’t have one). Feel free to pad it with extra words if they are not of the same length.

solution :





def translate_dna_to_protein(dna_sequence): # Codon table for translation codon_table = { 'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTS': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'TAA': '', 'TAG': '', 'TGA': '', # Stop codons 'TTA': 'L', 'TTC': 'F', 'TTT': 'F', 'TTG': 'L', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TAC': 'Y', 'TAT': 'Y', 'TAA': '', 'TAG': '',       'TGC': 'C', 'TGT': 'C', 'TGA': '', 'TGG': 'W',
    }
    
    protein = []
    # Translate DNA to protein
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        protein.append(codon_table.get(codon, ''))
    
    return ''.join(protein)


2. Function for Simulating Logistic Population Growth


import numpy as np import pandas as pd def logistic_growth_curve(K, r, t, lag_phase, exp_phase): # Initial population size P0 = 1 # Time array time = np.linspace(0, t, t) # Logistic growth formula P = K / (1 + (K/P0 - 1) * np.exp(-r * (time - lag_phase))) # Randomize the exponential phase growth = P * np.exp(r * (time - lag_phase) * (time < (lag_phase + exp_phase))) return time, growth

# Example usage: K = 100 # Carrying capacity r = 0.1 # Growth rate t = 100 # Time points lag_phase = np.random.randint(1, 10) exp_phase = np.random.randint(10, 30) time, population = logistic_growth_curve(K, r, t, lag_phase, exp_phase) # Generate a dataframe with 100 different growth curves growth_curves = [] for _ in range(100): lag_phase = np.random.randint(1, 10) exp_phase = np.random.randint(10, 30) time, population = logistic_growth_curve(K, r, t, lag_phase, exp_phase) growth_curves.append(population)
df_growth_curves = pd.DataFrame(growth_curves)

3. Function for Determining Time to Reach 80% of Maximum Growth

def time_to_reach_80_percent(K, r, lag_phase, exp_phase):
    target_population = 0.8 * K
    time, population = logistic_growth_curve(K, r, 100, lag_phase, exp_phase)
    time_to_target = time[np.where(population >= target_population)[0][0]]
    return time_to_target

4. Function for Calculating Hamming Distance

def hamming_distance(str1, str2):
    # Pad with extra words if needed
    max_length = max(len(str1), len(str2))
    str1 = str1.ljust(max_length)
    str2 = str2.ljust(max_length















Stage 1 personal.
def translate_dna_to_protein(dna_sequence):
    # Genetic code dictionary
    genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S…def hamming_distance(slack, twitter):
    #the shorter string will be padded with spaces so both have the same length
    if len(slack) < len(twitter):
        slack = slack + " " * (len(twitter) - len(slack))
    elif len(twitter) < len(slack):
        twitter = twitter + " " * (len(slack) - len(twitter))
    
    # the counter for hamming distance is initialized as zero
    distance = 0
    #next, the function is looped for each character
    for i in range(len(slack)):
        if slack[i] != twitter[i]:
            distance += 1
    return distance

# slack and twitter usernames
slack_username = "Oluwafisayo"
twitter_username = "hermeeenarh"# calculating and printing the final Hamming distance
result = hamming_distance(slack_username, twitter_username)
print("Hamming distance:", result)
def population(x , k = 2, x_mid = 5):
    solution = 1/(1 + math.exp(-k*(x-x_mid)))
    return solution
print(population(k = 2, x_mid = 5, x = 10))
od_600 = []
for i in range(0,24):
    od_600.append(population(x = i, k = 0.5, x_mid = 10))
print(od_600)


#the function to generate 100 different growth curves.
def generate_growth_curves(num_curves, time_points):
    curves = {}
    for i in range(num_curves):
        curve = simulate_logistic_curve(time_points)
        curves["Curve_" + str(i + 1)] = curve
    return curves


#the function to determine the time when the growth curve reaches 80% of the carrying capacity.
def time_to_reach_80_percent(time_points, growth_curve, K=1):
    threshold = 0.8 * K
    for i in range(len(time_points)):
        if growth_curve[i] >= threshold:
            return time_points[i]
    return None


