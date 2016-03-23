#!/usr/bin/python

#import handy python modules
import numpy as np
import pylab as plt

# returns True is A is prime
def is_prime(A):
    return all(A % i for i in xrange(2, A))

# My rule for accepting or rejecting a suggested Monte-Carlo move
# Return "True" to accpt move, or "False" to reject move
def PrimeRule(stateThere):
    if(is_prime(stateThere)):
        return(True)
    else:
        chance = np.random.rand()
        if(chance<(1.0/6.0)):
            return(True)
        else:
            return(False)

# HERE YOU SHOULD DEFINE AND LATER USE A NEW, BETTER RULE FOR BOLTZMANN STATISTICS
def BoltzmannRule(e0, e1, k, T):
    if np.random.rand() < np.exp( ( -(e1 - e0)) / (k * T) ): 
        return (True)
    else:
        return (False)
        




def simulate_rule(stateEnergies,Nsteps):
        # set physical constants
        k = 8.6e-5
        T = 100.0

        # check how many states were provided
        Nstates = len(stateEnergies)

        # Create an array to show were we were (which states we were in)
        # throughout the simulation
        statesOverTime = np.zeros(Nsteps)

        # Make the first state random by using the randint()-function
        statesOverTime[0] = np.random.randint(0,Nstates)

        # Now perform the simulation, by looping thourhg all the steps
        for i in np.arange(Nsteps-1):

            # Randomize a new state until we get one which is
            # different from the one we are allready in
            currentState = statesOverTime[i].astype(int)
            newState = currentState.astype(int)
            while newState == currentState:
                newState = np.random.randint(0,Nstates)

            # Make a random number between 0 and 1
            chance = np.random.rand()

            #should we move?
            move=False
            #move=PrimeRule(newState)
            move= BoltzmannRule(stateEnergies[currentState], stateEnergies[newState], k , T)

            if(move):
                # move to the proposed state
                statesOverTime[i+1] = newState
            else:
                # set the next state to be the same as the current
                statesOverTime[i+1] = currentState

        # Let's plot the resulting disitrbution in blue
        bins=np.arange(Nstates+1)+0.5
        y,dummy = plt.histogram(statesOverTime+1, bins=bins, normed=True)
        plt.bar(np.arange(len(energies))+0.9, y, width=0.4, fc='blue',alpha=0.5, label='sampledByRule')

        # Let's also return the array so that we can plot more things if we want to
        return statesOverTime

#  --------------- THE PROGAM STARTS HERE, AND USES THE ABOVE FUNCTION -----------------

# Set the temperature and boltzmann constant in [eV/K] & [K]
k = 8.6e-5
T = 100.0

# Create a number of energy-states, in this case three of them
energies=[-0.30, -0.35, -0.25, -0.30, -0.20 ]

#   FIRST, LET'S PLOT THE PREDICTED BOLTZMANN DISTRIBUTION

# Create an array to hold the probabilites that is
# the same length as the array with state-energies
predictedDistribution = np.zeros(len(energies))

PartitionFunction = 0.0
for i in np.arange(len(energies)):
    # calculate the boltzmann-factor
    predictedDistribution[i]  = np.exp( (-1 * energies[i]) / (k * T) ) 
    # And keep adding all boltzmann factors to the Partion Function
    PartitionFunction = PartitionFunction + predictedDistribution[i]

# Normalize by PartitionFunction "Z"
predictedDistribution = predictedDistribution / PartitionFunction

# MAke a bar-plot in red
plt.bar(np.arange(len(energies))+0.5, predictedDistribution, width=0.4, fc='red',alpha=0.5, label='prediction')

# Randomize numbers to see if it mathces the predicted distribution.
output = simulate_rule(energies,5000)


plt.legend()
plt.show()
