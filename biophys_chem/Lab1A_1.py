#!/usr/bin/python

#import handy python modules
import numpy as np
import pylab as plt

# This is a function which takes a list of energies, one for each state
# in a system, and samples it for the required number of observations.
def sampleStates(stateEnergies, numberOfObservations):
    # Set the temperature and boltzmann constant in [eV/K] & [K]
    k = 8.6e-5
    T = 1000.0

    # Count how many state-energies were provided, that is how many states to use
    Nstates = len(stateEnergies)

    # Create an array to hold the probability interval (ceiling) for each state
    Probs = np.zeros(Nstates)

    # Begin by setting the Partition Function to zero.
    PartitionFunction = 0.0

    # Go through all state energies in the input
    for i in np.arange(Nstates):
        # calculate the boltzmann factor
        Probs[i]            = np.exp( (-1 * stateEnergies[i]) / (k * T) ) 
        # And keep adding all boltzmann factors to the Partion Function
        PartitionFunction   = PartitionFunction + Probs[i]

    # When all Boltzmann factors have been calulated and the Partition Function is
    # complete, we finish the probabilities by dividing by the Partition Function,
    # as Boltzmann statistics instructs us to ("Z" in Eq. 1 in the lab instruction)
    Probs = Probs /  PartitionFunction

    # Now we see each probability as a part of the line between 0 and 1. If the
    # probability of the first state "A" is e.g. 0.23, then all numbers between
    # 0.00 and 0.23 belong to state A. If the next state "B" has probability 0.31,
    # then all numbers between 0.23 and 0.54 (0.23+0.31) belong to state B, and
    # so on. An easy way to find these thresholds is to use a cumulative sum, where
    # each element in the array is the sum of all previous elements.
    ProbThresholds=np.cumsum(Probs)
    # If you would like to see how this works out, try printing them out by
    # un-commenting the lines below:

    # print("The probabilites of the states are")
    # print(Probs)
    # print("And the thresholds of the states we will use are")
    # print(ProbThresholds)

    # Now we will make lots of random numbers and assign each to a state. Hopefully,
    # the probabilites we made will guide this assignement to look like we expect
    # Boltzmann statistics to look.

    # Make an array to fill with state assignments
    states = np.zeros(numberOfObservations)

    # Make as many observations as required
    for i in np.arange(numberOfObservations):
        # Make a random number this observation
        current = np.random.rand()

        # Test if the current number is above the threshold for a state.
        # as long as it IS, keep testing the next state.
        while current>ProbThresholds[states[i]]:
                states[i]+=1

    # We now have a number of state assignments, as many as required. Let's create a plot
    # of them in histogram, so that the first state is "1", then "2", and so on.
    bins=np.arange(Nstates+1)+0.5
    y,dummy = plt.histogram(states+1, bins=bins, normed=True)

    plt.bar(np.arange(len(energies))+0.9, y, width=0.4, fc='blue',alpha=0.5, label='sampled states')


#  --------------- THE PROGAM STARTS HERE, AND USES THE ABOVE FUNCTION -----------------

# Set the temperature and boltzmann constant in [eV/K] & [K]
k = 8.6e-5
T = 1000.0

# Create a number of energy-states, in this case three of them. I.e. there are 3 states here.
energies=[-0.30, -0.35, -0.25]



# Set the number of observations of the system
steps = 500

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
#predictedDistribution = np.exp( (-1 * predictedDistribution) / (k * T) ) / PartitionFunction

# MAke a bar-plot in red
plt.bar(np.arange(len(energies))+0.5, predictedDistribution, width=0.4, fc='red',alpha=0.5, label='predicted by Boltzmann')

# Randomize numbers to see if it mathces the predicted distribution.
sampleStates(energies,steps)

plt.legend()
plt.show()
