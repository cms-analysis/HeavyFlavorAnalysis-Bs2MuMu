#!/usr//bin/env python
import sys
import random

if len(sys.argv) < 3:
    print "Not enough arguments"
    print sys.argv[0] + " <output_file> <Nbr of configs>"
    sys.exit()

DEF = "!H:V:TestRate=10"
FILE = open(sys.argv[1],"w");
TOTAL = int(sys.argv[2])

print "Output file: '" + sys.argv[1] + "'"
print "Generate " + str(TOTAL) + " configurations"

# BOUNDARIES
VAR_TRANSFORM = ("Gauss","Deco")
ACT_FUNCTIONS = ("sigmoid","tanh","radial")
NUM_CYCLES = (600,2000)
LEARN_RATE = (0.01,0.5)
DECAY_RATE = (0.001,0.05)
FIRST_LAYER = (-5,10)
SECOND_LAYER = (-5,5)

print "Input Variable Transformation " + str(VAR_TRANSFORM)
print "Activation function " + str(ACT_FUNCTIONS)
print "NCycles in" + str(NUM_CYCLES)
print "Learn rate in " + str(LEARN_RATE)
print "Decay rate in " + str(DECAY_RATE)
print "First layer in " + str(FIRST_LAYER)
print "Second layer in " + str(SECOND_LAYER)

for nbr in range(1,TOTAL+1):
    # Start with default string
    OPTS = DEF
    # Input Variable Transformation
    OPTS = OPTS + ":VarTransform=Norm"
    for intrans in VAR_TRANSFORM:
        if random.random() >= 0.5: # Add Gaussianization
            OPTS = OPTS + "," + intrans
    # Activation function
    OPTS = OPTS + ":NeuronType=" + ACT_FUNCTIONS[int(random.random()*len(ACT_FUNCTIONS))]
    # Number of Cycles
    OPTS = OPTS  + ":NCycles=" + str(int(random.uniform(NUM_CYCLES[0],NUM_CYCLES[1])))
    # Hidden Layer Configuration
    OPTS = OPTS + ":HiddenLayers=N+" + str(int(round(random.uniform(FIRST_LAYER[0],FIRST_LAYER[1]))))
    if random.random() >= 0.5: # add a second layer
        OPTS = OPTS + ",N+" + str(int(round(random.uniform(SECOND_LAYER[0],SECOND_LAYER[1]))))
    # Learning Rate
    OPTS = OPTS + ":LearningRate=" + str(round(random.uniform(LEARN_RATE[0],LEARN_RATE[1]),3))
    OPTS = OPTS + ":DecayRate=" + str(round(random.uniform(DECAY_RATE[0],DECAY_RATE[1]),4))
    # Fix hidden layer issue
    OPTS = OPTS.replace("+-","-")
    FILE.write(OPTS + "\n")

FILE.close()
