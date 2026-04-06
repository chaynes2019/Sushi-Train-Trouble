import numpy as np
import matplotlib.pyplot as plt

from reaction_flask import ReactionFlask


rxnFlask = ReactionFlask(['A', 'B', 'C'], components = {"Main" : ['A', 'B', 'C']})
rxnFlask.setInitialCondition([1, 1, 0])

rxnFlask.addReaction("Test Rxn: A + B -> C", 1, reactantsDict={'A' : 1, 'B' : 1}, productsDict={'C' : 1})
rxnFlask.addReaction("Test Rxn: A + B -> A", 0, reactantsDict={'A' : 1, 'B' : 1}, productsDict={'A' : 1})

rxnFlask.runSystem(100)

timeVals = rxnFlask.latestSimulationOutput["t"]
cVals = rxnFlask.latestSimulationOutput["y"][2]
EightyPercent = np.where(cVals >= 0.8)[0]
print(timeVals[EightyPercent[0]])

rxnFlask.plotSystem()