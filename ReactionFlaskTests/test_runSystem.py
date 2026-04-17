from reaction_flaskTest import ReactionFlask
import pytest
import numpy as np

def test_runningUninitializedSystem():
    testRxnFlask = ReactionFlask(['A', 'B', 'C'], components = {"Main" : ['A', 'B', 'C']})
    testRxnFlask.addReaction("A + B -> C", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'C' : 1})
    
    with pytest.raises(AttributeError, match = "Reaction flask has not yet been initialized."):
        testRxnFlask.runSystem(100)

def test_negativeTimeEndpoint():
    testRxnFlask = ReactionFlask(['A', 'B', 'C'], components = {"Main" : ['A', 'B', 'C']})
    testRxnFlask.addReaction("A + B -> C", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'C' : 1})

    testRxnFlask.setInitialCondition([1, 1, 1])

    with pytest.raises(ValueError, match = "Specified time endpoint is negative!"):
        testRxnFlask.runSystem(-100)

def test_neutralRxnFlaskZeroInitCondition():
    neutralRxnFlask = ReactionFlask(['A', 'B', 'C'], components = {"Main" : ['A', 'B', 'C']})
    neutralRxnFlask.addReaction("A + B -> C", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'C' : 1})

    neutralRxnFlask.setInitialCondition([0, 0, 0])

    neutralRxnFlask.runSystem(100)

    finalSimulationOutput = neutralRxnFlask.latestSimulationOutput["y"]

    finalSimulationValues = [finalSimulationOutput[k][-1] for k in range(3)]

    #This works, without using .all(),
    #  because the finalSimulationOutput is np.float()
    assert finalSimulationValues == [0, 0, 0]

def test_neutralRxnFlaskOnesInitCondition():
    neutralRxnFlask = ReactionFlask(['A', 'B', 'C'], components = {"Main" : ['A', 'B', 'C']})
    neutralRxnFlask.addReaction("A + B -> A + B", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'A' : 1, 'B' : 1})

    neutralRxnFlask.setInitialCondition([1, 1, 1])

    neutralRxnFlask.runSystem(100)

    finalSimulationOutput = neutralRxnFlask.latestSimulationOutput["y"]

    finalSimulationValues = [finalSimulationOutput[k][-1] for k in range(3)]

    #This works, without using .all(),
    #  because the finalSimulationOutput is np.float()
    assert finalSimulationValues == [1, 1, 1]

def test_nominalRxnFlask1():
    nominalRxnFlask = ReactionFlask(['A', 'B', 'C'], components = {"Main" : ['A', 'B', 'C']})
    nominalRxnFlask.addReaction("A + B -> C", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'C' : 1})

    nominalRxnFlask.setInitialCondition([1, 1, 0])

    nominalRxnFlask.runSystem(1000)

    finalSimulationOutput = nominalRxnFlask.latestSimulationOutput["y"]

    finalSimulationValues = np.array([finalSimulationOutput[k][-1] for k in range(3)])

    assert (np.abs(np.array([0, 0, 1]) - finalSimulationValues) < 0.01).all()

def test_nominalRxnFlask2():
    nominalRxnFlask = ReactionFlask(['A', 'B', 'C'], components = {"Main" : ['A', 'B', 'C']})
    nominalRxnFlask.addReaction("A + B -> C", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'C' : 1})
    nominalRxnFlask.addReaction("C -> A + B", 1, reactantsDict = {'C' : 1}, productsDict = {'A' : 1, 'B' : 1})

    nominalRxnFlask.setInitialCondition([1, 1, 0])

    nominalRxnFlask.runSystem(1000)

    finalSimulationOutput = nominalRxnFlask.latestSimulationOutput["y"]

    finalSimulationValues = np.array([finalSimulationOutput[k][-1] for k in range(3)])

    #These can be computed by realizing that,
    # at equilbrium, 1 * A * B = 1 * C, and
    # A = B, and A = startingVal - C = 1 - C
    # Compute from there using algebra 
    finalAandBVal = (-1 / 2) + (np.sqrt(5) / 2)

    finalCVal = (1 / 4) * (6 - 2 * np.sqrt(5))

    assert (np.abs(np.array([finalAandBVal, finalAandBVal, finalCVal]) - finalSimulationValues) < 0.01).all()

