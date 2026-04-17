from reaction_flaskTest import ReactionFlask
import pytest

def test_negativeTime():
    testRxnFlask = ReactionFlask(['A', 'B', 'C'], components = {'Main' : ['A', 'B', 'C']})
    testRxnFlask.addReaction("A + B -> C", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'C' : 1})

    with pytest.raises(ValueError, match = "Input time has a negative value"):
        testRxnFlask.reactionDeriv(-20, [1, 1, 1])

def test_negativeYValues():
    testRxnFlask = ReactionFlask(['A', 'B', 'C'], components = {'Main' : ['A', 'B', 'C']})
    testRxnFlask.addReaction("A + B -> C", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'C' : 1})

    with pytest.raises(ValueError, match = "Element"):
        testRxnFlask.reactionDeriv(0, [1, -1, 1])

def test_nominalRxnFlaskOneRxn():
    testRxnFlask = ReactionFlask(['A', 'B', 'C'], components = {'Main' : ['A', 'B', 'C']})
    testRxnFlask.addReaction("A + B -> C", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'C' : 1})

    testRxnFlask.setInitialCondition([1, 1, 0])

    rxnDeriv = testRxnFlask.reactionDeriv(0, testRxnFlask._concentrations)

    assert (rxnDeriv == [-1, -1, 1]).all()
    
    rxnDeriv = testRxnFlask.reactionDeriv(0, [1, 2, 3])

    assert (rxnDeriv == [-2, -2, 2]).all()

def test_nominalRxnFlaskTwoRxn():
    testRxnFlask = ReactionFlask(['A', 'B', 'C'], components = {"Main" : ['A', 'B', 'C']})
    testRxnFlask.addReaction("A + B -> C", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'C' : 1})
    testRxnFlask.addReaction("C -> A + B", 1, reactantsDict = {'C' : 1}, productsDict = {'A' : 1, 'B' : 1})

    testRxnFlask.setInitialCondition([1, 1, 0])

    rxnDeriv = testRxnFlask.reactionDeriv(0, testRxnFlask._concentrations)

    assert (rxnDeriv == [-1 + 0, -1 + 0, 1 - 0]).all()

    rxnDeriv = testRxnFlask.reactionDeriv(0, [1, 2, 3])

    assert (rxnDeriv == [-2 + 3, -2 + 3, 2 - 3]).all()

