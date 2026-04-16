from reaction_flaskTest import ReactionFlask
import pytest

testRxnFlask = ReactionFlask(['A', 'B', 'C'], components = {"Main" : ['A', 'B', 'C']})

def test_resetUninitializedFlask():
    with pytest.raises(TypeError, match = "Initial condition has not yet been set"):
        testRxnFlask.resetInitialCondition()

def test_resetUninstantiatedFlask():
    ReactionFlask._initialCondition = [2, 2, 2]
    ReactionFlask.resetInitialCondition(ReactionFlask)

def test_resetInitializedFlask():
    testRxnFlask.setInitialCondition([1, 0, 0])

    testRxnFlask._concentrations = [0, 1, 0]
    testRxnFlask.resetInitialCondition()

    assert (testRxnFlask._concentrations == [1, 0, 0]).all()
    assert testRxnFlask._concentrationsInitialized == True

def test_resetFlaskAfterRun():
    secondFlask = ReactionFlask(['M', 'N', 'P'], components = {'Secondary' : ['M', 'N', 'P']})
    secondFlask.setInitialCondition([1, 0, 0])
    secondFlask.addReaction("M -> N + P", 1, reactantsDict = {'M' : 1}, productsDict = {'N' : 1, 'P' : 1})

    secondFlask.runSystem(10)

    assert secondFlask._concentrationsInitialized == False

    secondFlask.resetInitialCondition()

    assert (secondFlask._concentrations == [1, 0, 0]).all()
    assert secondFlask._concentrationsInitialized == True