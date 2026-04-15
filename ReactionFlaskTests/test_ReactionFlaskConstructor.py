import pytest
from reaction_flaskTest import ReactionFlask

def test_nullEntityList():
    testFlask = ReactionFlask([], components = [])

    assert len(testFlask._entityList) == 0
    assert (testFlask._concentrations == []).all()

    #These won't be repeated in later tests, because they are
    # imperatively set this way
    # in the constructor
    assert testFlask._concentrationsInitialized == False
    assert testFlask._initialCondition == None
    assert testFlask.latestSimulationOutput == None

def test_nominalEntityList():
    testFlask = ReactionFlask(['A', 'B', 'C'], components=[])

    assert len(testFlask._entityList) == 3
    assert (testFlask._concentrations == [0, 0, 0]).all()

def test_abnormalTypeEntityList():
    with pytest.raises(TypeError, match = "The entity list contains non-string elements"):
        testFlask = ReactionFlask([1, 'A', []], components=[])

def test_abnormalTypeComponents():
    with pytest.raises(TypeError, match = "A component key has a non-string value"):
        testFlask = ReactionFlask(['A', 'B', 'C'], components= {1 : ['A', 'B'], 2 : ['C']})

def test_nominalReactionFlask():
    testFlask = ReactionFlask(['A', 'B', 'C'], components = {'Main' : ['A', 'B'], 'Auxiliary' : ['C']})

    assert len(testFlask._entityList) == 3
    assert len(testFlask._components) == 2
    assert (testFlask._concentrations == [0, 0, 0]).all()
    