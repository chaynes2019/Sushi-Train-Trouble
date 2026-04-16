import pytest
from reaction_flaskTest import ReactionFlask

testReactionFlask = ReactionFlask(['A', 'B', 'C'], components = {"Main" : ['A', 'B', 'C']})

def test_setNoneInitialCondition():
    with pytest.raises(TypeError, match = "y0 has type None"):
        testReactionFlask.setInitialCondition(None)

def test_setNegativeInitialCondition():
    with pytest.raises(ValueError, match = "y0 has been initialized with negative values"):
        testReactionFlask.setInitialCondition([-1, -1, -1])

def test_setZeroInitialCondition():
    testReactionFlask.setInitialCondition([0, 0, 0])

    assert (testReactionFlask._initialCondition == [0, 0, 0]).all()

def test_positiveInitialCondition():
    testReactionFlask.setInitialCondition([1, 2, 3])

    assert (testReactionFlask._initialCondition == [1, 2, 3]).all()