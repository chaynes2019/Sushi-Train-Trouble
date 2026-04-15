import pytest
from reaction_flaskTest import ReactionFlask

testReactionFlask = ReactionFlask(['A', 'B', 'C'], components = {"Main" : ['A', 'B', 'C']})

def test_setNoneInitialCondition():
    with pytest.raises(TypeError, match = "y0 has type None"):
        testReactionFlask.setInitialCondition(None)