from reactionTest import Reaction
import pytest

def test_NoneEntityListInput():
    with pytest.raises(ValueError):
        testReaction = Reaction(None, "Test Rxn", 1, reactantsDict={'A' : 1, 'B' : 1}, productsDict={'C' : 1})