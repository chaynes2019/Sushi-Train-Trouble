from reactionTest import Reaction
import pytest

def test_NoneEntityListInput():
    with pytest.raises(TypeError):
        testReaction = Reaction(None, "Test Rxn", 1, reactantsDict={'A' : 1, 'B' : 1}, productsDict={'C' : 1})

def test_NoneNameInput():
    with pytest.raises(ValueError):
        testReaction = Reaction(['A', 'B', 'C'], None, 1, reactantsDict={'A' : 1, 'B' : 1}, productsDict={'C' : 1})

def test_NoneReactionRateInput():
    with pytest.raises(ValueError):
        testReaction = Reaction(['A', 'B', 'C'], "Test Rxn", None, reactantsDict={'A' : 1, 'B' : 1}, productsDict={'C' : 1})

def test_NoneReactantsDictInput():
    with pytest.raises(TypeError):
        testReaction = Reaction(['A', 'B', 'C'], "Test Rxn", 1, reactantsDict = None, productsDict={'C' : 1})

def test_NoneProductsDictInput():
    with pytest.raises(TypeError):
        testReaction = Reaction(['A', 'B', 'C'], "Test Rxn", 1, reactantsDict={'A' : 1, 'B' : 1}, productsDict=None)

def test_nullEntityListInput():
    testReaction = Reaction([], "Test Rxn", 1, reactantsDict = {}, productsDict = {})

def test_inputCharacterForMultiplicity():
    with pytest.raises(ValueError):
        testReaction = Reaction(['A', 'B', 'C'], "Test Rxn", 1, reactantsDict = {'A' : 'a', 'B' : 1}, productsDict = {'C' : 1})

def test_inputNegativeValueForMultiplicity():
    with pytest.raises(ValueError):
        testReaction = Reaction(['A', 'B', 'C'], "Test Rxn", 1, reactantsDict = {'A' : -1, 'B' : 1}, productsDict = {'C' : 1})

