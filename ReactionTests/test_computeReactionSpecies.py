import pytest
from reactionTest import Reaction

def test_nullReaction():
    _entityList = ['A', 'B', 'C']
    nullReaction = Reaction(_entityList, "Null Rxn", 1, reactantsDict = {}, productsDict = {})
    nullReaction.computeReactionSpecies()

    assert nullReaction._reactants == []
    assert nullReaction._products == []

def test_neutralReaction():
    _entityList = ['A', 'B', 'C']
    neutralReaction = Reaction(_entityList, "Neutral Rxn", 1, reactantsDict = {'A' : 0, 'B' : 0}, productsDict = {'C' : 0})
    neutralReaction.computeReactionSpecies()

    assert neutralReaction._reactants == ['A', 'B']
    assert neutralReaction._products == ['C']

def test_positiveReaction():
    _entityList = ['A', 'B', 'C']
    positiveReaction = Reaction(_entityList, "Positive Rxn", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'C' : 1})
    positiveReaction.computeReactionSpecies()

    assert positiveReaction._reactants == ['A', 'B']
    assert positiveReaction._products == ['C']

def test_onlyReactantsChanging():
    _entityList = ['A', 'B', 'C']
    reactantsChangingRxn = Reaction(_entityList, "Reactants Changing Rxn", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'C' : 0})
    reactantsChangingRxn.computeReactionSpecies()

    assert reactantsChangingRxn._reactants == ['A', 'B']
    assert reactantsChangingRxn._products == ['C']
