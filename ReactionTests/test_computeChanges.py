import pytest

from reactionTest import Reaction

def test_nullReaction():
    _entityList = ['A', 'B', 'C']
    nullReaction = Reaction(_entityList, "Null Rxn", 1, reactantsDict = {}, productsDict = {})
    nullReaction.computeChanges(_entityList)

    assert (nullReaction._changes[0] == [0, 0, 0]).all()
    assert (nullReaction._changes[1] == [0, 0, 0]).all()

def test_neutralReaction():
    _entityList = ['A', 'B', 'C']
    neutralReaction = Reaction(_entityList, "Neutral Rxn", 1, reactantsDict = {'A' : 0, 'B' : 0}, productsDict = {'C' : 0})
    neutralReaction.computeChanges(_entityList)

    assert (neutralReaction._changes[0] == [0, 0, 0]).all()
    assert (neutralReaction._changes[1] == [0, 0, 0]).all()

def test_onlyOneReactantChange():
    _entityList = ['A', 'B', 'C']
    reactantChangeRxn = Reaction(_entityList, "Reactant Change Rxn", 1, reactantsDict = {'A' : 1, 'B' : 0}, productsDict = {'C' : 0})
    reactantChangeRxn.computeChanges(_entityList)

    assert (reactantChangeRxn._changes[0] == [1, 0, 0]).all()
    assert (reactantChangeRxn._changes[1] == [0, 0, 0]).all()

def test_onlyReactantsChange():
    _entityList = ['A', 'B', 'C']
    reactantsChangingRxn = Reaction(_entityList, "Reactants Changing Rxn", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'C' : 0})
    reactantsChangingRxn.computeChanges(_entityList)

    assert (reactantsChangingRxn._changes[0] == [1, 1, 0]).all()
    assert (reactantsChangingRxn._changes[1] == [0, 0, 0]).all()

def test_allChange():
    _entityList = ['A', 'B', 'C']
    testReaction = Reaction(_entityList, "Test Rxn", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'C' : 1})
    testReaction.computeChanges(_entityList)

    assert (testReaction._changes[0] == [1, 1, 0]).all()
    assert (testReaction._changes[1] == [0, 0, 1]).all()

