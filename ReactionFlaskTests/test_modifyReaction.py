import pytest
from reaction_flaskTest import ReactionFlask

testReactionFlask = ReactionFlask(['A', 'B', 'C'], components={'Main' : ['A', 'B', 'C']})

testReactionFlask.setInitialCondition([1, 1, 0])

testReactionFlask.addReaction('Test Rxn', 1, reactantsDict={'A' : 1}, productsDict={'C' : 1})

def testNewReactionRate():
    testReactionFlask.modifyReaction('Test Rxn', newRxnRate=0)
    testReaction = testReactionFlask._reactions['Test Rxn']
    assert testReaction._rxnK == 0

def testNewReactantsDict():
    testReactionFlask.modifyReaction('Test Rxn', newRxnReactantsDict={'A' : 2})
    testReaction = testReactionFlask._reactions['Test Rxn']
    assert testReaction._reactantsDict == {'A' : 2}

def testNewProductsDict():
    testReactionFlask.modifyReaction('Test Rxn', newRxnProductsDict={'C' : 2})
    testReaction = testReactionFlask._reactions['Test Rxn']
    assert testReaction._productsDict == {'C' : 2}

def testNewProductsList():
    testReactionFlask.modifyReaction('Test Rxn', newRxnProductsDict={'B' : 1})
    testReaction = testReactionFlask._reactions['Test Rxn']
    assert testReaction._products == ['B']

def testNewReactantsList():
    testReactionFlask.modifyReaction('Test Rxn', newRxnReactantsDict={'C' : 1})
    testReaction = testReactionFlask._reactions['Test Rxn']
    assert testReaction._reactants == ['C']