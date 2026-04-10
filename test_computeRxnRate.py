import numpy as np
import pytest

from reaction import Reaction

_entityList = ['A', 'B', 'C']

def test_testReactionRate():
    testReaction = Reaction(_entityList, "A + B -> C", 1, reactantsDict={'A' : 1, 'B' : 1}, productsDict={'C' : 1})

    y = [1, 1, 0]
    assert testReaction.computeRxnRate(_entityList, y) == 1

def test_decimalReactionRate():
    testReaction = Reaction(_entityList, "A + B -> C", 0.4, reactantsDict={'A' : 1, 'B' : 1}, productsDict={'C' : 1})

    y = [1, 1, 0]
    assert testReaction.computeRxnRate(_entityList, y) == 0.4

def test_decimalMultiplicity():
    testReaction = Reaction(_entityList, "A + B -> C", 1, reactantsDict={'A' : 1, 'B' : 0.5}, productsDict={'C' : 1})
    
    y = [1, 1, 0]

    assert testReaction.computeRxnRate(_entityList, y) == 1

def test_decimalMultiplicityAndRxnRate():
    testReaction = Reaction(_entityList, "A + B -> C", 0.4, reactantsDict={'A' : 1, 'B' : 0.5}, productsDict={'C' : 1})
    
    y = [1, 1, 0]

    assert testReaction.computeRxnRate(_entityList, y) == 0.4

def test_decimalMultiplicityNonzeroConcentration():
    testReaction = Reaction(_entityList, "A + B -> C", 1, reactantsDict={'A' : 1, 'B' : 0.5}, productsDict={'C' : 1})

    y = [1, 2, 0]

    assert testReaction.computeRxnRate(_entityList, y) == 2**0.5