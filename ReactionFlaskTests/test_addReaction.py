import numpy as np
import pytest
import re

from reaction_flaskTest import ReactionFlask

testRxnFlask = ReactionFlask(['A', 'B', 'C'], components={"Main" : ['A', 'B', 'C']})

def test_addTestReaction():
    testRxnFlask.addReaction('A + B -> C', 1, reactantsDict={'A' : 1, 'B' : 1}, productsDict={'C' : 1})
    assert len(testRxnFlask._reactions) == 1
    reactionNames = list(testRxnFlask._reactions.keys())
    assert reactionNames[0] == "A + B -> C"

def test_add2ndReaction():
    testRxnFlask.addReaction('C -> A + B', 1, reactantsDict={'C' : 1}, productsDict={'A' : 1, 'B' : 1})
    assert len(testRxnFlask._reactions) == 2
    reactionNames = list(testRxnFlask._reactions.keys())
    assert reactionNames[1] == "C -> A + B"

def test_addReactantOutsideEntityList():
    with pytest.raises(ValueError):
        testRxnFlask.addReaction("A + D -> B", 1, reactantsDict={'A' : 1, 'D' : 1}, productsDict={'B' : 1})

def test_addProductOutsideEntityList():
    with pytest.raises(ValueError):
        testRxnFlask.addReaction("A + B -> D", 1, reactantsDict={'A' : 1, 'B' : 1}, productsDict={'D' : 1})

def test_addReactionOfSameName():
    name = 'A + B -> C'
    with pytest.raises(NameError):
        testRxnFlask.addReaction(name, 1, reactantsDict={'A' : 1, 'B' : 1}, productsDict={'C' : 0})