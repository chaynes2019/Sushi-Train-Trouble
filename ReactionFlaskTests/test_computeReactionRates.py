from reaction_flaskTest import ReactionFlask

def test_NullRxnFlask():
    nullRxnFlask = ReactionFlask(['A', 'B', 'C'], components = {'Main' : ['A', 'B', 'C']})

    rxnRates = nullRxnFlask.computeReactionRates([1, 1, 1])
    rxnRatesList = [rxnRates[k] for k in range(len(rxnRates))]

    assert len(rxnRates) == 0
    assert rxnRatesList == []

def test_neutralRxnFlask():
    neutralRxnFlask = ReactionFlask(['A', 'B', 'C'], components = {'Main' : ['A', 'B', 'C']})
    neutralRxnFlask.addReaction("A + B -> A + B", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'A' : 1, 'B' : 1})

    rxnRates = neutralRxnFlask.computeReactionRates([1, 1, 1])
    rxnRatesList = [float(rxnRates[k]) for k in range(len(rxnRates))]

    assert rxnRatesList == [1]

def test_nominalRxnFlask():
    nominalRxnFlask = ReactionFlask(['A', 'B', 'C'], components = {'Main' : ['A', 'B', 'C']})
    nominalRxnFlask.addReaction("A + B -> C", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'C' : 1})
    nominalRxnFlask.addReaction("2A + B -> 2C", 1, reactantsDict = {'A' : 2, 'B' : 1}, productsDict = {'C' : 2})

    noReactants = [0, 0, 0]
    assert (nominalRxnFlask.computeReactionRates(noReactants) == [0, 0]).all()

    onlyOneReactant = [1, 0, 0]
    assert (nominalRxnFlask.computeReactionRates(onlyOneReactant) == [0, 0]).all()

    bothReactantsOne = [1, 1, 0]
    assert (nominalRxnFlask.computeReactionRates(bothReactantsOne) == [1, 1]).all()

    firstReactantTwo = [2, 1, 0]
    assert (nominalRxnFlask.computeReactionRates(firstReactantTwo) == [2, 4]).all()

    secondReactantTwo = [1, 2, 0]
    assert (nominalRxnFlask.computeReactionRates(secondReactantTwo) == [2, 2]).all()

    productsHaveNoEffect = [1, 2, 1000000]
    assert (nominalRxnFlask.computeReactionRates(productsHaveNoEffect) == [2, 2]).all()
