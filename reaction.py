import numpy as np

class Reaction():
  def __init__(self, entityList, name, rxnK, reactantsDict = {}, productsDict = {}):
    self._reactantsDict = reactantsDict
    self._reactants = list(self._reactantsDict.keys())

    self._productsDict = productsDict
    self._products = list(self._productsDict.keys())

    self._netChanges = self.computeNetChanges(entityList)

    self._rxnK = rxnK

    self._name = name

  def computeNetChanges(self, entityList):
    #1: Set up vectors that are the length of the entity list
    # to catalog the number of the different reactants and
    # products present
    productsMultiplicities = np.zeros(len(entityList))
    reactantsMultiplicities = np.zeros(len(entityList))

    #2: Search through products and catalog the multiplicities thereof
    # in the appropriate place in the productsMultiplicities vector
    for product in self._products:
      entityIdx = entityList.index(product)

      #2b: Use the productsDict and the product key to look up the
      # multiplicity in the reaction
      productsMultiplicities[entityIdx] = self._productsDict[product]

    #3: Repeat step 2 for reactants
    for reactant in self._reactants:
      entityIdx = entityList.index(reactant)

      #3b: Use the reactantsDict and the reactant key to look up the
      # multiplicity in the reaction
      reactantsMultiplicities[entityIdx] = self._reactantsDict[reactant]

    #4: Now, their subtraction will yield the net difference in materials
    netChanges = productsMultiplicities - reactantsMultiplicities

    return netChanges

  def computeRxnRate(self, entityList, y):
    rxnRate = 1

    for reactant in self._reactantsDict:
      idx = entityList.index(reactant)

      if(self._name == "Reaction 30: M_empty + g * phi -> M_empty BAD"):
        reactantConcentration = 0
      else:
        reactantConcentration = y[idx]
      reactantMultiplicity = self._reactantsDict[reactant]

      rxnRate *= (reactantConcentration ** reactantMultiplicity)

      if(self._name == "Reaction 30: M_empty + g * phi -> M_empty BAD"):
        print(f"Reactant Concentration Product = {(reactantConcentration ** reactantMultiplicity)}")

    rxnRate *= self._rxnK

    if (self._name == "Reaction 30: M_empty + g * phi -> M_empty"):
      return int(rxnRate)
    
    else:
      return rxnRate