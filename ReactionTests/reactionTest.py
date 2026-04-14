import numpy as np

class Reaction():
  def __init__(self : Reaction, 
               entityList : list[str], 
               name : str, 
               rxnK : float, 
               reactantsDict : dict[str, float] = {}, 
               productsDict : dict[str, float] = {}) -> None:
    '''The Reaction class constructor takes in the entity list
     of its containing ReactionFlask in order to compute its
     changes to system variables in context. Then, it is
     parameterized with the information contained in a
     normal chemical reaction equation:

     reactants -> products'''
    
    #This section evaluates the inputs to ensure they
    # are of the expected type.

    for species in entityList:
      if type(species) != str:
        raise(ValueError("The Entity List has a non-string element"))
      
    if type(name) != str:
        raise(
          ValueError(
            "The reaction name is not a string!"
            )
          )
    
    if type(rxnK) != float and type(rxnK) != int:
      raise(
        ValueError(
          "The reaction rate (i.e., reaction constant) is not a float!"
        )
      )
    
    if type(reactantsDict) != dict[str, float] and type(reactantsDict) != dict[str, int]:
      raise(
        ValueError(
          "The reactants dictionary has keys which are not strings or values which are not floats"
        )
      )
    
    if type(productsDict) != dict[str, float] and type(productsDict) != dict[str, int]:
      raise(
        ValueError(
          "The products dictionary has keys which are not strings or values which are not floats"
        )
      )
    
    #The Reactants Dictionary stores the information
    # about the reactants involved and their relative
    # multiplicities on a species-by-species basis

    self._reactantsDict = reactantsDict

    #The Products Dictionary stores the information
    # about the products involved and their relative
    # multiplicities on a species-by-species basis

    self._productsDict = productsDict

    #The _changes attribute provides two vectors,
    # one that shows how the reactants change per
    # reaction and another that shows how the
    # products change per reaction.

    self._changes = None

    self._reactants = None
    self._products = None

    self.computeReactionSpecies()
    
    self.computeChanges(entityList)

    #The reaction constant provides information
    # about how spontaneous the reaction is and
    # its energetics. If the reactants are placed
    # together in space and they react readily,
    # the reaction is thermodynamically favorable,
    # and the reaction constant is accordingly high.
    # If, when placed together, the reactants react
    # slowly with one another, the reaction is
    # thermondynamically unfavorable, and the
    # reaction constant will be small.

    self._rxnK = rxnK

    #This namestring allows the specific reaction
    # to be obtained from ReactionFlasks _reactions
    # dictionary. This enables the precise use of
    # reaction information within ReactionFlask.

    self._name = name

  def computeChanges(self : Reaction,
                     entityList : list[str]) -> None:
    '''Once a reaction has been instantiated with 
    reactant and product dictionaries, this function
    computes and returns the reaction's changes to the 
    systemvariables in a format that is very useful for
    ReactionFlask's reactionDeriv() function'''

    #1: Set up vectors that are the length of the entity list
    # to catalog the number of the different reactants and
    # products present

    productsMultiplicities = np.zeros(len(entityList))
    reactantsMultiplicities = np.zeros(len(entityList))


    #2: Search through products and catalog the multiplicities thereof
    # in the appropriate place in the productsMultiplicities vector

    for product in self._products:
      entityIdx = entityList.index(product)

      #2a: Use the productsDict and the product key to look up the
      # multiplicity in the reaction
      productsMultiplicities[entityIdx] = self._productsDict[product]

    #3: Repeat step 2 for reactants

    for reactant in self._reactants:
      entityIdx = entityList.index(reactant)

      #3a: Use the reactantsDict and the reactant key to look up the
      # multiplicity in the reaction

      reactantsMultiplicities[entityIdx] = self._reactantsDict[reactant]

    self._changes = (reactantsMultiplicities, productsMultiplicities)
  
  def computeReactionSpecies(self : Reaction) -> None:
    '''This function parses the reactants and
    products dictionaries to provide lists of
    the namestrings for the reactants and
    products involved in the reaction.'''

    #1. List of the string name of reactants is generated
    #from the keys of the reactants dictionary

    self._reactants = list(self._reactantsDict.keys())

    #2. List of the string name of products is generated
    #from the keys of the products dictionary

    self._products = list(self._productsDict.keys())

  def computeRxnRate(self : Reaction, 
                     entityList : list[str], 
                     y : list[float]) -> float:
    '''Given a reaction aA + bB -> cC, the rate
      at which the reaction happens over time
      is given by mass-action kinetics, assuming
      that the system is well-mixed. In this case,
      the reaction rate is
      
      rate = k_rxn * [A]^a * [B]^b
      
      This function computes and returns the value
      of rate.'''
    
    #1. While we could also initialize the variable to
    # return, rxnRate, with the first concentration,
    # multiplicity exponentiation term in the product,
    # we find it easiest to understand if we start the
    # product computed in this function with 1, much
    # like one starts a sum with a summing variable
    # set equal to 0.

    rxnRate = 1

    #2. Compute reaction rate product (apart from rate
    # constant) by iterating through each reactant
    # and multiplying by its contribution 

    for reactant in self._reactantsDict:
      #2a. Retrieve index of reactant in entity list b/c
      # entityList is a list of reactant names

      idx = entityList.index(reactant)

      #2b. Use this index to get the concentration from y,
      # because y has the same shape as entityList

      reactantConcentration = y[idx]

      #2c. Get the multiplicity of the reactant in the reaction
      # using the value given in reactantsDict

      reactantMultiplicity = self._reactantsDict[reactant]

      #2d. For aA + bB -> P, Reaction Rate = k_rxn * [A]^a * [B]^b
      # So, first form product of [A]^a, [B]^b, etc.

      rxnRate *= (reactantConcentration ** reactantMultiplicity)

      '''if(self._name == "Reaction 30: M_empty + g * phi -> M_ePhi or M_empty"):
        print(f"Reactant Concentration Product = {(reactantConcentration ** reactantMultiplicity)}")'''

    #3. Now, multiply resultant product by k_rxn, the reaction constant

    rxnRate *= self._rxnK

    return rxnRate