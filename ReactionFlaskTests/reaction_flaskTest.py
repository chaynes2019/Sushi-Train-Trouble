import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import re

from reactionTest import Reaction

type SystemState = list[float]

class ReactionFlask():
  '''Implementing biochemical reaction network
  models can be tedious and error-prone.
  The loss of one negative sign or parameter could
  lead to dramatically different output, but the
  models are so complex that noticing such an
  error is almost out of the question.
  
  This class provides a solution to this problem
  by abstracting the implementation of ODE BRNs.
  Instead of explicitly writing down the differential
  equations, the user inputs the species involved in
  the network, the reactions occurring, and their
  related rate parameters.

  Then, a flexible framework for simulating and
  plotting the network's output is given. This
  facilitates rapid development of a BRN model
  and allows the user to evaluate its output over
  a range of parameters with ease.

  For more information on the member attributes
  of the ReactionFlask class, see the docstring for
  the class constructor. Enjoy!'''

  def __init__(self : ReactionFlask, 
               entityList : list[str], 
               components: dict[str : list[str]]) -> None:
    '''This constructor requires two inputs to get a
    Reaction Flask started.
    
    1. The entity list, which
    maintains a record of what entities or reaction
    species are actually in the reaction network.


    2. Components are merely subsets of the entity list
    that the user desires to have graphed together.
    PlotSystem() uses the components dictionary to
    graph subsets of the network on the same plot,
    but there is no restriction on what subsets can
    be specified.
    
    For example, every reaction species could be 
    plotted together or individually, or a user 
    could even partition the entity list into subsets
    that are mostly distinct but have one or two
    overlapping species in order to make it easier
    to compare the rise and fall of associated
    compounds.
    
    While the components feature makes graphing the
    output of the ReactionFlask much easier, it has
    no bearing on the simulation of the system itself.'''

    #The following commands check the entity list and
    # components inputs for any erroneous data types
    # before they are used for downstream operations
    self._entityList = entityList
    for entity in self._entityList:
      if type(entity) != str:
        raise(
          TypeError(
            "The entity list contains non-string elements"
          )
        )

    self._components = components

    for component in self._components:
      if type(component) != str:
        raise(
          TypeError(
            "A component key has a non-string value"
          )
        )

    #The concentrations list can act as a record of the
    # amount of each reaction species at certain points
    # in time
    self._concentrations = np.zeros(len(self._entityList))

    #The initial condition list provides a (mutable) vessel
    # for the user to input an initial value for the system
    self._initialCondition = None

    #The reactions dictionary keeps track of the Reaction
    # objects that have been added to the system. A user
    # can easily get a Reaction object by inputting the
    # reaction's name.
    self._reactions = {}

    #This flag prevents the system from being run before
    # it is properly initialized with an initial value
    self._concentrationsInitialized = False

    #This variable stores the output from scipy.integrate's
    # solveIVP function when it is run in reactionDeriv()
    self.latestSimulationOutput = None

  def setInitialCondition(self : ReactionFlask, 
                          y0 : SystemState) -> None:
    '''This function takes in an initial condition for
    the biochemical reaction network in the form of a
    list of values, each one corresponding to a species
    (that is, an entity) in the reaction network and
    fixes that as a value for runSystem() to automatically
    use upon being called.
    
    Note: this could be called multiple times on the same
    ReactionFlask within a loop; for example, if one was
    interested in seeing how changing the initial value
    of one of the system's species changed the outcome,
    one could call this function with the new initial
    condition within each loop, without having to create
    a new ReactionFlask object.'''

    #Prevent y0 from having type None
    if y0 == None:
      raise(
        TypeError(
          "y0 has type None"
        )
      )

    #Because the ReactionFlask's concentrations list is
    # initialized based on the length of the entity list
    # a mismatch between the lengths of these vectors
    # implies that an entity/reaction species has been
    # forgotten
    if len(y0) != len(self._concentrations):
      raise ValueError("y0 and the entity list have different dimensions")

    self._initialCondition = np.array(y0)

    #While the structure of the differential equations
    # negative values, the biological interpretation of
    # the biochemical reaction network favors enforcing
    # an initial condition with non-negative values, to
    # help avoid error
    for val in y0:
      if val == None:
        raise(
          TypeError(
            "y0 has a None value"
          )
        )

      if val < 0:
        raise(
          ValueError(
            "y0 has been initialized with negative values"
          )
        )
    
    #Map y0 values into the concentrations list:
    # this prepares the ReactionFlask to be run w/
    # runSystem()
    self._concentrations = [y0[k] for k in range(len(y0))]

    #This flag ensures that the system is not run
    # before being initialized
    self._concentrationsInitialized = True

  def resetInitialCondition(self : ReactionFlask) -> None:
    '''This function remaps the initial condition to the
    ReactionFlask's concentration list in order to
    reinitialize the system. This is useful to call when
    one is running the system repeatedly with the same
    initial conditions, perhaps changing the value of a
    parameter on each iteration.'''

    try:
      if self._initialCondition == None:
        raise(
          TypeError(
            "Initial condition has not yet been set"
          )
        )
    except(ValueError):

      #Remap y0 values, stored in _initialCondition,
      # to the system's concentrations list
      self._concentrations = self._initialCondition

      #Reset this flag so the system can be rerun
      self._concentrationsInitialized = True

  def addReaction(self : Reaction, 
                  name : str, 
                  rxnK : float, 
                  reactantsDict : dict[str : float], 
                  productsDict : dict[str : float]) -> None:
    '''This function adds a new Reaction object
    to the ReactionFlask's reactions dictionary.
    It first performs input checking before
    instantiating the new object.'''
    
    #The list of reactant names is generated
    # from the keys of the reactants dictionary
    reactants = list(reactantsDict.keys())

    #The list of product names is generated
    # from the keys of the products dictionary
    products = list(productsDict.keys())

    #Check the reactant and product names to
    # make sure that they are within the scope 
    # of the reaction network
    for reactant in reactants:
      if reactant not in self._entityList:
        raise(ValueError(f"There is a reaction with {reactants} as reactants, but {reactant} is not in the entity list"))

    for product in products:
      if product not in self._entityList:
        raise(ValueError(f"There is a reaction with {products} as products, but {products} is not in the entity list"))
    
    #Check to see if Reaction name is 
    # already used by an existing reaction
    if name in list(self._reactions.keys()):
      raise(
        NameError(
          f"There is already a reaction of the name \"{name}\""
        )
      )

    #Add the new Reaction object
    self._reactions[name] = Reaction(self._entityList, name, rxnK, reactantsDict, productsDict)

  def computeReactionRates(self : ReactionFlask, 
                           y : SystemState) -> list[float]:
    '''Each Reaction object includes a function to
    compute its reaction rate according to mass-action
    kinetics. Now, this function will iterate through
    the Reaction objects of the ReactionFlask and call
    this function for each one. computeReactionRates()
    is called at each timestep during reactionDeriv(),
    and it forms an essential step of simulating the
    reaction network.'''

    #To iterate through the reaction objects themselves
    # one must use the reactions dictionary values
    reactionObjects = self._reactions.values()

    #Initialize an empty numpy array to hold the
    # computed reaction rates 
    reactionRates = np.zeros(len(reactionObjects))

    #Call each Reaction object's computeRxnRate function
    # and store the result. By using enumerate, we ensure
    #that we iterate through the reactionObjects in a 
    # repeatable order
    for k, rxn in enumerate(reactionObjects):
      reactionRates[k] = rxn.computeRxnRate(self._entityList, y)
    
    return reactionRates

  def reactionDeriv(self : ReactionFlask, 
                    t : float, 
                    y : SystemState):
    '''This function forms the heart of the
    ReactionFlask class. It uses the abstracted
    Reaction objects to form, piece-by-piece,
    the derivative of the system. It computes
    all of the reaction rates and then uses
    the changes in reactants and products reported
    by each reaction to compute the contribution
    of that reaction to the derivative vector.
    This function can then be called by an
    ODE solver, like scipy.integrate's solveIVP().'''

    #0: Check inputs for non-negativity
    if t < 0:
      raise(
        ValueError(
          "Input time has a negative value"
        )
      )
    
    for k, yVal in enumerate(y):
      #If the value is negative and large enough
      # that it is unlikely to have originated
      # from float arithmetic errors, then there's
      # an issue  
      if yVal < 0 and abs(yVal) > 0.01:
        raise(
          ValueError(
            f"Element {k} of the y vector is negative = {yVal}"
          )
        )

    #1: Set up receptacle for derivative values
    derivative = np.zeros(len(self._concentrations))

    #2: Compute the 1/s rates of each reaction,
    #based on current concentrations
    rxnRates = self.computeReactionRates(y)

    #3: For each reaction, go to each reactant
    #  and add the change (# / s) to its respective derivative

    #To iterate through the reaction objects themselves
    # one must use the reactions dictionary values
    reactionObjects = self._reactions.values()

    #Loop through the Reaction objects and compute 
    # its contribution to the reaction derivative
    for rxnIdx, rxn in enumerate(reactionObjects):
      #Obtain the changes in reactants and products
      reactantChanges = rxn._changes[0]
      productChanges = rxn._changes[1]

      #Loop through each reactant and compute the
      # reaction's mass-action term in the reactant's
      # differential equation
      for reactant in rxn._reactants:
        idx = self._entityList.index(reactant)

        #The mass-action term in the reactant's differential
        # equation will be of the form reactantChange * reactionRate 
        derivativeChange = -reactantChanges[idx] * rxnRates[rxnIdx]

        #By using the idx of the reactant, the mass-action term
        # can be added to the differential equation for the 
        # reactant, specifically 
        derivative[idx] = derivative[idx] + derivativeChange

      #Loop through each product and compute the
      # reaction's mass-action term in the product's
      # differential equation
      for product in rxn._products:
        idx = self._entityList.index(product)

        #The mass-action term in the product's differential
        # equation will be of the form productChange * reactionRate 
        derivativeChange = productChanges[idx] * rxnRates[rxnIdx]

        #By using the idx of the product, the mass-action term
        # can be added to the differential equation for the 
        # product, specifically 
        derivative[idx] = derivative[idx] + derivativeChange

    return derivative

  def runSystem(self : ReactionFlask, 
                timeEndpoint : list[float]) -> None:
    '''This function gets the simulation of the
    reaction network underway. After ensuring that
    the system has been properly initialized, it calls
    scipy.integrate's solveIVP function to simulate
    the BRN forward in time, calling reactionDeriv()
    as the system derivative, until a specified time
    endpoint.
    
    At this point, after assembling the ReactionFlask
    by using the abstracted tools provided above, calling
    this function should be easy!'''

    #Ensure timeEndpoint >= 0
    if timeEndpoint < 0:
      raise(
        ValueError(
          "Specified time endpoint is negative!"
        )
      )

    #Check to see if system has been properly initialized
    if self._concentrationsInitialized == True:
      #1: Reset flag to show that concentrations are no longer
      #at initial values
      self._concentrationsInitialized = False

      #2: Run initial value problem for specified time length
      self.latestSimulationOutput = solve_ivp(self.reactionDeriv, [0, timeEndpoint], self._concentrations)

    else:
      raise(
        AttributeError(
          "Reaction flask has not yet been initialized. Please initialize before running."
        )
      )

  def plotSystem(self, widthSpacing = 0.5, heightSpacing = 0.5, leftEdgeOfPlots = None, rightEdgeOfPlots = None, numberRows = 1):
    #0: Grab time values from simulation
    tVals = self.latestSimulationOutput["t"]

    #There's no list of axes if there's only one component!
    if len(self._components) > 1:
      [figure, axes] = plt.subplots(1, len(self._components))

    else:
      fig, axes = plt.subplots()
      axes = [axes]


    #1. Loop through entity strings, getting their official index
    # with enumerate for use in getting timeseries data
    for k, entity in enumerate(self._entityList):

      entityVals = self.latestSimulationOutput["y"][k]

      for j, component in enumerate(self._components.values()):
        if entity in component:
          axes[j].plot(tVals, entityVals, label = entity)

      #1a. Use the index from above to procure timeseries data
      # for the entity


      #1b. Plot the entityVals vs. tVals, providing the entity
      #name/string as a plot label
      #plt.plot(tVals, entityVals, label = entity)

    for j, component in enumerate(self._components):
      #2a. Title = "System vs. Time"
      axes[j].set_title(f"{component} vs. Time")

      #2b. Let the x-axis be "Time"
      axes[j].set_xlabel("Time")

      #2c. Let the y-axis be "System Concentration Values"
      axes[j].set_ylabel(f"{component} Concentration Values")

      #2c 1/2: Set bottom of y-axis to 0.1, just to help make sure my comparisons are correct
      axes[j].set_ylim(bottom = -0.1)

      #2d: Activate Legend
      axes[j].legend()

    #3a. Give a little wiggle room around the component plots.
    plt.subplots_adjust(hspace = heightSpacing, wspace = widthSpacing, left = leftEdgeOfPlots, right = rightEdgeOfPlots)

    #3b. Display plot.
    plt.show()

  def getFinalValueOfVariable(self, variable):
    #1. Retrieve Index of variable in the simulation output
    # through its index in the entity list
    variableIdx = self._entityList.index(variable)

    #2. Retrieve timeseries in question from simulation output
    timeseries = self.latestSimulationOutput["y"][variableIdx]

    #3. Return final value
    return timeseries[-1]

  def modifyReaction(self, rxnName, newRxnRate = None, newRxnReactantsDict = None, newRxnProductsDict = None):
    reactionInQuestion = self._reactions[rxnName]

    if newRxnRate != None:
      reactionInQuestion._rxnK = newRxnRate

    if newRxnReactantsDict != None:
      reactionInQuestion._reactantsDict = newRxnReactantsDict
      reactionInQuestion.computeReactionSpecies()
      reactionInQuestion.computeChanges(self._entityList)

    if newRxnProductsDict != None:
      reactionInQuestion._productsDict = newRxnProductsDict
      reactionInQuestion.computeReactionSpecies()
      reactionInQuestion.computeChanges(self._entityList)
