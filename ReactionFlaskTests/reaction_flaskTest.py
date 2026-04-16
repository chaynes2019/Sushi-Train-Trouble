import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from reactionTest import Reaction


class ReactionFlask():
  def __init__(self, entityList, components):
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

    self._concentrations = np.zeros(len(self._entityList))
    self._initialCondition = None

    #Changing this to be a dictionary to make it more readily searchable
    self._reactions = {}

    self._concentrationsInitialized = False

    self.latestSimulationOutput = None

  def setInitialCondition(self : ReactionFlask, 
                          y0 : list[float]) -> None:
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

    #Remap y0 values, stored in _initialCondition,
    # to the system's concentrations list
    self._concentrations = self._initialCondition

    #Reset this flag so the system can be rerun
    self._concentrationsInitialized = True

  def addReaction(self, name, rxnK, reactantsDict, productsDict):
    reactants = list(reactantsDict.keys())
    products = list(productsDict.keys())

    for reactant in reactants:
      if reactant not in self._entityList:
        raise(ValueError(f"There is a reaction with {reactants} as reactants, but {reactant} is not in the entity list"))

    for product in products:
      if product not in self._entityList:
        raise(ValueError(f"There is a reaction with {products} as products, but {products} is not in the entity list"))
      
    self._reactions[name] = Reaction(self._entityList, name, rxnK, reactantsDict, productsDict)

  def computeReactionRates(self, y):
    #Because the reactions member is now a dictionary,
    #to iterate through the reactions, one now must
    #use the values
    reactionObjects = self._reactions.values()

    reactionRates = np.zeros(len(reactionObjects))
    for k, rxn in enumerate(reactionObjects):
      reactionRates[k] = rxn.computeRxnRate(self._entityList, y)

    return reactionRates

  def reactionDeriv(self, t, y):
    #1: Set up receptacle for derivative values
    derivative = np.zeros(len(self._concentrations))

    #2: Compute the 1/s rates of each reaction,
    #based on current concentrations
    rxnRates = self.computeReactionRates(y)

    #3: For each reaction, go to each reactant
    #  and add the change (# / s) to its respective derivative

    #Because the reactions member is a dictionary,
    #to iterate through the reactions, one now must
    #use the values
    reactionObjects = self._reactions.values()

    for rxnIdx, rxn in enumerate(reactionObjects):
      if(rxn._name == "Reaction 30: M_empty + g * phi -> M_ePhi"):
        print(f"Reaction Rate = {rxn.computeRxnRate(self._entityList, y)}")

      reactantChanges = rxn._changes[0]
      productChanges = rxn._changes[1]

      for reactant in rxn._reactants:
        idx = self._entityList.index(reactant)

        #The change in reactants reported in the first element of the changes vector,
        # which can then be used, through the index obtained from the entityList,
        # to have the right reaction rate in the concentration
        derivativeChange = -reactantChanges[idx] * rxnRates[rxnIdx]
        if(rxn._name == "Test Rxn: A + B -> C"):
          print(f"Reactant Changes = {-reactantChanges[idx]}")
        derivative[idx] = derivative[idx] + derivativeChange

      for product in rxn._products:
        idx = self._entityList.index(product)

        #The net change is reported in the netChanges vector, which can then
        # be used, through the index obtained from the entityList, to have
        # the right reaction rate in the concentration
        derivativeChange = productChanges[idx] * rxnRates[rxnIdx]
        if(rxn._name == "Test Rxn: A + B -> C"):
          print(f"Product Change = {productChanges[idx]}")
        derivative[idx] = derivative[idx] + derivativeChange

    return derivative

  def runSystem(self, timeEndpoint, fromSteadyState = False):
    if fromSteadyState:
      if self._concentrationsInitialized == True:
        #1: Reset flag to show that concentrations are no longer
        #at initial values
        self._concentrationsInitialized = False

        #2: Run initial value problem for specified time length
        #First, get the wound size for later, and then remove the
        #wound from the initial condition
        perturbationValue = self._concentrations[0]

        self._concentrations[0] -= perturbationValue

        simulationOutput = solve_ivp(self.reactionDeriv, [0, 1000], self._concentrations)

        newInitialCondition = [simulationOutput["y"][k][-1] for k in range(len(self._entityList))]

        if(abs(newInitialCondition[0]) > 0.00001):
          print("Something's wrong with the 0-indexed variable: it can't be a perturbation, because the system is changing it from zero")
        else:
          newInitialCondition[0] = perturbationValue

        newSimulationOutput = solve_ivp(self.reactionDeriv, [0, timeEndpoint], newInitialCondition)

        newSimulationOutput["t"] += 1000

        self.latestSimulationOutput = {}

        self.latestSimulationOutput["t"] = [simulationOutput["t"][i] for i in range(len(simulationOutput["t"]) - 1)] + [newSimulationOutput["t"][i] for i in range(len(newSimulationOutput["t"]))]

        timeSeriesOutputs = [[simulationOutput["y"][k][i] for i in range(len(simulationOutput["y"][k]) - 1)] + [newSimulationOutput["y"][k][i] for i in range(len(newSimulationOutput["y"][k]))] for k in range(len(self._entityList))]

        self.latestSimulationOutput["y"] = timeSeriesOutputs

    else:
      if self._concentrationsInitialized == True:
        #1: Reset flag to show that concentrations are no longer
        #at initial values
        self._concentrationsInitialized = False

        #2: Run initial value problem for specified time length
        self.latestSimulationOutput = solve_ivp(self.reactionDeriv, [0, timeEndpoint], self._concentrations)

      else:
        print("Please initialize reaction concentrations and try again.")

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
