import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from reaction import Reaction

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
    no bearing on the simulation of the system itself.
    
    Inputs:

    self (ReactionFlask): the ReactionFlask that
    is being constructed

    entityList (list of strings): a list of the
    names of the reaction species as strings.
    Generally, abbreviated forms of the names
    of reaction species or single letters will
    work best.

    components (dictionary of lists): a dictionary
    mapping names of system "components" to the
    subsets of reaction species forming those
    components. These will then be used to define
    which subsets of reaction species are graphed
    together when plotSystem() is run

    Outputs:

    None

    '''

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
    a new ReactionFlask object.

    Inputs:
    
    self (ReactionFlask) : the ReactionFlask object whose
    initial condition is being set

    y0 (list of floats/SystemState) : the initial state of
    the system, formatted perfectly for use in the solveIVP
    function of scipy.integrate in runSystem(). Each element
    is the initial concentration of the corresponding
    reaction species. Units are not specified here, so these
    values can take on units based on model parametrization.

    Outputs:

    None
    '''

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
    parameter on each iteration.
    
    Inputs:

    self (ReactionFlask) : the ReactionFlask whose
    initial condition is being reset

    Outputs:

    None
    '''

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
    instantiating the new object.
    
    Inputs:

    self (ReactionFlask) : the ReactionFlask that the
    Reaction object is being added to

    name (str) : the namestring of the reaction. It
    is often helpful to write out the reaction equation
    for the reaction being modeled, and then, one uses 
    that reaction equation as the reaction name.
    For example, nameRxn = "A + B -> C"

    rxnK (float) : the reaction rate constant. This
    value defines how fast the reaction proceeds for
    every given concentration tuple of the reactants.
    In the mass-action kinetics term corresponding to
    the reaction "A + B -> ?", k_rxn * [A]^a * [B]^b,
    the reaction rate constant is k_rxn. The value of
    this constant is related to the thermodynamic
    properties of the reaction and how spontaneously
    it occurs. If the change in Gibbs free energy between
    the reactants and products is known, this value can
    be explicitly calculated as 
    
    k_rxn = e^(-(deltaGibbs) / RT)


    reactantsDict (dictionary, str -> int/float) : this
    dictionary contains the namestrings of the reactants
    as keys and the multiplicities of those reactants as
    values. Thus, for example, "2A + B -> 3C" would become

    {'A' : 2, 'B' : 1}


    productsDict (dictionary, str -> int/float) : this
    dictionary contains the namestrings of the products
    as keys and the multiplicities of those products as
    values. Thus, for example, "2A + B -> 3C" would become

    {'C' : 3}

    Outputs:

    None

    '''
    
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
    reaction network.
    
    Inputs:

    self (ReactionFlask) : the ReactionFlask that is
    having its Reaction objects compute their reaction
    rates

    y (list of floats/SystemState) : this is a list of
    the concentrations of reaction species, maintained
    in the same order as the entity list and initial
    condition. This provides an input for reactionDeriv()
    and enables the Reaction objects to accurately
    compute reaction rates based on the species
    concentration values at the present moment.

    Outputs:

    reactionRates (list of floats) : the reaction rates
    computed from each of the system's reactions, in the
    same order that the Reaction objects are stored in
    the reactants dictionary
    
    '''

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
                    y : SystemState) -> list[float]:
    '''This function forms the heart of the
    ReactionFlask class. It uses the abstracted
    Reaction objects to form, piece-by-piece,
    the derivative of the system. It computes
    all of the reaction rates and then uses
    the changes in reactants and products reported
    by each reaction to compute the contribution
    of that reaction to the derivative vector.
    This function can then be called by an
    ODE solver, like scipy.integrate's solveIVP().
    
    Inputs:

    self (ReactionFlask) : the ReactionFlask that
    is computing its reaction derivative


    t (float) : the current time elapsed during
    the simulation

    y (list of floats/SystemState) : this is a list of
    the concentrations of reaction species, maintained
    in the same order as the entity list and initial
    condition. This provides an input for reactionDeriv()
    and enables the Reaction objects to accurately
    compute reaction rates based on the species
    concentration values at the present moment.

    Outputs:

    derivative (list of floats) : the derivative of the
    biochemical reaction network, expressed as an
    Ordinary Differential Equation. This value can be
    used in an ODE Solver, such as scipy.integrate's
    solveIVP function, to simulate the system forward
    in time.
    
    '''

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
                timeEndpoint : float) -> None:
    '''This function gets the simulation of the
    reaction network underway. After ensuring that
    the system has been properly initialized, it calls
    scipy.integrate's solveIVP function to simulate
    the BRN forward in time, calling reactionDeriv()
    as the system derivative, until a specified time
    endpoint.
    
    At this point, after assembling the ReactionFlask
    by using the abstracted tools provided above, calling
    this function should be easy!

    Inputs:
    
    self (ReactionFlask) : the ReactionFlask that is
    running a simulation of its biochemical reaction
    network

    timeEndpoint (float) : the time elapsed during
    simulation when the ODE solver should stop
    integrating the system forward and the simulation
    ends

    Outputs:

    None

    '''

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

  def plotSystem(self : ReactionFlask, 
                 widthSpacing : float = 0.5, 
                 heightSpacing : float = 0.5, 
                 leftEdgeOfPlots : float = None, 
                 rightEdgeOfPlots : float = None) -> None:
    '''After running the system, the user can call
    this function to automatically plot the results.
    It will generate #components plots, each with time
    on the x-axis and concentration on the y-axis. By
    default, these are sized so that several can be fit
    into one row, but the sizing of each can be edited
    using the parameters.

    The entities listed together in each component will
    be plotted together, with a legend automatically
    generated.

    Inputs:

    self (ReactionFlask): the ReactionFlask that
    has simulated the system with the reaction
    species and contains data on its temporal
    dynamics
    
    widthSpacing (float): changes the horizontal spacing
    between consecutive plots
    
    heightSpacing (float): changes the vertical spacing
    between consecutive plots
    
    leftEdgeOfPlots (float): defines the position on the
    left-right axis of the canvas where the plot array
    will begin
    
    rightEdgeOfPlots (float): defines the position on the
    left-right axis of the canvas where the plot array
    will end

    Outputs:

    None
    '''

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

  def getFinalValueOfSpecies(self : ReactionFlask, 
                             species : str) -> float:
    '''This function retrieves the final value
    for any system reaction species: the user
    simply inputs the namestring of the species
    and the function does the rest.
    
    This function can be particularly helpful for
    steady state analyses.

    Inputs:
    
    self (ReactionFlask): the ReactionFlask that
    has simulated the system with the reaction
    species and contains data on its temporal
    dynamics

    species (str): the namestring of the reaction
    species whose final value is being accessed

    Outputs:

    finalValue (float) : the value of the specified
    reaction species at the end of the simulation
    '''

    #Ensure that system has been run
    
    if self.latestSimulationOutput == None:
      raise(
        AttributeError(
          "Reaction flask has not yet been run"
        )
      )

    #1. Retrieve Index of species in the simulation output
    # through its index in the entity list
    try:
      speciesIdx = self._entityList.index(species)
    except(ValueError):
      raise(
        ValueError(
          "Requested reaction species is not in the entity list"
        )
      )

    #2. Retrieve timeseries in question from simulation output
    timeseries = self.latestSimulationOutput["y"][speciesIdx]

    #3. Return final value
    finalValue = timeseries[-1]
    return finalValue

  def modifyReaction(self : ReactionFlask, 
                     rxnName : str, 
                     newRxnRate : float = None, 
                     newRxnReactantsDict : dict[str : float] = None, 
                     newRxnProductsDict : dict[str : float] = None) -> None:
    '''After adding a Reaction object to a ReactionFlask,
    it is sometimes necessary to modify the Reaction object
    at runtime. This is particularly the case when one is
    simulating the ReactionFlask for many different values
    of a parameter, such as a reaction rate. In this case,
    modifyReaction() would be called once per iteration to
    update the reaction to have the new parameter value.
    
    It can also be used to update the reactants and products
    that are involved in a reaction as well as their
    multiplicities. It is hoped that this offers a flexible
    set of options for experimentation.
    
    Inputs:

    self (ReactionFlask) : the ReactionFlask containing the 
    Reaction object that is being modified

    rxnName (str) : the namestring of the Reaction object and
    the key for that Reaction object in the _reactions dictionary

    newRxnRate (float) : the new reaction rate that is desired for
    the Reaction in question. If none is specified, the reaction
    rate of the Reaction will be left unchanged.

    newRxnReactantsDict (dictionary, str -> float) : the new
    reactants dictionary that is desired for the Reaction
    in question. It can be different than the original in
    reactants involved, their multiplicities, or both. If
    none is specified, the reactants dictionary of the
    Reaction will be left unchanged.

    newRxnProductsDict (dictionary, str -> float) : the new
    products dictionary that is desired for the Reaction
    in question. It can be different than the original in
    products involved, their multiplicities, or both. If
    none is specified, the products dictionary of the
    Reaction will be left unchanged.

    Outputs:

    None
    
    '''

    #Obtain the Reaction object that is to be
    # modified by using its namestring as a key
    # in the _reactions dictionary 
    reactionInQuestion = self._reactions[rxnName]


    #If no new reaction rate has been specified,
    # leave the original. If one has been given, 
    # replace the original with the new. 
    if newRxnRate != None:
      reactionInQuestion._rxnK = newRxnRate

    #If no new reactants dictionary has been specified,
    # leave the original. If one has been given,
    # replace the original with the new 
    if newRxnReactantsDict != None:
      reactionInQuestion._reactantsDict = newRxnReactantsDict

      #Because the reaction species and changes in
      # reactants involves the reactantsDict, these
      # need to be recomputed before the ReactionFlask
      # can be run again 
      reactionInQuestion.computeReactionSpecies()
      reactionInQuestion.computeChanges(self._entityList)

    #If no new products dictionary has been specified,
    # leave the original. If one has been given,
    # replace the original with the new 
    if newRxnProductsDict != None:
      reactionInQuestion._productsDict = newRxnProductsDict

      #Because the reaction species and changes in
      # products involves the productsDict, these
      # need to be recomputed before the ReactionFlask
      # can be run again 
      reactionInQuestion.computeReactionSpecies()
      reactionInQuestion.computeChanges(self._entityList)
