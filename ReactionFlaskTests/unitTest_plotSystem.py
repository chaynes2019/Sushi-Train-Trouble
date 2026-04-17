from reaction_flaskTest import ReactionFlask

#Instantiate test flask with test reaction
testFlask = ReactionFlask(['A', 'B', 'C'], components = {"Main" : ['A', 'B', 'C']})
testFlask.addReaction("A + B -> C", 1, reactantsDict = {'A' : 1, 'B' : 1}, productsDict = {'C' : 1})

'''Test null initial condition'''
testFlask.setInitialCondition([0, 0, 0])

testFlask.runSystem(50)

testFlask.plotSystem()

'''Test neutral: only one needed reactant'''
testFlask.resetInitialCondition()
testFlask.setInitialCondition([1, 0, 0])

testFlask.runSystem(50)

testFlask.plotSystem()

'''Test neutral: only product'''
testFlask.resetInitialCondition()
testFlask.setInitialCondition([0, 0, 1])

testFlask.runSystem(50)

testFlask.plotSystem()

'''Test nominal: both needed reactants'''
testFlask.resetInitialCondition()
testFlask.setInitialCondition([1, 1, 0])

testFlask.runSystem(50)

testFlask.plotSystem()

'''Test nominal: reactants start at different values'''
testFlask.resetInitialCondition()
testFlask.setInitialCondition([2, 1, 0])

testFlask.runSystem(50)

testFlask.plotSystem()