from VirusSpreadingSimulation.VirusSpredingSimulation import VirusSpreadingSimulation

pop = 500
xlim = 100.
ylim = 100.
radius_infection = 1.

if __name__ == "__main__":

    #simulation = VirusSpreadingSimulation(population=pop, xlim=xlim, ylim=ylim)
    simulation = VirusSpreadingSimulation()
    #simulation += 50
    simulation.get_plot(save=True)
    # simulation.get_plot()
