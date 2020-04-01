from VirusSpreadingSimulation.VirusSpredingSimulation import VirusSpreadingSimulation

xlim = 100.
ylim = 100.

if __name__ == "__main__":

    simulation = VirusSpreadingSimulation(population=500, xlim=xlim, ylim=ylim, radius_infection=1.)

    simulation.get_plot()
