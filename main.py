from VirusSpreadingSimulation.VirusSpredingSimulation import VirusSpreadingSimulation

xlim = 100.
ylim = 100.

if __name__ == "__main__":

    simulation = VirusSpreadingSimulation(population=50, xlim=xlim, ylim=ylim)

    while simulation.day < 100:
        simulation.step()
