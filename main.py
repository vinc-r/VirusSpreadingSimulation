from VirusSpreadingSimulation.VirusSpredingSimulation import VirusSpreadingSimulation

xlim = 4.
ylim = 4.

if __name__ == "__main__":

    simulation = VirusSpreadingSimulation(population=10, xlim=xlim, ylim=ylim)

    while simulation.day < 100:
        simulation.step()
