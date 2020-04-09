from VirusSpreadingSimulation.VirusSpredingSimulation import VirusSpreadingSimulation

pop = 500
xlim = 100.
ylim = 100.
radius_infection = 1.

if __name__ == "__main__":

    simulation = VirusSpreadingSimulation()
    simulation.get_plot(save=False, save_name='quarantine_and_containment_result')

    simulation = VirusSpreadingSimulation(isolation_threshold=.99)
    simulation.get_plot(save=False, save_name='only_quarantine_result')

    simulation = VirusSpreadingSimulation(isolation_threshold=.99, quarantine_zone=False)
    simulation.get_plot(save=False, save_name='no_measure_result')
