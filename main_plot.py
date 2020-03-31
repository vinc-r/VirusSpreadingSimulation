import matplotlib.pyplot as plt
import matplotlib.animation as animation
from VirusSpreadingSimulation.VirusSpredingSimulation import VirusSpreadingSimulation

xlim = 100.
ylim = 100
density = 0.1
pop = int(density * xlim * ylim)

if __name__ == "__main__":
    # set up initial state
    simulation = VirusSpreadingSimulation(population=500, xlim=xlim, ylim=ylim, speed_avg=15, p_infection=.9,
                                          radius_infection=4.)

    # set up figure and animation
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                         xlim=(-xlim * 0.005, xlim * 1.005), ylim=(-ylim * 0.005, ylim * 1.005))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    # particles holds the locations of the particles
    h_points, = ax.plot(simulation.pop.loc[simulation.pop["situation"] == "not_infected", "x"],
                        simulation.pop.loc[simulation.pop["situation"] == "not_infected", "y"], 'bo', ms=3)
    i_points, = ax.plot(simulation.pop.loc[simulation.pop["situation"] == "infected", "x"],
                        simulation.pop.loc[simulation.pop["situation"] == "infected", "y"], 'o', color="orange", ms=3)
    s_points, = ax.plot(simulation.pop.loc[simulation.pop["situation"] == "sick", "x"],
                        simulation.pop.loc[simulation.pop["situation"] == "sick", "y"], 'ro', ms=3)
    r_points, = ax.plot(simulation.pop.loc[simulation.pop["situation"] == "recovered", "x"],
                        simulation.pop.loc[simulation.pop["situation"] == "recovered", "y"], 'go', ms=3)
    d_points, = ax.plot(simulation.pop.loc[simulation.pop["situation"] == "dead", "x"],
                        simulation.pop.loc[simulation.pop["situation"] == "dead", "y"], 'ko', ms=3)

    # rect is the box edge
    rect = plt.Rectangle((0, 0), xlim, ylim, ec='none', lw=2, fc='none')
    ax.add_patch(rect)

    def animate(i):
        """perform animation step"""
        global simulation, rect, ax, fig
        simulation.step()

        # update pieces of the animation
        rect.set_edgecolor('k')
        h_points.set_data(simulation.pop.loc[simulation.pop["situation"] == "not_infected", "x"],
                          simulation.pop.loc[simulation.pop["situation"] == "not_infected", "y"])
        i_points.set_data(simulation.pop.loc[simulation.pop["situation"] == "infected", "x"],
                          simulation.pop.loc[simulation.pop["situation"] == "infected", "y"])
        s_points.set_data(simulation.pop.loc[simulation.pop["situation"] == "sick", "x"],
                          simulation.pop.loc[simulation.pop["situation"] == "sick", "y"])
        r_points.set_data(simulation.pop.loc[simulation.pop["situation"] == "recovered", "x"],
                          simulation.pop.loc[simulation.pop["situation"] == "recovered", "y"])
        d_points.set_data(simulation.pop.loc[simulation.pop["situation"] == "dead", "x"],
                          simulation.pop.loc[simulation.pop["situation"] == "dead", "y"])
        return h_points, i_points, s_points, r_points, d_points, rect


    ani = animation.FuncAnimation(fig, animate, frames=600, interval=10, blit=True)

    ani.save('animation.gif', writer='PillowWriter', fps=20)

    plt.show()
