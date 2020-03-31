import matplotlib.pyplot as plt
import matplotlib.animation as animation
from VirusSpreadingSimulation.VirusSpredingSimulation import VirusSpreadingSimulation

xlim = 100.
ylim = 100.
density = 0.1
pop = int(density * xlim * ylim)

if __name__ == "__main__":

    # set up initial state
    simulation = VirusSpreadingSimulation(population=pop, xlim=xlim, ylim=ylim)

    # set up figure and animation
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                         xlim=(-xlim*0.005, xlim*1.005), ylim=(-ylim*0.005, ylim*1.005))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    # particles holds the locations of the particles
    points, = ax.plot(simulation.pop["x"], simulation.pop["y"], 'bo', ms=3)

    # rect is the box edge
    rect = plt.Rectangle((0, 0), xlim, ylim, ec='none', lw=2, fc='none')
    ax.add_patch(rect)


    def init():
        """initialize animation"""
        global rect
        points.set_data([], [])
        rect.set_edgecolor('none')
        return points, rect

    def animate(i):
        """perform animation step"""
        global simulation, rect, ax, fig
        simulation.step()

        # update pieces of the animation
        rect.set_edgecolor('k')
        points.set_data(simulation.pop["x"], simulation.pop["y"])
        return points, rect


    ani = animation.FuncAnimation(fig, animate, frames=600, interval=10, blit=True, init_func=init)

    plt.show()
