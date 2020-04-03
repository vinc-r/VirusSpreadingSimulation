import os
import numpy as np
import pandas as pd
from math import pi
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.animation as animation

MAX_DAY_KEEP_RAD = 11
FPS_SAVE = 30

np.random.seed(123)


class VirusSpreadingSimulation:
    """
    @title : VirusSpreadingSimulation (VSS)
    @author : Vincent Rosset
    @version : 0.0.11
    """

    def __init__(self, population=100, xlim=5., ylim=5.,
                 radius_infection="default", p_infection=0.1, incubation_period=14,
                 healing_duration=14, death_rate=0.1, without_symptom_rate=0.3,
                 quarantine_zone=True, isolation_threshold=0.05, pct_pop_isolated=0.9,
                 immunity_spreading=True, speed_avg="default", mpd=5, fpm=2):

        """
        __init__

        Parameters :
        ------------
            - population : int ∈ N+
                number of point in population
            - xlim : float ∈ R+*
                x width
            - ylim : float ∈ R+*
                y width
            - radius_infection : "default", float ∈ R+*
                contagion distance
                if default radius_infection is set up at : (xlim + ylim) / 100
            - p_infection : float ∈ [0,1]
                probability to be infected by a point in contagion distance
            - incubation_period : int ∈ N+
                number of days to have symptoms
            - healing_duration : int ∈ N+
                number of days to heal
            - death_rate : float ∈ [0,1]
                probability to die if infected
            - without_symptom_rate : float ∈ [0,1]
                probability to not show any symptoms
            - quarantine_zone : bool
                create a quarantine area (separation between sick and healthy points)
            - isolation_threshold : float ∈ [0,1]
                threshold of sick points from which movements will stop
            - pct_pop_isolated : float ∈ [0,1]
                pct of population who will really stop moving
            - speed_avg : 'default', float ∈ R+*
                average speed of a point
                if default speed_avg is set up at : (xlim + ylim) / 20
            - mpd : int ∈ N+*
                number of movement per day (possible to be infect at each movement)
            - fpm : int ∈ N+*
                number of frame per movement (actualisation of point positions/color on DataViz)
        """

        self.__name__ = 'VirusSpreadingSimulation'
        self.__nickname__ = 'VSS'
        self.__version__ = '0.0.11'
        self.__author__ = 'Vincent ROSSET'
        self.__email__ = 'vinc.rosset@gmail.com'
        self.__status__ = 'Production'
        self.__coding__ = 'utf-8'

        # set up default arg
        if speed_avg == "default":
            speed_avg = (xlim + ylim) / 10
        if radius_infection == "default":
            radius_infection = (xlim + ylim) / 50

        # Set population
        self.pop_size = population
        self.xlim = xlim
        self.ylim = ylim
        self.pop = pd.DataFrame({
            "x": np.random.uniform(0, xlim, population),
            "y": np.random.uniform(0, ylim, population),
            "speed": np.random.chisquare(speed_avg, population),
            "healthy": [True] * population,
            "infected_day": [-1] * population,
            "quarantine_zone": [False] * population,
            "stay_confined": np.random.choice([True, False], population, p=[pct_pop_isolated, 1 - pct_pop_isolated]),
            "recovered": [False] * population,
            "dead": [False] * population,
            # day of death (will not die if "inf")
            "death_day": np.where(np.random.uniform(0, 1, population) > death_rate, float("inf"),
                                  np.random.choice(
                                      range(incubation_period, incubation_period + healing_duration), population)),
            "radian": np.random.uniform(-pi, pi, population),
            "days_keeping_radian": np.random.randint(1, MAX_DAY_KEEP_RAD, population),
            # Possible situations ["not_infected", "infected", "sick", "recovered", "dead"]
            "situation": ["not_infected"] * population,
            "without_symptoms": np.random.choice([True, False], population,
                                                 p=[without_symptom_rate, 1 - without_symptom_rate])
        })
        # first sick person in middle of area => patient zero is average person
        self.pop.iloc[0, [0, 1, 2, 3, 4, 12, 13]] = [xlim / 2, ylim / 2, speed_avg, False, 0, "infected", False]

        # parameters about how virus is spreading and how population react
        self.radius_infection = radius_infection
        self.p_infection = p_infection
        self.incubation_period = incubation_period
        self.healing_duration = healing_duration
        self.death_rate = death_rate
        self.without_symptom_rate = without_symptom_rate
        self.quarantine_zone = quarantine_zone
        self.isolation_threshold = isolation_threshold
        self.pct_pop_isolated = pct_pop_isolated
        self.immunity_spreading = immunity_spreading

        # no isolation at the beginning
        self.isolation = False

        # counting time
        self.mpd = mpd
        self.fpm = fpm
        self.day = 0
        self.mvt = 0
        self.mvt_in_day = 0
        self.nb_step = 0

        # init spreading counters
        # MAJ each day
        self.nb_not_infected = population - 1
        self.nb_infected = 1
        self.nb_sick = 0
        self.nb_recovered = 0
        self.nb_dead = 0
        self.nb_immune = 0

        # init spreading statistics
        self.stats = pd.DataFrame({
            "nb_not_infected": [self.nb_not_infected],
            "nb_infected": [self.nb_infected],
            "nb_sick": [self.nb_sick],
            "nb_recovered": [self.nb_recovered],
            "nb_dead": [self.nb_dead],
            "nb_immune": [self.nb_immune]},
            index=['Day_' + str(0)]
        )

    def __repr__(self):
        """print class"""
        return '<ViruseSpredingSimulation class>'

    def __str__(self):
        """print class"""
        return '<ViruseSpredingSimulation class>'

    def __len__(self):
        """number of days"""
        return self.day

    def __add__(self, nb):
        """increase simulation of nb day"""
        try:
            nb = int(nb)
        except:
            raise ValueError
        if nb < 0:
            raise ValueError
        start = self.day
        while start + nb != self.day:
            self.frame(verbose=0)

    def frame(self, verbose=1):
        """every frame update positions"""
        self.nb_step += 1

        # update positions
        immobile_points = (self.pop["situation"].isin(["dead", "sick"])) | (self.isolation & self.pop["stay_confined"])
        self.pop["x"] = \
            self.pop["x"] + self.pop["speed"] * (1 - immobile_points) / self.mpd / self.fpm * np.cos(self.pop["radian"])
        self.pop["y"] = \
            self.pop["y"] + self.pop["speed"] * (1 - immobile_points) / self.mpd / self.fpm * np.sin(self.pop["radian"])

        # check for crossing boundary
        crossed_x1 = (self.pop["x"] < 0)
        crossed_x2 = (self.pop["x"] > self.xlim)
        crossed_y1 = (self.pop["y"] < 0)
        crossed_y2 = (self.pop["y"] > self.xlim)

        # Calculate new position and new radian if boundary crossed
        # new position on border, new radian calculate for each point (on available directions)
        # -> crossing line
        self.pop.loc[crossed_x1, ["x", "radian"]] = [0, np.random.uniform(-pi / 4, pi / 4, sum(crossed_x1))]
        self.pop.loc[crossed_x2, ["x", "radian"]] = [self.xlim,
                                                     np.random.uniform(3 * pi / 4, 5 * pi / 4, sum(crossed_x2))]
        self.pop.loc[crossed_y1, ["y", "radian"]] = [0, np.random.uniform(pi / 4, 3 * pi / 4, sum(crossed_y1))]
        self.pop.loc[crossed_y2, ["y", "radian"]] = [self.ylim,
                                                     np.random.uniform(5 * pi / 4, 7 * pi / 4, sum(crossed_y2))]
        # -> crossing corner
        self.pop.loc[crossed_x1 & crossed_y1, ["x", "y", "radian"]] = \
            [0, 0, np.random.uniform(pi / 8, 3 * pi / 8, sum(crossed_x1 & crossed_y1))]
        self.pop.loc[crossed_x1 & crossed_y2, ["x", "y", "radian"]] = \
            [0, self.ylim, np.random.uniform(13 * pi / 8, 15 * pi / 8, sum(crossed_x1 & crossed_y2))]
        self.pop.loc[crossed_x2 & crossed_y2, ["x", "y", "radian"]] = \
            [self.xlim, self.ylim, np.random.uniform(9 * pi / 8, 11 * pi / 8, sum(crossed_x2 & crossed_y2))]
        self.pop.loc[crossed_x2 & crossed_y1, ["x", "y", "radian"]] = \
            [self.xlim, 0, np.random.uniform(5 * pi / 8, 7 * pi / 8, sum(crossed_x2 & crossed_y1))]

        if self.nb_step == self.fpm:
            self.nb_step = 0
            self.movement(verbose=verbose)

    def movement(self, verbose=0):
        """
        Movement is represent in a second on the plot, with fpm frame per movement
        At each movement, can be contaminate
        """
        # Contaminate with probability
        self.contamination()
        # immunity spreading
        if self.immunity_spreading:
            self.immunity_spreading_function()

        self.mvt += 1
        self.mvt_in_day += 1
        if self.mvt_in_day == self.mpd:
            self.mvt_in_day = 0
            self.day += 1
            self.new_day(verbose=verbose)

    def contamination(self):
        """Contamination depending p_infection"""
        # identify caring points
        if not self.quarantine_zone:
            caring_points = self.pop.loc[(self.pop["situation"].isin(["infected", "sick"]) &
                                          (self.pop["infected_day"] != self.day)), ["x", "y"]]
        else:
            # if quarantine, sick points can't infect
            caring_points = self.pop.loc[(self.pop["situation"] == "infected") &
                                         (self.pop["infected_day"] != self.day), ["x", "y"]]
        if len(caring_points) == 0:
            return
        # contaminate close points with probability
        for index, row in self.pop.iterrows():
            if row["situation"] == "not_infected":
                # probability depending of the number of points around
                # each point have probability self.p_infection to infect
                first_filter = (caring_points["x"].between(row["x"] - self.radius_infection,
                                                           row["x"] + self.radius_infection, inclusive=False)) & \
                               (caring_points["y"].between(row["y"] - self.radius_infection,
                                                           row["y"] + self.radius_infection, inclusive=False))
                if sum(first_filter) == 0:
                    continue
                p_stay_safe = (1 - self.p_infection) ** self.get_nb_caring_points_around((row["x"], row["y"]),
                                                                                         caring_points.loc[
                                                                                             first_filter])
                if p_stay_safe < np.random.uniform():
                    self.pop.loc[index, ["healthy", "infected_day", "situation"]] = [False, self.day, "infected"]

    def immunity_spreading_function(self):
        """Contamination depending p_infection"""
        # identify caring points of immunity
        caring_points = self.pop.loc[self.pop["situation"].isin(["immune", "recovered"]), ["x", "y"]]
        if len(caring_points) == 0:
            return
        # contaminate close points with probability
        for index, row in self.pop.iterrows():
            # only not_infected can receive immunity
            if row["situation"] == "not_infected":
                # probability depending of the number of points around
                # each point have probability self.p_infection to infect
                first_filter = (caring_points["x"].between(row["x"] - self.radius_infection,
                                                           row["x"] + self.radius_infection, inclusive=False)) & \
                               (caring_points["y"].between(row["y"] - self.radius_infection,
                                                           row["y"] + self.radius_infection, inclusive=False))
                if sum(first_filter) == 0:
                    continue
                p_stay_usual = (1 - self.p_infection) ** self.get_nb_caring_points_around((row["x"], row["y"]),
                                                                                          caring_points.loc[first_filter])
                if p_stay_usual < np.random.uniform():
                    self.pop.loc[index, ["situation"]] = "immune"

    def get_nb_caring_points_around(self, pos, caring_points):
        """
        get_nb_person_around

        Parameters :
        ------------
            - pos : tuple
                coordinates of a point (x, y)
            - caring_points : DataFrame
                DataFrame of caring points with columns ["x", "y"]

        Returns :
        ---------
            - int : number of infected point around point position
        """
        # don't need to subtract himself because point is'nt infected
        return caring_points.loc[
            (np.sqrt((caring_points["x"] - pos[0]) ** 2 + (caring_points["y"] - pos[1]) ** 2) < self.radius_infection),
            ["x", "y"]].shape[0]

    def new_day(self, verbose=0):

        # update points situation
        self.update_dead_points()
        self.update_infected_to_sick_points()
        self.update_recovered_points()
        self.update_radian_points()

        # update counters
        self.update_spreading_counters()

        # update isolation of pop
        sick_prop = self.nb_sick / self.pop_size
        if not self.isolation and self.isolation_threshold <= sick_prop:
            self.isolation = True
        elif self.isolation and self.isolation_threshold > sick_prop:
            self.isolation = False

        if verbose > 0:
            self.print_situation()

    def update_spreading_counters(self):
        """update all spreading counters"""
        value_count = self.pop["situation"].value_counts()
        for situation in ["not_infected", "infected", "sick", "recovered", "dead", "immune"]:
            try:
                exec('self.nb_' + situation + ' = value_count["' + situation + '"]')
            except KeyError:
                exec('self.nb_' + situation + ' = 0')

        # update spreading statistics
        self.stats = self.stats.append(pd.DataFrame({
            "nb_not_infected": [self.nb_not_infected],
            "nb_infected": [self.nb_infected],
            "nb_sick": [self.nb_sick],
            "nb_recovered": [self.nb_recovered],
            "nb_dead": [self.nb_dead],
            "nb_immune": [self.nb_immune]},
            index=['Day_' + str(self.day)]
        ))

    def update_dead_points(self):
        """update dead points"""
        self.pop.loc[(self.pop["death_day"] == self.day - self.pop["infected_day"]) &
                     (self.pop["infected_day"] != -1) & (~self.pop["without_symptoms"]),
                     ["quarantine_zone", "stay_confined", "dead", "days_keeping_radian", "situation"]] = \
            [True, True, True, 999999, "dead"]

    def update_infected_to_sick_points(self):
        """points at the end of incubation period passing from infected to sick"""
        self.pop.loc[(self.pop["infected_day"] + self.incubation_period == self.day) & (~self.pop["without_symptoms"]) &
                     (self.pop["situation"] == "infected"), ["quarantine_zone", "situation"]] = [True, "sick"]

    def update_recovered_points(self):
        """points at the end of healing period passing from sick to recovered"""
        self.pop.loc[(self.day - self.pop["infected_day"] == self.healing_duration + self.incubation_period) &
                     (self.pop["situation"].isin(["sick", "infected"])),
                     ["healthy", "quarantine_zone", "stay_confined", "recovered", "situation"]] = \
            [True, False, False, True, "recovered"]

    def update_radian_points(self):
        """decrease days_keeping_radian and change radian if needed"""
        self.pop["days_keeping_radian"] = self.pop["days_keeping_radian"] - 1
        change_radian = (self.pop["days_keeping_radian"] == 0)
        if sum(change_radian) > 1:
            self.pop.loc[change_radian, ["radian", "days_keeping_radian"]] = \
                np.array([np.random.uniform(-pi, pi, sum(change_radian)),
                          np.random.randint(1, MAX_DAY_KEEP_RAD, sum(change_radian)).astype(int)]).T
        elif sum(change_radian) == 1:
            self.pop.loc[change_radian, ["radian", "days_keeping_radian"]] = \
                [np.random.uniform(-pi, pi, sum(change_radian))[0],
                 int(np.random.randint(1, MAX_DAY_KEEP_RAD, sum(change_radian))[0])]

    def get_stats(self):
        return self.stats

    def get_plot(self, show=True, save=False, save_name=None):

        fig = plt.figure(constrained_layout=False, figsize=(9, 6))

        # Write suptitle with Latex
        fig.suptitle(r'$Virus~Spreading~Simulation$' + '\n' +
                     r'$_{population:' + str(self.pop_size) + '~~~}$' +
                     r'$_{radius~infection:' + str(self.radius_infection) + '~~~}$' +
                     r'$_{\mathbb{P}_{infection}:' + str(self.p_infection) + '~~~}$' +
                     r'$_{incubation~period:' + str(self.incubation_period) + '}$' +
                     r'$_{~~~healing~duration:' + str(self.healing_duration) + '~~~}$' +
                     r'$_{death~rate:' + str(self.death_rate) + '}$' +
                     '\n' + r'$_{\%~without~symptoms:' + str(self.without_symptom_rate * 100) + '~~~}$' +
                     r'$_{containment~threshold:' + str(self.isolation_threshold) + '~~~}$' +
                     r'$_{\%~population~contain:' + str(self.pct_pop_isolated * 100) + '~~~}$' +
                     r'$_{immunity~spreading:' + str(self.immunity_spreading) + '}$')

        # init axis
        grid = fig.add_gridspec(nrows=2 + self.quarantine_zone, ncols=3)
        ax1 = fig.add_subplot(grid[:, :-1])
        ax2 = fig.add_subplot(grid[0, -1])
        ax3 = fig.add_subplot(grid[1, -1])
        if self.quarantine_zone:
            ax4 = fig.add_subplot(grid[2, -1])

        # ax1.axis('off')
        ax1.autoscale(enable=True, axis='both', tight=True)
        ax1.add_patch(mpatches.Rectangle((0, 0), self.xlim, self.ylim, color="black", fill=False))
        if self.quarantine_zone:
            ax1.set_xlabel('Global area')
        ax2.axis('off')
        ax3.set_title('Evolution by days')
        if self.quarantine_zone:
            ax4.axis('off')
            ax4.add_patch(mpatches.Rectangle((0, 0), self.xlim, self.ylim, color="red", fill=False))
            ax4.set_title('Quarantine area', color="red", y=-0.12)

        # plot points (except sick and dead)
        h_points, = ax1.plot(self.pop.loc[self.pop["situation"] == "not_infected", "x"],
                             self.pop.loc[self.pop["situation"] == "not_infected", "y"], 'bo', ms=3)
        i_points, = ax1.plot(self.pop.loc[self.pop["situation"] == "infected", "x"],
                             self.pop.loc[self.pop["situation"] == "infected", "y"], 'o', color="orange", ms=3)
        r_points, = ax1.plot(self.pop.loc[self.pop["situation"] == "recovered", "x"],
                             self.pop.loc[self.pop["situation"] == "recovered", "y"],
                             marker='$☺$', color='g', linestyle='None', ms=10)
        im_points, = ax1.plot(self.pop.loc[self.pop["situation"] == "immune", "x"],
                              self.pop.loc[self.pop["situation"] == "immune", "y"], 'co', ms=3)

        # plot sick points
        if self.quarantine_zone:
            # plot on quarantine area
            s_points, = ax4.plot(self.pop.loc[self.pop["situation"] == "sick", "x"],
                                 self.pop.loc[self.pop["situation"] == "sick", "y"], 'ro', ms=3)
            d_points, = ax4.plot(self.pop.loc[self.pop["situation"] == "dead", "x"],
                                 self.pop.loc[self.pop["situation"] == "dead", "y"],
                                 marker='$†$', color="k", linestyle='None', ms=8)
        else:
            # plot on main area
            s_points, = ax1.plot(self.pop.loc[self.pop["situation"] == "sick", "x"],
                                 self.pop.loc[self.pop["situation"] == "sick", "y"], 'ro', ms=3)
            d_points, = ax1.plot(self.pop.loc[self.pop["situation"] == "dead", "x"],
                                 self.pop.loc[self.pop["situation"] == "dead", "y"],
                                 marker='$†$', color="k", linestyle='None', ms=8)

        # plot cumulative area
        ax3.fill_between(range(self.day + 1), self.stats.cumsum(axis=1)["nb_immune"], facecolor='cyan')
        ax3.fill_between(range(self.day + 1), self.stats.cumsum(axis=1)["nb_dead"], facecolor='black')
        ax3.fill_between(range(self.day + 1), self.stats.cumsum(axis=1)["nb_recovered"], facecolor='green')
        ax3.fill_between(range(self.day + 1), self.stats.cumsum(axis=1)["nb_sick"], facecolor='red')
        ax3.fill_between(range(self.day + 1), self.stats.cumsum(axis=1)["nb_infected"], facecolor='orange')
        ax3.fill_between(range(self.day + 1), self.stats.cumsum(axis=1)["nb_not_infected"],
                         facecolor='navy')

        h_patch = mpatches.Patch(color='navy', label='Not Infected : ' + str(self.nb_not_infected))
        i_patch = mpatches.Patch(color='orange', label='Infected : ' + str(self.nb_infected))
        s_patch = mpatches.Patch(color='red', label='Sick : ' + str(self.nb_sick))
        r_patch = mpatches.Patch(color='green', label='Recovered : ' + str(self.nb_recovered))
        d_patch = mpatches.Patch(color='black', label='Dead : ' + str(self.nb_dead))
        im_patch = mpatches.Patch(color='cyan', label='Immune : ' + str(self.nb_immune))
        # hide immune legend if immunity can't spread
        if self.immunity_spreading:
            handles_legend = [h_patch, i_patch, s_patch, r_patch, d_patch, im_patch]
        else:
            handles_legend = [h_patch, i_patch, s_patch, r_patch, d_patch]
        prop_legend = {'size': 9 + 3 * (1 - self.quarantine_zone)}
        ax2.legend(handles=handles_legend, loc="center", title="Day 0", prop=prop_legend)

        def animate(i):
            """perform animation frame"""

            self.frame()

            # update pieces of the animation
            h_points.set_data(self.pop.loc[self.pop["situation"] == "not_infected", "x"],
                              self.pop.loc[self.pop["situation"] == "not_infected", "y"])
            i_points.set_data(self.pop.loc[self.pop["situation"] == "infected", "x"],
                              self.pop.loc[self.pop["situation"] == "infected", "y"])
            s_points.set_data(self.pop.loc[self.pop["situation"] == "sick", "x"],
                              self.pop.loc[self.pop["situation"] == "sick", "y"])
            r_points.set_data(self.pop.loc[self.pop["situation"] == "recovered", "x"],
                              self.pop.loc[self.pop["situation"] == "recovered", "y"])
            d_points.set_data(self.pop.loc[self.pop["situation"] == "dead", "x"],
                              self.pop.loc[self.pop["situation"] == "dead", "y"])
            im_points.set_data(self.pop.loc[self.pop["situation"] == "immune", "x"],
                               self.pop.loc[self.pop["situation"] == "immune", "y"])

            h_patch.set_label('Not Infected : ' + str(self.nb_not_infected))
            i_patch.set_label('Infected : ' + str(self.nb_infected))
            s_patch.set_label('Sick : ' + str(self.nb_sick))
            r_patch.set_label('Recovered : ' + str(self.nb_recovered))
            d_patch.set_label('Dead : ' + str(self.nb_dead))
            im_patch.set_label('Immune : ' + str(self.nb_immune))
            ax2.legend(handles=handles_legend, loc="center", title="Day " + str(self.day), prop=prop_legend)

            # plot cumulative area
            ax3.fill_between(range(self.day + 1), self.stats.cumsum(axis=1)["nb_immune"], facecolor='cyan')
            ax3.fill_between(range(self.day + 1), self.stats.cumsum(axis=1)["nb_dead"], facecolor='black')
            ax3.fill_between(range(self.day + 1), self.stats.cumsum(axis=1)["nb_recovered"],
                             facecolor='green')
            ax3.fill_between(range(self.day + 1), self.stats.cumsum(axis=1)["nb_sick"], facecolor='red')
            ax3.fill_between(range(self.day + 1), self.stats.cumsum(axis=1)["nb_infected"],
                             facecolor='orange')
            ax3.fill_between(range(self.day + 1), self.stats.cumsum(axis=1)["nb_not_infected"],
                             facecolor='navy')

            if self.quarantine_zone:
                return ax1, ax2, ax3, ax4
            else:
                return ax1, ax2, ax3

        ani = animation.FuncAnimation(fig, animate, frames=600)

        if save:
            if os.path.isfile(save_name + '.gif'):
                i = 1
                while os.path.isfile(save_name + str(i) + '.gif'):
                    i += 1
                ani.save(save_name + str(i) + '.gif', writer='PillowWriter', fps=FPS_SAVE, dpi=1)
            else:
                ani.save(save_name + '.gif', writer='PillowWriter', fps=FPS_SAVE, dpi=2)

        if show:
            plt.show()

    def compute_n_days(self, n, get_summary=True, verbose=0):
        """increase simulation of nb day"""
        try:
            n = int(n)
        except:
            raise ValueError
        if n < 0:
            raise ValueError

        for i in range(n):
            self.__add__(1)
            if verbose > 0:
                self.print_situation()

        if get_summary:
            return self.stats

    def print_situation(self, pct=True):
        print(">>>-- DAY", self.day, "--<<<")
        if pct:
            print("  -> pct not_infected :", round(self.nb_not_infected / self.pop_size * 100, 2), "%",
                  self.nb_not_infected)
            print("  -> pct infected :", round(self.nb_infected / self.pop_size * 100, 2), "% \t")
            print("  -> pct sick :", round(self.nb_sick / self.pop_size * 100, 2), "%")
            print("  -> pct recovered :", round(self.nb_recovered / self.pop_size * 100, 2), "%")
            print("  -> pct dead :", round(self.nb_dead / self.pop_size * 100, 2), "%")
            print("  -> pct immune :", round(self.nb_immune / self.pop_size * 100, 2), "%\n")
        else:
            print("  -> nb not_infected :", self.nb_not_infected)
            print("  -> nb infected :", self.nb_infected)
            print("  -> nb sick :", self.nb_sick)
            print("  -> nb recovered :", self.nb_recovered)
            print("  -> nb dead :", self.nb_dead)
            print("  -> nb immune :", self.nb_immune, "\n")
