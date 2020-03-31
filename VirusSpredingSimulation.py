import numpy as np
import pandas as pd
from math import pi

np.random.seed(123)

# 30 frame per second
fps = 30


class VirusSpreadingSimulation:

    def __init__(self, population=20, xlim=20., ylim=20.,
                 radius_infection=1., p_infection=0.2, incubation_period=14,
                 healing_duration=14, death_rate=0.1, without_symptom_rate=0.3,
                 quarantine_zone=False,
                 isolation_threshold=0.2, pct_pop_isolated=0.8,
                 speed_avg=10., step_per_day=5):

        """
        __init__
        :param population: int => number of point in population
        :param xlim: float => x width
        :param ylim: float => y width
        :param radius_infection: float => contagion distance
        :param p_infection: float ∈ [0,1] => probability to be infected by a point in contagion distance
        :param incubation_period: int => number of days to have symptoms
        :param healing_duration: int => number of days to heal
        :param death_rate: float ∈ [0,1] => probability to die if infected
        :param without_symptom_rate: float ∈ [0,1] => probability to not show any symptoms
        :param quarantine_zone: bool => create a quarantine area (separation between sick and healthy points)
        :param isolation_threshold: float ∈ [0,1] => threshold from which movements will stop
        :param pct_pop_insolated: float ∈ [0,1] => pct of population who will really stop moving
        :param speed_avg: float ∈ R+* => average speed of a point
        :param step_per_day: int => number of step per day (possible to be infect at each step)
        """

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
            "confined": [False] * population,
            "recovered": [False] * population,
            "dead": [False] * population,
            # day of death (will not die if "inf")
            "death_day": np.where(np.random.uniform(0, 1, population) > death_rate, float("inf"),
                                  np.random.choice(
                                      range(incubation_period, incubation_period + healing_duration), population)),
            "radian": np.random.uniform(-pi, pi, population),
            "days_keeping_radian": np.random.randint(1, 5, population),
            # Possible situations ["not_infected", "infected", "sick", "recovered", "dead"]
            "situation": ["not_infected"] * population
        })
        # first sick person in middle of area => patient zero is average person
        self.pop.iloc[0, [0, 1, 2, 3, 4, 12]] = [xlim / 2, ylim / 2, speed_avg, False, 0, "infected"]

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

        # counting time
        self.step_per_day = step_per_day
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

        # init spreading statistics
        self.stats = pd.DataFrame({
            "nb_not_infected": [self.nb_not_infected],
            "nb_infected": [self.nb_infected],
            "nb_sick": [self.nb_sick],
            "nb_recovered": [self.nb_recovered],
            "nb_dead": [self.nb_dead]},
            index=['Day_' + str(0)]
        )

    def step(self):
        """every frame update positions"""
        self.nb_step += 1

        # update positions
        self.pop["x"] = self.pop["x"] + self.pop["speed"] / self.step_per_day / fps * np.cos(self.pop["radian"])
        self.pop["y"] = self.pop["y"] + self.pop["speed"] / self.step_per_day / fps * np.sin(self.pop["radian"])

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

        if self.nb_step == fps:
            self.nb_step = 0
            self.movement()

    def movement(self):
        """
        Movement is represent in a second on the plot, with fps frame per seconds
        At each movement, can be contaminate
        """
        # Contaminate with probability
        self.contamination()

        self.mvt += 1
        self.mvt_in_day += 1
        if self.mvt_in_day == self.step_per_day:
            self.mvt_in_day = 0
            self.day += 1
            self.new_day()

    def contamination(self):
        """Contamination depending p_infection"""
        # identify caring points
        caring_points = self.pop.loc[(self.pop["situation"].isin(["infected", "sick"]) &
                                      (self.pop["infected_day"] != self.day)),
                                     ["x", "y"]]
        print(caring_points)
        # contaminate close points with probability
        for index, row in self.pop.iterrows():
            if row["situation"] == "not_infected":
                # probability depending of the number of points around
                # each point have probability self.p_infection to infect
                p_stay_safe = (1 - self.p_infection) ** \
                              self.get_nb_caring_points_around((row["x"], row["y"]), caring_points)
                if p_stay_safe < np.random.uniform():
                    self.pop.loc[index, ["healthy", "infected_day", "situation"]] = [False, self.day, "infected"]

    def get_nb_caring_points_around(self, pos, caring_points):
        """
        get_nb_person_around
        :param pos: tuple => coordinates of a point (x, y)
        :param caring_points: DataFrame => DataFrame of caring points with columns ["x", "y"]
        :return: int => number of infected point around point position
        """
        # don't need to subtract himself because point is'nt infected
        return caring_points.loc[
            (np.sqrt((caring_points["x"] - pos[0]) ** 2 + (caring_points["y"] - pos[1]) ** 2) < self.radius_infection),
            ["x", "y"]].shape[0]

    def get_nb_caring_points_around_previous(self, pos, ):
        """
        get_nb_person_around
        :param pos: tuple => coordinates of a point (x, y)
        :return: int => number of infected point around point position
        """
        # don't need to subtract himself because point is'nt infected
        return self.pop.loc[
            (self.pop["infected_day"] > 0) & (~self.pop["dead"]) & (~self.pop["recovered"]) &
            (np.sqrt((self.pop["x"] - pos[0]) ** 2 + (self.pop["y"] - pos[1]) ** 2) < self.radius_infection),
            ["x", "y"]].shape[0]

    def new_day(self):

        # update points situation
        self.update_dead_points()
        self.update_infected_to_sick_points()
        self.update_recovered_points()
        self.update_radian_points()

        # update counters
        self.update_spreading_counters()

        # TODO
        # if self.quarantine_zone:
        # ----> Remove sick people of area

        # TODO
        # Confined points (if self.isolation_threshold > self.nb_sick / self.pop_size)
        # ----> Points stay immobile if isolation_threshold is exceeded

        print("-- END OF DAY", self.day, "--")
        print("  -> pct not_infected :", self.nb_not_infected / self.pop_size * 100, "% \t", self.nb_not_infected)
        print("  -> pct infected :", self.nb_infected / self.pop_size * 100, "% \t", self.nb_infected)
        print("  -> pct sick :", self.nb_sick / self.pop_size * 100, "% \t", self.nb_sick)
        print("  -> pct recovered :", self.nb_recovered / self.pop_size * 100, "% \t", self.nb_recovered)
        print("  -> pct dead :", self.nb_dead / self.pop_size * 100, "% \t", self.nb_dead, "\n")

    def update_spreading_counters(self):
        """update all spreading counters"""
        value_count = self.pop["situation"].value_counts()
        for situation in ["not_infected", "infected", "sick", "recovered", "dead"]:
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
            "nb_dead": [self.nb_dead]},
            index=['Day_' + str(self.day)]
        ))

    def update_dead_points(self):
        """update dead points"""
        self.pop.loc[(self.pop["death_day"] == self.day - self.pop["infected_day"]) & (self.pop["infected_day"] != -1),
                     ["quarantine_zone", "confined", "dead", "days_keeping_radian", "situation"]] = \
            [False, True, True, 999999, "dead"]

    def update_infected_to_sick_points(self):
        """points at the end of incubation period passing from infected to sick"""
        self.pop.loc[(self.pop["infected_day"] + self.incubation_period == self.day) & (self.pop["infected_day"] != -1),
                     ["quarantine_zone", "situation"]] = [True, "sick"]

    def update_recovered_points(self):
        """points at the end of healing period passing from sick to recovered"""
        self.pop.loc[(self.day - self.pop["infected_day"] == self.healing_duration + self.incubation_period) &
                     (self.pop["infected_day"] != -1),
                     ["healthy", "quarantine_zone", "confined", "recovered", "situation"]] = \
            [True, False, False, True, "recovered"]

    def update_radian_points(self):
        """decrease days_keeping_radian and change radian if needed"""
        self.pop["days_keeping_radian"] = self.pop["days_keeping_radian"] - 1
        change_radian = (self.pop["days_keeping_radian"] == 0)
        if sum(change_radian) > 1:
            self.pop.loc[change_radian, ["radian", "days_keeping_radian"]] = \
                np.array([np.random.uniform(-pi, pi, sum(change_radian)),
                          np.random.randint(1, 6, sum(change_radian)).astype(int)]).T
        elif sum(change_radian) == 1:
            self.pop.loc[change_radian, ["radian", "days_keeping_radian"]] = \
                [np.random.uniform(-pi, pi, sum(change_radian))[0],
                 int(np.random.randint(1, 6, sum(change_radian))[0])]

    def get_stats(self):
        return self.stats
