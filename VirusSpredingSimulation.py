import numpy as np
import pandas as pd
import random as rd
from math import pi

rd.seed(123)
np.random.seed(123)

# 30 frame per second
fps = 30


class VirusSpreadingSimulation:

    def __init__(self, population=1000, xlim=100., ylim=100.,
                 radius_infection=1, p_infection=0.2, incubation_period=14,
                 healing_duration=14, death_rate=0.1, without_symptom_rate=0.3,
                 quarantine_zone=False,
                 isolation_threshold=0.2, pct_pop_insolated=0.8,
                 speed_avg=1, step_per_day=5):

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
        :param speed_avg: float => average speed of a point
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
            "incubation_period": [-1] * population,
            "quarantine_zone": [False] * population,
            "confined": [False] * population,
            "recovered": [False] * population,
            "dead": [False] * population,
            # day of death (will not die if "inf")
            "death_day": np.where(np.random.uniform(0, 1, population) > self.death_rate, float("inf"),
                                  np.random.choice(range(incubation_period, healing_duration+1), population)),
            "radian": np.random.uniform(-pi, pi, population),
            "days_keeping_radian": np.random.randint(1, 11, population),
            # Possible situations ["not_infected", "infected", "sick", "recovered", "dead"]
            "situation": ["not_infected"] * population
        })
        # first sick person in middle of area => patient zero is average person
        self.pop.iloc[0, [0, 1, 2, 3, 4, 12]] = [xlim/2, ylim/2, speed_avg, False, 0, "infected"]

        # parameters about how virus is spreading and how population react
        self.radius_infection = radius_infection
        self.p_infection = p_infection
        self.incubation_period = incubation_period
        self.healing_duration = healing_duration
        self.death_rate = death_rate
        self.without_symptom_rate = without_symptom_rate
        self.quarantine_zone = quarantine_zone
        self.isolation_threshold = isolation_threshold
        self.pct_pop_insolated = pct_pop_insolated

        # counting time
        self.step_per_day = step_per_day
        self.day = 0
        self.mvt = 0
        self.mvt_in_day = 0
        self.nb_step = 0

    def step(self):
        """
        every frame update positions
        """
        self.nb_step += 1

        # update positions
        self.pop["x"] = self.pop["x"] + self.pop["speed"] * fps * np.cos(self.pop["radian"])
        self.pop["y"] = self.pop["y"] + self.pop["speed"] * fps * np.sin(self.pop["radian"])

        # check for crossing boundary
        crossed_x1 = (self.pop["x"] < 0)
        crossed_x2 = (self.pop["x"] > self.xlim)
        crossed_y1 = (self.pop["y"] < 0)
        crossed_y2 = (self.pop["y"] > self.xlim)

        # Calculate new position and new radian
        # new position on border, new radian calculate for each point (on available directions)
        # -> crossing corner
        self.pop.loc[crossed_x1 & crossed_y1,
                     ["x", "y", "radian"]] = [0, 0, rd.uniform(0, pi / 2, sum(crossed_x1 & crossed_y1))]
        self.pop.loc[crossed_x1 & crossed_y2,
                     ["x", "y", "radian"]] = [0, self.ylim, rd.uniform(pi / 2, pi, sum(crossed_x1 & crossed_y2))]
        self.pop.loc[crossed_x2 & crossed_y2,
                     ["x", "y", "radian"]] = [self.xlim, self.ylim,
                                              rd.uniform(pi, 3 * pi / 2, sum(crossed_x2 & crossed_y2))]
        self.pop.loc[crossed_x2 & crossed_y1,
                     ["x", "y", "radian"]] = [self.xlim, 0, rd.uniform(3 * pi / 2, 2 * pi),
                                              sum(crossed_x2 & crossed_y1)]
        # -> crossing line
        self.pop.loc[crossed_x1, ["x", "radian"]] = [0, rd.uniform(0, pi, sum(crossed_x1))]
        self.pop.loc[crossed_x2, ["x", "radian"]] = [self.xlim, rd.uniform(pi, 2 * pi, sum(crossed_x2))]
        self.pop.loc[crossed_y1, ["y", "radian"]] = [0, rd.uniform(-pi / 2, pi / 2, sum(crossed_y1))]
        self.pop.loc[crossed_y2, ["y", "radian"]] = [self.ylim, rd.uniform(pi / 2, 3 * pi / 2, sum(crossed_y2))]

        if self.nb_step == fps:
            self.nb_step = 0
            self.movement()

    def movement(self):
        # Contaminate with probability
        self.contamination()

        self.mvt += 1
        self.mvt_in_day += 1
        if self.mvt_in_day == self.step_per_day:
            self.mvt_in_day = 0
            self.day += 1
            self.new_day()

    def contamination(self):
        # contaminate close points with probability
        for index, row in self.pop.iterrows():
            # probability depending of the number of points around
            # each point have probability self.p_infection to infect
            p_stay_safe = (1 - self.p_infection)**self.get_nb_person_around((row["x"], row["y"]))
            if p_stay_safe < rd.uniform(0, 1):
                self.pop.loc[index, ["healthy", "incubation_period", "situation"]] = [False, 0, "infected"]

    def get_nb_person_around(self, pos):
        # number of person around point (minus himself)
        return sum(np.sqrt((self.pop["x"] - pos[0])**2 + (self.pop["y"] - pos[1])**2) < self.radius_infection) - 1

    def new_day(self):
        # decrease days_keeping_radian and change radian if needed
        self.pop["days_keeping_radian"] = self.pop["days_keeping_radian"] - 1
        change_radian = (self.pop["days_keeping_radian"] == 0)
        self.pop[change_radian, ["radian", "days_keeping_radian"]] = \
            [np.random.uniform(-pi, pi, sum(change_radian)), np.random.randint(1, 11, sum(change_radian))]

        # Vérifie quanrentine and sick people
        # TODO



print("FINISH")
