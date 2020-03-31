import unittest
import pandas as pd
from VirusSpreadingSimulation.VirusSpredingSimulation import VirusSpreadingSimulation

# GLOBAL VARIABLES
if True:
    pop = 10
    xlim = 5.
    ylim = 6.
    min_days_keeping_radian = 1
    max_days_keeping_radian = 5
    speed_avg = 1.5
    radius_infection = 1.
    p_infection = .5
    incubation_period = 14
    healing_duration = 10
    death_rate = .1
    without_symptom_rate = .2
    quarantine_zone = True
    isolation_threshold = .333
    pct_pop_isolated = .75
    step_per_day = 5


class TestVSS(unittest.TestCase):
    """Test case for class VirusSpreadingSimulation"""

    def setUp(self):
        """Initialization using global variables"""
        self.VSS = VirusSpreadingSimulation(population=pop, xlim=xlim, ylim=ylim, speed_avg=speed_avg,
                                            radius_infection=radius_infection, p_infection=p_infection,
                                            incubation_period=incubation_period, healing_duration=healing_duration,
                                            death_rate=death_rate, without_symptom_rate=without_symptom_rate,
                                            quarantine_zone=quarantine_zone, isolation_threshold=isolation_threshold,
                                            pct_pop_isolated=pct_pop_isolated, step_per_day=step_per_day)

    def test_init(self):

        # check parameters
        self.assertEqual(self.VSS.pop_size, pop)
        self.assertEqual(self.VSS.xlim, xlim)
        self.assertEqual(self.VSS.ylim, ylim)
        self.assertEqual(self.VSS.radius_infection, radius_infection)
        self.assertEqual(self.VSS.p_infection, p_infection)
        self.assertEqual(self.VSS.incubation_period, incubation_period)
        self.assertEqual(self.VSS.healing_duration, healing_duration)
        self.assertEqual(self.VSS.death_rate, death_rate)
        self.assertEqual(self.VSS.without_symptom_rate, without_symptom_rate)
        self.assertEqual(self.VSS.quarantine_zone, quarantine_zone)
        self.assertEqual(self.VSS.isolation_threshold, isolation_threshold)
        self.assertEqual(self.VSS.pct_pop_isolated, pct_pop_isolated)
        self.assertEqual(self.VSS.step_per_day, step_per_day)

        # check created variables
        self.assertEqual(self.VSS.day, 0)
        self.assertEqual(self.VSS.mvt, 0)
        self.assertEqual(self.VSS.mvt_in_day, 0)
        self.assertEqual(self.VSS.nb_step, 0)
        self.assertEqual(self.VSS.nb_not_infected, pop - 1)
        self.assertEqual(self.VSS.nb_infected, 1)
        self.assertEqual(self.VSS.nb_sick, 0)
        self.assertEqual(self.VSS.nb_recovered, 0)
        self.assertEqual(self.VSS.nb_dead, 0)

        # check pop
        self.assertEqual(self.VSS.pop.shape, (pop, 13))
        self.assertListEqual(list(self.VSS.pop.columns),
                             ['x', 'y', 'speed', 'healthy', 'infected_day',
                              'quarantine_zone', 'confined', 'recovered', 'dead',
                              'death_day', 'radian', 'days_keeping_radian', 'situation'])
        self.assertTrue(self.VSS.pop['x'].min() >= 0)
        self.assertTrue(self.VSS.pop['x'].max() <= xlim)
        self.assertTrue(self.VSS.pop['y'].min() >= 0)
        self.assertTrue(self.VSS.pop['y'].max() <= ylim)
        self.assertTrue(self.VSS.pop['speed'].min() > 0)
        self.assertTrue(self.VSS.pop['healthy'].sum() == pop - 1)
        self.assertTrue((self.VSS.pop['infected_day'].value_counts() == [pop - 1, 1]).all())
        self.assertListEqual(list(self.VSS.pop['infected_day'].value_counts().index), [-1, 0])
        self.assertTrue(self.VSS.pop['quarantine_zone'].sum() == 0)
        self.assertTrue(self.VSS.pop['confined'].sum() == 0)
        self.assertTrue(self.VSS.pop['recovered'].sum() == 0)
        self.assertTrue(self.VSS.pop['dead'].sum() == 0)
        self.assertTrue(self.VSS.pop['days_keeping_radian'].min() >= min_days_keeping_radian)
        self.assertTrue(self.VSS.pop['days_keeping_radian'].max() <= max_days_keeping_radian)
        self.assertTrue((self.VSS.pop['situation'].value_counts() == [pop - 1, 1]).all())
        self.assertListEqual(list(self.VSS.pop['situation'].value_counts().index), ["not_infected", "infected"])
        self.assertListEqual(list(self.VSS.pop.loc[0])[:9],
                             [xlim / 2, ylim / 2, speed_avg, False, 0, False, False, False, False])
        self.assertEqual(self.VSS.pop.loc[0, 'situation'], "infected")

        # check stats
        self.assertEqual(self.VSS.stats.shape, (1, 5))
        self.assertTrue(self.VSS.stats.index == ["Day_0"])
        self.assertListEqual(list(self.VSS.stats.columns),
                             ["nb_not_infected", "nb_infected", "nb_sick", "nb_recovered", "nb_dead"])

    def test_step(self):
        # call function
        self.VSS.step()

        # check incrementation nb_step
        self.assertEqual(self.VSS.nb_step, 1)

        # check positions
        self.assertGreaterEqual(self.VSS.pop['x'].min(), 0)
        self.assertLessEqual(self.VSS.pop['x'].max(), xlim)
        self.assertGreaterEqual(self.VSS.pop['y'].min(), 0)
        self.assertLessEqual(self.VSS.pop['y'].max(), ylim)

    def test_movement(self):
        """new_day() and contamination() test in respective function test"""
        # call function
        self.VSS.movement()

        # check incrementation
        self.assertEqual(self.VSS.mvt, 1)
        self.assertEqual(self.VSS.mvt_in_day, 1)

        # call function having VSS.mvt_in_day == VSS.step_per_day
        # minus 1 because step_per_day will be increase firstly
        self.VSS.mvt_in_day = step_per_day - 1
        self.VSS.movement()

        # check incrementation
        self.assertEqual(self.VSS.mvt_in_day, 0)
        self.assertEqual(self.VSS.day, 1)

    def test_contamination(self):
        """
        Check how contamination work well.
        Create a fake pop with choose specific positions

        create specific pop to test contamination function
        pop having 2 infected points :
        ----> 1 isolated point (will not contaminate)
        ----> 1 points having neighbours
              ----> 2 points not infected
              ----> 1 points not infected at distance of 1 (no contagion risk)
              ----> 1 point dead
              ----> 1 point recovered (can't be infected again)
        -----> 3 other isolated points (not_infected * 2, recovered)
        """
        self.VSS.pop = pd.DataFrame({
            "x": [0, .5, .5, 1, 1.5, 1.5, 2, 2, 3, 4],
            "y": [0, .25, 1.25, 1, 1, 1.5, 1, 2, 2, 1],
            "healthy": [True, True, False, False, True, True, True, True, True, False],
            "infected_day": [-1, -1, 2, 25, -1, 4, -1, 6, -1, 18],
            "dead": [False, False, True, False, False, False, False, False, False, False],
            "recovered": [False, False, False, False, False, True, False, True, False, False],
            "situation": ["not_infected", "not_infected", "dead", "infected", "not_infected",
                          "recovered", "not_infected", "recovered", "not_infected", "sick"]
        })
        # call function
        # self.VSS.contamination()

    def test_get_nb_caring_points_around(self):
        # TODO
        # self.VSS.get_nb_caring_points_around(pos, caring_points)
        pass

    def test_new_day(self):
        # TODO
        self.VSS.new_day()

    def test_update_spreading_counters(self):
        # TODO
        self.VSS.update_spreading_counters()

    def test_update_dead_points(self):
        # TODO
        self.VSS.update_dead_points()

    def test_update_infected_to_sick_points(self):
        # TODO
        self.VSS.update_infected_to_sick_points()

    def test_update_recovered_points(self):
        # TODO
        self.VSS.update_recovered_points()

    def test_update_radian_points(self):
        # TODO
        self.VSS.update_radian_points()

    def test_get_stats(self):
        # TODO
        self.VSS.get_stats()
