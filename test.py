import unittest
import pandas as pd
from pandas.util.testing import assert_frame_equal
from VirusSpreadingSimulation.VirusSpredingSimulation import VirusSpreadingSimulation

# GLOBAL VARIABLES
if True:
    pop = 10
    xlim = 5.
    ylim = 6.
    min_days_keeping_radian = 1
    max_days_keeping_radian = 11
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
    mpd = 5
    fpm = 5


class TestVSS(unittest.TestCase):
    """Test case for class VirusSpreadingSimulation"""

    def setUp(self):
        """Initialization using global variables"""
        self.VSS = VirusSpreadingSimulation(population=pop, xlim=xlim, ylim=ylim, speed_avg=speed_avg,
                                            radius_infection=radius_infection, p_infection=p_infection,
                                            incubation_period=incubation_period, healing_duration=healing_duration,
                                            death_rate=death_rate, without_symptom_rate=without_symptom_rate,
                                            quarantine_zone=quarantine_zone, isolation_threshold=isolation_threshold,
                                            pct_pop_isolated=pct_pop_isolated, mpd=mpd, fpm=fpm)

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
        self.assertEqual(self.VSS.mpd, mpd)
        self.assertEqual(self.VSS.fpm, fpm)

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
        self.assertEqual(self.VSS.pop.shape, (pop, 14))
        self.assertListEqual(list(self.VSS.pop.columns),
                             ['x', 'y', 'speed', 'healthy', 'infected_day', 'quarantine_zone',
                              'stay_confined', 'recovered', 'dead', 'death_day',
                              'radian', 'days_keeping_radian', 'situation', 'without_symptoms'])
        self.assertTrue(self.VSS.pop['x'].min() >= 0)
        self.assertTrue(self.VSS.pop['x'].max() <= xlim)
        self.assertTrue(self.VSS.pop['y'].min() >= 0)
        self.assertTrue(self.VSS.pop['y'].max() <= ylim)
        self.assertTrue(self.VSS.pop['speed'].min() > 0)
        self.assertTrue(self.VSS.pop['healthy'].sum() == pop - 1)
        self.assertTrue((self.VSS.pop['infected_day'].value_counts() == [pop - 1, 1]).all())
        self.assertListEqual(list(self.VSS.pop['infected_day'].value_counts().index), [-1, 0])
        self.assertTrue(self.VSS.pop['quarantine_zone'].sum() == 0)
        self.assertTrue(self.VSS.pop['stay_confined'].isin([False, True]).sum() == pop)
        self.assertTrue(self.VSS.pop['recovered'].sum() == 0)
        self.assertTrue(self.VSS.pop['dead'].sum() == 0)
        self.assertTrue(self.VSS.pop['days_keeping_radian'].min() >= min_days_keeping_radian)
        self.assertTrue(self.VSS.pop['days_keeping_radian'].max() <= max_days_keeping_radian)
        self.assertTrue((self.VSS.pop['situation'].value_counts() == [pop - 1, 1]).all())
        self.assertListEqual(list(self.VSS.pop['situation'].value_counts().index), ["not_infected", "infected"])
        self.assertListEqual(list(self.VSS.pop.loc[0])[:6], [xlim / 2, ylim / 2, speed_avg, False, 0, False])
        self.assertListEqual(list(self.VSS.pop.loc[0])[7:9], [False, False])
        self.assertEqual(self.VSS.pop.loc[0, 'situation'], "infected")
        self.assertEqual(self.VSS.pop.loc[0, 'without_symptoms'], False)
        self.assertTrue(self.VSS.pop['without_symptoms'].isin([False, True]).sum() == pop)

        # check stats
        self.assertEqual(self.VSS.stats.shape, (1, 6))
        self.assertTrue(self.VSS.stats.index == ["Day_0"])
        self.assertListEqual(list(self.VSS.stats.columns),
                             ["nb_not_infected", "nb_infected", "nb_sick", "nb_recovered", "nb_dead", "nb_immune"])

    def test_repr(self):
        """test __repr__ / test print"""
        self.assertEqual(str(self.VSS), '<ViruseSpredingSimulation class>')

    def test_frame(self):
        # call function
        self.VSS.frame()

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

        # call function having VSS.mvt_in_day == VSS.mpd
        # minus 1 because mpd will be increase firstly
        self.VSS.mvt_in_day = mpd - 1
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
        # create fake pop
        pop = pd.DataFrame({
            "x": [0, .5, .5, 1, 1.5, 1.5, 2, 2, 3, 4],
            "y": [0, .25, 1.25, 1, 1, 1.5, 1, 2, 2, 1],
            "healthy": [True, True, False, False, True, True, True, True, True, False],
            "infected_day": [-1, -1, 2, 25, -1, 4, -1, 6, -1, 18],
            "dead": [False, False, True, False, False, False, False, False, False, False],
            "recovered": [False, False, False, False, False, True, False, True, False, False],
            "situation": ["not_infected", "not_infected", "dead", "infected", "not_infected",
                          "recovered", "not_infected", "recovered", "not_infected", "sick"]
        })
        # check contamination if p_infection is 0 => pop have to stay the same
        self.VSS.pop = pop.copy()
        self.VSS.p_infection = 0
        self.VSS.contamination()
        assert_frame_equal(self.VSS.pop, pop)

        # check contamination if p_infection is 1 => points close to infected points are infected
        self.VSS.p_infection = 1
        self.VSS.contamination()
        self.assertFalse(self.VSS.pop.equals(pop))
        pop2 = pop.copy()
        pop2.loc[[1, 4], ["healthy", "infected_day", "situation"]] = [False, self.VSS.day, "infected"]
        assert_frame_equal(self.VSS.pop, pop2)
        self.assertFalse(self.VSS.pop.equals(pop))

        # check contamination if p_infection is 0.5 => points close to infected points may be infected
        self.VSS.pop = pop.copy()
        self.VSS.p_infection = .4
        self.VSS.contamination()
        # check not close points still the same
        assert_frame_equal(self.VSS.pop.loc[~self.VSS.pop.index.isin([1, 4])], pop.loc[~pop.index.isin([1, 4])])
        # check if 2 points are infected/not_infected
        self.assertTrue((list(self.VSS.pop.loc[1]) == [0.5, 0.25, True, -1, False, False, 'not_infected']) |
                        (list(self.VSS.pop.loc[1]) == [0.5, 0.25, False, self.VSS.day, False, False, 'infected']))
        self.assertTrue((list(self.VSS.pop.loc[4]) == [1.5, 1, True, -1, False, False, 'not_infected']) |
                        (list(self.VSS.pop.loc[4]) == [1.5, 1, False, self.VSS.day, False, False, 'infected']))

    def test_get_nb_caring_points_around(self):
        """
        Check how get_nb_caring_points_around_count.
        Create same specific pop as test_contamination but point [1.5, 1.5] is modify to a caring point
        """
        # create fake pop
        pop = pd.DataFrame({
            "x": [0, .5, .5, 1, 1.5, 1.5, 2, 2, 3, 4],
            "y": [0, .25, 1.25, 1, 1, 1.5, 1, 2, 2, 1],
            "healthy": [True, True, False, False, True, False, True, True, True, False],
            "infected_day": [-1, -1, 2, 25, -1, 7, -1, 6, -1, 18],
            "dead": [False, False, True, False, False, False, False, False, False, False],
            "recovered": [False, False, False, False, False, False, False, True, False, False],
            "situation": ["not_infected", "not_infected", "dead", "infected", "not_infected",
                          "sick", "not_infected", "recovered", "not_infected", "sick"]
        })
        caring_points = pop.loc[pop["situation"].isin(["infected", "sick"]), ["x", "y"]]
        self.assertEqual(self.VSS.get_nb_caring_points_around((1.5, 1), caring_points), 2)
        self.assertEqual(self.VSS.get_nb_caring_points_around((2, 1), caring_points), 1)
        self.assertEqual(self.VSS.get_nb_caring_points_around((3, 2), caring_points), 0)
        self.assertEqual(self.VSS.get_nb_caring_points_around((.5, .25), caring_points), 1)

        # with no caring points
        caring_points = pop.loc[pop["situation"] == "no_caring_point_match", ["x", "y"]]
        self.assertEqual(self.VSS.get_nb_caring_points_around((1.5, 1), caring_points), 0)
        self.assertEqual(self.VSS.get_nb_caring_points_around((2, 1), caring_points), 0)
        self.assertEqual(self.VSS.get_nb_caring_points_around((3, 2), caring_points), 0)
        self.assertEqual(self.VSS.get_nb_caring_points_around((.5, .25), caring_points), 0)

        # with only one caring points (only one point have situation "infected" at [1, 1])
        caring_points = pop.loc[pop["situation"] == "infected", ["x", "y"]]
        self.assertEqual(self.VSS.get_nb_caring_points_around((1.5, 1), caring_points), 1)
        self.assertEqual(self.VSS.get_nb_caring_points_around((2, 1), caring_points), 0)
        self.assertEqual(self.VSS.get_nb_caring_points_around((3, 2), caring_points), 0)
        self.assertEqual(self.VSS.get_nb_caring_points_around((.5, .25), caring_points), 1)

    def test_new_day(self):
        # for the moment only call function already tested

        # TODO TEST
        # if self.quarantine_zone:
        # ----> Remove sick people of area

        # TODO TEST
        # Confined points (if self.isolation_threshold > self.nb_sick / self.pop_size)
        # ----> Points stay immobile if isolation_threshold is exceeded
        pass

    def test_update_spreading_counters(self):
        self.assertEqual(self.VSS.nb_not_infected, pop - 1)
        self.assertEqual(self.VSS.nb_infected, 1)
        self.assertEqual(self.VSS.nb_sick, 0)
        self.assertEqual(self.VSS.nb_recovered, 0)
        self.assertEqual(self.VSS.nb_dead, 0)
        self.assertEqual(self.VSS.nb_immune, 0)
        self.VSS.pop = pd.DataFrame({
            "x": [0] * 18,
            "y": [0] * 18,
            "situation": ["not_infected", "not_infected", "dead", "infected", "not_infected",
                          "sick", "not_infected", "recovered", "not_infected", "sick", "dead",
                          "dead", "infected", "infected", "infected", "infected", "infected", "immune"]
        })
        self.VSS.day += 1
        self.VSS.update_spreading_counters()
        self.assertEqual(self.VSS.nb_not_infected, 5)
        self.assertEqual(self.VSS.nb_infected, 6)
        self.assertEqual(self.VSS.nb_sick, 2)
        self.assertEqual(self.VSS.nb_recovered, 1)
        self.assertEqual(self.VSS.nb_dead, 3)
        self.assertEqual(self.VSS.nb_immune, 1)
        assert_frame_equal(self.VSS.stats,
                           pd.DataFrame([[pop - 1, 1, 0, 0, 0, 0], [5, 6, 2, 1, 3, 1]], index=['Day_0', 'Day_1'],
                                        columns=['nb_not_infected', 'nb_infected',
                                                 'nb_sick', 'nb_recovered', 'nb_dead', 'nb_immune']))

    def test_update_dead_points(self):
        """Check points correctly change to dead"""
        # create fake pop
        pop = pd.DataFrame({
            "infected_day": [-1, -1, 2, 5, -1, 4, -1, 10, 20, 17],
            "quarantine_zone": [False, False, True, False, False, False, False, False, False, True],
            "stay_confined": [False] * 10,
            "dead": [False, False, True, False, False, False, False, False, False, False],
            "days_keeping_radian": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            "situation": ["not_infected", "not_infected", "dead", "infected", "not_infected",
                          "recovered", "not_infected", "infected", "infected", "sick"],
            "death_day": [1, 2, 3, 4, 2, 5, 1, 12, 1, 4],
            "without_symptoms": [False] * 10,
        })

        # check a situation
        self.VSS.day = 20
        self.VSS.pop = pop.copy()
        self.VSS.update_dead_points()
        assert_frame_equal(self.VSS.pop, pop)

        # check an other situation
        self.VSS.day = 21
        self.VSS.pop = pop.copy()
        self.VSS.update_dead_points()
        pop_test1 = pop.copy()
        pop_test1.loc[[8, 9], ["quarantine_zone", "stay_confined", "dead", "days_keeping_radian",
                               "situation"]] = [True, True, True, 999999, "dead"]
        assert_frame_equal(self.VSS.pop, pop_test1)

        # check a 3rd situation
        self.VSS.day = 22
        self.VSS.pop = pop.copy()
        self.VSS.update_dead_points()
        pop_test2 = pop.copy()
        pop_test2.loc[7, ["quarantine_zone", "stay_confined", "dead", "days_keeping_radian",
                          "situation"]] = [True, True, True, 999999, "dead"]
        assert_frame_equal(self.VSS.pop, pop_test2)

        # check without symptoms (can't die)
        pop["without_symptoms"] = [True] * pop.shape[0]
        self.VSS.day = 21
        self.VSS.pop = pop.copy()
        self.VSS.update_dead_points()
        assert_frame_equal(self.VSS.pop, pop)

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
        assert_frame_equal(self.VSS.get_stats(), self.VSS.stats)

    def test_get_plot(self):
        # TODO
        pass

    def test_print_situation(self):
        # TODO
        pass
