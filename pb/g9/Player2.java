package pb.g9;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.Random;
import java.util.*;

public class Player implements pb.sim.Player {

	int iteration =1;

	// used to pick asteroid and velocity boost randomly
	private Random random = new Random();

	// current time, time limit
	private long time = -1;
	private long time_limit = -1;

	// time until next push
	private long time_of_push = 0;

	// number of retries
	private int retries_per_turn = 1;
	private int turns_per_retry = 3;

	public Map<Asteroid, Integer> asteroidToMassRank = null;
	public Map<Asteroid, Integer> asteroidToDistRank = null;
	public Map<Asteroid, Double> asteroidToWeight = null;
	public Map<Asteroid, Integer> asteroidToIndex = null;

	// print orbital information
	public void init(Asteroid[] asteroids, long time_limit)
	{
		if (Orbit.dt() != 24 * 60 * 60)
			throw new IllegalStateException("Time quantum is not a day");
		this.time_limit = time_limit;
		storeIndexOfAsteroid(asteroids);
	}

	public void storeIndexOfAsteroid(Asteroid[] asteroids) {
		asteroidToIndex = new HashMap<Asteroid, Integer>();
		int i = 0;
		for (Asteroid a: asteroids) {
			asteroidToIndex.put(a, i);
			i++;
		}
	}

	public void rankAsteroidsByMass(Asteroid[] asteroids) {

		Arrays.sort(asteroids, new Comparator<Asteroid>() {
			public int compare(Asteroid a, Asteroid b) {
				return Double.valueOf(b.mass).compareTo(Double.valueOf(a.mass));
			}
		});

		asteroidToMassRank = new HashMap<Asteroid, Integer>();
		for (int i = 0; i < asteroids.length; i++) {
			asteroidToMassRank.put(asteroids[i], i+1);
		}
	}

	public void rankAsteroidsByDist(Asteroid[] asteroids) {

		Arrays.sort(asteroids, new Comparator<Asteroid>() {
			public int compare(Asteroid a, Asteroid b) {
				Point a_center,  b_center;
				a_center = new Point();
				b_center = new Point();
				a.orbit.positionAt(time - a.epoch, a_center);
				b.orbit.positionAt(time - b.epoch, b_center);

				return Double.valueOf(Point.distance(a_center, new Point(0,0))).compareTo(Double.valueOf(Point.distance(b_center, new Point(0,0))));
			}
		});

		asteroidToDistRank = new HashMap<Asteroid, Integer>();
		for (int i = 0; i< asteroids.length; i++) {
			asteroidToDistRank.put(asteroids[i], i+1);
		}
	}

	public void assignWeights() {
		asteroidToWeight = new HashMap<Asteroid, Double>();
		for (Map.Entry<Asteroid, Integer> entry: asteroidToMassRank.entrySet()) {
			asteroidToWeight.put(entry.getKey(), ((0.5) * entry.getValue() + (0.5) * asteroidToDistRank.get(entry.getKey())));
		}
	}

	public Asteroid getHighestWeightAsteroid() {
		Double max = Double.MIN_VALUE;
		Asteroid highest = null;
		for (Map.Entry<Asteroid, Double> entry: asteroidToWeight.entrySet()) {
			if (entry.getValue() > max) {
				max = entry.getValue();
				highest = entry.getKey();
			}
		}
		return highest;
	}

	public Asteroid getLowestWeightAsteroid() {
		Double min = Double.MAX_VALUE;
		Asteroid lowest = null;
		for (Map.Entry<Asteroid, Double> entry: asteroidToWeight.entrySet()) {
			if (entry.getValue() < min) {
				min = entry.getValue();
				lowest = entry.getKey();
			}
		}
		return lowest;
	}

	// try to push asteroid
	public void play(Asteroid[] asteroids,
	                 double[] energy, double[] direction)
	{

		Asteroid asteroidToPush = null;
		// if not yet time to push do nothing
		if (++time <= time_of_push) return;
		// System.out.println("Year: " + (1 + time / 365));
		// System.out.println("Day: "  + (1 + time % 365));
		storeIndexOfAsteroid(asteroids);
		rankAsteroidsByMass(asteroids);
		rankAsteroidsByDist(asteroids);
		assignWeights();

		for (int retry = 1 ; retry <= retries_per_turn ; ++retry) {

			asteroidToPush = getHighestWeightAsteroid();
			int i = asteroidToIndex.get(asteroidToPush);

			Point v = asteroids[i].orbit.velocityAt(time);
			// add 5-50% of current velocity in magnitude
			// System.out.println("Try: " + retry + " / " + retries_per_turn);
			double v1 = Math.sqrt(v.x * v.x + v.y * v.y);

			double v2 = v1 * (random.nextDouble() * 0.45 + 0.05);
			// double v2 = v1 * (random.nextDouble() * 0.15 + 0.05);

			// System.out.println("  Speed: " + v1 + " +/- " + v2);
			// apply push at -π/8 to π/8 of current angle
			double d1 = Math.atan2(v.y, v.x);
			double d2 = d1 + (random.nextDouble() * 0.5) * Math.PI * 0.25;
			// System.out.println("  Angle: " + d1 + " -> " + d2);
			// compute energy
			double E = 0.5 * asteroids[i].mass * v2 * v2;
			// try to push asteroid
			Asteroid a1 = null;
			try {
				a1 = Asteroid.push(asteroids[i], time, E, d2);
			} catch (InvalidOrbitException e) {
				System.out.println("  Invalid orbit: " + e.getMessage());
				continue;
			}
			// avoid allocating a new Point object for every position
			Point p1 = v, p2 = new Point();

			if (iteration == 1) {
				for (int j = 0 ; j != asteroids.length ; ++j) {
					if (i == j) continue;
					Asteroid a2 = asteroids[j];
					double r = a1.radius() + a2.radius();
					// look 10 years in the future for collision
					for (long ft = 0 ; ft != 3650 ; ++ft) {
						long t = time + ft;
						if (t >= time_limit) break;
						a1.orbit.positionAt(t - a1.epoch, p1);
						a2.orbit.positionAt(t - a2.epoch, p2);
						// if collision, return push to the simulator
						if (Point.distance(p1, p2) < r) {
							energy[i] = E;
							direction[i] = d2;
							iteration++;
							// do not push again until collision happens
							time_of_push = t + 1;
							
							System.out.println("  Collision prediction !");
							System.out.println("  Year: " + (1 + t / 365));
							System.out.println("  Day: "  + (1 + t % 365));
							return;
						}
					}
				}
				System.out.println("  No collision ...");
			}
			else {
				Asteroid a2 = getLowestWeightAsteroid();
				double r = a1.radius() + a2.radius();
				// look 10 years in the future for collision
				for (long ft = 0 ; ft != 3650 ; ++ft) {
					long t = time + ft;
					if (t >= time_limit) break;
					a1.orbit.positionAt(t - a1.epoch, p1);
					a2.orbit.positionAt(t - a2.epoch, p2);
					// if collision, return push to the simulator
					if (Point.distance(p1, p2) < r) {
						energy[i] = E;
						direction[i] = d2;
						// do not push again until collision happens
						time_of_push = t + 1;
						System.out.println("  Collision prediction !");
						System.out.println("  Year: " + (1 + t / 365));
						System.out.println("  Day: "  + (1 + t % 365));
						return;
					}
				}
			}
			time_of_push = time + turns_per_retry;
		}
	}
}
