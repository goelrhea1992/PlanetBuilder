package pb.g9;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.*;

public class Player implements pb.sim.Player {

	//iteration number
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

	private Point origin = new Point(0,0);

	private ArrayList<Point> find_intersection(Asteroid a, Asteroid b){
		ArrayList<Point> intersection_list = new ArrayList<Point>();
		double r = a.radius() + b.radius();
		long period = a.orbit.period();
		for (long ft = 0 ; ft != period ; ++ft) {
			long t = time + ft;
			Point p1 = new Point();
			Point c = new Point();
			a.orbit.positionAt(t - a.epoch, p1);
			b.orbit.center(c);
			Point foci = new Point(c.x*2,c.y*2);
			double dist = Point.distance(p1,foci)+Point.distance(p1,foci);
			if (Math.abs(dist - b.orbit.a*2) < r){
				intersection_list.add(p1);
			}
				
		}
		return intersection_list;
	}
	// print orbital information
	public void init(Asteroid[] asteroids, long time_limit)
	{
		if (Orbit.dt() != 24 * 60 * 60)
			throw new IllegalStateException("Time quantum is not a day");
		this.time_limit = time_limit;
	}
	
	// try to push asteroid
	public void play(Asteroid[] asteroids,
	                 double[] energy, double[] direction)
	{

		// if not yet time to push do nothing
		if (++time <= time_of_push) return;
		System.out.println("Year: " + (1 + time / 365));
		System.out.println("Day: "  + (1 + time % 365));
		for (int retry = 1 ; retry <= retries_per_turn ; ++retry) {
			// pick a random asteroid and get its velocity
			// int i = random.nextInt(asteroids.length);

			int i = getHighestWeightAsteroid(asteroids);

			Point v = asteroids[i].orbit.velocityAt(time);
			// add 5-50% of current velocity in magnitude
			System.out.println("Try: " + retry + " / " + retries_per_turn);
			double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
			for (double k=0; k< 5; k=k+0.1)
			{
				double v2 = v1 * (k * 0.45 + 0.05);
				System.out.println("  Speed: " + v1 + " +/- " + v2);
				
				// apply push at -π/8 to π/8 of current angle
				double d1 = Math.atan2(v.y, v.x);
				double d2 = d1 + (random.nextDouble() - 0.5) * Math.PI * 0.25;
				System.out.println("  Angle: " + d1 + " -> " + d2);
				
				// compute energy
				double E = 0.5 * asteroids[i].mass * v2 * v2;
				
				// try to push asteroid
				Asteroid a1 = null;
				try 
				{
					a1 = Asteroid.push(asteroids[i], time, E, d2);
				} catch (InvalidOrbitException e) 
				{
					System.out.println("  Invalid orbit: " + e.getMessage());
					continue;
				}
				
				// avoid allocating a new Point object for every position
				Point p1 = v, p2 = new Point();
				
				if (iteration ==1)
				{
				// search for collision with any other asteroids
				for (int j = 0 ; j != asteroids.length ; ++j) 
				{
					if (i == j) continue;
					Asteroid a2 = asteroids[j];
					ArrayList<Point> intersections = find_intersection(a1,a2);
					if (intersections.size() == 0) continue;
					double r = a1.radius() + a2.radius();
					// look 10 years in the future for collision
					for (long ft = 0 ; ft != 3650 ; ++ft) 
					{
						long t = time + ft;
						if (t >= time_limit) break;
						a1.orbit.positionAt(t - a1.epoch, p1);
						a2.orbit.positionAt(t - a2.epoch, p2);
						
						// if collision, return push to the simulator
						if (Point.distance(p1, p2) < r) 
						{
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
				if(iteration > 1)
				{
					int j = getLeastWeightAsteroid(asteroids);
					Asteroid a2 = asteroids[j];
					ArrayList<Point> intersections = find_intersection(a1,a2);
					if (intersections.size() == 0) continue;
					double r = a1.radius() + a2.radius();
					
					// look 10 years in the future for collision
					for (long ft = 0 ; ft != 3650 ; ++ft) 
					{
						long t = time + ft;
						if (t >= time_limit) break;
						a1.orbit.positionAt(t - a1.epoch, p1);
						a2.orbit.positionAt(t - a2.epoch, p2);
						
						// if collision, return push to the simulator
						if (Point.distance(p1, p2) < r) 
						{
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
				System.out.println("  No collision ...");
			}
		}
		time_of_push = time + turns_per_retry;
	}
	private int getLeastWeightAsteroid(Asteroid asteroids[]) 
	{
		int min =0;
		double distance = 0.0;
		double mass = 0.0;
		double weight = 0.0;
		double minWeight = 0.0;
		
		for(int i =0; i<asteroids.length; i++)
		{
			Asteroid a = asteroids[i];
			Point a_center = new Point();
			a.orbit.positionAt(time - a.epoch, a_center);
			distance = Double.valueOf(Point.distance(a_center, new Point(0,0)));
			mass = asteroids[i].mass;
			weight = distance/mass;
			if(weight <= minWeight) {
				min = i;
				minWeight = weight;
			}
		}
		return min;
	}

	private int getHighestWeightAsteroid(Asteroid asteroids[]) 
	{
		int max =0;
		double distance = 0.0;
		double mass = 0.0;
		double weight = 0.0;
		double maxWeight = 0.0;

		for(int i =0; i<asteroids.length; i++)
		{
			Asteroid a = asteroids[i];
			Point a_center = new Point();
			a.orbit.positionAt(time - a.epoch, a_center);
			distance = Double.valueOf(Point.distance(a_center, new Point(0,0)));
			mass = asteroids[i].mass;
			weight = distance/mass;
			if(weight > maxWeight) {
				max = i;
				maxWeight = weight;
			}
		}
		return max;
	}
}
