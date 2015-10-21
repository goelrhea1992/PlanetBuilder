package pb.g9;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.*;

public class Player implements pb.sim.Player {

	// iteration number
	int iteration = 1;

	// used to pick asteroid and velocity boost randomly
	private Random random = new Random();

	// current time, time limit
	private long time = -1;
	private long time_limit = -1;

	// time until next push
	private long time_of_push = 0;

	// number of retries
	private int retries_per_turn = 2;
	private int turns_per_retry = 0;

	private Point origin = new Point(0, 0);

	private onePush bestpush;
	private double Total_mass = 0;
	private double target_mass = 0;
	private double max_mass = 0;
	private long max_period = 0;
	private long wait_time = 100;
	private int num_asteroids = 0;

	int sink = -1;
	int asteroidsToConsider = -1;

	/**
	* Returns an array list of all intersection points between the orbits of asteroids a and b.
	*/
	private ArrayList<Point> find_intersection(Asteroid a, Asteroid b, HashSet<Long> timelist) {
		ArrayList<Point> intersection_list = new ArrayList<Point>();
		long period = a.orbit.period();

		if (period > 3650)
			period = 3650;
		double r = a.radius() + b.radius();
		for (long ft = 0 ; ft <= period ; ++ft) {
			long t = time + wait_time + ft;
			Point p1 = new Point();
			Point c = new Point();
			a.orbit.positionAt(t - a.epoch, p1);
			b.orbit.center(c);
			Point foci = new Point(c.x * 2, c.y * 2);
			double dist = Point.distance(p1, foci) + Point.distance(p1, new Point(0, 0));
			if (Math.abs(dist - b.orbit.a * 2) < r) {
				intersection_list.add(p1);
				timelist.add(t);
			}
		}
		return intersection_list;
	}

	/**
	* Returns an array list of all intersection points between the orbits of asteroids a and b.
	*
	* Faster version of the algorithm.
	*/
	private ArrayList<Point> fast_find_intersection(Asteroid a, Asteroid b, HashSet<Long> timelist) {
		ArrayList<Point> intersection_list = new ArrayList<Point>();
		long period = a.orbit.period();
		long origin_period = a.orbit.period();

		if (period  > max_period * 2) {
			return intersection_list;
		}

		long time_left = time_limit - time + wait_time;
		long time_for_finding_collision = (long) (time_left / ((target_mass - a.mass)/b.mass));
		time_for_finding_collision = time_left < time_for_finding_collision? time_left : time_for_finding_collision;
		long multiple = (long) (0.1 * time_left/time_for_finding_collision);
		if (multiple < 1) 
			multiple = 1;
		if (period > multiple * time_for_finding_collision)
			period = multiple * time_for_finding_collision;
		double r = a.radius() + b.radius();
		Point c = new Point();
		b.orbit.center(c);
		Point foci = new Point(c.x * 2, c.y * 2);
		long step_length = Math.round(Math.sqrt(period));
		double[] diff = new double[3];
		
		long t = time + wait_time;
		Point p1 = new Point();
		a.orbit.positionAt(t - a.epoch, p1);

		double dist = Point.distance(p1, foci) + Point.distance(p1, new Point(0, 0));
		diff[1] = dist - b.orbit.a * 2;
		diff[0] = diff[1];
		long ft;
		for (ft = step_length ; ft <= period ; ft += step_length) {
			t = time + wait_time + ft;
			p1 = new Point();
			a.orbit.positionAt(t - a.epoch, p1);
			dist = Point.distance(p1, foci)+Point.distance(p1,new Point(0, 0));
			diff[2] = dist - b.orbit.a * 2;
			if (diff[1] >= 0 && diff[2] <0 || diff[1] <= 0 && diff[2] >0 ) {
				for(long ft1 = ft-step_length;ft1 < ft;ft1++) {
					t = time + wait_time + ft1;
					p1 = new Point();
					a.orbit.positionAt(t - a.epoch, p1);
					dist = Point.distance(p1, foci) + Point.distance(p1, new Point(0, 0));
					if (Math.abs(dist - b.orbit.a * 2) < r) {
						intersection_list.add(p1);
						timelist.add(t % origin_period);
						break;
					}
				}
			}
			else if (diff[1] <= diff[2] && diff[1] <= diff[0]) {
				for (long ft1 = ft - 2 * step_length; ft1 < ft; ft1++) {
					t = time + wait_time + ft1;
					p1 = new Point();
					a.orbit.positionAt(t - a.epoch, p1);
					dist = Point.distance(p1, foci) + Point.distance(p1, new Point(0, 0));
					if (Math.abs(dist - b.orbit.a * 2) < r) {
						intersection_list.add(p1);
						timelist.add(t % origin_period);
						break;
					}
				}
			}
			diff[0] = diff[1];
			diff[1] = diff[2];
		}

		t = time + wait_time + period;
		p1 = new Point();
		a.orbit.positionAt(t - a.epoch, p1);		
		dist = Point.distance(p1, foci) + Point.distance(p1, new Point(0, 0));
		diff[2] = dist - b.orbit.a * 2;
		if (diff[1] >= 0 && diff[2] < 0 || diff[1] <= 0 && diff[2] > 0) {
			for (long ft1 = ft - step_length; ft1 < period; ft1++) {
				t = time + wait_time + ft1;
				p1 = new Point();
				a.orbit.positionAt(t - a.epoch, p1);
				dist = Point.distance(p1, foci) + Point.distance(p1, new Point(0, 0));
				if (Math.abs(dist - b.orbit.a * 2) < r) {
					intersection_list.add(p1);
					timelist.add(t % origin_period);
				}
			}
		}
		else if (diff[1] <= diff[2] && diff[1] <= diff[0]) {
			for (long ft1 = ft - 2 * step_length; ft1 < period; ft1++) {
				t = time + wait_time + ft1;
				p1 = new Point();
				a.orbit.positionAt(t - a.epoch, p1);
				dist = Point.distance(p1, foci) + Point.distance(p1, new Point(0, 0));
				if (Math.abs(dist - b.orbit.a * 2) < r) {
					intersection_list.add(p1);
					timelist.add(t % origin_period);
					break;
				}
			}
		}
		return intersection_list;
	}

	/**
	* Returns an array of 30 points - (Energy, Angle) based on which two asteroids to push, 
	* and if the player wants a bigger search space or a smaller search space.
	* The boolean variable bigger will usually be set when the player is not looking ahead much into the future, 
	* and is willing to make high energy pushes as there is limited time left.
	*/
	private double[][] find_search_space(Asteroid[] asteroids, int i, int j, boolean bigger) {
		double p = 1.5;
		if (bigger)
			p = 3;
		int num = 30;
		double[][] space = new double[num][2];
		Point v = asteroids[i].orbit.velocityAt(time + wait_time
				- asteroids[i].epoch);
		double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
		double d1 = Math.atan2(v.y, v.x);
		double E_p = asteroids[i].mass * Orbit.GM / (asteroids[i].orbit.a);
		double e_j= Math.sqrt(asteroids[j].orbit.a * asteroids[j].orbit.a
				- asteroids[j].orbit.b * asteroids[j].orbit.b) / asteroids[j].orbit.a;
		double d2;
		if ((1+e_j) * asteroids[j].orbit.a <= asteroids[i].orbit.a) {
			d2 = Math.PI + d1;
			double E_j1, E_j2;
			E_j1 = E_p - Math.abs(asteroids[i].mass * Orbit.GM / (asteroids[i].orbit.a + (1+e_j) * asteroids[j].orbit.a));
			E_j2 = E_p - Math.abs(asteroids[i].mass * Orbit.GM / (asteroids[i].orbit.a + (1-e_j) * asteroids[j].orbit.a / p));
			double v_j1 = Math.sqrt(2 * E_j1 / asteroids[i].mass);
			v_j1 = Math.abs(v_j1 - v1);
			double v_j2 = Math.sqrt(2 * E_j2 / asteroids[i].mass);
			v_j2 = Math.abs(v_j2 - v1);
			E_j1 = 0.5 * asteroids[i].mass * v_j1 * v_j1;
			E_j2 = 0.5 * asteroids[i].mass * v_j2 * v_j2;
			double E_diff = E_j2-E_j1;
			for(int k = 0; k < num; k++) {
				space[k][0] = E_j1 + k * E_diff/num;
				space[k][1] = d2;
			}
			return space;
		}
		else if ((1-e_j) * asteroids[j].orbit.a >= asteroids[i].orbit.a) {
			d2 = d1;
			double E_j1, E_j2;
			E_j1 = E_p - Math.abs(asteroids[i].mass * Orbit.GM / (asteroids[i].orbit.a + (1-e_j) * asteroids[j].orbit.a));
			E_j2 = E_p - Math.abs(asteroids[i].mass * Orbit.GM / (asteroids[i].orbit.a + (1+e_j) * asteroids[j].orbit.a * p));
			double v_j1 = Math.sqrt(2 * E_j1 / asteroids[i].mass);
			v_j1 = Math.abs(v_j1 - v1);
			double v_j2 = Math.sqrt(2 * E_j2 / asteroids[i].mass);
			v_j2 = Math.abs(v_j2 - v1);
			E_j1 = 0.5 * asteroids[i].mass * v_j1 * v_j1;
			E_j2 = 0.5 * asteroids[i].mass * v_j2 * v_j2;
			double E_diff = E_j2-E_j1;
			for(int k = 0; k < num; k++) {
				space[k][0] = E_j1 + k * E_diff/num;
				space[k][1] = d2;
			}
			return space;
		}
		else {
			d2 = d1;
			double E_j = E_p - Math.abs(asteroids[i].mass * Orbit.GM / (asteroids[i].orbit.a + (1+e_j) * asteroids[j].orbit.a * p));
			double v_j = Math.sqrt(2*E_j/asteroids[i].mass);
			v_j = Math.abs(v_j - v1);
			E_j = 0.5 * asteroids[i].mass * v_j * v_j;
			for(int k = 0; k < num/2; k++) {
				space[k][0] = k * E_j/num;
				space[k][1] = d2;
			}
			
			d2 = Math.PI + d1;
			E_j= E_p - Math.abs(asteroids[i].mass * Orbit.GM / (asteroids[i].orbit.a + (1-e_j) * asteroids[j].orbit.a / p));
			v_j = Math.sqrt(2 * E_j/asteroids[i].mass);
			v_j = Math.abs(v_j - v1);
			E_j = 0.5 * asteroids[i].mass * v_j * v_j;
			for(int k = 0; k < num/2; k++) {
				space[k][0] = k * E_j/num;
				space[k][1] = d2;
			}
			return space;
		}
	}

	/**
	* Returns the index of the chosen sink asteroid based on the following:
	* 
	* Finds the shortest range [i, j] of orbits that constitute at least 50% mass, 
	* and has minimum Sum (over k) (1/r_k - 1/r_median), where k = i..j, and median = (j - 1) / 2.
	* 
	* Returns one of the indices in the optimal range (based on a similar metric) to serve as the sink.
	*/
	private int findTargetForSink(Asteroid[] asteroids, double target_mass) {
		Comparator<Map.Entry<Asteroid, Integer>> comparator = new Comparator<Map.Entry<Asteroid, Integer>>() {
			public int compare(Map.Entry<Asteroid, Integer> o1, Map.Entry<Asteroid, Integer> o2) {
				return (new Double(o1.getKey().orbit.a)).compareTo( o2.getKey().orbit.a );
			}
		};

		HashMap<Asteroid, Integer> asteroidSorted = new HashMap<Asteroid, Integer>();
		for (int i = 0; i < asteroids.length; i++) {
			asteroidSorted.put(asteroids[i], i);
		}
		ArrayList<Map.Entry<Asteroid, Integer>>  asteroidlist =	new ArrayList<Map.Entry<Asteroid, Integer>>(asteroidSorted.entrySet());
		Collections.sort( asteroidlist, comparator);
		
		double sum_mass = 0;
		double min_value = Double.MAX_VALUE;
		int i = 0, j = 0;
		int target_index = -1;

		long max_period = 0;
		long time_left = time_limit;
		long time_for_finding_collision = (long) (time_left / (asteroids.length * target_mass / Total_mass));
		long multiple = (long) (0.1 * time_left/time_for_finding_collision);
		if (multiple < 1) 
			multiple = 1;
		max_period = multiple*time_for_finding_collision;		
		
		for(;i<asteroids.length;) {
			for (;sum_mass < target_mass && j < asteroids.length; j++) {
				sum_mass += asteroidlist.get(j).getKey().mass;
			}
			if (sum_mass < target_mass)
				break;
			for (int k = i; k < j; k++) {
				if (asteroidlist.get(k).getKey().orbit.period() > max_period && target_index != -1)
					break;
				double sum_value = 0;
				for(int kk = i; kk < j; kk++) {
					sum_value += asteroidlist.get(kk).getKey().mass*
							Math.abs(1/asteroidlist.get(kk).getKey().orbit.a
									-1/asteroidlist.get(k).getKey().orbit.a);
				}
				if (sum_value < min_value) {
					target_index = asteroidlist.get(k).getValue();
					min_value = sum_value;
				}
			}
			sum_mass -=asteroidlist.get(i).getKey().mass;
			i++;
		}
		return target_index;
	}
	
	public void init(Asteroid[] asteroids, long time_limit) {

		if (Orbit.dt() != 24 * 60 * 60)
			throw new IllegalStateException("Time quantum is not a day");
		this.time_limit = time_limit;
		
		for (int i=0;i<asteroids.length;i++) {
			Total_mass += asteroids[i].mass;
			if (asteroids[i].mass > max_mass)
				max_mass = asteroids[i].mass;
			if (asteroids[i].orbit.period() > max_period)
				max_period = asteroids[i].orbit.period();
		}
		if (max_mass > .5*Total_mass)
			target_mass = Total_mass;
		else
			target_mass = .5*Total_mass;
		num_asteroids = asteroids.length;
		bestpush = new onePush(-1, Double.MAX_VALUE, 1.0, 0, 0, 0, 0, 0);
	}

	/**
	* Returns the index of the sink asteroid, assuming the latest collided asteroid has the maximum id.
	*/
	public int getSinkIndex(Asteroid[] asteroids) {
		long maxID = Long.MIN_VALUE;
		int maxIndex = -1;

		for (int i = 0; i < asteroids.length; i++) {
			if (asteroids[i].id > maxID) {
				maxID = asteroids[i].id;
				maxIndex = i;
			}
		}
		return maxIndex;
	}

	public void play(Asteroid[] asteroids, double[] energy, double[] direction) {

		if (asteroids.length < num_asteroids) {
			bestpush = new onePush(-1, Double.MAX_VALUE, 1.0, 0, 0, 0, 0, 0);
		}
		num_asteroids = asteroids.length;
		time++;
		if (time == bestpush.time) {
		 	int i = bestpush.i;
		 	energy[i] = bestpush.energy;
		 	direction[i] = bestpush.direction;
		 	
		 	// do not push again until collision happens
		 	time_of_push = bestpush.collision_time + 1;
		 	bestpush = new onePush(-1, Double.MAX_VALUE, 1.0, 0, 0, 0, 0, 0);
		 	iteration++;
		 	return;
		}
			
		// if not yet time to push do nothing
		if (time <= time_of_push)
			return;

		boolean found = false;
		int j;

		// first iteration
		if (iteration == 1) {
			if (sink == -1) {
				j = findTargetForSink(asteroids, target_mass);
				sink = j;
			}
			else
				j = sink;
		}
		else
			j = getSinkIndex(asteroids);
		
		if (asteroids.length < 10)
			// if fewer than 10 asteroids, consider all.
			asteroidsToConsider = asteroids.length - 1;
		else
			// if more than 10 asteorids, consider 20%
			asteroidsToConsider = asteroids.length/5;

		for (int retry = 1; retry <= retries_per_turn; ++retry) {

			List<Integer> favorableAsteroidsOrbitDistance = getKHighestWeightOrbitDistance(asteroids, asteroidsToConsider, asteroids[j]);
			Set<Integer> set = new HashSet<Integer>();
        	set.addAll(favorableAsteroidsOrbitDistance);
        	List<Integer> favorableAsteroids = new ArrayList<Integer>(set);

			for (int i : favorableAsteroids) {
				if (i == j)
					continue;
				double[][] search_space;
				if (retry > 1) {
					// Using bigger search space
					search_space = find_search_space(asteroids, i, j, true);
				}
				else {
					// Using smaller search space
					search_space = find_search_space(asteroids, i, j, false);
				}

				Asteroid a2 = asteroids[j];
				long time_left = time_limit - time + wait_time;
				long time_for_finding_collision = (long) (time_left / ((target_mass - a2.mass)/asteroids[i].mass));
				time_for_finding_collision = time_left < time_for_finding_collision? time_left : time_for_finding_collision;
				long multiple = (long) (0.1 * time_left/time_for_finding_collision);

				if (multiple < 1) 
					multiple = 1;
				
				Point p1 = new Point();
				Point p2 = new Point();

				
				for (int k = 0; k <search_space.length; k ++) {
					double E = search_space[k][0];
					double d2 = search_space[k][1];
					// try to push asteroid
					Asteroid a1 = null;
					try {
						a1 = Asteroid.push(asteroids[i], time + wait_time
								- asteroids[i].epoch, E, d2);
					} catch (InvalidOrbitException e) {
						continue;
					}

					HashSet<Long> timelist = new HashSet<Long>();
					ArrayList<Point> intersections = fast_find_intersection(a1,
							a2, timelist);
					if (intersections.size() == 0)
						continue;
					double r = a1.radius() + a2.radius();
				
					for (long ft = 0; ft < multiple * time_for_finding_collision; ++ft) {

						long t = time + wait_time + ft;
						if (t >= time_limit)
							break;
						if (!timelist.contains(t%a1.orbit.period()))
							continue;
						a1.orbit.positionAt(t - a1.epoch, p1);
						a2.orbit.positionAt(t - a2.epoch, p2);

						// if collision, compare it with the best push
						if (Point.distance(p1, p2) < r) {
							found = true;
							onePush push = new onePush(time + wait_time, E, asteroids[i].mass,
									d2, i, a1.id, t, a2.id);
							if (E / push.mass < bestpush.energy/bestpush.mass) {
								if (CheckIncidentalCollisions(push, asteroids, a1)) {
									bestpush = push;
								}
								System.out.println("  Collision prediction !");
								System.out.println("  Year: " + (1 + t / 365)+"  Day: " + (1 + t % 365));
							}
						}
					}
				}
			}
			if (found) 
				break;
		}
		time_of_push = time + turns_per_retry;
	}

	/**
	* Checks if there is any incidental collisions (other than the one that is desired) resulting from push.
	*/
	private boolean CheckIncidentalCollisions(onePush push, Asteroid[] asteroids, Asteroid a1) {
		Asteroid newAst = a1;
		
		Point newAstPos;
		Point otherAstPos;
		
		for(int i = 0; i < asteroids.length; i++) {
			long pushTime = push.time;
			if(asteroids[i].id == newAst.id || asteroids[i].id == push.asteroidCollidedId) {
				continue;
			}
			HashSet<Long> timelist = new HashSet<Long>();
			fast_find_intersection(newAst, asteroids[i], timelist);
			for(;pushTime <= push.collision_time; pushTime++) {

				if (!timelist.contains(pushTime % a1.orbit.period()))
					continue;
				newAstPos   = newAst.orbit.positionAt(pushTime - newAst.epoch) ;
				otherAstPos = asteroids[i].orbit.positionAt(pushTime - asteroids[i].epoch); 
				if(Point.distance(newAstPos, otherAstPos) < newAst.radius() + asteroids[i].radius()) {
					System.out.println("Bad Push: Incidental Collision Detected with Asteroid ");
					return false;
				}
			}
		}
		return true;
	}

	/** 
	* Returns the indices of the K asteroids with the maximum distance to mass ratio, multiplied by 1/distance between orbits.
	* (Higher the distance mass ratio and lesser the orbit distance, more "pushable" the asteroid - costs less energy)
	**/
	private List<Integer> getKHighestWeightOrbitDistance(Asteroid[] asteroids, int k, Asteroid target) {
		
		double[] weights = new double[asteroids.length];
		double distance = 0.0;
		double distanceToTarget = 0.0;
		double diistanceBetweenOrbits = 0.0;
		double mass = 0.0;

		for (int i = 0; i < asteroids.length; i++) {
			Asteroid a = asteroids[i];
			Point a_center = new Point();
			a.orbit.positionAt(time - a.epoch, a_center);
			Point target_center = new Point();
			target.orbit.positionAt(time - target.epoch, target_center);
			distance = Double
					.valueOf(Point.distance(a_center, new Point(0, 0)));
			distanceToTarget = Double
					.valueOf(Point.distance(a_center, target_center));
			diistanceBetweenOrbits = Math.abs(1/a.orbit.a - 1/target.orbit.a);
			mass = asteroids[i].mass;
			// weights[i] = distance / (mass * distanceToTarget);
			weights[i] = distance / (mass * diistanceBetweenOrbits);
		}

		List<Integer> highestKIndices;
		double[] temp = new double[k];
		double minWeight = Integer.MAX_VALUE;
		int minIndex = -1;
		Integer[] tempIndexInWeights = new Integer[k];

		// put k weights in temp, and keep track of max seen so far
		for (int i = 0; i < k; i++) {
			temp[i] = weights[i];
			tempIndexInWeights[i] = i;
			if (temp[i] < minWeight) {
				minWeight = temp[i];
				minIndex = i;
			}
		}

		// check remaining weights and maintain lowest k in temp
		for (int i = k; i < weights.length; i++) {
			if (weights[i] > minWeight) {
				// replace the max so far with this value
				temp[minIndex] = weights[i];
				tempIndexInWeights[minIndex] = i;

				//find new max
				minWeight = Double.MAX_VALUE;
				minIndex = -1;

				for (int j = 0; j< k; j++) {
					if (temp[j] < minWeight) {
						minWeight = temp[j];
						minIndex = j;
					}
				}
			}
		}
		highestKIndices = new ArrayList<Integer>(Arrays.asList(tempIndexInWeights));
		return highestKIndices;
	}
}

class onePush {

	// the time the push is made
	long time;

	// energy used to make the push
	double energy;
	
	// mass of the asteroids to push
	double mass;

	// the direction of the push
	double direction;

	// the time at which the push will cause a collision
	long collision_time;
	
	// index of the asteroid to push in the array asteroids[]
	int i;
	
	// id of the asteroid to push
	long asteroidPushId;
	
	// id of the asteroid to push towards
    long asteroidCollidedId;
	
	public onePush(long time, double energy, double mass, double direction, int i, long pushId,
			long collision_time, long collidedId) {
		this.time = time;
		this.energy = energy;
		this.mass = mass;
		this.direction = direction;
		this.i = i;
		this.asteroidPushId = pushId;
		this.collision_time = collision_time;
		this.asteroidCollidedId = collidedId;
		
	}
}


