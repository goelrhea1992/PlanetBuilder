package pb.g9;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.*;


public class Utils {


	public static List<Integer> getKHighestWeightEuclideanDistance(Asteroid[] asteroids, int k, Asteroid target) {
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
			diistanceBetweenOrbits = Math.abs(a.orbit.a - target.orbit.a);
			mass = asteroids[i].mass;
			weights[i] = distance / (mass * distanceToTarget);
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

	public static List<Integer> getKLowestWeight(Asteroid[] asteroids, int k) {
		/** 
		* Returns the indices of the K asteroids with the minimum distance to mass ratio.
		* (Higher the distance mass ratio, more "pushable" the asteroid - costs less energy)
		* 
		* Hence, it is not good to push the one with the lowest distance to mass ratio.
		* because it is too close to the sun, and too heavy.
		**/

		double[] weights = new double[asteroids.length];
		double distance = 0.0;
		double mass = 0.0;

		for (int i = 0; i < asteroids.length; i++) {
			Asteroid a = asteroids[i];
			Point a_center = new Point();
			a.orbit.positionAt(time - a.epoch, a_center);
			distance = Double
					.valueOf(Point.distance(a_center, new Point(0, 0)));
			mass = asteroids[i].mass;
			weights[i] = distance / mass;
		}

		List<Integer> lowestKIndices;
		double[] temp = new double[k];
		double maxWeight = Integer.MIN_VALUE;
		int maxIndex = -1;
		Integer[] tempIndexInWeights = new Integer[k];

		// put k weights in temp, and keep track of max seen so far
		for (int i = 0; i < k; i++) {
			temp[i] = weights[i];
			tempIndexInWeights[i] = i;
			if (temp[i] > maxWeight) {
				maxWeight = temp[i];
				maxIndex = i;
			}
		}

		// check remaining weights and maintain lowest k in temp
		for (int i = k; i < weights.length; i++) {
			if (weights[i] < maxWeight) {
				// replace the max so far with this value
				temp[maxIndex] = weights[i];
				tempIndexInWeights[maxIndex] = i;

				//find new max
				maxWeight = Double.MIN_VALUE;
				maxIndex = -1;

				for (int j = 0; j< k; j++) {
					if (temp[j] > maxWeight) {
						maxWeight = temp[j];
						maxIndex = j;
					}
				}
			}
		}
		lowestKIndices = new ArrayList<Integer>(Arrays.asList(tempIndexInWeights));
		return lowestKIndices;
	}

	public static int getHeaviestAsteroidAmong(Asteroid asteroids[], List<Integer> desiredOrbits) {
		/**
		* Returns the heaviest asteroid among the asteroids on desiredOrbits.
		*/

		int max = 0;
		double mass = 0.0;
		double maxMass = Double.MIN_VALUE;

		for (int i: desiredOrbits) {
			mass = asteroids[i].mass;
			if (mass > maxMass) {
				max = i;
				maxMass = mass;
			}
		}
		return max;
	}

	public static int getHeaviestAsteroid(Asteroid asteroids[]) {
		/**
		* Returns the heaviest asteroid among all asteroids.
		*/

		int max = 0;
		double mass = 0.0;
		double maxMass = Double.MIN_VALUE;

		for (int i = 0; i < asteroids.length; i++) {
			mass = asteroids[i].mass;
			if (mass > maxMass) {
				max = i;
				maxMass = mass;
			}
		}
		return max;
	}

	public static int findCandidateOrbitsForSink(Asteroid[] asteroids, double target_mass) {
		/**
		* Finds the shortest range [1, j] of orbits that constitute at least 50% mass, 
		* and have minimum Sum (over k) (1/r_k - 1/r_median), where k = i..j, and median = (j - 1) / 2.
		* 
		* Returns the median of the optimal range to serve as the sink.
		*/

		Map<Integer, Double> asteroidToRadius = new HashMap<Integer, Double>();
		Map<Integer, Double> asteroidToRadiusSorted = new LinkedHashMap<>();
		for (int i = 0; i < asteroids.length; i++) {
			asteroidToRadius.put(i, asteroids[i].orbit.a);
		}

		asteroidToRadiusSorted = sortByValue(asteroidToRadius);

		Map<Integer, Integer> radiusIndexToAsteroid = new HashMap<Integer, Integer>();
		int i = 0;
		double[] radii = new double[asteroids.length];

		for (Map.Entry<Integer, Double> entry: asteroidToRadiusSorted.entrySet()) {
			radiusIndexToAsteroid.put(i, entry.getKey());
			radii[i] = entry.getValue();
			i++;
			
		}

		double minMetricSum = Double.MAX_VALUE;
		int bestI = 0;
		int bestJ = 0;
		int j, k;
		boolean found;

		for (i = 0; i < radii.length; i++) {
			double massSum = 0.0;
			found = false;
			for (j = i; j < radii.length; j++) {
				massSum += asteroids[j].mass;
				if (massSum >= target_mass) {
					found = true;
					break;
				}
			}
			if (found) {
				int median = (j - i)/2;
				double thisMetricSum = 0.0;

				for (k = i; k <= j; k++) {
					thisMetricSum += Math.abs((1/radii[k]) - (1/radii[median]));
				}
				if (thisMetricSum <= minMetricSum) {
					minMetricSum = thisMetricSum;
					bestI = i;
					bestJ = j;
				}
			}
		}

		int bestMedian = (bestJ - bestI) / 2;
		return radiusIndexToAsteroid.get(bestMedian);
	}

	public static <Integer, Double extends Comparable<Double>> Map<Integer, Double> sortByValue( Map<Integer, Double> map )
	{
		/**
		* Sorts a hashmap by value.
		*/

	    List<Map.Entry<Integer, Double>> list =
	        new LinkedList<>( map.entrySet() );
	    Collections.sort( list, new Comparator<Map.Entry<Integer, Double>>()
	    {
	        @Override
	        public int compare( Map.Entry<Integer, Double> o1, Map.Entry<Integer, Double> o2 )
	        {
	            return (o1.getValue()).compareTo( o2.getValue() );
	        }
	    } );

	    Map<Integer, Double> result = new LinkedHashMap<>();
	    for (Map.Entry<Integer, Double> entry : list)
	    {
	        result.put( entry.getKey(), entry.getValue() );
	    }
	    return result;
	}

}
	