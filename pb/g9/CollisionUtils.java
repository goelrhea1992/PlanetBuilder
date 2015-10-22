package pb.g9;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.*;

public class CollisionUtils {

	

	public static int getMinDistanceTime(Asteroid a, Asteroid  b) {
		/**
		* Returns the time at which the Euclidean distance between asteroids is minimized.
		**/

		double minDistance = -1;
		int minDistTime = Integer.MAX_VALUE;
		for(int t = 0; t < 20*365; t++){
			double distance = Point.distance(a.orbit.positionAt((long)t - a.epoch),b.orbit.positionAt((long)t - b.epoch));
			if(distance < minDistance){
				minDistance = distance;
				minDistTime = t;
			}			
		}
		return minDistTime;
	} 

	public static void attemptCollisions(Asteroid a, Asteroid b, int minDistTime){
		// Calculate which needs to be pushed
		Point posA = a.orbit.positionAt((long)minDistTime - a.epoch - 365);
		Point posB = b.orbit.positionAt((long)minDistTime - b.epoch);		
		// For now, just push "a"
		boolean didPush = false;
		// Create vector that draws a line from location of push to where I 
		// want the intersection to take place
		Point dirVector = new Point(posB.x - posA.x, posB.y - posA.y);
		double dir = dirVector.direction();
		Point v = a.orbit.velocityAt((long)minDistTime - a.epoch - 365);
		double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
		double r = a.radius() + b.radius();
		Point newPos;
		for(int i = 0; i < 500; i++){
			double v2 = v1 * ((1+i/100f));
			double E = 0.5 * a.mass * v2 * v2;
			Asteroid testAst = pushTest(a, minDistTime - 365, E, dir);
			newPos = testAst.orbit.positionAt((long)minDistTime - testAst.epoch);
			if (Point.distance(posB, newPos) < r){
				// System.out.println("Good push!");
				Asteroid.push(a, minDistTime - 365, E, dir);
				didPush = true;
				break;
			} 
			
		}
		if(didPush == false){
			System.out.println("Could not find good push!");			
		}
	}

	public static Asteroid pushTest(Asteroid asteroid, long time,
            double energy, double direction)
	{
		if (Double.isNaN(energy) || Double.isInfinite(energy)
		             || energy < 0.0)
		throw new IllegalArgumentException("Invalid energy");
		if (Double.isNaN(direction) || Double.isInfinite(direction))
		throw new IllegalArgumentException("Invalid direction");
		// find current position and velocity of asteroid
		long t = time - asteroid.epoch;
		Point r = asteroid.orbit.positionAt(t);
		Point v = asteroid.orbit.velocityAt(t);
		// translate push energy to velocity and combine
		double magnitude = Math.sqrt(2.0 * energy / asteroid.mass);
		v.x += magnitude * Math.cos(direction);
		v.y += magnitude * Math.sin(direction);
		// return (new object of) the "same" asteroid with new orbit
		return new Asteroid(new Orbit(r, v), asteroid.mass, time);
	}

}