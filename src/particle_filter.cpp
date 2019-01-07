/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 300;  // TODO: Set the number of particles

 	// Set standard deviations for x, y, and theta. (for better readability)
   double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

   // prepare random generator
   std::default_random_engine gen;

	// Create normal (gaussian) distributions for x, y and theta.
	std::normal_distribution<double> dist_x(x, std_x);
	std::normal_distribution<double> dist_y(y, std_y);
	std::normal_distribution<double> dist_theta(theta, std_theta);

   for(int i=0; i<num_particles; i++) {
      Particle p;
      p.id = i;
      p.x = dist_x(gen);
      p.y = dist_y(gen);
      p.theta = dist_theta(gen);
      p.weight = 1.0;

      particles.push_back(p);
      weights.push_back(p.weight);
   }

   is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

 	// Set standard deviations for x, y, and theta. (for better readability)
   double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];

   // prepare random generator
   std::default_random_engine gen;

   for(unsigned int i=0; i<particles.size(); i++) {
      Particle p = particles[i];

      double new_x, new_y, new_theta;

      if(yaw_rate == 0) {
         new_x = p.x + velocity * delta_t * cos(p.theta);
         new_y = p.y + velocity * delta_t * sin(p.theta);
         new_theta = p.theta;
      } else {
         new_x = p.x + velocity / yaw_rate * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta));
         new_y = p.y + velocity / yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate * delta_t));
         new_theta = p.theta + yaw_rate * delta_t;
      }

      // Create normal (gaussian) distributions for x, y and theta.
      std::normal_distribution<double> dist_x(new_x, std_x);
      std::normal_distribution<double> dist_y(new_y, std_y);
      std::normal_distribution<double> dist_theta(new_theta, std_theta);

      p.x = dist_x(gen);
      p.y = dist_y(gen);
      p.theta = dist_theta(gen);
      
      particles[i] = p;
   }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

   // Vector of weights of all particles
   std::vector<double> new_weights;

   // Set standard deviations for x and y (for better readability)
   double std_x = std_landmark[0];
   double std_y = std_landmark[1];

   //  go through all particles
   for(unsigned int index=0; index<particles.size(); index++) {
      // pick according particle data
      Particle p = particles[index];

      // prepare memory for associations
      std::vector<int> associations;
      std::vector<double> sense_x, sense_y;

      // transform observation coordinates from car-coord-system to map-system
      for(unsigned int index_obs=0; index_obs<observations.size(); index_obs++) {
         // get data for current observation
         LandmarkObs obs = observations[index_obs];

         // do the actual transformation
         double map_x = p.x + (cos(p.theta) * obs.x) - (sin(p.theta) * obs.y);
         double map_y = p.y + (sin(p.theta) * obs.x) + (cos(p.theta) * obs.y);

         // find closest landmarks and associate

         // if measurment is not within sensor range, skip...
         if(dist(p.x, p.y, map_x, map_y) > sensor_range)
            continue;

         double closest_dist = sensor_range;
         double current_dist = sensor_range;
         Map::single_landmark_s closest_lm = {-1, 0, 0};
         for(unsigned int index_lm=0; index_lm<map_landmarks.landmark_list.size(); index_lm++) {
            Map::single_landmark_s lm = map_landmarks.landmark_list[index_lm];
            current_dist = dist(map_x, map_y, lm.x_f, lm.y_f);
            if(current_dist < closest_dist) {
               closest_lm = lm;
               closest_dist = current_dist;
            }
         }
         if(closest_lm.id_i > -1) {
            associations.push_back(closest_lm.id_i);
            sense_x.push_back(map_x);
            sense_y.push_back(map_y);
         }
      }
      // set associations
      ParticleFilter::SetAssociations(p, associations, sense_x, sense_y);
         
      // update particle weight
      double weight = 1.0;
      for(unsigned int i=0; i<p.associations.size(); i++) {
         double map_x = p.sense_x[i];
         double map_y = p.sense_y[i];
         int lm_id = associations[i];

         Map::single_landmark_s lm;
         // pick the correct landmark
         for(unsigned int lm_i=0; lm_i<map_landmarks.landmark_list.size(); lm_i++) {
            if(map_landmarks.landmark_list[lm_i].id_i == lm_id) {
               lm = map_landmarks.landmark_list[lm_i];
               break;
            }
         }

         //double lX = map_landmarks.landmark_list[lm_id].x_f;
         //double lY = map_landmarks.landmark_list[lm_id].y_f;
         double lX = lm.x_f;
         double lY = lm.y_f;

         weight *= (1.0 / (2.0 * M_PI * std_x * std_y)) * exp(-1.0 * ( ((map_x-lX)*(map_x-lX))/(2*std_x*std_x) + ((map_y-lY)*(map_y-lY))/(2*std_y*std_y) ));
      }
      p.weight = weight;
      new_weights.push_back(weight);

      particles[index] = p;
   }
   
   weights = new_weights;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

   std::default_random_engine gen;
	std::discrete_distribution<int> dist(weights.begin(), weights.end());

	std::vector<Particle> resampled_particles;
	for(unsigned int i = 0; i < particles.size(); i++){
      resampled_particles.push_back(particles[dist(gen)]);
      // fix particle ids
      resampled_particles[i].id = i;
	}
	particles = resampled_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
