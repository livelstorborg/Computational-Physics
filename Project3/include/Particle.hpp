#ifndef __Particle_hpp__
#define __Particle_hpp__

class Particle
{

	public:

		double charge;
		double mass;
		arma::vec position;
		arma::vec velocity;

		//constructor that assigns values to the member variables
		Particle(double charge_in, double mass_in, arma::vec position_in, arma::vec velocity_in);

		// Function that returns the charge 
		double return_charge();


		// Function that returns the mass 
		double return_mass(); 


		// Function that returns the position 
		arma::vec return_position();


		// Function that returns the velocity 
		arma::vec return_velocity();

		//Function that sets a new position to Particle
		void set_position(const arma::vec& new_position);

		//Funciton that sets a new velocity to Particle.
    	void set_velocity(const arma::vec& new_velocity);
};

#endif