#include <iostream>
#include <armadillo>
#include <vector>

#include "Particle.hpp"

// Constructor 
Particle::Particle(double charge_in, double mass_in, arma::vec position_in, arma::vec velocity_in)
{
	charge = charge_in;
	mass = mass_in;
	position = position_in;
	velocity = velocity_in;

}

// Function that returns the charge 
double Particle::return_charge()
{
	return charge;
}


// Function that returns the mass 
double Particle::return_mass()
{
	return mass;
}


// Function that returns the position 
arma::vec Particle::return_position()
{
	return position;
}


// Function that returns the velocity 
arma::vec Particle::return_velocity()
{
	return velocity;
}

//setter for position
void Particle::set_position(const arma::vec& new_position)
{
    position = new_position;
}

// Setter for velocity
void Particle::set_velocity(const arma::vec& new_velocity)
{
    velocity = new_velocity;
}