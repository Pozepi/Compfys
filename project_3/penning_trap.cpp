#include "particle.hpp"
#include "penning_trap.hpp"


PenningTrap::PenningTrap(double magnetic_field, double potential, double dimension)
{
    magnetic_field_ = magnetic_field;
    potential_ = potential;
    dimension_ = dimension;
}

void PenningTrap::add_particle(Particle particle_in)
{
    particles.push_back(particle_in);
}

void PenningTrap::add_n_random_particles(int n, double charge, double mass)
{
    /*
    Adds n particles to the box with random position and velocity following a distribution
    Args:
        n       (int)       : amount of particles to add
        charge  (double)    : the charge of all the n particles
        mass    (double)    : mass of all the n particles
    */
    for(int i=0; i<n; i++)
    {
        arma::arma_rng::set_seed(i);
        arma::vec r = arma::vec(3).randn()*0.1*dimension_;
        arma::vec v = arma::vec(3).randn()*0.1*dimension_;
        Particle new_particle(charge, mass, r, v);
        add_particle(new_particle);
    }
}

double PenningTrap::particle_count()
{
    /*
    Find the current number of particles in the penning trap
    Returns: 
        particle.size()     (double)    : current number of particles in the penning trap
    */
    return particles.size();
}

bool PenningTrap::particle_outside_trap_check(arma::vec r)
{
    /*
    Checks if the particle is outside the trap
    Args:
        r               (arma::vec) : the position of the particle
    Returns: 
        outside_trap    (bool)      : true if particle is outside the trap
    */
    bool outside_trap = false;
    if(std::abs(r(0)) > dimension_ || std::abs(r(1)) > dimension_ || std::abs(r(2)) > dimension_)
    {
        outside_trap = true;
    }

    return outside_trap;
}

arma::vec PenningTrap::external_E_field(arma::vec r, double time, double f, double wv)
{
    /*
    Calculates the external electric field. If particle is outside of PenningTrap dimensions, E is set to zero. 
    Args:
        r       (arma::vec) : position of the particle
        time    (double)    : the current time of the simulation
        f       (double)    : amplitue of the time-dependent part of the time-dependent potential
        wv      (double)    : frequency of the time-dependent potential
    Returns: 
        E       (arma::vec) : the electric field
    */
    double pot = potential_+potential_*f*std::cos(wv*time);

    double x_der = 2*r[0]*pot/(2*dimension_*dimension_);
    double y_der = 2*r[1]*pot/(2*dimension_*dimension_);
    double z_der = -4*r[2]*pot/(2*dimension_*dimension_);
    
    arma::vec E = {x_der,y_der,z_der};

    if(PenningTrap::particle_outside_trap_check(r) == true)
        E = {0,0,0};

    return E;
}

arma::vec PenningTrap::external_B_field(arma::vec r)
{
    /*
    Finds the external magnetic field. If the particle is outside the box gives {0, 0, 0}
    Args:
        r   (arma::vec) : the position of the particle
    Returns: 
        B   (arma::vec) : the external magnetic field
    */
    arma::vec B = {0,0,magnetic_field_};
    if(PenningTrap::particle_outside_trap_check(r) == true)
        B = {0,0,0};

    return B;
}
// Force on particle i from particle j
arma::vec PenningTrap::force_particle(int i, int j)
{
    /*
    Calculates the force acting on particle i from particle j
    Args:
        i   (int)   : our particle index
        j   (int)   : index of particle acting on our particle
    Returns: 
        F   (arma::vec) : the force acting on particle i from particle j
    */
    double ke = 1.38935333e+05;
    Particle particle_i = particles[i];
    Particle particle_j = particles[j];

    double rxi = particle_i.position()(0);
    double ryi = particle_i.position()(1);
    double rzi = particle_i.position()(2);

    double rxj = particle_j.position()(0);
    double ryj = particle_j.position()(1);
    double rzj = particle_j.position()(2);
    
    double A = ke*particle_i.charge()*particle_j.charge();
    arma::vec R = {rxi-rxj,ryi-ryj,rzi-rzj};

    double value = std::sqrt(R(0)*R(0) + R(1)*R(1) + R(2)*R(2));
    double value3 = value*value*value;
    
    double Fx = A*R(0)/value3;
    double Fy = A*R(1)/value3;
    double Fz = A*R(2)/value3;
    
    arma::vec F = {Fx,Fy,Fz};

    return F;
}
// Total force on particle i from external fields
arma::vec PenningTrap::total_force_external(int i, double time, double f, double wv)
{  
    /*
    Sum of the force due to the external magnetic and electric fields
    Args:
        i       (int)       : particle index i to calculate force of
        time    (double)    : time used to calculate time-dependent potential 
        f       (double)    : amplitude of the time-dependent part of the potential
        wv      (double)    : frequency of the time-dependent potential
    Returns:
        F       (arma::vec) : total force due to the magnetic and electric fields
    */
    Particle particle_i = particles[i];
    arma::vec v = particle_i.velocity();
    arma::vec B = PenningTrap::external_B_field(particle_i.position());
    arma::vec E = PenningTrap::external_E_field(particle_i.position(), time, f, wv);

    double Fx = particle_i.charge()*E(0)+particle_i.charge()*(v(1)*B(2)-v(2)*B(1));
    double Fy = particle_i.charge()*E(1)-particle_i.charge()*(v(0)*B(2)-v(2)*B(0));
    double Fz = particle_i.charge()*E(2)+particle_i.charge()*(v(0)*B(1)-v(1)*B(0));
    arma::vec F = {Fx, Fy, Fz};
    return F;
}

arma::vec PenningTrap::total_force_particles(int i)
{
    /*
    Calculates the force acting on particle i due to all other particles
    Args:
        i   (int)       : particle index to calcualte forces of
    Returns:
        F   (arma::vec) : Force acting on particle i 
    */
    double Fx = 0;
    double Fy = 0;
    double Fz = 0;
    for (int j=0;j<PenningTrap::particle_count();j++)
    {
        if (i!=j)
        {
            Fx += PenningTrap::force_particle(i,j)(0);
            Fy += PenningTrap::force_particle(i,j)(1);
            Fz += PenningTrap::force_particle(i,j)(2);
        }
    }
    arma::vec F = {Fx, Fy, Fz};
    return F;
}
arma::vec PenningTrap::total_force(int i, bool interaction, double time, double f, double wv)
{
    /*
    Total force on particle i from both external fields and other particles
    Args:
        i           (int)       : particle index
        interaction (bool)      : true if particle interaction is enabled
        time        (double)    : current time of the simulation
        f           (double)    : amplitude of the time-dependent part of the potential
        wv          (double)    : frequency of the time-dependent potential
    Returns: 
        F           (arma::vec) : the total force acting on particle i
    */
    arma::vec F;
    if(interaction)
    {
        double Fx = PenningTrap::total_force_external(i, time, f, wv)(0)+PenningTrap::total_force_particles(i)(0);
        double Fy = PenningTrap::total_force_external(i, time, f, wv)(1)+PenningTrap::total_force_particles(i)(1);
        double Fz = PenningTrap::total_force_external(i, time, f, wv)(2)+PenningTrap::total_force_particles(i)(2);
        F = {Fx,Fy,Fz};
    }
    if(!interaction)
    {  
        F = PenningTrap::total_force_external(i, time, f, wv);
    }
    return F;
}

void PenningTrap::evolve_RK4(double dt, double time_stop, bool interaction, double f, double wv, bool makefile, std::string filename)
{
    /*
    Integrates the position and velocity of the particle(s) using the Runge Kutta 4 method
    Args:
        dt          (double)        : timestep
        time_stop   (double)        : tells the system when to stop simulating (microseconds)
        interation  (bool)          : false turns off particle interactions
        f           (double)        : amplitude of the time dependent part of the time dependent potential
        wv          (double)        : frequency of the time dependent potential
        makefile    (bool)          : if true writes results to a file with filename given by paramtere filename
        filename    (std::string)   : name of the file to write results to
    */
    int size = PenningTrap::particle_count()*3;
    int N = time_stop/dt;
    arma::mat v(N, size);
    arma::mat pos(N, size);
    arma::vec times(N);
    
    int jump = 0;
    // setting initial values
    times(0) = 0;
    for(int k=0; k<PenningTrap::particle_count(); k++)
    {
        Particle particle = particles[k];
        v(0,jump) = particle.velocity()(0);
        v(0,jump+1) = particle.velocity()(1);
        v(0,jump+2) = particle.velocity()(2);
        pos(0,jump) = particle.position()(0);
        pos(0,jump+1) = particle.position()(1);
        pos(0,jump+2) = particle.position()(2);
        jump += 3;
    }
    
    for(int i=0; i<N-1; i++)
    {   
        // solving for the position and velocity of the particles
        //
        // h = t[i+1] - t[i] -> dt
        // k1 = f(y[i], t[i], *args)
        // k2 = f(y[i] + k1 * h / 2., t[i] + h / 2., *args)
        // k3 = f(y[i] + k2 * h / 2., t[i] + h / 2., *args)
        // k4 = f(y[i] + k3 * h, t[i] + h, *args)
        // y[i+1] = y[i] + (h / 6.) * (k1 + 2*k2 + 2*k3 + k4)
        arma::mat k1v(3, PenningTrap::particle_count());
        arma::mat k1p(3, PenningTrap::particle_count());

        arma::mat k2v(3, PenningTrap::particle_count());
        arma::mat k2p(3, PenningTrap::particle_count());

        arma::mat k3v(3, PenningTrap::particle_count());
        arma::mat k3p(3, PenningTrap::particle_count());

        arma::mat k4v(3, PenningTrap::particle_count());
        arma::mat k4p(3, PenningTrap::particle_count());


        jump=0;
        for(int j=0; j<PenningTrap::particle_count(); j++)
        { 
            double charge = particles[j].charge();
            double mass = particles[j].mass();
            arma::vec tmp_pos = particles[j].position();
            arma::vec tmp_vel = particles[j].velocity();
            arma::vec F = PenningTrap::total_force(j, interaction, times(i), f, wv);

            k1v(0,j) = F(0)/mass;
            k1v(1,j) = F(1)/mass;
            k1v(2,j) = F(2)/mass;
            
            k1p(0,j) = tmp_vel(0);
            k1p(1,j) = tmp_vel(1);
            k1p(2,j) = tmp_vel(2);
            
            tmp_pos = {pos(i,jump) + k1p(0,j)*dt/2, pos(i,jump+1) + k1p(1,j)*dt/2, pos(i,jump+2) + k1p(2,j)*dt/2};
            tmp_vel = {v(i,jump) + k1v(0,j)*dt/2, v(i,jump+1) + k1v(1,j)*dt/2, v(i,jump+2) + k1v(2,j)*dt/2};
            particles[j] = {charge, mass, tmp_pos, tmp_vel};
            jump += 3;
        }
        
        jump = 0;
        for(int j=0; j<PenningTrap::particle_count(); j++)
        {
            double charge = particles[j].charge();
            double mass = particles[j].mass();
            arma::vec tmp_pos = particles[j].position();
            arma::vec tmp_vel = particles[j].velocity();
            arma::vec F = PenningTrap::total_force(j, interaction, times(i)+dt/2, f, wv);

            k2v(0,j) = F(0)/mass;
            k2v(1,j) = F(1)/mass;
            k2v(2,j) = F(2)/mass;
            
            k2p(0,j) = tmp_vel(0);
            k2p(1,j) = tmp_vel(1);
            k2p(2,j) = tmp_vel(2);
            
            tmp_pos = {pos(i,jump) + k2p(0,j)*dt/2, pos(i,jump+1) + k2p(1,j)*dt/2, pos(i,jump+2) + k2p(2,j)*dt/2};
            tmp_vel = {v(i,jump) + k2v(0,j)*dt/2, v(i,jump+1) + k2v(1,j)*dt/2, v(i,jump+2) + k2v(2,j)*dt/2};
            particles[j] = {charge, mass, tmp_pos, tmp_vel};  
            jump+=3;
        }

        jump=0;
        for(int j=0; j<PenningTrap::particle_count(); j++)
        {  
            double charge = particles[j].charge();
            double mass = particles[j].mass();
            arma::vec tmp_pos = particles[j].position();
            arma::vec tmp_vel = particles[j].velocity();
            arma::vec F = PenningTrap::total_force(j, interaction, times(i)+dt/2, f, wv);

            k3v(0,j) = F(0)/mass;
            k3v(1,j) = F(1)/mass;
            k3v(2,j) = F(2)/mass;
            
            k3p(0,j) = tmp_vel(0);
            k3p(1,j) = tmp_vel(1);
            k3p(2,j) = tmp_vel(2);
            
            tmp_pos = {pos(i,jump) + k3p(0,j)*dt/2, pos(i,jump+1) + k3p(1,j)*dt/2, pos(i,jump+2) + k3p(2,j)*dt/2};
            tmp_vel = {v(i,jump) + k3v(0,j)*dt/2, v(i,jump+1) + k3v(1,j)*dt/2, v(i,jump+2) + k3v(2,j)*dt/2};
            particles[j] = {charge, mass, tmp_pos, tmp_vel};
            jump+=3;  
        }

        for(int j=0; j<PenningTrap::particle_count(); j++)
        {
            double charge = particles[j].charge();
            double mass = particles[j].mass();
            arma::vec tmp_pos = particles[j].position();
            arma::vec tmp_vel = particles[j].velocity();
            arma::vec F = PenningTrap::total_force(j, interaction, times(i)+dt, f, wv);

            k4v(0,j) = F(0)/mass;
            k4v(1,j) = F(1)/mass;
            k4v(2,j) = F(2)/mass;
            
            k4p(0,j) = tmp_vel(0);
            k4p(1,j) = tmp_vel(1);
            k4p(2,j) = tmp_vel(2);
        }
        
        // k4 = f(y[i] + k3 * h, t[i] + h, *args)
        // same as above but with k3 and with dt not dt/2

        // FINALLY update velocity and position using RK4
        jump = 0;
        for(int j=0; j<PenningTrap::particle_count(); j++)
        {
            v(i+1,jump) = v(i,jump) + (dt / 6) * (k1v(0,j) + 2*k2v(0,j) + 2*k3v(0,j) + k4v(0,j));
            v(i+1,jump+1) = v(i,jump+1) + (dt / 6) * (k1v(1,j) + 2*k2v(1,j) + 2*k3v(1,j) + k4v(1,j));
            v(i+1,jump+2) = v(i,jump+2) + (dt / 6) * (k1v(2,j) + 2*k2v(2,j) + 2*k3v(2,j) + k4v(2,j));
            
            pos(i+1,jump) = pos(i,jump) + (dt / 6) * (k1p(0,j) + 2*k2p(0,j) + 2*k3p(0,j) + k4p(0,j));
            pos(i+1,jump+1) = pos(i,jump+1) + (dt / 6) * (k1p(1,j) + 2*k2p(1,j) + 2*k3p(1,j) + k4p(1,j));
            pos(i+1,jump+2) = pos(i,jump+2) + (dt / 6) * (k1p(2,j) + 2*k2p(2,j) + 2*k3p(2,j) + k4p(2,j));
            
            jump += 3;
        }
        jump = 0;
        // updating particles position and velocity for next iteration
        for(int k=0; k<PenningTrap::particle_count(); k++)
        {
            Particle particle = particles[0];
            double mass = particle.mass();
            double charge = particle.charge();

            arma::vec update_vel = {v(i+1,jump), v(i+1,jump+1), v(i+1,jump+2)};
            arma::vec update_pos = {pos(i+1,jump), pos(i+1,jump+1), pos(i+1,jump+2)};
            Particle update_particle(charge, mass, update_pos, update_vel);
            particles.erase(particles.begin());
            PenningTrap::add_particle(update_particle);

            jump += 3;
        }
        times(i+1) = times(i)+dt;
    }
    if(makefile)
    {
        mkdir("output_files",0777);
        mkdir(("output_files//"+filename).c_str(),0777);
        times.save("output_files//"+filename+"//"+filename+"_t.txt");
        v.save("output_files//"+filename+"//"+filename+"_v.txt");
        pos.save("output_files//"+filename+"//"+filename+"_pos.txt");
    }
}

void PenningTrap::evolve_forward_Euler(double dt, double time_stop, bool interaction, double f, double wv, bool makefile, std::string filename)
{
    /*
    Integrates the position and velocity of the particle(s) using the Euler-Cromer method
    Args:
        dt          (double)        : timestep
        time_stop   (double)        : tells the system when to stop simulating (microseconds)
        interation  (bool)          : false turns off particle interactions
        f           (double)        : amplitude of the time dependent part of the time dependent potential
        wv          (double)        : frequency of the time dependent potential
        makefile    (bool)          : if true writes results to a file with filename given by paramtere filename
        filename    (std::string)   : name of the file to write results to
    */
    int size = PenningTrap::particle_count()*3;
    int N = time_stop/dt;
    arma::mat v(N, size);
    arma::mat pos(N, size);
    arma::vec times(N);;
    
    int jump = 0;
    // setting initial values
    times(0) = 0;
    for(int k=0; k<PenningTrap::particle_count(); k++)
    {
        Particle particle = particles[k];
        v(0,jump) = particle.velocity()(0);
        v(0,jump+1) = particle.velocity()(1);
        v(0,jump+2) = particle.velocity()(2);
        pos(0,jump) = particle.position()(0);
        pos(0,jump+1) = particle.position()(1);
        pos(0,jump+2) = particle.position()(2);
        jump += 3;
    }

    for(int i=0; i<N-1; i++)
    {   
        // solving for the position and velocity of the particles
        jump = 0;
        for(int j=0; j<PenningTrap::particle_count(); j++)
        {
            Particle particle_j = particles[j];

            arma::vec F = PenningTrap::total_force(j, interaction, times(i), f, wv);
            arma::vec a = {F(0)/particle_j.mass(), F(1)/particle_j.mass(), F(2)/particle_j.mass()};
            v(i+1,jump) = v(i,jump) + a(0)*dt;
            v(i+1,jump+1) = v(i,jump+1) + a(1)*dt;
            v(i+1,jump+2) = v(i,jump+2) + a(2)*dt;
            pos(i+1,jump) = pos(i,jump)+v(i+1,jump)*dt;
            pos(i+1,jump+1) = pos(i,jump+1)+v(i+1,jump+1)*dt;
            pos(i+1,jump+2) = pos(i,jump+2)+v(i+1,jump+2)*dt;
            
            jump += 3;
        }
        // updating particles position and velocity for next iteration
        jump = 0;
        for(int k=0; k<PenningTrap::particle_count(); k++)
        {
            Particle particle = particles[0];
            double mass = particle.mass();
            double charge = particle.charge();

            arma::vec update_vel = {v(i+1,jump), v(i+1,jump+1), v(i+1,jump+2)};
            arma::vec update_pos = {pos(i+1,jump), pos(i+1,jump+1), pos(i+1,jump+2)};
            Particle update_particle(charge, mass, update_pos, update_vel);
            particles.erase(particles.begin());
            PenningTrap::add_particle(update_particle);

            jump += 3;
        }
        times(i+1) = times(i)+dt;
    }
    if(makefile)
    {
        mkdir("output_files",0777);
        mkdir(("output_files//"+filename).c_str(),0777);
        times.save("output_files//"+filename+"//"+filename+"_t.txt");
        v.save("output_files//"+filename+"//"+filename+"_v.txt");
        pos.save("output_files//"+filename+"//"+filename+"_pos.txt");
        /*
        std::fstream file;
        file.open(filename+".txt", std::ios::out);
        file << times << '\n';
        file << v << '\n';
        file << pos << '\n';
        */
    }
}

double PenningTrap::particles_inside_trap_count()
{
    /*
    Counts the number of particles within the trap 
    Returns:
        count   (int)    :  number of particles within the trap
    */
    double N = PenningTrap::particle_count();
    int count = 0;
    for(int i=0; i<N; i++)
    {
        Particle particle_i = particles[i];
        arma::vec pos = particle_i.position();
        if(PenningTrap::particle_outside_trap_check(pos) == false)
            count += 1;
    }
    
    return count;
}

