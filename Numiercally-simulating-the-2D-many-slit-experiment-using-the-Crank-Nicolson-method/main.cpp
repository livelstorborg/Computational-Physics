#include "PDEModel.hpp"


void save_matrix_to_csv(const arma::mat& matrix, const std::string& filename) 
{
    std::ofstream file(filename);
    for (size_t i = 0; i < matrix.n_rows; ++i) {
        for (size_t j = 0; j < matrix.n_cols; ++j) {
            file << matrix(i, j);
            if (j < matrix.n_cols - 1) {
                file << ","; 
            }
        }
        file << "\n";  
    }
    file.close();
}



void save_vector_to_csv(const arma::vec& vec, const std::string& filename) 
{
    std::ofstream file(filename);

    file << std::setprecision(16) << std::scientific;

    for (size_t i = 0; i < vec.n_elem; ++i) {
        file << vec[i];
        if (i < vec.n_elem - 1) {
            file << ","; 
        }
    }
    file << "\n"; 
    file.close();
}






int main()
{





    // ---------- Setting up simulation parameters ----------
    double dt = 2.5e-5;
    double T = 0.002;
    arma::vec t = arma::regspace(0, dt, T);
    int timesteps = t.n_elem;                   
    double dx = 0.005;
    double dy = 0.005;
    int M = 1 / dx;
    double x_c = 0.25;
    double y_c = 0.5;
    double p_x = 200.0;
    double p_y = 0.0;
    double sigma_x = 0.05;
    double sigma_y = 0.2;

    PDEModel model(dt, dx, dy, x_c, y_c, sigma_x, sigma_y, p_x, p_y); 




    // ---------- Setting up the potential matrix V ----------
    arma::mat V(M-2, M-2, arma::fill::zeros);               
    double v_0 = 1e10;
    int nr_slits = 2;
    double thickness = 0.02;
    double centre = 0.5;
    double middle_wall = 0.05;
    double opening = 0.05;
    model.construct_potential(V, v_0, M, nr_slits, thickness, centre, middle_wall, opening);

    //Fixing boundary conditions before plotting
    arma::mat V_full(M, M, arma::fill::zeros);             
    V_full.submat(1, 1, M-2, M-2) = V; 

    save_matrix_to_csv(V_full, "Potential_2.csv");        






    // ---------- Construction A and B matrices ----------
    auto result = model.construct_A_B(M, dx, dt, V);        
    arma::sp_cx_mat A = std::get<0>(result);
    arma::sp_cx_mat B = std::get<1>(result);

    arma::cx_vec u = model.initial_state(M);  
    u = model.normalised_initial_state(u);     






    // ---------- Looping over time steps and storing each u ----------

    arma::cx_mat U((M-2) * (M-2), timesteps);   //Will store all the time iterations of u
    U.col(0) = u;                               //Saving first u


    for (size_t i = 1; i < timesteps; ++i)                         //Finding the u of all time steps
    {
        arma::cx_vec u_1 = model.crank_nicolson(A, B, U.col(i-1)); //Finding next u
        U.col(i) = u_1;                                            //Saving the next u
        std::cout << i << std::endl;
    }


    arma::cx_mat U_conj = arma::conj(U);
    arma::mat P = arma::real(U_conj % U);       //This P matrix is a 2D one representing a 3D one

    arma::mat U_real = arma::real(U);
    arma::mat U_imag = arma::imag(U);



    // ---------- Saving data to produce plots and animations ----------
    save_matrix_to_csv(P, "P_2.csv");
    save_matrix_to_csv(U_real, "U_real_2.csv");
    save_matrix_to_csv(U_imag, "U_imag_2.csv");

    arma::vec p(timesteps);
    for(size_t i = 0; i < timesteps; ++i){
        p(i) = arma::accu(P.col(i));        //Summing up all the probabilites of each time step
    }
    save_vector_to_csv(p, "Prob_vec_2.csv");


    //Plotting code to see structure of A or B matrix
    // std::cout << "Shape of Matrices A and B" << std::endl;
    // model.print_sp_matrix_structure(A);
    // model.print_sp_matrix_structure(B);



	return 0;
}