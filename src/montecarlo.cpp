#include "montecarlo.hpp" 


FishMCMC::FishMCMC(size_t n, double beta, double si, double sh, double dx)
    : n_{n}, beta_{beta}, si_{si}, sh_{sh}, dx_{dx} {
        positions_ = tank_.get_random_positions(n);
}

void FishMCMC::sweep(){
    Coord3D rand_move{3, n_};
    rand_move.setRandom();  // ~ U(-1, 1)
    rand_move *= dx_;  // ~ U(-dx, dx)
    Property rand_acc{n_};
    rand_acc.setRandom();  // ~ U(-1, -1)
    rand_acc = (rand_acc + 1) * 0.5;  // ~U(0, 1)
    Vec3D p_new, p_old;
    double ndE;  // -(E_new - E_old)
    for (size_t i = 0; i < n_; i++){
        p_new = positions_.col(i) + rand_move.col(i);
        if (tank_.is_inside(p_new)){
            p_old = positions_.col(i);
            ndE = this->get_dE(p_old, p_new, i);
            if (ndE > 0) {
                positions_.col(i) = p_new;
            } else if (exp(beta_ * ndE) > rand_acc(i)){
                positions_.col(i) = p_new;
            }
        }
    }
}

void FishMCMC::evolve(size_t steps){
    for (size_t s = 0; s < steps; s++){
        this->sweep();
    }
}

double FishMCMC::get_dE(Vec3D p_new, Vec3D p_old, size_t i){
    double E_new = get_energy(p_new, i);
    double E_old = get_energy(p_old, i);
    return E_new - E_old;
}

double FishMCMC::get_energy(Vec3D pos, size_t i){
    // pairwise interaction
    double e_pair = 0, dij = 0;
    for (size_t j = 0; j < n_; j++){
        if (j != i){
            dij = (pos - positions_.col(j)).norm();
            e_pair += -log(h_) + log(2) * pow(
                log(2.0 * a_ * (dij - c_) / w_ + 1) / a_, 2
            );
            //tmp = -pow(log(1 + (2 * a_ * (dij - c_) / w_)), 2) / (a_ * a_);
            //e_pair -= log(pow(2, tmp) * h_);
        }
    }
    // gravity
    double z = pos(2);
    // holes
    double r = sqrt(pos(0)*pos(0) + pos(1)*pos(1));
    double e_holes = 0;
    for (auto rh : r_holes){
        e_holes += pow(
            pow(r - rh, 2) + pow(z - c_ * rh*rh, 2),
            -2
        ) / rh;
    }
    e_holes = sh_ * e_holes;
    return z + e_holes + e_pair * si_;
}


double FishMCMC::get_total_energy(){
    double energy = 0;
    for (size_t i = 0; i < n_; i++){
        energy += this->get_energy(positions_.col(i), i);
    }
    return energy;
}

