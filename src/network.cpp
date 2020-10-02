#include "network.h"

random_device rd;
mt19937 g(rd());


Network3D::Network3D(int n, int k, double eta)
    : Vicsek3D{n, 1.0, eta, 1.0} {
        for (int i = 0; i < n; i++){
            indices_.push_back(i);
        }
    }


void Network3D::random_align(){
    Coord3D new_velocities{3, n_};
    new_velocities.setZero();

    for (int i1 = 0; i1 < n_; i1++){
        bool self_included{false};
        shuffle(indices_.begin(), indices_.end(), g);
        for (int j = 0; j < k_ - 1; j++){
            int i2 = indices_[j];
            if (i2 == i1){ self_included = true; }
            new_velocities.col(i1) += velocities_.col(i2);
        }
        if (self_included){
            new_velocities.col(i1) += velocities_.col(indices_[k_]);
        } else {
            new_velocities.col(i1) += velocities_.col(indices_[i1]);
        }
    }
    velocities_ = new_velocities / k_;
}


void Network3D::move(bool rebuild){
    random_align();
    add_noise();
}
