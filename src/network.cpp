#include "network.hpp"

std::random_device rd;
std::mt19937 g(rd());


Network3D::Network3D(int n, int k, double eta)
    : AVS3D{n, eta, 1.0}, k_(k) {
        for (int i = 0; i < n; i++){
            indices_.push_back(i);
        }
    }


void Network3D::random_align(){
    Coord3D new_velocities{3, n_};
    new_velocities.setZero();

    for (int i1 = 0; i1 < n_; i1++){
        shuffle(indices_.begin(), indices_.end(), g);

        bool self_included = false;

        for (int j = 0; j < k_ - 1; j++){
            int i2 = indices_[j];
            if (i2 == i1){ self_included = true; }
            new_velocities.col(i1) += velocities_.col(i2);
        }

        if (self_included){
            new_velocities.col(i1) += velocities_.col(indices_[k_ - 1]);
        } else {
            new_velocities.col(i1) += velocities_.col(i1);
        }
    }

    velocities_ = new_velocities / k_;
}


void Network3D::move(){
    random_align();
    normalise(velocities_);
    add_noise();
}

void Network3D::evolve(int steps){
    for (int i = 0; i < steps; i++){
        this->move();
    }
}
