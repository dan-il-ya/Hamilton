#include <iostream>
#include <vector>

struct Solver{
public:
    float dt = 1e-4;
    std::vector<float> q;
    std::vector<float> p;
    uint16_t n;
    Solver(const std::vector<float>& q_, std::vector<float> p_);
    virtual float Hamiltonian(std::vector<float> Q,std::vector<float> P) = 0;
    float Hamiltonian(){
        return Hamiltonian(q,p);
    }
    void update(){
        std::vector<float> p_(p);
        std::vector<float> q_(q);
        std::vector<float> pdp(p);
        std::vector<float> qdq(q);
        for(int i = 0;i < n;i++){
            qdq[i] += dt;
            pdp[i] += dt;
            q_[i] += Hamiltonian(q,pdp) - Hamiltonian();
            p_[i] += Hamiltonian() - Hamiltonian(qdq,p);
            qdq[i] -= dt;
            pdp[i] -= dt;
        }
        q = std::move(q_);
        p = std::move(p_);
    };

    void updateRK4(){
        std::vector<float> p_(p);
        std::vector<float> q_(q);
        std::vector<float> pdp(p);
        std::vector<float> qdq(q);
        for(int i = 0;i < n;i++){
            qdq[i] += dt;
            pdp[i] += dt;
            q_[i] += (Hamiltonian(q,pdp) - Hamiltonian())/2;
            p_[i] += (Hamiltonian() - Hamiltonian(qdq,p))/2;
            qdq[i] -= dt;
            pdp[i] -= dt;
        }
        pdp = p_;
        qdq = q_;
        q = q_;
        p = p_;
        for(int i = 0;i < n;i++){
            qdq[i] += dt;
            pdp[i] += dt;
            q_[i] += (Hamiltonian(q,pdp) - Hamiltonian())/2;
            p_[i] += (Hamiltonian() - Hamiltonian(qdq,p))/2;
            qdq[i] -= dt;
            pdp[i] -= dt;
        }
        q = std::move(q_);
        p = std::move(p_);
    };
    void update(float  t){
        for(int i = 0;i < (int)(t/dt);i++) update();
    }
};

Solver::Solver(const std::vector<float>& q_, std::vector<float> p_) : q(q_), p(std::move(p_)), n(q_.size()) {}

struct Oscillator: public Solver{
    Oscillator(std::vector<float> q_, std::vector<float> p_);
    float Hamiltonian(std::vector<float> Q,std::vector<float> P) override {
        return 39.4784176*Q[0]*Q[0]/2 + P[0]*P[0]/2 +
                39.4784176*Q[1]*Q[1]*2 + P[1]*P[1]/2;
    };
};

Oscillator::Oscillator(std::vector<float> q_, std::vector<float> p_) : Solver(std::move(q_),std::move(p_)) {}

int main() {
    std::vector<float> q({1, 1});
    std::vector<float> p({0, 0});

    Oscillator O(q,p);

    for(int i = 0;i <= 10;i++){
        std::cout << O.q[0] << ' ' << O.q[1] << '\n';
        O.update(0.1);
    }

    return 0;
}
