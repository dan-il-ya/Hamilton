#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#define PI 3.14159265358979323846

struct Solver{
public:
    float dt = 1e-4;
    std::vector<float> q;
    std::vector<float> p;
    std::vector<float> dHdq;
    std::vector<float> dHdp;
    uint16_t n;
    Solver(const std::vector<float>& q_, std::vector<float> p_, float dt);
    virtual float Hamiltonian(std::vector<float> Q,std::vector<float> P) = 0;
    virtual void updateDer() = 0;

    float Hamiltonian(){
        return Hamiltonian(q,p);
    }
    void update(){
        updateDer();
        for(int i = 0;i < n;i++){
            q[i] += dt * dHdp[i];
            p[i] -= dt * dHdq[i];
        }
    };

    void updateEK2(){
        std::vector<float> p_(p);
        std::vector<float> q_(q);
        std::vector<float> dHdp_(n,0);
        std::vector<float> dHdq_(n,0);
        updateDer();
        for(int i = 0;i < n;i++){
            dHdp_[i] = dHdp[i];
            dHdq_[i] = dHdq[i];
            q[i] += dt * dHdp[i];
            p[i] += - dt * dHdq[i];
        }
        updateDer();
        for(int i = 0;i < n;i++){
            q_[i] += dt * (dHdp[i]+dHdp_[i])/2;
            p_[i] += - dt * (dHdq[i]+dHdq_[i])/2;
        }
        q = std::move(q_);
        p = std::move(p_);
    };

    void update(float  t){
        for(int i = 0;i < (int)(t/dt);i++) updateEK2();
    }
};

Solver::Solver(const std::vector<float>& q_, std::vector<float> p_, float dt = 1e-4) : q(q_), p(std::move(p_)), n(q_.size()), dt(dt) {}



struct Writer{
    std::ofstream file;
    int num = 0;
    Solver *S;
    Writer(Solver &s) : S(&s){};

    void snapshot(){
        file.open ("snp"+ std::to_string(num) + ".csv");
        file << S->q[0] << ',' << S->p[0];
        file.close();
        num++;
    }

};



struct Oscillator: public Solver{
    Oscillator(const std::vector<float>& q_, std::vector<float> p_);
    float Hamiltonian(std::vector<float> Q,std::vector<float> P) override {
        return 39.4784176*Q[0]*Q[0]/2 + P[0]*P[0]/2 +
               39.4784176*Q[1]*Q[1]*2 + P[1]*P[1]/2;
    };

    void updateDer(){
        dHdp.clear();
        dHdq.clear();
        /*
        dHdp.push_back(p[0]);
        dHdp.push_back(p[1]);
        dHdq.push_back(39.4784176*q[0]);
        dHdq.push_back(9*39.4784176*q[1]);
         */
        float den = (16-9* cos(q[0]-q[1])*cos(q[0]-q[1]));
        float num1 = (2*p[1]-3*p[0]*cos(q[0]-q[1]));

        dHdp.push_back(6*num1/den);
        dHdp.push_back(6*(6*p[1]+num1)/den);
        float prod = dHdp[0]*dHdp[1]*sin(q[0]-q[1]);

        dHdq.push_back((prod+sin(q[0]))/2);
        dHdq.push_back((-prod+sin(q[1]))/2);
    };
};

Oscillator::Oscillator(const std::vector<float>& q_, std::vector<float> p_) : Solver(std::move(q_),std::move(p_)) {}




int main() {
    std::vector<float> q({1, 1});
    std::vector<float> p({0, 0});

    Oscillator O(q,p);
/*
    for(int i = 0;i < 10;i++){
        std::cout << O.q[0] << ' ' << O.q[1] << ' ' << O.Solver::Hamiltonian() << '\n';
        O.update(0.1);
    }
*/
    Writer F(O);

    F.snapshot();

    int N = 50;
    float x_min = -PI,x_max = PI,y_min = -PI,y_max = PI;

    std::vector<std::vector<Oscillator*>> OO(N+1);
    for(int i = 0;i <= N;i++){
        OO[i].resize(N+1);
        for(int j = 0;j <= N;j++){
            float x = x_min + (float)i/(float)(N-1)*(x_max - x_min);
            float y = y_min + (float)j/(float)(N-1)*(y_max - y_min);
            OO[i][j] = new Oscillator(std::vector<float>({x,y}),std::vector<float>({0,0}));
        }
    }


    float  dt = 0.01;
    std::ofstream file;

    for(int i = 0;i < 10;i++){
        file.open ("snp"+ std::to_string(i) + ".csv");
        for(auto oo: OO){
            for(auto o:oo){
                o->update(dt);
                file << std::sin(o->q[0]) * std::cos(o->q[1]) << ", ";
            }
            file << '\n';
        }

        file.close();
    }



    return 0;
}
