#include <iostream>
#include <pwd.h>
#include <cmath>
#include <random>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>


namespace pt = boost::property_tree;

class Taub_Miles {
private:
    double _v;
    double _t;
    double _n;
    double _m;
    double _h;
    double _z;
    double _eta;
    double _xi;
    double _xi_eta;

    double dt = 0.0001;

    float D_v;
    float D_eta;
    float tau_eta = 10;
    float C = 1.;
    float g_na = 100;
    float g_k = 80;
    float g_l = 0.1;
    float g_z = 5.;
    float E_na = 50;
    float E_l = -67;
    float E_k = -100;

    float _mu;

    std::random_device rd;
    std::mt19937 _generator;
    std::normal_distribution<double> dist;

public:
    Taub_Miles(double v0, double n0, double m0, double h0, double z0, float I, float Dv, float Dn) : _generator(rd()),
                                                                                 dist(std::normal_distribution<double>(
                                                                                         0, 1)) {
        _v = v0;
        _n = n0;
        _m = m0;
        _h = h0;
        _z = z0;
        _mu = I;
        D_v = Dv;
        D_eta = Dn;
    }

    void set_current(float I) {
        _mu = I;
    }

    double get_time() const {
        return _t;
    }

    double get_voltage() const {
        return _v;
    }

    double get_adaptation() const {
        return _z;
    }

    double get_cnoise() const{
        return _eta;
    }

    double get_I_adaptation() const {
        return I_adap();
    }

    double get_n() const {
        return _n;
    }

    double get_m() const {
        return _m;
    }

    double get_h() const {
        return _h;
    }

    void integrate_with_adap() {
        _xi = dist(_generator);
        update_h();
        update_m();
        update_n();
        update_z();
        _v += (_mu - I_ion() - I_adap() + I_noise())/C * dt + sqrt(2 * D_v * dt)/C * _xi;
        _t += dt;
    }

    void integrate_without_adap() {
        _xi = dist(_generator);
        update_h();
        update_m();
        update_n();
        update_z();
        _v += (_mu - I_ion() + I_noise())/C * dt + sqrt(2 * D_v * dt)/C * _xi;
        _t += dt;
    }

    void kick_voltage(double e) {
        _v += e;
    }

    double I_ion() const {
        return g_na * _h * pow(_m, 3) * (_v - E_na) + g_l * (_v - E_l) + g_k * pow(_n, 4) * (_v - E_k);
    }

    void update_h() {
        _h += (a_h() * (1. - _h) - b_h() * _h) * dt;
    }

    double a_h() const {
        return 0.128 * exp(-(50. + _v) / 18.);
    }

    double b_h() const {
        return 4. / (1. + exp(-(_v + 27.) / 5.));
    }

    void update_m() {
        _m += (a_m() * (1. - _m) - b_m() * _m) * dt;
    }

    double a_m() const {
        return 0.32 * (_v + 54) / (1. - exp(-(_v + 54.) / 4.));
    }

    double b_m() const {
        return 0.28 * (_v + 27.) / (exp((_v + 27.) / 5.) - 1.);
    }

    void update_n() {
        _n += (a_n() * (1. - _n) - b_n() * _n) * dt;
    }

    double a_n() const {
        return 0.032 * (_v + 52.) / (1. - exp(-(_v + 52.) / 5.));
    }

    double b_n() const {
        return 0.5 * exp(-(_v + 57.) / 40.);
    }

    double I_adap() const {
        return g_z * _z * (_v - E_k);
    }

    void update_z() {
        _z += 0.01 * (1. / (1. + exp(-(_v + 20.) / 5.)) - _z) * dt;
    }

    double I_noise() {
        _xi_eta = dist(_generator);
        _eta += (-_eta / tau_eta) * dt + sqrt(2 * D_eta * dt) * _xi_eta / tau_eta;
        return _eta;
    }
};


int main(int argc, char *argv[]) {
    float t_kick = atof(argv[1]);
    float e_kick = atof(argv[2]);
    bool kicked = false;

    float Dv = 1.0;
    float Dn = 0.0;
    double v0 = 0.;
    double n0 = 0.;
    double m0 = 0.;
    double h0 = 0.;
    double z0 = 0.;

    /* I = 10 with adap?
   double v0 = 44.8812;
   double n0 = 0.305847;
   double m0 = 0.93392;
   double h0 = 0.639683;
   double z0 = 0.0344068;
   */

    // I = 5 with adap
    /*
    double v0 = 45.0041;
    double n0 = 0.340271;
    double m0 = 0.947488;
    double h0 = 0.616439;
    double z0 = 0.0221737;
    */

    // I = 5 no adap
    /*double v0 = 45.0028;
    double n0 = 0.326704;
    double m0 = 0.94069;
    double h0 = 0.628454;
    double z0 = 0.0222072;
    */

    //float I = 3.6365; 10
    float I = 0; // 1.4167;

    struct passwd *pw = getpwuid(getuid());
    const char *homedir = pw->pw_dir;

    std::string path;
    //path = "../out/PRC_cnoise/";
    path = "../out/";
    char parameters[200];
    //std::sprintf(parameters, "_I%.1f_t%.2f_e%.2f", I, t_kick, e_kick);


    std::sprintf(parameters, "_I%.1f_Dv%.2f_Dn%.2f", I, Dv, Dn);


    std::string out_file;
    //out_file = path + "taubs_miles_cnoise_PRC" + parameters + ".dat";
    out_file = path + "taubs_miles_adap_cnoise" + parameters + ".dat";
    std::ofstream file;
    file.open(out_file);
    if (!file.is_open()) {
        std::cout << "Could not open file at: " << out_file << std::endl;
        std::cout << "This is where I am: " << std::string(homedir) << std::endl;
        return 1;
    }

    std::string out_file_spikes;
    //out_file_spikes = path + "spikes_taubs_miles_cnoise_PRC" + parameters + ".dat";
    out_file_spikes = path + "spikes_taubs_miles_adap_cnoise" + parameters + ".dat";
    std::ofstream file_spikes;
    file_spikes.open(out_file_spikes);
    if (!file_spikes.is_open()) {
        std::cout << "Could not open file at: " << out_file_spikes << std::endl;
        std::cout << "This is where I am: " << std::string(homedir) << std::endl;
        return 1;
    }

    Taub_Miles adap_taub_miles(v0, n0, m0, h0, z0, I, Dv, Dn);
    double v, a, eta, n, m , h, Iadap, t;
    double v_tmp;

    double t_since_spike = 0;
    double t_since_print = 0;
    int num_spikes = 0;

    //(while(true){
    while (num_spikes < 100){
        /*if(v_tmp < 45. and v > 45. and t > 2.){
            break;
        }*/

        /*if (t >= t_kick and not kicked){
            adap_taub_miles.kick_voltage(e_kick);
            kicked = true;
        }*/


        v_tmp = adap_taub_miles.get_voltage();

        adap_taub_miles.integrate_with_adap();

        v = adap_taub_miles.get_voltage();
        a = adap_taub_miles.get_adaptation();
        eta = adap_taub_miles.get_cnoise();
        n = adap_taub_miles.get_n();
        m = adap_taub_miles.get_m();
        h = adap_taub_miles.get_h();

        Iadap = adap_taub_miles.get_I_adaptation();
        t = adap_taub_miles.get_time();
        t_since_print += 0.0001;
        t_since_spike += 0.0001;

        if (t >= 25){
            adap_taub_miles.set_current(5.);
        }

        if (v_tmp < 45. and v > 45. and t_since_spike > 2.) {
            num_spikes += 1;
            file_spikes << t_since_spike << "\n";
            t_since_spike = 0;
        }

        if (t_since_print > 0.001 and t < 1000) {
            file << v << " " << a << " " << eta << " " << n << " " << m << " " << h << " " << Iadap << " " << t << "\n";
            t_since_print = 0;
        }
    }
}
