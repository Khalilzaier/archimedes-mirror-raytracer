// main.cpp
// Miroirs d'Archimède : simulation par lancer de rayons (Monte-Carlo).
//
// Modes :
//   ./archimede map      [cle=valeur ...]  carte de flux sur la coque pour un scénario
//   ./archimede validate [cle=valeur ...]  un miroir parfait, comparaison avec la théorie
//   ./archimede sweep    [cle=valeur ...]  miroirs nécessaires selon distance et erreur de visée
//   ./archimede day      [cle=valeur ...]  miroirs nécessaires selon l'heure de la journée
//
// Les résultats (CSV) sont écrits dans le dossier results/.

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "config.hpp"
#include "simulation.hpp"
#include "soleil.hpp"

namespace fs = std::filesystem;

// Ajoute des miroirs un par un et renvoie le pic de flux après chaque ajout.
// peaks[k] = pic avec k+1 miroirs.
static std::vector<double> courbePic(const Config& c, const Vec3& s, double dni, unsigned graine) {
    Simulation sim(c, s, dni, graine);
    std::vector<double> peaks;
    peaks.reserve(c.n_max);
    for (int k = 0; k < c.n_max; ++k) {
        sim.ajouterMiroir();
        peaks.push_back(sim.grille().pic());
    }
    return peaks;
}

// Plus petit nombre de miroirs qui atteint le seuil, ou -1 si jamais atteint.
static int miroirsRequis(const std::vector<double>& peaks, double seuil) {
    for (size_t k = 0; k < peaks.size(); ++k)
        if (peaks[k] >= seuil) return static_cast<int>(k) + 1;
    return -1;
}

static void afficherSoleil(const Vec3& s, double dni) {
    double el = deg(std::asin(s.z));
    double az = deg(std::atan2(s.x, s.y));
    if (az < 0) az += 360.0;
    std::printf("Soleil : elevation %.1f deg, azimut %.1f deg, DNI %.0f W/m2\n", el, az, dni);
}

// ---------------------------------------------------------------------
static int modeMap(const Config& c) {
    Vec3 s = directionSoleil(c);
    double dni = eclairementDirect(c, s);
    afficherSoleil(s, dni);

    Simulation sim(c, s, dni, c.graine);
    for (int k = 0; k < c.n_miroirs; ++k) sim.ajouterMiroir();

    const GrilleFlux& g = sim.grille();
    std::ofstream f("results/flux_map.csv");
    f << "# demi_largeur=" << g.demiLargeur() << " pas=" << g.pas() << " N=" << c.n_miroirs
      << " D=" << c.distance << " aim_deg=" << c.sigma_visee_deg << "\n";
    for (int j = 0; j < g.taille(); ++j) {
        for (int i = 0; i < g.taille(); ++i) f << g.flux(i, j) << (i + 1 < g.taille() ? "," : "");
        f << "\n";
    }

    double recu = g.puissanceTotale();
    double renvoye = c.reflectivite * sim.puissanceInterceptee();
    std::printf("Miroirs : %d de %.2f m2, distance %.0f m, erreur de visee %.2f deg\n", c.n_miroirs,
                c.cote_miroir * c.cote_miroir, c.distance, c.sigma_visee_deg);
    std::printf("Angle d'incidence moyen sur les miroirs : %.1f deg\n", deg(sim.angleIncidenceMoyen()));
    std::printf("Puissance renvoyee par les miroirs : %.1f kW\n", renvoye / 1000.0);
    std::printf("Puissance tombant sur la grille de %.0f x %.0f m : %.1f kW (%.0f %%)\n",
                2 * g.demiLargeur(), 2 * g.demiLargeur(), recu / 1000.0, 100.0 * recu / renvoye);
    std::printf("Pic de flux (moyenne sur %.0f cm x %.0f cm) : %.1f kW/m2\n", 300 * g.pas(), 300 * g.pas(),
                g.pic() / 1000.0);
    std::printf("Surface au-dessus de %.0f kW/m2 : %.2f m2\n", c.q_crit / 1000, g.surfaceAuDessus(c.q_crit));
    std::printf("Surface au-dessus de %.0f kW/m2 : %.2f m2\n", c.q_incandescence / 1000,
                g.surfaceAuDessus(c.q_incandescence));
    std::printf("-> results/flux_map.csv\n");
    return 0;
}

// ---------------------------------------------------------------------
// Validation : un seul miroir parfait (aucune erreur de visée ni de polissage), soleil
// placé dans le plan vertical qui contient le miroir et la cible. Le flux au centre de la
// tache doit valoir rho * DNI * cos(gamma) * F, où F est la fraction du disque solaire
// visible à travers le miroir depuis la cible (calcul analytique dans le script Python).
static int modeValidate(Config c) {
    c.colonnes = 1;
    c.sigma_visee_deg = 0.0;
    c.sigma_pente_mrad = 0.0;
    if (c.soleil_el_deg < -90) c.soleil_el_deg = 40.0;
    if (c.soleil_az_deg < -900) c.soleil_az_deg = 90.0;
    if (c.dni <= 0) c.dni = 900.0;
    if (c.rayons == 10000) c.rayons = 10000000;
    c.pas_grille = 0.01;  // fenêtre de lissage de 3 cm : on veut le flux en un point

    Vec3 s = directionSoleil(c);
    afficherSoleil(s, c.dni);

    std::ofstream f("results/validation.csv");
    f << "D,flux_centre,theta_inc_rad,cos_gamma,cote,alpha,rho,dni\n";
    for (double D : {5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 120.0, 140.0,
                     160.0, 200.0, 250.0, 300.0}) {
        c.distance = D;
        Simulation sim(c, s, c.dni, c.graine);
        sim.ajouterMiroir();
        const GrilleFlux& g = sim.grille();
        int m = g.taille() / 2;
        // Le point visé est au coin commun de 4 mailles : on moyenne leurs valeurs lissées.
        double centre = 0.25 * (g.fluxLisse(m - 1, m - 1) + g.fluxLisse(m, m - 1) +
                                g.fluxLisse(m - 1, m) + g.fluxLisse(m, m));
        const Miroir& mi = sim.miroirs()[0];
        Vec3 versCible = normaliser(sim.pointVise() - mi.centre);
        double theta = std::acos(dot(mi.normale, s));
        f << D << "," << centre << "," << theta << "," << versCible.x << "," << c.cote_miroir << ","
          << c.alpha_soleil << "," << c.reflectivite << "," << c.dni << "\n";
        std::printf("D = %5.0f m   flux au centre = %7.1f W/m2\n", D, centre);
    }
    std::printf("-> results/validation.csv\n");
    return 0;
}

// ---------------------------------------------------------------------
static int modeSweep(Config c) {
    Vec3 s = directionSoleil(c);
    double dni = eclairementDirect(c, s);
    afficherSoleil(s, dni);

    std::vector<double> visees = {0.0, 0.25, 0.5, 1.0};
    std::vector<double> distances;
    for (double D = 20.0; D <= 150.0 + 1e-9; D += 10.0) distances.push_back(D);

    struct Tache { double aim, D; int graine; };
    std::vector<Tache> taches;
    for (double a : visees)
        for (double D : distances)
            for (int g = 0; g < c.graines; ++g) taches.push_back({a, D, g});

    std::vector<std::vector<double>> resultats(taches.size());

    // Chaque tâche est indépendante : on peut les répartir sur tous les coeurs.
#pragma omp parallel for schedule(dynamic)
    for (int t = 0; t < static_cast<int>(taches.size()); ++t) {
        Config ct = c;
        ct.sigma_visee_deg = taches[t].aim;
        ct.distance = taches[t].D;
        resultats[t] = courbePic(ct, s, dni, 1000u + 7919u * taches[t].graine);
    }

    std::ofstream fr("results/sweep_required.csv");
    std::ofstream fp("results/sweep_peak.csv");
    fr << "aim_deg,D,seed,N_req_flamme,N_req_incandescence\n";
    fp << "aim_deg,D,seed,N,peak\n";
    for (size_t t = 0; t < taches.size(); ++t) {
        const auto& pk = resultats[t];
        fr << taches[t].aim << "," << taches[t].D << "," << taches[t].graine << ","
           << miroirsRequis(pk, c.q_crit) << "," << miroirsRequis(pk, c.q_incandescence) << "\n";
        for (size_t k = 0; k < pk.size(); ++k)
            fp << taches[t].aim << "," << taches[t].D << "," << taches[t].graine << "," << k + 1 << ","
               << pk[k] << "\n";
    }
    std::printf("%zu simulations, jusqu'a %d miroirs chacune -> results/sweep_*.csv\n", taches.size(),
                c.n_max);
    return 0;
}

// ---------------------------------------------------------------------
static int modeDay(Config c) {
    std::vector<double> heures;
    for (double h = 5.5; h <= 18.5 + 1e-9; h += 0.5) heures.push_back(h);

    struct Tache { double h; int graine; };
    std::vector<Tache> taches;
    for (double h : heures)
        for (int g = 0; g < c.graines; ++g) taches.push_back({h, g});

    struct Res { double el, az, dni, theta; int nreq; double pic100; };
    std::vector<Res> res(taches.size());

#pragma omp parallel for schedule(dynamic)
    for (int t = 0; t < static_cast<int>(taches.size()); ++t) {
        Config ct = c;
        ct.heure_solaire = taches[t].h;
        Vec3 s = directionSoleil(ct);
        double dni = eclairementDirect(ct, s);
        Res r{deg(std::asin(s.z)), deg(std::atan2(s.x, s.y)), dni, 0.0, -1, 0.0};
        if (r.az < 0) r.az += 360.0;
        if (r.el > 3.0) {
            Simulation sim(ct, s, dni, 1000u + 7919u * taches[t].graine);
            for (int k = 0; k < ct.n_max; ++k) {
                sim.ajouterMiroir();
                double p = sim.grille().pic();
                if (r.nreq < 0 && p >= ct.q_crit) r.nreq = k + 1;
                if (k + 1 == 100) { r.pic100 = p; r.theta = deg(sim.angleIncidenceMoyen()); }
            }
        }
        res[t] = r;
    }

    std::ofstream f("results/day.csv");
    f << "heure,seed,elevation_deg,azimut_deg,dni,theta_moyen_deg,N_req_flamme,pic_N100\n";
    for (size_t t = 0; t < taches.size(); ++t) {
        const Res& r = res[t];
        f << taches[t].h << "," << taches[t].graine << "," << r.el << "," << r.az << "," << r.dni << ","
          << r.theta << "," << r.nreq << "," << r.pic100 << "\n";
    }
    std::printf("Journee simulee (D = %.0f m, visee %.2f deg) -> results/day.csv\n", c.distance,
                c.sigma_visee_deg);
    return 0;
}

// ---------------------------------------------------------------------
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage : " << argv[0] << " map|validate|sweep|day [cle=valeur ...]\n";
        return 1;
    }
    std::string mode = argv[1];
    Config c;
    if (!lireArguments(argc, argv, 2, c)) return 1;
    fs::create_directories("results");

    if (mode == "map") return modeMap(c);
    if (mode == "validate") return modeValidate(c);
    if (mode == "sweep") return modeSweep(c);
    if (mode == "day") return modeDay(c);
    std::cerr << "Mode inconnu : " << mode << "\n";
    return 1;
}
