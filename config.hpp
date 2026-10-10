// config.hpp
// Tous les paramètres physiques et numériques de la simulation, avec leurs valeurs
// par défaut. Chacun peut être modifié en ligne de commande : ./archimede map D=60 aim=0.5
#pragma once

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

struct Config {
    // ----- Soleil -----
    double latitude_deg = 37.07;     // Syracuse
    double declinaison_deg = 21.0;   // début juillet (été)
    double heure_solaire = 10.0;     // heure solaire vraie (12 = midi solaire)
    double soleil_el_deg = -999.0;   // si renseignés, remplacent le modèle d'éphémérides
    double soleil_az_deg = -999.0;   // azimut compté depuis le nord, vers l'est
    double dni = -1.0;               // éclairement direct (W/m2). <= 0 : calculé (modèle de Meinel)
    double alpha_soleil = 4.65e-3;   // demi-angle apparent du disque solaire (rad)

    // ----- Miroirs (boucliers de bronze poli) -----
    double cote_miroir = 0.8;        // m, miroir carré de 0,64 m2 (ordre de grandeur d'un bouclier)
    double reflectivite = 0.60;      // bronze poli, spectre solaire (hypothèse à discuter)
    double sigma_pente_mrad = 2.0;   // défauts de planéité du polissage (écart-type par axe)
    double sigma_visee_deg = 0.25;   // erreur d'orientation à la main de chaque soldat (par axe)

    // ----- Disposition sur le rivage -----
    int colonnes = 25;               // miroirs par rangée
    double pas_lateral = 1.0;        // m entre centres dans une rangée
    double recul_rangee = 1.2;       // m de recul entre rangées
    double montee_rangee = 0.6;      // m de surélévation entre rangées (gradins, remparts)
    double hauteur_quai = 3.0;       // m, altitude de la première rangée

    // ----- Cible (flanc d'une trière) -----
    double distance = 40.0;          // m entre la première rangée et la coque
    double hauteur_cible = 2.0;      // m, point visé sur la coque

    // ----- Numérique -----
    int n_miroirs = 60;              // nombre de miroirs (mode map)
    int n_max = 500;                 // nombre maximal de miroirs testés (balayages)
    int rayons = 10000;              // rayons par miroir
    double demi_largeur = 5.0;       // m, demi-côté de la grille de flux sur la coque
    double pas_grille = 0.05;        // m, taille d'une maille
    unsigned graine = 1;             // graine du générateur aléatoire
    int graines = 3;                 // nombre de tirages indépendants (balayages)

    // ----- Seuils d'inflammation du bois (flux incident) -----
    double q_crit = 20000.0;         // W/m2, auto-inflammation avec flamme
    double q_incandescence = 10000.0;// W/m2, inflammation par incandescence

    double sigma_visee_rad() const { return sigma_visee_deg * 3.14159265358979323846 / 180.0; }
    double sigma_pente_rad() const { return sigma_pente_mrad * 1e-3; }
};

// Lecture des arguments "cle=valeur". Renvoie false si une clé est inconnue.
inline bool lireArguments(int argc, char** argv, int debut, Config& c) {
    std::map<std::string, double*> reels = {
        {"lat", &c.latitude_deg},      {"decl", &c.declinaison_deg}, {"heure", &c.heure_solaire},
        {"el", &c.soleil_el_deg},      {"az", &c.soleil_az_deg},     {"dni", &c.dni},
        {"cote", &c.cote_miroir},      {"rho", &c.reflectivite},     {"pente", &c.sigma_pente_mrad},
        {"aim", &c.sigma_visee_deg},   {"pas", &c.pas_lateral},      {"recul", &c.recul_rangee},
        {"montee", &c.montee_rangee},  {"quai", &c.hauteur_quai},    {"D", &c.distance},
        {"zc", &c.hauteur_cible},      {"grille", &c.demi_largeur},  {"maille", &c.pas_grille},
        {"qcrit", &c.q_crit}};
    std::map<std::string, int*> entiers = {
        {"N", &c.n_miroirs}, {"Nmax", &c.n_max}, {"rayons", &c.rayons},
        {"colonnes", &c.colonnes}, {"graines", &c.graines}};

    for (int i = debut; i < argc; ++i) {
        std::string arg = argv[i];
        auto pos = arg.find('=');
        if (pos == std::string::npos) {
            std::cerr << "Argument ignore (format cle=valeur attendu) : " << arg << "\n";
            return false;
        }
        std::string cle = arg.substr(0, pos);
        std::string val = arg.substr(pos + 1);
        if (reels.count(cle)) *reels[cle] = std::atof(val.c_str());
        else if (entiers.count(cle)) *entiers[cle] = std::atoi(val.c_str());
        else if (cle == "seed") c.graine = static_cast<unsigned>(std::atoi(val.c_str()));
        else {
            std::cerr << "Cle inconnue : " << cle << "\n";
            return false;
        }
    }
    return true;
}
