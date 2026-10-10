// simulation.hpp
// Le coeur du projet : un champ de miroirs plans, chacun orienté (imparfaitement)
// vers un point de la coque, et une grille qui compte l'énergie reçue.
#pragma once

#include <random>
#include <vector>

#include "config.hpp"
#include "vec3.hpp"

// Grille de flux posée sur le flanc du navire (plan vertical face au rivage).
// On stocke la puissance reçue par maille, et en parallèle une version "lissée"
// sur 3x3 mailles : le pic de flux est lu sur cette version, ce qui évite qu'une
// seule maille chanceuse (bruit Monte-Carlo) fasse croire à un point chaud.
class GrilleFlux {
public:
    GrilleFlux(double demiLargeur, double pas);

    // Ajoute une puissance (W) au point (y, z) exprimé dans le repère de la cible.
    // Renvoie false si le rayon tombe hors de la grille.
    bool deposer(double y, double z, double puissance);

    double flux(int i, int j) const;        // W/m2, maille brute
    double fluxLisse(int i, int j) const;   // W/m2, moyenne 3x3
    double pic() const { return picLisse_ / (9.0 * pas_ * pas_); }  // W/m2
    double surfaceAuDessus(double seuil) const;  // m2 où le flux lissé dépasse le seuil

    int taille() const { return n_; }
    double pas() const { return pas_; }
    double demiLargeur() const { return demi_; }
    double puissanceTotale() const { return total_; }

private:
    double demi_, pas_;
    int n_;
    std::vector<double> brut_, lisse_;
    double picLisse_ = 0.0;
    double total_ = 0.0;
    int index(int i, int j) const { return j * n_ + i; }
};

struct Miroir {
    Vec3 centre;
    Vec3 normale;  // normale réelle (avec l'erreur de visée du soldat)
    Vec3 u, v;     // bords du miroir
};

class Simulation {
public:
    Simulation(const Config& c, const Vec3& soleil, double dni, unsigned graine);

    // Place le miroir suivant sur le rivage, l'oriente et lance ses rayons.
    void ajouterMiroir();

    int nombreMiroirs() const { return static_cast<int>(miroirs_.size()); }
    const GrilleFlux& grille() const { return grille_; }
    const std::vector<Miroir>& miroirs() const { return miroirs_; }

    double puissanceInterceptee() const { return puissanceMiroirs_; }  // W reçus par les miroirs
    double angleIncidenceMoyen() const;                                 // rad

    // Position du k-ième miroir : rangées remplies de gauche à droite, puis en reculant.
    Vec3 positionMiroir(int k) const;

    Vec3 pointVise() const { return {cfg_.distance, 0.0, cfg_.hauteur_cible}; }

private:
    Config cfg_;
    Vec3 soleil_;
    double dni_;
    GrilleFlux grille_;
    std::vector<Miroir> miroirs_;
    std::mt19937_64 rng_;
    std::uniform_real_distribution<double> uniforme_{0.0, 1.0};
    std::normal_distribution<double> gauss_{0.0, 1.0};
    double puissanceMiroirs_ = 0.0;
    double sommeIncidence_ = 0.0;

    Vec3 tirerDirectionSoleil(const Vec3& u, const Vec3& v);
    Vec3 perturber(const Vec3& n, const Vec3& u, const Vec3& v, double sigma);
};
